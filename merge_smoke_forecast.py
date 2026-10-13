#!/usr/bin/env python3
"""Join NOAA smoke altitude, column and aerosol time series with CAMS context.

Reads ONLY existing JSON under docs/. No NOAA downloads, new dependencies, or
changes to the source products. AGL values are never interpolated in time.
A strict UTC valid-time + NOAA model-cycle match is required to join data.
"""
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from zoneinfo import ZoneInfo

UTC = timezone.utc
LOCAL = ZoneInfo("America/Chicago")
ROOT = Path(__file__).resolve().parent
DOCS = ROOT / "docs"
OUTPUT = DOCS / "smoke_forecast.json"


def iso_utc(s: str) -> datetime:
    if not isinstance(s, str):
        raise ValueError("timestamp must be a string")
    dt = datetime.fromisoformat(s.replace("Z", "+00:00"))
    if dt.tzinfo is None:
        raise ValueError("timestamp must have timezone")
    return dt.astimezone(UTC)


def z(t: datetime) -> str:
    return t.astimezone(UTC).isoformat().replace("+00:00", "Z")


def numeric(v) -> bool:
    return (isinstance(v, (float, int)) and not isinstance(v, bool)
            and math.isfinite(v) and v >= 0)


def read_json(name: str) -> tuple[dict | None, str | None]:
    try:
        with (DOCS / name).open(encoding="utf-8") as file:
            obj = json.load(file)
        if not isinstance(obj, dict):
            raise ValueError("JSON root is not an object")
        return obj, None
    except (OSError, ValueError) as exc:
        return None, f"{name}: {type(exc).__name__}: {exc}"[:300]


def valid_sample(sample: dict) -> bool:
    """Confirm the reported local clock represents the actual UTC valid time."""
    if not isinstance(sample, dict):
        return False
    try:
        return (iso_utc(sample["valid_utc"]) ==
                iso_utc(sample["valid_local"]))
    except (KeyError, ValueError, TypeError):
        return False


def matching_site(a: dict, b: dict, max_degrees: float = 0.03) -> bool:
    try:
        pa, pb = a["location"], b["location"]
        return all(abs(float(pa[k]) - float(pb[k])) < max_degrees
                   for k in ("latitude", "longitude"))
    except (KeyError, TypeError, ValueError):
        return False


def cycle_agrees(column: dict, vertical: dict) -> bool:
    try:
        return (iso_utc(column["model_cycle_utc"]) ==
                iso_utc(vertical["model_cycle_utc"]))
    except (KeyError, TypeError, ValueError):
        return False


def forecast_matches(sample: dict, cycle_string: str) -> bool:
    try:
        lead = (iso_utc(sample["valid_utc"]) - iso_utc(cycle_string)).total_seconds() / 3600
        hour = sample["forecast_hour"]
        return lead.is_integer() and 1 <= lead <= 48 and int(hour) == int(lead)
    except (KeyError, TypeError, ValueError):
        return False


def regional_values(column: dict, vertical: dict | None) -> list[dict]:
    v_by_label = {r.get("label"): r for r in (vertical or {}).get("region", [])
                  if isinstance(r, dict)}
    rows = []
    for r in column.get("region", []):
        if not isinstance(r, dict):
            continue
        vr = v_by_label.get(r.get("label"), {})
        coords_ok = (vr and abs(float(vr.get("latitude", 999)) - float(r.get("latitude", 0))) < 0.03
                     and abs(float(vr.get("longitude", 999)) - float(r.get("longitude", 0))) < 0.03)
        rows.append({
            "label": r.get("label"),
            "latitude": r.get("latitude"),
            "longitude": r.get("longitude"),
            "smoke_column_mg_m2": r.get("smoke_column_mg_m2") if numeric(r.get("smoke_column_mg_m2")) else None,
            "hrrr_aod": r.get("hrrr_aod") if numeric(r.get("hrrr_aod")) else None,
            "smoke_1000ft_agl_ug_m3": vr.get("smoke_1000ft_agl_ug_m3") if coords_ok and numeric(vr.get("smoke_1000ft_agl_ug_m3")) else None,
            "smoke_6000ft_agl_ug_m3": vr.get("smoke_6000ft_agl_ug_m3") if coords_ok and numeric(vr.get("smoke_6000ft_agl_ug_m3")) else None,
        })
    return rows


def merge(column: dict | None, vertical: dict | None, cams: dict | None,
          now: datetime | None = None) -> dict:
    now = now or datetime.now(UTC)
    column = column or {}
    vertical = vertical or {}
    cams = cams or {}
    warnings = []
    matching_noaa = (cycle_agrees(column, vertical)
                     and matching_site(column, vertical))
    if vertical and not matching_noaa:
        warnings.append("NOAA vertical and column model cycles or locations differ; vertical results were not joined")

    v_by_time = {}
    if matching_noaa:
        for s in vertical.get("samples", []):
            if (valid_sample(s) and s.get("status") == "available"
                    and forecast_matches(s, vertical["model_cycle_utc"])):
                v_by_time[s["valid_utc"]] = s
    cam_by_time = {}
    if cams and matching_site(column, cams):
        for h in cams.get("hours", []):
            if valid_sample(h) and numeric(h.get("site", {}).get("aod")):
                cam_by_time[z(iso_utc(h["valid_utc"]))] = h
    elif cams:
        warnings.append("CAMS location does not match NOAA site; CAMS values were not joined")

    data = []
    for c in column.get("samples", []):
        if not valid_sample(c) or not forecast_matches(c, column.get("model_cycle_utc", "")):
            warnings.append("Skipped a NOAA column sample with invalid timestamps or forecast lead")
            continue
        utc_str = z(iso_utc(c["valid_utc"]))
        v = v_by_time.get(utc_str)
        if v and int(v.get("forecast_hour", -1)) != int(c.get("forecast_hour", -2)):
            v = None
            warnings.append(f"NOAA forecast-hour mismatch at {utc_str}")
        valid_c = c.get("status") == "available"
        a = v.get("smoke_1000ft_agl_ug_m3") if v else None
        b = v.get("smoke_6000ft_agl_ug_m3") if v else None
        column_smoke = c.get("smoke_column_mg_m2") if valid_c else None
        hrrr_aod = c.get("hrrr_aod") if valid_c else None
        cam = cam_by_time.get(utc_str)
        layer_ok = v is not None and numeric(a) and numeric(b)
        statuses = {
            "smoke_column_mg_m2": "available" if numeric(column_smoke) else "unavailable",
            "hrrr_aod": "available" if numeric(hrrr_aod) else "unavailable",
            "smoke_1000ft_agl_ug_m3": "available" if layer_ok else "not_sampled_or_unavailable",
            "smoke_6000ft_agl_ug_m3": "available" if layer_ok else "not_sampled_or_unavailable",
            "cams_total_aod": "available" if cam else "unavailable",
        }
        data.append({
            "valid_utc": utc_str,
            "valid_local": iso_utc(utc_str).astimezone(LOCAL).isoformat(),
            "forecast_hour": c.get("forecast_hour"),
            "measurement_status": statuses,
            "smoke_1000ft_agl_ug_m3": a if layer_ok else None,
            "smoke_6000ft_agl_ug_m3": b if layer_ok else None,
            "smoke_column_mg_m2": column_smoke if numeric(column_smoke) else None,
            "hrrr_aod": hrrr_aod if numeric(hrrr_aod) else None,
            "cams_total_aod": cam["site"]["aod"] if cam else None,
            "region": regional_values(c, v if layer_ok else None),
        })
    data.sort(key=lambda s: s["valid_utc"])
    n_vertical = sum(s["measurement_status"]["smoke_1000ft_agl_ug_m3"] == "available" for s in data)
    n_column = sum(s["measurement_status"]["smoke_column_mg_m2"] == "available" for s in data)
    status = ("available" if data and n_vertical == n_column == len(data)
              else "partial" if n_column or n_vertical else "unavailable")
    return {
        "schema_version": "1.0",
        "generated_utc": z(now),
        "status": status,
        "location": column.get("location", vertical.get("location")),
        "timezone": "America/Chicago",
        "noaa_model_cycle_utc": column.get("model_cycle_utc"),
        "source_generated_utc": {
            "noaa_column": column.get("generated_utc"),
            "noaa_vertical": vertical.get("generated_utc"),
            "cams": cams.get("generated_utc"),
        },
        "source_files": ["noaa_timeseries.json", "noaa_vertical_timeseries.json", "latest.json"],
        "summary": {
            "noaa_column_valid_hours": n_column,
            "noaa_vertical_valid_hours": n_vertical,
            "noaa_column_sample_count": len(data),
            "cams_aod_valid_hours": sum(s["measurement_status"]["cams_total_aod"] == "available" for s in data),
            "all_four_noaa_fields_available_at": [s["valid_local"] for s in data
                                                   if all(s["measurement_status"][k] == "available" for k in (
                                                       "smoke_column_mg_m2", "hrrr_aod",
                                                       "smoke_1000ft_agl_ug_m3", "smoke_6000ft_agl_ug_m3"))],
        },
        "notes": [
            "Three native-level vertical times are joined to hourly column smoke and AOTK; no time interpolation for missing AGL hours.",
            "CAMS AOD represents all aerosol species and is NOT interchangeable with NOAA HRRR AOTK; NOAA AOTK wavelength remains unverified.",
            "NOAA MASSDEN kg/m^3 was converted to micrograms/m^3 by 1e9 in source extractor.",
            "The JSON may become outdated; downstream consumers must check valid times and source model cycle before interpreting it as tonight's forecast.",
            "Status partial means four NOAA fields are not available every hour; it does not invalidate completed 3-time profiles.",
        ],
        "warnings": warnings,
        "samples": data,
    }


def main() -> None:
    column, err1 = read_json("noaa_timeseries.json")
    vertical, err2 = read_json("noaa_vertical_timeseries.json")
    cams, err3 = read_json("latest.json")
    result = merge(column, vertical, cams)
    result["warnings"].extend(x for x in (err1, err2, err3) if x)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    s = result["summary"]
    print("Merged smoke forecast:", result["status"])
    print("Model cycle:", result["noaa_model_cycle_utc"])
    print("Valid NOAA column hours:", s["noaa_column_valid_hours"])
    print("Valid NOAA vertical hours:", s["noaa_vertical_valid_hours"])
    print("CAMS AOD hours:", s["cams_aod_valid_hours"])
    for sample in result["samples"]:
        if sample["measurement_status"]["smoke_1000ft_agl_ug_m3"] == "available":
            print(sample["valid_local"],
                  "1000ft=", sample["smoke_1000ft_agl_ug_m3"],
                  "6000ft=", sample["smoke_6000ft_agl_ug_m3"],
                  "column=", sample["smoke_column_mg_m2"],
                  "HRRR AOD=", sample["hrrr_aod"])
    if result["warnings"]:
        print("Warnings:", "; ".join(result["warnings"]))
    print("Saved", OUTPUT)


if __name__ == "__main__":
    main()
