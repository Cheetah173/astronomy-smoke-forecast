#!/usr/bin/env python3
"""Three-time regional NOAA HRRR smoke/AOD diagnostic.

Uses the *validated* GRIB field extraction in noaa_probe.py; leaves latest.json
and noaa_probe.json unchanged. NOAA NOMADS requests are paced 12 seconds apart.
Results are forecasts, not measured extinction or altitude-resolved smoke.
"""
from __future__ import annotations

import json
import os
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from zoneinfo import ZoneInfo

import requests

from noaa_probe import URL, build_params, decode, select_cycle

UTC = timezone.utc
LOCAL = ZoneInfo("America/Chicago")
OUTPUT = Path(__file__).resolve().parent / "docs" / "noaa_timeseries.json"
PAUSE_SECONDS = 12
MAX_DOWNLOAD_BYTES = 15_000_000


def to_utc_z(moment: datetime) -> str:
    return moment.astimezone(UTC).isoformat().replace("+00:00", "Z")


def target_hours(now: datetime) -> list[datetime]:
    """Sample the observing night containing now or the next coming evening.

    At 00:00–05:59 local, the observing night started the prior evening.
    Otherwise, evaluate the evening of the local calendar day.
    """
    local = now.astimezone(LOCAL)
    day = local.date() - timedelta(days=1 if local.hour < 6 else 0)
    base = datetime(day.year, day.month, day.day, 20, tzinfo=LOCAL)
    return [(base + timedelta(hours=offset)).astimezone(UTC)
            for offset in (0, 5, 8)]  # 8 PM, 1 AM, 4 AM CDT/CST


def field_matches_cycle(candidate_fields: list[dict], cycle: datetime, fxx: int) -> bool:
    """Do not accept a response carrying unexpected run/hour metadata."""
    if len(candidate_fields) != 2:
        return False
    for item in candidate_fields:
        meta = item.get("metadata", {})
        try:
            date = int(meta["dataDate"])
            hm = int(meta["dataTime"])
            fc = int(meta["forecastTime"])
            step = int(meta["stepUnits"])
        except (ValueError, TypeError, KeyError):
            return False
        if date != int(cycle.strftime("%Y%m%d")) or hm != int(cycle.strftime("%H%M")):
            return False
        if fc != fxx or step != 1:
            return False
    return True


def site_from_region(items: list[dict], field: str) -> float:
    return float(next(p[field] for p in items if p["label"] == "SITE"))


def extract_one(session: requests.Session, cycle: datetime, valid: datetime) -> dict:
    forecast_hours = (valid - cycle).total_seconds() / 3600
    result = {
        "valid_utc": to_utc_z(valid),
        "valid_local": valid.astimezone(LOCAL).isoformat(),
        "forecast_hour": forecast_hours,
        "status": "unavailable",
        "smoke_column_mg_m2": None,
        "hrrr_aod": None,
        "region": [],
    }
    if not forecast_hours.is_integer() or not 1 <= forecast_hours <= 48:
        result["reason"] = "Forecast valid time falls outside supported HRRR window"
        return result
    fxx = int(forecast_hours)
    result["forecast_hour"] = fxx
    try:
        response = session.get(URL, params=build_params(cycle, fxx), timeout=90)
        response.raise_for_status()
        raw = response.content
        print(f"F{fxx:02d} NOAA response bytes: {len(raw)}", flush=True)
        if not raw.startswith(b"GRIB") or len(raw) > MAX_DOWNLOAD_BYTES:
            raise ValueError("NOAA response not valid GRIB2 or outside download size limit")
        parsed = decode(raw)
        if not field_matches_cycle(parsed.get("candidate_fields", []), cycle, fxx):
            raise ValueError("NOAA GRIB metadata did not match requested cycle/hour")
        if "column_mg_m2" not in parsed or "hrrr_aod" not in parsed:
            raise ValueError("Smoke column or AOD missing/ambiguous in decoded NOAA data")
        columns = {p["label"]: p for p in parsed["column_mg_m2"]}
        aods = {p["label"]: p for p in parsed["hrrr_aod"]}
        if set(columns) != set(aods) or len(columns) != 9:
            raise ValueError("Incomplete or inconsistent NOAA regional grid sampling")
        result["region"] = [{
            "label": label,
            "latitude": columns[label]["latitude"],
            "longitude": columns[label]["longitude"],
            "smoke_column_mg_m2": columns[label]["smoke_column_mg_m2"],
            "hrrr_aod": aods[label]["aod"],
        } for label in ("NW", "N", "NE", "W", "SITE", "E", "SW", "S", "SE")]
        result["smoke_column_mg_m2"] = site_from_region(result["region"], "smoke_column_mg_m2")
        result["hrrr_aod"] = site_from_region(result["region"], "hrrr_aod")
        result["status"] = "available"
    except Exception as exc:
        result["reason"] = f"{type(exc).__name__}: {exc}"[:400]
        print(f"F{fxx:02d} unavailable: {result['reason']}", flush=True)
    return result


def build_report(now: datetime, cycle: datetime, samples: list[dict]) -> dict:
    return {
        "schema_version": "1.0",
        "generated_utc": to_utc_z(now),
        "timezone": "America/Chicago",
        "location": {"name": "Robinson TX approximate", "latitude": 31.45, "longitude": -97.14},
        "source": "NOAA NCEP HRRR CONUS NOMADS filtered GRIB2",
        "source_url": URL,
        "model_cycle_utc": to_utc_z(cycle),
        "measurement_status": {
            "smoke_column_mg_m2": "NOAA HRRR COLMD - see each sample status",
            "aerosol_optical_depth": "NOAA HRRR AOTK - wavelength unverified",
            "smoke_1000ft_agl_ug_m3": "unavailable - altitude decoder not connected",
            "smoke_6000ft_agl_ug_m3": "unavailable - altitude decoder not connected",
        },
        "notes": [
            "Only three diagnostic forecast times are sampled (8 PM, 1 AM, 4 AM).",
            "Column smoke is from NOAA HRRR documented COLMD product; GRIB record lacks smoke species tag.",
            "AOD is total aerosol extinction, not wildfire-smoke-only; HRRR wavelength unverified.",
            "Nine regional sample points approximate trends but are not a continuous smoke map.",
            "Missing or forecast-hour mismatched values remain null, never assumed zero.",
        ],
        "samples": samples,
    }


def main() -> None:
    now = datetime.now(UTC)
    cycle, _ = select_cycle(now)
    times = target_hours(now)
    report = {"status": "not-run"}
    try:
        session = requests.Session()
        session.headers["User-Agent"] = "AstronomySmokeForecast/0.3 (personal research)"
        samples = []
        requested = 0
        for valid in times:
            fxx = int((valid - cycle).total_seconds() / 3600)
            if 1 <= fxx <= 48:
                if requested:
                    time.sleep(PAUSE_SECONDS)
                requested += 1
            sample = extract_one(session, cycle, valid)
            samples.append(sample)
            print(f"{sample['valid_local']} -> {sample['status']}: smoke={sample['smoke_column_mg_m2']} mg/m2; AOD={sample['hrrr_aod']}", flush=True)
        report = build_report(now, cycle, samples)
    except Exception as exc:
        report = build_report(now, cycle, [])
        report["error"] = f"{type(exc).__name__}: {exc}"[:400]
        print("NOAA timeseries process error:", report["error"])
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("Saved", OUTPUT)


if __name__ == "__main__":
    main()
