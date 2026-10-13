#!/usr/bin/env python3
"""NOAA HRRR native-level smoke at 1000/6000 ft AGL, three observing-night times.

Diagnostic companion to noaa_vertical_probe.py. It uses the SAME model cycle
selection and observing-night timestamps as noaa_timeseries.py. It does not
modify the existing NOAA column-smoke or Open-Meteo results.

The model forecast values are NOT measured aerosol extinction. Native HRRR
hybrid-level MASSDEN (kg/m^3) is converted to ug/m^3 by multiplying by 1e9.
For each site, HGT at the surface is subtracted from hybrid-level HGT before
linear interpolation to 304.8 and 1828.8 metres above local terrain.
"""
from __future__ import annotations

import json
import math
import tempfile
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from zoneinfo import ZoneInfo

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from noaa_probe import POINTS, select_cycle
from noaa_vertical_probe import (
    COARSE, MAX_INDEX, MAX_MESSAGE, TARGETS, URL_PREFIX,
    bracket, find_entry, index_entries, interp, utc,
)

OUTPUT = Path(__file__).resolve().parent / "docs" / "noaa_vertical_timeseries.json"
LOCAL = ZoneInfo("America/Chicago")
PAUSE_BETWEEN_FORECAST_HOURS = 10
TIMES_OF_NIGHT = (0, 5, 8)  # 8 PM, 1 AM, 4 AM local.


def observing_times(now: datetime) -> list[datetime]:
    """Match the established column-smoke calendar-night selection rule."""
    local = now.astimezone(LOCAL)
    day = local.date() - timedelta(days=1 if local.hour < 6 else 0)
    base = datetime(day.year, day.month, day.day, 20, tzinfo=LOCAL)
    return [(base + timedelta(hours=hour)).astimezone(timezone.utc)
            for hour in TIMES_OF_NIGHT]


def make_session() -> requests.Session:
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=Retry(
        total=3, connect=3, read=3, status=3,
        backoff_factor=1.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=frozenset({"GET"}),
        respect_retry_after_header=True,
    )))
    session.headers["User-Agent"] = "AstronomySmokeForecast/0.6 (noncommercial vertical profile)"
    return session


def sample_one(session: requests.Session, cycle: datetime, valid: datetime,
               grid_reference: dict) -> dict:
    from eccodes import (
        codes_get, codes_get_double_elements, codes_grib_find_nearest_multiple,
        codes_grib_new_from_file, codes_release,
    )

    lead_float = (valid - cycle).total_seconds() / 3600.0
    sample = {
        "valid_utc": utc(valid),
        "valid_local": valid.astimezone(LOCAL).isoformat(),
        "forecast_hour": int(lead_float) if lead_float.is_integer() else lead_float,
        "status": "unavailable",
        "smoke_1000ft_agl_ug_m3": None,
        "smoke_6000ft_agl_ug_m3": None,
        "region": [],
    }
    if not lead_float.is_integer() or not (1 <= lead_float <= 48):
        sample["reason"] = "Outside native HRRR supported 1–48 forecast-hour range"
        return sample
    lead = int(lead_float)
    name = f"hrrr.t{cycle:%H}z.wrfnatf{lead:02d}.grib2"
    file_url = f"{URL_PREFIX}/hrrr.{cycle:%Y%m%d}/conus/{name}"
    sample["source_url"] = file_url
    bytes_read = 0
    downloaded_messages = 0
    file_cache = {}

    def fetch(entry: dict) -> dict[str, float]:
        nonlocal bytes_read, downloaded_messages
        key = (entry["name"], entry["level"], entry["surface"])
        if key in file_cache:
            return file_cache[key]
        start, end = entry["start"], entry["end"]
        length = end - start + 1
        if not 0 < length <= MAX_MESSAGE:
            raise ValueError(f"Unexpected GRIB message size for {key}: {length}")
        response = session.get(
            file_url, headers={"Range": f"bytes={start}-{end}"}, timeout=90)
        response.raise_for_status()
        if (response.status_code != 206 or len(response.content) != length
                or not response.content.startswith(b"GRIB")):
            raise ValueError(f"Invalid range response for {key}: "
                             f"status={response.status_code}, size={len(response.content)}")
        with tempfile.TemporaryFile(mode="w+b") as tmp:
            tmp.write(response.content)
            tmp.seek(0)
            gid = codes_grib_new_from_file(tmp)
            if gid is None:
                raise ValueError(f"GRIB decoder returned no field for {key}")
            try:
                if (int(codes_get(gid, "dataDate")) != int(cycle.strftime("%Y%m%d"))
                        or int(codes_get(gid, "dataTime")) != int(cycle.strftime("%H%M"))
                        or int(codes_get(gid, "forecastTime")) != lead):
                    raise ValueError(f"Model timestamp mismatch for {key}")
                signature = str(codes_get(gid, "md5GridSection"))
                if not signature or signature in ("missing", "undef"):
                    raise ValueError(f"Missing grid checksum for {key}")
                if grid_reference["signature"] is None:
                    nearest = codes_grib_find_nearest_multiple(
                        gid, False,
                        [lat for _, lat, _ in POINTS],
                        [lon for _, _, lon in POINTS],
                    )
                    if len(nearest) != len(POINTS):
                        raise ValueError("Incomplete first regional grid lookup")
                    indices = [int(hit.index) for hit in nearest]
                    values = codes_get_double_elements(gid, "values", indices)
                    for hit, v in zip(nearest, values):
                        if not math.isclose(float(hit.value), float(v),
                                            rel_tol=1e-7, abs_tol=1e-12):
                            raise ValueError("Cached GRIB cell index verification failed")
                    grid_reference["signature"] = signature
                    grid_reference["indices"] = indices
                    print("Grid indexes cached and verified", flush=True)
                else:
                    if signature != grid_reference["signature"]:
                        raise ValueError("NOAA grid geometry changed; refusing cached index reuse")
                    values = codes_get_double_elements(
                        gid, "values", grid_reference["indices"])
                if len(values) != len(POINTS):
                    raise ValueError("Returned value count differs from location count")
                result = {}
                for (label, _, _), value in zip(POINTS, values):
                    v = float(value)
                    if not math.isfinite(v):
                        raise ValueError(f"Non-finite decoded field at {label}: {key}")
                    result[label] = v
            finally:
                codes_release(gid)
        file_cache[key] = result
        bytes_read += length
        downloaded_messages += 1
        return result

    try:
        index_response = session.get(file_url + ".idx", timeout=45)
        index_response.raise_for_status()
        if len(index_response.content) > MAX_INDEX:
            raise ValueError("Unexpectedly large native-level index")
        entries = index_entries(index_response.text)
        hgt_levels = {e["level"] for e in entries
                      if e["name"] == "HGT" and e["level"] is not None}
        mass_levels = {e["level"] for e in entries
                       if e["name"] == "MASSDEN" and e["level"] is not None}
        levels = hgt_levels & mass_levels
        if len(levels) < 45:
            raise ValueError(f"Only {len(levels)} matching HGT/MASSDEN hybrid levels")
        terrain = fetch(find_entry(entries, "HGT", surface=True))
        coarse = [level for level in COARSE if level in levels]
        for level in coarse:
            fetch(find_entry(entries, "HGT", level=level))
        refined = set()
        for label, _, _ in POINTS:
            heights = {level: fetch(find_entry(entries, "HGT", level=level))[label]
                       - terrain[label] for level in coarse}
            for target in TARGETS.values():
                pair = bracket(heights, target)
                if pair is None:
                    raise ValueError(f"Coarse levels do not bracket {target}m AGL at {label}")
                refined.update(level for level in levels
                               if min(pair) <= level <= max(pair))
        for level in sorted(refined):
            fetch(find_entry(entries, "HGT", level=level))
            fetch(find_entry(entries, "MASSDEN", level=level))
        for label, lat, lon in POINTS:
            heights = {
                level: fetch(find_entry(entries, "HGT", level=level))[label]
                - terrain[label] for level in refined
            }
            concentrations = {
                level: fetch(find_entry(entries, "MASSDEN", level=level))[label] * 1e9
                for level in refined
            }
            row = {
                "label": label, "latitude": lat, "longitude": lon,
                "surface_height_m": round(terrain[label], 2),
            }
            for field, target in TARGETS.items():
                pair = bracket(heights, target)
                if pair is None:
                    row[field] = None
                    row[field + "_reason"] = "No bracketing hybrid heights"
                    continue
                lo, hi = sorted(pair, key=lambda i: heights[i])
                h0, h1 = heights[lo], heights[hi]
                v0, v1 = concentrations[lo], concentrations[hi]
                if (abs(hi - lo) != 1 or h1 - h0 > 750 or h0 < -100
                        or min(v0, v1) < -0.001):
                    row[field] = None
                    row[field + "_reason"] = "Missing/implausible adjacent-level data"
                    continue
                row[field] = round(interp((h0, max(v0, 0.0)),
                                          (h1, max(v1, 0.0)), target), 4)
                row[field + "_levels"] = [lo, hi]
                row[field + "_height_brackets_m_agl"] = [round(h0, 1), round(h1, 1)]
            sample["region"].append(row)
        site = next(x for x in sample["region"] if x["label"] == "SITE")
        for field in TARGETS:
            sample[field] = site[field]
        sample["status"] = "available" if all(site[field] is not None
                                              for field in TARGETS) else "partial"
        sample["native_levels_used"] = sorted(refined)
    except Exception as exc:
        sample["reason"] = f"{type(exc).__name__}: {exc}"[:500]
    finally:
        sample["download_stats"] = {"messages": downloaded_messages, "bytes": bytes_read}
    return sample


def report(now: datetime, cycle: datetime, samples: list[dict]) -> dict:
    return {
        "schema_version": "1.0",
        "generated_utc": utc(now),
        "timezone": "America/Chicago",
        "location": {"name": "Robinson TX approximate", "latitude": 31.45,
                     "longitude": -97.14},
        "source": "NOAA NCEP HRRR native-level wrfnat MASSDEN and HGT",
        "model_cycle_utc": utc(cycle),
        "measurement_status": {
            "smoke_1000ft_agl_ug_m3": "See each hourly sample status (NOAA hybrid MASSDEN)",
            "smoke_6000ft_agl_ug_m3": "See each hourly sample status (NOAA hybrid MASSDEN)",
        },
        "notes": [
            "Three diagnostic overnight local times: 20:00, 01:00, 04:00; NOT hourly coverage.",
            "Smoke concentration = native HRRR MASSDEN kg/m3 * 1e9 to ug/m3.",
            "AGL elevations are model hybrid HGT minus model surface HGT.",
            "Forecast model validity is per-sample; do not use a past run as tonight's forecast.",
            "Grid-cell indices verified and cached across files only if grid checksum unchanged.",
            "Column smoke and AOTK remain in independent noaa_timeseries.json.",
        ],
        "samples": samples,
    }


def main() -> None:
    now = datetime.now(timezone.utc)
    cycle, _ = select_cycle(now)
    times = observing_times(now)
    session = make_session()
    grid_reference = {"signature": None, "indices": None}
    samples = []
    for i, valid in enumerate(times):
        if i:
            time.sleep(PAUSE_BETWEEN_FORECAST_HOURS)
        started = time.monotonic()
        try:
            item = sample_one(session, cycle, valid, grid_reference)
        except Exception as exc:
            item = {
                "valid_utc": utc(valid), "valid_local": valid.astimezone(LOCAL).isoformat(),
                "status": "unavailable", "smoke_1000ft_agl_ug_m3": None,
                "smoke_6000ft_agl_ug_m3": None, "region": [],
                "reason": f"{type(exc).__name__}: {exc}"[:500],
            }
        samples.append(item)
        print(
            f"{item['valid_local']} => {item['status']} "
            f"1000ft={item['smoke_1000ft_agl_ug_m3']} ug/m3 "
            f"6000ft={item['smoke_6000ft_agl_ug_m3']} ug/m3 "
            f"elapsed={time.monotonic()-started:.1f}s "
            f"messages={item.get('download_stats',{}).get('messages',0)}",
            flush=True,
        )
        if item.get("reason"):
            print("  Reason:", item["reason"], flush=True)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(report(now, cycle, samples), indent=2) + "\n",
                      encoding="utf-8")
    print("Saved", OUTPUT, flush=True)
    print(f"Valid vertical samples: {sum(s['status']=='available' for s in samples)}/{len(samples)}",
          flush=True)


if __name__ == "__main__":
    main()
