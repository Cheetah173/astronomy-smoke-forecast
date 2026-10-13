#!/usr/bin/env python3
"""Read-only HRRR native-level smoke profile diagnostic, one forecast hour.

Reads GRIB2 index and byte-range selected HGT surface, hybrid HGT and hybrid
MASSDEN records from the public NOAA HRRR object archive. NOAA documents
hybrid MASSDEN as smoke, in kg/m^3 (converted to ug/m^3 by 1e9).
Does not modify the existing hourly smoke-column or Open-Meteo products.
"""
from __future__ import annotations

import time
import json
import math
import re
import tempfile
from datetime import datetime, timedelta, timezone
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from noaa_probe import POINTS, select_cycle

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "docs" / "noaa_vertical_probe.json"
URL_PREFIX = "https://noaa-hrrr-bdp-pds.s3.amazonaws.com"
COARSE = [4, 5, 6, 11, 12]
TARGETS = {"smoke_1000ft_agl_ug_m3": 304.8,
           "smoke_6000ft_agl_ug_m3": 1828.8}
MAX_INDEX = 250_000
MAX_MESSAGE = 6_000_000


def utc(dt):
    return dt.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")


def index_entries(content):
    out = []
    for line in content.splitlines():
        parts = line.split(":")
        if len(parts) < 6 or not parts[1].isdigit():
            continue
        level = parts[4]
        match = re.fullmatch(r"(\d+) hybrid level", level)
        out.append({"start": int(parts[1]), "name": parts[3],
                    "level": int(match.group(1)) if match else None,
                    "surface": level == "surface", "text_level": level})
    for i, item in enumerate(out):
        item["end"] = out[i + 1]["start"] - 1 if i + 1 < len(out) else None
    return out


def find_entry(entries, name, level=None, surface=False):
    matches = [e for e in entries if e["name"] == name and
               (e["surface"] if surface else e["level"] == level)]
    if len(matches) != 1 or matches[0]["end"] is None:
        raise ValueError(f"Expected exactly one bounded {name} at {level or 'surface'}, found {len(matches)}")
    return matches[0]


def bracket(level_to_height, target):
    ordered = sorted(level_to_height.items(), key=lambda x: x[1])
    for (l0, h0), (l1, h1) in zip(ordered, ordered[1:]):
        if h0 <= target <= h1 and h1 > h0:
            return l0, l1
    return None


def interp(a0, a1, target):
    h0, v0 = a0
    h1, v1 = a1
    if not (math.isfinite(h0) and math.isfinite(h1) and
            math.isfinite(v0) and math.isfinite(v1) and h0 < h1):
        raise ValueError("Invalid vertical interpolation brackets")
    return v0 + (target - h0) / (h1 - h0) * (v1 - v0)


def extract():
    from eccodes import codes_get, codes_grib_find_nearest_multiple, codes_grib_new_from_file, codes_release

    now = datetime.now(timezone.utc)
    cycle, _ = select_cycle(now)
    # Use a published +6h lead (one forecast only, not a time series).
    lead = 6
    name = f"hrrr.t{cycle:%H}z.wrfnatf{lead:02d}.grib2"
    file_url = f"{URL_PREFIX}/hrrr.{cycle:%Y%m%d}/conus/{name}"
    report = {"generated_utc": utc(now), "cycle_utc": utc(cycle),
              "valid_utc": utc(cycle + timedelta(hours=lead)),
              "forecast_hour": lead, "source_url": file_url,
              "status": "unavailable", "units": "ug/m^3",
              "model_field": "HRRR wrfnat hybrid-level MASSDEN (kg/m^3)",
              "method": "linear interpolation over HGT minus surface HGT (AGL)",
              "results": []}
    session = requests.Session()
    retry_policy = Retry(
        total=3,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=frozenset(["GET"]),
    )
    session.mount(
        "https://",
        HTTPAdapter(max_retries=retry_policy)
    )
    session.headers["User-Agent"] = "AstronomySmokeForecast/0.5 (research profile diagnostic)"
    cache = {}
    stats = {"downloaded_messages": 0, "bytes": 0}

    def fetch(e):
        key = (e["name"], e["level"], e["surface"])
        if key in cache:
            return cache[key]
        start, end = e["start"], e["end"]
        length = end - start + 1
        if not 0 < length <= MAX_MESSAGE:
            raise ValueError(f"Unexpected GRIB message length {length} for {key}")
        started = time.monotonic()
        response = session.get(file_url, headers={"Range": f"bytes={start}-{end}"}, timeout=90)
        response.raise_for_status()
        download_seconds = time.monotonic() - started
        if response.status_code != 206 or len(response.content) != length or not response.content.startswith(b"GRIB"):
            raise ValueError(f"HTTP partial content failed for {key}: status={response.status_code}, bytes={len(response.content)}, wanted={length}")
        with tempfile.TemporaryFile(mode="w+b") as tmp:
            tmp.write(response.content)
            tmp.seek(0)
            gid = codes_grib_new_from_file(tmp)
            if gid is None:
                raise ValueError(f"GRIB decoder found no message for {key}")
            try:
                gd = int(codes_get(gid, "dataDate"))
                gh = int(codes_get(gid, "dataTime"))
                ft = int(codes_get(gid, "forecastTime"))
                if gd != int(cycle.strftime("%Y%m%d")) or gh != int(cycle.strftime("%H%M")) or ft != lead:
                    raise ValueError(f"Unexpected metadata date/time/hour for {key}")
                nearest = codes_grib_find_nearest_multiple(
                    gid,
                    False,
                    [lat for _, lat, _ in POINTS],
                    [lon for _, _, lon in POINTS]
                )
                if len(nearest) != len(POINTS):
                    raise ValueError("Incomplete regional grid lookup")

                values = {}
                for (label, _, _), closest in zip(POINTS, nearest):
                    v = float(closest.value)
                    if not math.isfinite(v):
                        raise ValueError(f"Non-finite {key} for {label}")
                    values[label] = v
            finally:
                codes_release(gid)
        cache[key] = values
        stats["downloaded_messages"] += 1
        stats["bytes"] += length
        
        print(
            f"{e['name']} level {e['level']}: "
            f"download={download_seconds:.1f}s, "
            f"total={time.monotonic()-started:.1f}s",
            flush=True
        )

        return values

    try:
        idx = session.get(file_url + ".idx", timeout=45)
        idx.raise_for_status()
        if len(idx.content) > MAX_INDEX:
            raise ValueError("Index unexpectedly large")
        entries = index_entries(idx.text)
        hybrid_heights = {e["level"] for e in entries if e["name"] == "HGT" and e["level"] is not None}
        hybrid_smoke = {e["level"] for e in entries if e["name"] == "MASSDEN" and e["level"] is not None}
        levels = hybrid_heights & hybrid_smoke
        if len(levels) < 45:
            raise ValueError(f"Too few common hybrid levels: {len(levels)}")
        terrain = fetch(find_entry(entries, "HGT", surface=True))
        sampled = [x for x in COARSE if x in levels]
        for level in sampled:
            fetch(find_entry(entries, "HGT", level=level))
        # Identify coarse bracketing intervals for both AGL targets, every regional point.
        finer = set()
        for label, _, _ in POINTS:
            heights = {level: fetch(find_entry(entries, "HGT", level=level))[label] - terrain[label]
                       for level in sampled}
            for target in TARGETS.values():
                pair = bracket(heights, target)
                if pair is None:
                    raise ValueError(f"Coarse hybrid height levels do not bracket {target} m at {label}")
                lo, hi = sorted(pair)
                finer.update(level for level in levels if lo <= level <= hi)
        for level in sorted(finer):
            fetch(find_entry(entries, "HGT", level=level))
            fetch(find_entry(entries, "MASSDEN", level=level))
        for label, lat, lon in POINTS:
            profile = {}
            for level in sorted(finer):
                h = fetch(find_entry(entries, "HGT", level=level))[label] - terrain[label]
                concentration = fetch(find_entry(entries, "MASSDEN", level=level))[label] * 1e9
                if h < -100 or concentration < -0.001:
                    raise ValueError(f"Implausible profile at {label} L{level}: {h}, {concentration}")
                profile[level] = (h, max(0.0, concentration))
            row = {"label": label, "latitude": lat, "longitude": lon,
                   "surface_height_m": round(terrain[label], 2)}
            for key, target in TARGETS.items():
                pair = bracket({level: values[0] for level, values in profile.items()}, target)
                if pair is None:
                    row[key] = None
                    row[key + "_reason"] = "No exact adjacent hybrid heights bracketing target"
                    continue
                l0, l1 = sorted(pair, key=lambda level: profile[level][0])
                h0, v0 = profile[l0]
                h1, v1 = profile[l1]
                if abs(l0 - l1) != 1 or (h1 - h0) > 750:
                    row[key] = None
                    row[key + "_reason"] = "Vertical samples too widely separated"
                    continue
                row[key] = round(interp((h0, v0), (h1, v1), target), 4)
                row[key + "_levels"] = [l0, l1]
                row[key + "_height_brackets_m_agl"] = [round(h0, 1), round(h1, 1)]
            report["results"].append(row)
        site = next(row for row in report["results"] if row["label"] == "SITE")
        report["status"] = "available" if all(site[key] is not None for key in TARGETS) else "partial"
        report["native_levels_used"] = sorted(finer)
    except Exception as exc:
        report["error"] = f"{type(exc).__name__}: {exc}"[:650]
    report["download_stats"] = stats
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("Vertical probe status:", report["status"])
    print("NOAA run:", report["cycle_utc"], "forecast hour:", lead)
    print("Downloaded messages:", stats["downloaded_messages"], "bytes:", stats["bytes"])
    if report.get("error"):
        print("Diagnostic error:", report["error"])
    if report.get("results"):
        print("SITE:", json.dumps(next(x for x in report["results"] if x["label"] == "SITE")))
    print("Saved", OUTPUT)
    return report


if __name__ == "__main__":
    extract()
