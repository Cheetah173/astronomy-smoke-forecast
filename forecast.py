#!/usr/bin/env python3
"""Regional AOD forecast extractor; stage-one foundation for NOAA HRRR Smoke integration.

This program *does not* claim to extract 1000ft/6000ft NOAA smoke. It intentionally
marks those unavailable until that numerical pipeline is implemented and validated.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from zoneinfo import ZoneInfo

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "docs" / "latest.json"
TZ = ZoneInfo("America/Chicago")
# Approximate central Robinson location; not a precise private observing address.
CENTER_LAT, CENTER_LON = 31.45, -97.14
# 9 regional points, roughly 70mi N/S and E/W, to detect spatial haze gradients.
LOCATIONS = [
    ("NW", 32.45, -98.34), ("N", 32.45, -97.14), ("NE", 32.45, -95.94),
    ("W", 31.45, -98.34), ("SITE", CENTER_LAT, CENTER_LON), ("E", 31.45, -95.94),
    ("SW", 30.45, -98.34), ("S", 30.45, -97.14), ("SE", 30.45, -95.94),
]
API = "https://air-quality-api.open-meteo.com/v1/air-quality"


def _round_number(value, digits=3):
    return round(float(value), digits) if value is not None else None


def format_utc_and_local(iso_utc: str) -> tuple[str, str]:
    dt = datetime.fromisoformat(iso_utc.replace("Z", "+00:00"))
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    dt = dt.astimezone(timezone.utc)
    return dt.isoformat().replace("+00:00", "Z"), dt.astimezone(TZ).isoformat()


def build_report(responses: list[dict], generated_utc: str) -> dict:
    """Build output from one response per regional point (also used by tests)."""
    if len(responses) != len(LOCATIONS):
        raise ValueError(f"Expected {len(LOCATIONS)} locations; received {len(responses)}")
    times = responses[0]["hourly"]["time"]
    if not times:
        raise ValueError("API returned zero hourly forecast entries")
    for entry in responses:
        if entry.get("hourly", {}).get("time") != times:
            raise ValueError("Locations have misaligned hourly forecast times")
        for field in ("aerosol_optical_depth", "pm2_5", "us_aqi"):
            values = entry["hourly"].get(field)
            if not isinstance(values, list) or len(values) != len(times):
                raise ValueError(f"Missing or invalid {field} hourly series")

    hours = []
    for i, t in enumerate(times):
        valid_utc, local = format_utc_and_local(t)
        regional = []
        for (label, lat, lon), entry in zip(LOCATIONS, responses):
            hourly = entry["hourly"]
            regional.append({
                "label": label, "latitude": lat, "longitude": lon,
                "aod": _round_number(hourly["aerosol_optical_depth"][i]),
                "pm2_5_ug_m3": _round_number(hourly["pm2_5"][i], 1),
                "us_aqi": _round_number(hourly["us_aqi"][i], 0),
            })
        site = next(p for p in regional if p["label"] == "SITE")
        hours.append({"valid_utc": valid_utc, "valid_local": local,
                      "site": site, "region": regional})

    return {
        "schema_version": "1.0", "generated_utc": generated_utc,
        "timezone": "America/Chicago",
        "location": {"name": "Robinson, Texas (approximate)",
                     "latitude": CENTER_LAT, "longitude": CENTER_LON},
        "sources": [{"name": "Open-Meteo Air Quality / CAMS",
                     "url": "https://open-meteo.com/en/docs/air-quality-api",
                     "model": "cams_global"}],
        "measurement_status": {
            "smoke_1000ft_agl_ug_m3": "unavailable: NOAA AGL decoder not yet connected",
            "smoke_6000ft_agl_ug_m3": "unavailable: NOAA AGL decoder not yet connected",
            "smoke_column_mg_m2": "unavailable: NOAA column decoder not yet connected",
            "aerosol_optical_depth": "available: Open-Meteo CAMS global (total aerosol, NOT smoke-specific)",
        },
        "notes": [
            "AOD is aerosol haze, not specifically wildfire smoke.",
            "The 9-point regional AOD grid can indicate where haze is forecast, but cannot independently identify a smoke plume.",
            "All values are modeled forecasts rather than measured transparency at the observing site.",
            "NOAA 1000ft, 6000ft, and column smoke remain unavailable until separately verified.",
        ],
        "hours": hours,
    }


def fetch_open_meteo() -> list[dict]:
    params = {
        "latitude": ",".join(str(lat) for _, lat, _ in LOCATIONS),
        "longitude": ",".join(str(lon) for _, _, lon in LOCATIONS),
        "hourly": "aerosol_optical_depth,pm2_5,us_aqi",
        "timezone": "UTC", "forecast_days": 4, "domains": "cams_global",
    }
    retry = Retry(total=3, backoff_factor=2, status_forcelist=[429, 500, 502, 503, 504],
                  allowed_methods=["GET"])
    session = requests.Session()
    session.headers["User-Agent"] = "AstronomySmokeForecast/0.1 (GitHub Actions; research)"
    session.mount("https://", HTTPAdapter(max_retries=retry))
    response = session.get(API, params=params, timeout=45)
    response.raise_for_status()
    data = response.json()
    if isinstance(data, dict):
        if data.get("error"):
            raise RuntimeError(f"Open-Meteo error: {data.get('reason', 'unknown')}")
        data = [data]
    if not isinstance(data, list):
        raise ValueError("Unexpected Open-Meteo response structure")
    return data


def main():
    created = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    data = build_report(fetch_open_meteo(), created)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    temp_file = OUTPUT.with_suffix(".json.tmp")
    temp_file.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    temp_file.replace(OUTPUT)
    available = sum(1 for h in data["hours"] if h["site"]["aod"] is not None)
    print(f"Forecast generated: {created}; {len(data['hours'])} hours, {available} site AOD values")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
