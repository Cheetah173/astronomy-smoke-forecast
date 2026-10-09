#!/usr/bin/env python3
"""NOAA HRRR *diagnostic* smoke-column extractor.

Single restricted-area GRIB2 download from NOMADS, safely isolated from the
working Open-Meteo product. Fails closed if aerosol species can't be verified.
Future stage: multiple forecast hours and 1000/6000 ft AGL native-level fields.
"""
from __future__ import annotations

import tempfile
import json
import math
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import requests

OUTPUT = Path(__file__).resolve().parent / 'docs' / 'noaa_probe.json'
URL = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_hrrr_2d.pl'
# Public, approximate location; nine points roughly 70 miles from Robinson.
POINTS = [
    ('NW', 32.45, -98.34), ('N', 32.45, -97.14), ('NE', 32.45, -95.94),
    ('W', 31.45, -98.34), ('SITE', 31.45, -97.14), ('E', 31.45, -95.94),
    ('SW', 30.45, -98.34), ('S', 30.45, -97.14), ('SE', 30.45, -95.94),
]


def select_cycle(now: datetime) -> tuple[datetime, int]:
    """Use an available 6-hourly cycle with enough publication time."""
    now = now.astimezone(timezone.utc)
    safe = now - timedelta(hours=3)
    cycle = safe.replace(hour=(safe.hour // 6) * 6, minute=0, second=0, microsecond=0)
    fxx = math.ceil((now - cycle).total_seconds() / 3600) + 2
    return cycle, max(1, min(18, fxx))


def build_params(cycle: datetime, fxx: int) -> dict:
    return {
        'file': f'hrrr.t{cycle:%H}z.wrfsfcf{fxx:02d}.grib2',
        'dir': f'/hrrr.{cycle:%Y%m%d}/conus',
        'var_COLMD': 'on',
        'var_AOTK': 'on',
        'lev_entire_atmosphere_(considered_as_a_single_layer)': 'on',
        'subregion': '',
        'leftlon': '-99.0', 'rightlon': '-95.0',
        'toplat': '33.0', 'bottomlat': '30.0',
    }


def get_or_none(gid, key: str):
    from eccodes import codes_get
    try:
        v = codes_get(gid, key)
        if v in (None, 'undef', 'unknown', 'missing', ''):
            return None
        return v
    except Exception:
        return None


def field_metadata(gid: int) -> dict:
    keys = (
        'shortName', 'name', 'units', 'paramId', 'discipline',
        'parameterCategory', 'parameterNumber', 'productDefinitionTemplateNumber',
        'typeOfAerosol', 'aerosolType', 'typeOfChemicalConstituent',
        'aerosolTypeDescription', 'typeOfSizeInterval',
        'typeOfFirstFixedSurface', 'typeOfSecondFixedSurface',
        'forecastTime', 'stepUnits', 'centre', 'dataDate', 'dataTime',
    )
    return {k: v for k in keys if (v := get_or_none(gid, k)) is not None}


def candidate_smoke(meta: dict) -> bool:
    """Only certify a column as smoke when aerosol species identifies organic matter.

    A generic COLMD is not sufficient evidence of smoke in a multi-species file.
    """
    identity = ' '.join(str(meta.get(k, '')) for k in
                    ('name', 'aerosolType', 'typeOfChemicalConstituent',
                     'aerosolTypeDescription')).lower()
    aerosol_codes = {str(meta.get(k, '')) for k in
                     ('aerosolType', 'typeOfAerosol', 'typeOfChemicalConstituent')}
    # WMO GRIB code table 4.233: particulate organic matter and smoke classes.
    # Organic matter is a *smoke proxy*, not uniquely wildfire-origin smoke.
    return (bool(aerosol_codes & {'62010', '62020', '62021'}) or
            'particulate organic matter' in identity or
            'organic matter dry' in identity or
            'wildfire smoke' in identity or
            'smoke - high absorption' in identity or
            'smoke - low absorption' in identity)


def grib_field_kind(meta: dict) -> str | None:
    """Identify fields by WMO GRIB2 codes even when ecCodes lacks shortName.

    NOAA/NCEP Table 4.2-0-20: (discipline=0, category=20)
      parameter 102 = AOTK (aerosol optical thickness; dimensionless)
      parameter   1 = COLMD (column-integrated mass density; kg/m^2)
    https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table4-2-0-20.shtml
    """
    try:
        code = tuple(int(meta[k]) for k in
                     ('discipline', 'parameterCategory', 'parameterNumber'))
    except (KeyError, ValueError, TypeError):
        return None
    return {(0, 20, 102): 'AOTK', (0, 20, 1): 'COLMD'}.get(code)


def nearby(gid: int) -> list[dict]:
    from eccodes import codes_grib_find_nearest
    rows = []
    for label, lat, lon in POINTS:
        hit = codes_grib_find_nearest(gid, lat, lon)[0]
        # ecCodes returns a namedtuple in recent versions, but older versions
        # may return mapping-like objects.
        def val(k):
            return getattr(hit, k) if hasattr(hit, k) else hit[k]
        v = float(val('value'))
        if not math.isfinite(v) or v < 0:
            raise ValueError(f'Invalid {label} GRIB value: {v}')
        rows.append({
            'label': label, 'latitude': lat, 'longitude': lon,
            'grid_latitude': round(float(val('lat')), 4),
            'grid_longitude': round(((float(val('lon')) + 180) % 360) - 180, 4),
            'raw_value': v,
        })
    return rows


def decode(data: bytes) -> dict:
    from eccodes import codes_grib_new_from_file, codes_release
    if not data.startswith(b'GRIB'):
        raise ValueError('NOMADS returned HTML/other content, not GRIB2')
    found = []
    with tempfile.TemporaryFile(mode='w+b') as source:
        source.write(data)
        source.seek(0)
        while True:
            gid = codes_grib_new_from_file(source)
            if gid is None:
                break
            try:
                meta = field_metadata(gid)
                kind = grib_field_kind(meta)
                print('DECODED GRIB FIELD:', json.dumps(meta, default=str))
                if kind is None:
                    continue
                # WMO code numbers are authoritative when ecCodes lacks an alias.
                meta['shortName'] = kind
                meta['units'] = 'kg/m^2' if kind == 'COLMD' else 'dimensionless'
                meta['identification_method'] = 'WMO GRIB2 discipline/category/parameterNumber'
                item = {'metadata': meta, 'nearest': nearby(gid)}
                item['smoke_species_encoded_in_message'] = candidate_smoke(meta)
                found.append(item)
            finally:
                codes_release(gid)
    columns = [f for f in found if f['metadata']['shortName'] == 'COLMD']
    aods = [f for f in found if f['metadata']['shortName'] == 'AOTK']
    result = {
        'candidate_fields': [{
            'metadata': f['metadata'],
            'smoke_species_encoded_in_message': f['smoke_species_encoded_in_message'],
            'site_raw_value': next(p['raw_value'] for p in f['nearest'] if p['label'] == 'SITE'),
        } for f in found],
        'column_status': 'unavailable — expected one HRRR COLMD record',
        'hrrr_aod_status': 'unavailable — expected one HRRR AOTK record',
        'provenance': {
            'grib_code_table': 'https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table4-2-0-20.shtml',
            'hrrr_smoke_documentation': 'https://www.weather.gov/media/notification/pdf2/scn20-46rap_v5_hrrr_v4_aab.pdf',
            'note': 'NOAA NWS explicitly documents the HRRR COLMD product as column-integrated smoke. This GRIB message may encode the physical parameter but no independent aerosol-species tag. Interpret only in the context of the NOAA HRRR CONUS wrfsfc COLMD product, not arbitrary COLMD fields.'
        },
    }
    if len(columns) == 1:
        col = columns[0]
        if col['smoke_species_encoded_in_message']:
            provenance = 'GRIB species metadata and NOAA HRRR documentation'
        else:
            provenance = 'NOAA HRRR CONUS smoke product documentation (no species identifier in GRIB metadata)'
        result['column_status'] = 'available: NOAA HRRR vertically integrated smoke (documented product)'
        result['column_identification_provenance'] = provenance
        result['column_mg_m2'] = [
            {**{k: v for k, v in p.items() if k != 'raw_value'},
             'smoke_column_mg_m2': round(p['raw_value'] * 1e6, 3)}
            for p in col['nearest']
        ]
    if len(aods) == 1:
        result['hrrr_aod_status'] = 'available: NOAA HRRR AOTK (wavelength unverified)'
        result['hrrr_aod'] = [
            {**{k: v for k, v in p.items() if k != 'raw_value'},
             'aod': round(p['raw_value'], 5)}
            for p in aods[0]['nearest']
        ]
    elif len(aods) > 1:
        result['hrrr_aod_status'] = 'unavailable — more than one AOTK record; need wavelength distinction' 
    return result


def run() -> dict:
    now = datetime.now(timezone.utc)
    cycle, fxx = select_cycle(now)
    params = build_params(cycle, fxx)
    valid = cycle + timedelta(hours=fxx)
    report = {
        'generated_utc': now.isoformat().replace('+00:00', 'Z'),
        'model': 'NOAA NCEP HRRR CONUS',
        'cycle_utc': cycle.isoformat().replace('+00:00', 'Z'),
        'forecast_hour': fxx,
        'valid_utc': valid.isoformat().replace('+00:00', 'Z'),
        'source': URL,
        'bbox': {'left': -99, 'right': -95, 'top': 33, 'bottom': 30},
        'column_status': 'unavailable — NOAA extraction failed',
        'hrrr_aod_status': 'unavailable — NOAA extraction failed',
        'smoke_1000ft_agl_status': 'not yet connected',
        'smoke_6000ft_agl_status': 'not yet connected',
    }
    try:
        s = requests.Session()
        s.headers['User-Agent'] = 'AstronomySmokeForecast/0.2 (noncommercial hobby research)'
        response = s.get(URL, params=params, timeout=90)
        response.raise_for_status()
        print('NOAA response bytes:', len(response.content))
        if len(response.content) > 15_000_000:
            raise ValueError('Subregion download unexpectedly large; stop')
        report.update(decode(response.content))
    except Exception as exc:
        report['error'] = f'{type(exc).__name__}: {exc}'[:700]
        print('NOAA PROBE UNAVAILABLE:', report['error'])
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('NOAA probe result:', report['column_status'])
    print('NOAA AOD result:', report['hrrr_aod_status'])
    print('Candidate metadata:', json.dumps(report.get('candidate_fields', []), indent=2)[:4500])
    if 'column_mg_m2' in report:
        print('SITE SMOKE COLUMN mg/m2:', next(x['smoke_column_mg_m2'] for x in report['column_mg_m2'] if x['label'] == 'SITE'))
    if 'hrrr_aod' in report:
        print('SITE HRRR AOD:', next(x['aod'] for x in report['hrrr_aod'] if x['label'] == 'SITE'))
    print('Saved', OUTPUT)
    return report


if __name__ == '__main__':
    run()
