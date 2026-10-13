import requests
from datetime import datetime, timedelta, timezone

now = datetime.now(timezone.utc)
safe = now - timedelta(hours=3)
cycle = safe.replace(
    hour=(safe.hour // 6) * 6,
    minute=0, second=0, microsecond=0
)
fxx = 6

url = (
    "https://nomads.ncep.noaa.gov/pub/data/nccf/"
    "com/hrrr/prod/"
    f"hrrr.{cycle:%Y%m%d}/conus/"
    f"hrrr.t{cycle:%H}z.wrfnatf{fxx:02d}.grib2.idx"
)

response = requests.get(url, timeout=45)
response.raise_for_status()

fields = ("MASSDEN", "HGT", "PMTF", "PMTC")
records = []

for line in response.text.splitlines():
    parts = line.split(":")
    if len(parts) >= 5 and parts[3] in fields:
        records.append((parts[3], parts[4]))

print("Model cycle:", cycle.isoformat())
print("Forecast hour:", fxx)
print("Index bytes:", len(response.content))

for field in fields:
    levels = [level for name, level in records
              if name == field]
    print(field, "records:", len(levels))
    print("First levels:", levels[:8])
