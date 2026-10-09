import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forecast import build_report, format_utc_and_local, LOCATIONS


class ForecastTests(unittest.TestCase):
    def test_timezone_conversion(self):
        utc, local = format_utc_and_local("2026-10-10T06:00")
        self.assertEqual(utc, "2026-10-10T06:00:00Z")
        self.assertEqual(local, "2026-10-10T01:00:00-05:00")

    def test_report_grid(self):
        inputs = [{"hourly": {
            "time": ["2026-10-10T06:00", "2026-10-10T07:00"],
            "aerosol_optical_depth": [0.14, 0.18],
            "pm2_5": [5, 6], "us_aqi": [20, 25],
        }} for _ in LOCATIONS]
        r = build_report(inputs, "2026-10-09T20:00:00Z")
        self.assertEqual(len(r["hours"]), 2)
        self.assertEqual(len(r["hours"][0]["region"]), 9)
        self.assertEqual(r["hours"][0]["site"]["aod"], 0.14)
        self.assertIn("unavailable", r["measurement_status"]["smoke_6000ft_agl_ug_m3"])


if __name__ == "__main__":
    unittest.main()
