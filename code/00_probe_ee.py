"""Check Earth Engine access and the daily MOD11A1 collection without exporting data."""

import argparse
import sys

sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")   # 只改 stdout 时 stderr 仍走 GBK，中文报错会乱码


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", default="evident-ocean-510002-t0")
    args = parser.parse_args()

    try:
        import ee

        ee.Initialize(project=args.project)
        image = (
            ee.ImageCollection("MODIS/061/MOD11A1")
            .filterDate("2020-01-01", "2020-01-03")
            .first()
        )
        bands = image.bandNames().getInfo()
        needed = {"LST_Day_1km", "LST_Night_1km", "QC_Day", "QC_Night"}
        print(f"EE access: OK; daily MOD11A1 bands: {len(bands)}")
        print(f"required bands present: {needed.issubset(set(bands))}")
        return 0 if needed.issubset(set(bands)) else 2
    except Exception as exc:
        # Avoid printing credentials or full request URLs from an auth exception.
        detail = str(exc).lower()
        if "authenticate" in detail or "credentials" in detail:
            reason = "credentials unavailable"
        elif "permission" in detail or "not registered" in detail:
            reason = "project permission unavailable"
        elif "connection" in detail or "network" in detail:
            reason = "network unavailable"
        else:
            reason = "unclassified"
        print(f"EE access: FAILED ({type(exc).__name__}: {reason})", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
