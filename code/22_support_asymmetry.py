"""Descriptive structure of the urban–rural date-set asymmetry (optional T8).

For every city–year–month the script asks how large a share of the observation days is
single-sided (exactly one zone passes the 20-pixel screen), then summarises that share by
overpass, by calendar month and by city latitude. This is a description of the sampling
asymmetry that drives the A−B step; it is not a causal statement, and it does not replace
any existing estimate.

Writes only new files: `outputs/daily/support_asymmetry_summary.csv` (per mode × calendar
month) and `tables/tableS5_00_data_products.md` (the product table, used by the supplement).
Existing analysis products are never overwritten.

⚠ 文件名故意用 `tableS5_` 而不是 `table5_`：`build_draft.py` 用 `glob("table5_*.md")` 取
正文表 5，若同目录再放一个 `table5_*` 文件，`sorted(...)[0]` 会按字母序选错表
（2026-10-04 实跑踩到：正文表 5 被渲染成产品表）。
"""

from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parents[2]
DATA = BASE / "data"
OUT = BASE / "outputs" / "daily"
TAB = BASE / "tables"
MIN_PX = 20

PRODUCTS = [
    ("MOD11A1.061 LST (Terra, daily)", "061", "1 km", "2001–2025",
     "Paired-date and independent-date estimators A and B"),
    ("MOD11A2.061 LST (Terra, 8-day)", "061", "1 km", "2001–2025",
     "Native product comparison D"),
    ("MOD13A2.061 NDVI", "061", "1 km", "2001–2025",
     "Rural-reference screen (regional 75th percentile)"),
    ("GHSL built-up surface, P2023A", "P2023A", "100 m, aggregated to 1 km", "2000–2020 epochs",
     "Urban mask (built fraction above 0.30)"),
    ("SRTMGL1 elevation", "v3", "30 m, aggregated to 1 km", "static",
     "Rural-reference elevation screen (within 50 m)"),
    ("JRC Global Surface Water", "v1.4", "30 m, aggregated to 1 km", "static",
     "Water exclusion (occurrence above 50%)"),
    ("FAO GAUL level 2 boundaries", "2025", "vector", "2025",
     "Administrative units and city identifiers"),
]


def summarise(mode: str) -> pd.DataFrame:
    rec = pd.read_csv(DATA / "daily" / f"daily_records_{mode}.csv",
                      usecols=["city_id", "year", "month", "u_n", "r_n"])
    rec["both"] = (rec.u_n >= MIN_PX) & (rec.r_n >= MIN_PX)
    rec["one"] = (rec.u_n >= MIN_PX) ^ (rec.r_n >= MIN_PX)
    rec["any"] = (rec.u_n >= MIN_PX) | (rec.r_n >= MIN_PX)
    g = rec.groupby(["city_id", "year", "month"]).agg(
        n_any=("any", "sum"), n_one=("one", "sum"), n_both=("both", "sum")).reset_index()
    g = g[g.n_any > 0].copy()
    g["share_one"] = 100.0 * g.n_one / g.n_any
    g["share_both"] = 100.0 * g.n_both / g.n_any
    return g


def main() -> None:
    cov = pd.read_csv(DATA / "city_covariates.csv")
    lat = cov[cov.year == 2020].set_index("city_id").lat

    rows = []
    print(f"单一侧有效日的占比（每城市—年—月一个观测，门槛 {MIN_PX} 像元）")
    for mode in ("night", "day"):
        g = summarise(mode)
        g["lat"] = g.city_id.map(lat)
        print(f"  {mode:5s} 城市—年—月 {len(g):,} 个 · 单侧占比 中位 {g.share_one.median():.1f}% "
              f"· 两侧均有效占比 中位 {g.share_both.median():.1f}% "
              f"· 与纬度相关 r={g.share_one.corr(g.lat):+.3f}")
        for m in range(1, 13):
            s = g[g["month"] == m]
            rows.append({"mode": mode, "month": m, "city_months": len(s),
                         "median_share_single_sided_pct": round(s.share_one.median(), 2),
                         "median_share_paired_pct": round(s.share_both.median(), 2)})

    summ = pd.DataFrame(rows)
    summ.to_csv(OUT / "support_asymmetry_summary.csv", index=False)
    print("  已写 outputs/daily/support_asymmetry_summary.csv")
    for mode in ("night", "day"):
        s = summ[summ["mode"] == mode]
        hi = s.loc[s.median_share_single_sided_pct.idxmax()]
        lo = s.loc[s.median_share_single_sided_pct.idxmin()]
        print(f"  {mode:5s} 单侧占比最高月 {int(hi['month'])} 月 {hi.median_share_single_sided_pct:.1f}%"
              f" · 最低月 {int(lo['month'])} 月 {lo.median_share_single_sided_pct:.1f}%")

    lines = ["**Table 5.** Input products and derived layers.",
             "", "| Product or layer | Version | Native resolution | Period used | Role in this analysis |",
             "|---|---|---|---|---|"]
    lines += [f"| {a} | {b} | {c} | {d} | {e} |" for a, b, c, d, e in PRODUCTS]
    lines += ["", "Note. All layers were accessed through the Google Earth Engine Data Catalog and "
                  "aggregated to the 1-km analysis grid described in Section 2.1. No layer was "
                  "resampled from a coarser product.", ""]
    (TAB / "tableS5_00_data_products.md").write_text("\n".join(lines), encoding="utf-8")
    print("  已写 tables/tableS5_00_data_products.md")


if __name__ == "__main__":
    main()
