"""Support diagnostics requested in review: zone size, monthly support, single-sided dates,
missing units and the candidate-city-year denominator.

Everything is read from frozen files (`daily_records_{day,night}.csv`,
`same_source_decomposition_city_year.csv`, `city_covariates.csv`, and the archived
`qc_tier_sensitivity_*.csv`). Nothing is re-estimated.

Outputs
    tables/table5_00_zone_size.md              main-text Table 5  (zone size vs the 20-pixel screen)
    （表 6「逐年严格质量」由 23_build_support_tables.py 产生，编号已改名，数值未动）
    tables/table9_00_monthly_support.md        main-text Table 9  (monthly support + single-sided split)
    outputs/daily/support_zone_size.csv        diagnostics
    outputs/daily/support_single_sided.csv     diagnostics
    outputs/daily/support_missing_units.csv    diagnostics
"""

from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parents[2]
DATA = BASE / "data"
OUT = BASE / "outputs" / "daily"
TAB = BASE / "tables"
ARCH = TAB / "_archive_pre_major_revision"
MIN_PX, MIN_DAYS = 20, 6
NEAR = 2 * MIN_PX          # “within a factor of two of the screen”


def read_records(mode: str, cols=None) -> pd.DataFrame:
    return pd.read_csv(DATA / "daily" / f"daily_records_{mode}.csv",
                       usecols=cols or ["city_id", "year", "month", "u_n", "r_n"])


def zone_size() -> tuple:
    rows, zh = [], []
    for mode in ("day", "night"):
        rec = read_records(mode)
        for zone, col in (("Urban", "u_n"), ("Rural", "r_n")):
            v = rec.loc[rec[col] > 0, col]
            rows.append({"mode": mode, "zone": zone, "n_zone_days": int(len(v)),
                         "median_px": float(np.median(v)),
                         "p10_px": float(np.percentile(v, 10)),
                         "p90_px": float(np.percentile(v, 90)),
                         "share_below_40_px": float((v < NEAR).mean())})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "support_zone_size.csv", index=False)
    body = ["| Zone | Overpass | Median valid pixels | 10th percentile | 90th percentile | Share of zone-days below 40 pixels |",
            "|---|---|---:|---:|---:|---:|"]
    zh_body = ["| 盘区 | 过境 | 有效像元中位数 | 第 10 百分位 | 第 90 百分位 | 不足 40 像元的 zone-day 占比 |",
               "|---|---|---:|---:|---:|---:|"]
    label = {("Urban", "day"): ("Urban", "Day"), ("Urban", "night"): ("Urban", "Night"),
             ("Rural", "day"): ("Rural", "Day"), ("Rural", "night"): ("Rural", "Night")}
    zh_label = {("Urban", "day"): ("城区", "白天"), ("Urban", "night"): ("城区", "夜间"),
                ("Rural", "day"): ("乡村", "白天"), ("Rural", "night"): ("乡村", "夜间")}
    for _, r in df.iterrows():
        a, b = label[(r.zone, r["mode"])]
        body.append(f"| {a} | {b} | {r.median_px:.0f} | {r.p10_px:.0f} | {r.p90_px:.0f} | {100*r.share_below_40_px:.1f}% |")
        za, zb = zh_label[(r.zone, r["mode"])]
        zh_body.append(f"| {za} | {zb} | {r.median_px:.0f} | {r.p10_px:.0f} | {r.p90_px:.0f} | {100*r.share_below_40_px:.1f}% |")
    note = ("Note. A zone-day is one city-day on which that zone has at least one valid 1-km pixel. "
            "On the 1-km grid 20 pixels correspond to about 20 km², so the 40-pixel column marks the "
            "zone-days that sit within a factor of two of the primary screen.")
    (TAB / "table5_00_zone_size.md").write_text(
        "\n".join(["**Table 5.** Zone size against the 20-pixel screen."] + [""] + body + ["", note, ""]),
        encoding="utf-8")
    print("已写 tables/table5_00_zone_size.md · 行", len(body) - 2)
    print("=== 中文表 5 ===")
    print("\n".join(zh_body))
    return df


def monthly_support() -> tuple:
    per_mode = {}
    rows = []
    for mode in ("day", "night"):
        rec = read_records(mode)
        rec["pair"] = (rec.u_n >= MIN_PX) & (rec.r_n >= MIN_PX)
        rec["u_ok"] = rec.u_n >= MIN_PX
        rec["r_ok"] = rec.r_n >= MIN_PX
        rec["any"] = rec.u_ok | rec.r_ok
        g = rec.groupby(["city_id", "year", "month"]).agg(
            n_pair=("pair", "sum"), n_u=("u_ok", "sum"), n_r=("r_ok", "sum"),
            n_any=("any", "sum")).reset_index()
        # 单侧日：恰好一侧通过
        rec["urban_only"] = rec.u_ok & ~rec.r_ok
        rec["rural_only"] = rec.r_ok & ~rec.u_ok
        g2 = rec.groupby(["city_id", "year", "month"]).agg(
            n_urban_only=("urban_only", "sum"), n_rural_only=("rural_only", "sum")).reset_index()
        g = g.merge(g2, on=["city_id", "year", "month"])
        g["estimable"] = (g.n_pair >= MIN_DAYS) & (g.n_u >= MIN_DAYS) & (g.n_r >= MIN_DAYS)
        g["urban_only_pct"] = 100 * g.n_urban_only / g.n_any.replace(0, np.nan)
        g["rural_only_pct"] = 100 * g.n_rural_only / g.n_any.replace(0, np.nan)
        g["mode"] = mode
        per_mode[mode] = g
        s = g.groupby("month").agg(median_pair=("n_pair", "median"),
                                   share_estimable=("estimable", "mean"),
                                   urban_only=("urban_only_pct", "median"),
                                   rural_only=("rural_only_pct", "median"),
                                   n_city_months=("estimable", "size")).reset_index()
        s["mode"] = mode
        rows.append(s)
        print(f"{mode}: 城市—月 {len(g)}；总体单侧占比中位 仅城区 "
              f"{g.urban_only_pct.median():.1f}% · 仅乡村 {g.rural_only_pct.median():.1f}%")
    allm = pd.concat(rows, ignore_index=True).sort_values(["mode", "month"])
    allm.to_csv(OUT / "support_single_sided.csv", index=False)
    body = ["| Month | Night: median paired days | Night: estimable share | Night: urban-only share | Night: rural-only share | "
            "Day: median paired days | Day: estimable share | Day: urban-only share | Day: rural-only share |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    zh = ["| 日历月 | 夜间：同日有效日中位数 | 夜间：可估计占比 | 夜间：仅城区占比 | 夜间：仅乡村占比 | "
          "白天：同日有效日中位数 | 白天：可估计占比 | 白天：仅城区占比 | 白天：仅乡村占比 |",
          "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for mth in range(1, 13):
        n = allm[(allm["mode"] == "night") & (allm["month"] == mth)].iloc[0]
        d = allm[(allm["mode"] == "day") & (allm["month"] == mth)].iloc[0]
        body.append(f"| {mth} | {n.median_pair:.0f} | {100*n.share_estimable:.1f}% | {n.urban_only:.1f}% | {n.rural_only:.1f}% | "
                    f"{d.median_pair:.0f} | {100*d.share_estimable:.1f}% | {d.urban_only:.1f}% | {d.rural_only:.1f}% |")
        zh.append(f"| {mth} | {n.median_pair:.0f} | {100*n.share_estimable:.1f}% | {n.urban_only:.1f}% | {n.rural_only:.1f}% | "
                  f"{d.median_pair:.0f} | {100*d.share_estimable:.1f}% | {d.urban_only:.1f}% | {d.rural_only:.1f}% |")
    note = ("Note. A city-month is estimable when the paired-day count and the valid-day count of each "
            "zone all reach six. The urban-only and rural-only shares are medians, over the city-months "
            "of that calendar month, of the share of observation days on which exactly one zone passes "
            "the 20-pixel screen; the denominator is the days on which at least one zone passes.")
    (TAB / "table9_00_monthly_support.md").write_text(
        "\n".join(["**Table 9.** Monthly support and the composition of single-sided days."] + [""] + body + ["", note, ""]),
        encoding="utf-8")
    print("已写 tables/table9_00_monthly_support.md · 行", len(body) - 2)
    print("=== 中文表 9 ===")
    print("\n".join(zh))
    return allm


def missing_units() -> pd.DataFrame:
    annual = pd.read_csv(OUT / "same_source_decomposition_city_year.csv",
                         usecols=["basis", "mode", "city_id", "year", "pixel_rule", "min_daily_days",
                                  "min_8day_periods", "common_months_all_methods",
                                  "delta_A_to_B_common_months_K"])
    q = annual.query("basis == 'dynamic' and pixel_rule == '20px' and "
                     "min_daily_days == 6 and min_8day_periods == 2")
    cov = pd.read_csv(DATA / "city_covariates.csv")
    cov = cov[cov.year == 2020][["city_id", "city"]]
    out = []
    for mode in ("day", "night"):
        rec = read_records(mode)
        g = q[q["mode"] == mode]
        usable = set(g.loc[g.delta_A_to_B_common_months_K.notna(), "city_id"])
        for _, row in cov.iterrows():
            sub = rec[rec.city_id == row.city_id]
            out.append({"mode": mode, "city_id": int(row.city_id), "city": row.city,
                        "in_records": bool(len(sub)),
                        "max_urban_px": int(sub.u_n.max()) if len(sub) else 0,
                        "median_urban_px": float(sub.u_n.median()) if len(sub) else np.nan,
                        "paired_days": int(((sub.u_n >= MIN_PX) & (sub.r_n >= MIN_PX)).sum()) if len(sub) else 0,
                        "estimable": bool(row.city_id in usable)})
    df = pd.DataFrame(out)
    df.to_csv(OUT / "support_missing_units.csv", index=False)
    n_units = cov.city_id.nunique()
    print(f"\n单元总数 {n_units} · 25 年 → 可能的城市—年 {n_units * 25}")
    for mode in ("day", "night"):
        rec = read_records(mode)
        present = rec.groupby(["city_id", "year"]).ngroups
        usable = int(q[(q["mode"] == mode)].delta_A_to_B_common_months_K.notna().sum())
        print(f"  {mode}: daily_records 中城市—年 {present} · 可估计（≥6 共同月份）{usable} · "
              f"缺口 {n_units*25 - present}（无合格逐日记录）+ {present - usable}（共同月份不足）")
    miss = df[~df.estimable]
    print("\n未进入可估计面板的单元：")
    for _, r in miss.iterrows():
        print(f"  {r['mode']:5s} {r.city:22s} 有无逐日记录={r.in_records} "
              f"城区像元 max={r.max_urban_px} 中位={r.median_urban_px} 同日有效日={r.paired_days}")
    return df


def main() -> None:
    zone_size()
    monthly_support()
    missing_units()


if __name__ == "__main__":
    main()
