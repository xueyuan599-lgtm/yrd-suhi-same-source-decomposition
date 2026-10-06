"""Sensitivity to assigning only wholly within-month 8-day windows to a month.

Uses the already exported daily-rebuilt/native 8-day city-period summaries. It
does not reopen raster inputs. The alternative estimate is matched to the
primary estimate on identical city-years and calendar months.
"""
from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[2]
OUT = BASE / "outputs" / "daily"
MONTHS = OUT / "same_source_decomposition_monthly.csv"
PERIODS = OUT / "same_source_decomposition_8day_city_period.csv"
PRIMARY = dict(basis="dynamic", pixel_rule="20px", min_daily_days=6, min_8day_periods=2)


def main():
    monthly = pd.read_csv(MONTHS, parse_dates=[])
    periods = pd.read_csv(PERIODS, parse_dates=[])
    periods["date"] = pd.to_datetime(periods["date"])
    end = periods["date"] + pd.Timedelta(days=7)
    periods = periods[(periods.date.dt.year == end.dt.year) & (periods.date.dt.month == end.dt.month)].copy()
    periods["month"] = periods.date.dt.month
    for x in ["c_u_n", "c_r_n", "d_u_n", "d_r_n"]:
        periods[x] = pd.to_numeric(periods[x], errors="coerce")
    periods["c_valid"] = periods.c_u_n.ge(20) & periods.c_r_n.ge(20)
    periods["d_valid"] = periods.d_u_n.ge(20) & periods.d_r_n.ge(20)
    periods["cd_valid"] = periods.c_valid & periods.d_valid
    periods["c_value"] = periods.c_u_mean - periods.c_r_mean
    periods["d_value"] = periods.d_u_mean - periods.d_r_mean
    periods["cd_value"] = periods.c_value - periods.d_value
    keys = ["basis", "year", "mode", "city_id", "month"]
    # Construct one monthly record from the C/D windows that fit wholly inside
    # the calendar month, then retain only the exact primary city-month support.
    alt_rows = []
    for key, g in periods.groupby(keys, sort=False):
        d = dict(zip(keys, key))
        for label, valid, value in [("c", "c_valid", "c_value"), ("d", "d_valid", "d_value"), ("cd", "cd_valid", "cd_value")]:
            s = g[g[valid]]
            d[label + "_n"] = len(s)
            d[label + "_k"] = float(s[value].mean()) if len(s) >= 2 else np.nan
        alt_rows.append(d)
    alternative = pd.DataFrame(alt_rows)
    q = monthly.query("basis == @PRIMARY['basis'] and pixel_rule == @PRIMARY['pixel_rule'] and min_daily_days == @PRIMARY['min_daily_days'] and min_8day_periods == @PRIMARY['min_8day_periods']").copy()
    q = q.merge(alternative, on=keys, how="left")
    q = q[q[["A_paired_daily_K", "B_independent_dates_K", "C_daily_rebuilt_8day_K", "D_native_MOD11A2_K", "C_minus_D_paired_K"]].notna().all(axis=1)]
    q = q[q[["c_k", "d_k", "cd_k"]].notna().all(axis=1)]
    rows = []
    for (mode, year, city), g in q.groupby(["mode", "year", "city_id"]):
        if len(g) < 6:
            continue
        r = {"mode": mode, "year": int(year), "city_id": int(city), "months": len(g),
             "primary_B_minus_C_K": (g.B_independent_dates_K - g.C_daily_rebuilt_8day_K).mean(),
             "within_month_B_minus_C_K": (g.B_independent_dates_K - g.c_k).mean(),
             "primary_C_minus_D_K": g.C_minus_D_paired_K.mean(),
             "within_month_C_minus_D_K": g.cd_k.mean(),
             "primary_paired_periods": int(g.CD_n_paired_periods.sum()),
             "within_month_C_periods": int(g.c_n.sum()), "within_month_D_periods": int(g.d_n.sum()),
             "within_month_paired_periods": int(g.cd_n.sum())}
        rows.append(r)
    city_year = pd.DataFrame(rows)
    city_year.to_csv(OUT / "same_source_same_month_window_city_year.csv", index=False)
    rng = np.random.default_rng(20261004)
    output = []
    for mode, g in city_year.groupby("mode"):
        for metric in ["primary_B_minus_C_K", "within_month_B_minus_C_K", "primary_C_minus_D_K", "within_month_C_minus_D_K"]:
            city = g.groupby("city_id")[metric].mean().dropna().to_numpy()
            boots = np.asarray([rng.choice(city, size=len(city), replace=True).mean() for _ in range(10000)]) if len(city) else np.array([])
            output.append({"mode": mode, "metric": metric, "city_years": len(g), "cities": g.city_id.nunique(),
                           "mean_common_months": g.months.mean(), "mean_K": city.mean() if len(city) else np.nan,
                           "ci95_lo_K": np.quantile(boots, .025) if len(boots) else np.nan,
                           "ci95_hi_K": np.quantile(boots, .975) if len(boots) else np.nan,
                           "within_month_C_periods": int(g.within_month_C_periods.sum()),
                           "within_month_D_periods": int(g.within_month_D_periods.sum()),
                           "primary_paired_periods": int(g.primary_paired_periods.sum()),
                           "within_month_paired_periods": int(g.within_month_paired_periods.sum())})
    summary = pd.DataFrame(output)
    summary.to_csv(OUT / "same_source_same_month_window_sensitivity.csv", index=False)
    labels = {"primary_B_minus_C_K": "B−C, start-month windows",
              "within_month_B_minus_C_K": "B−C, within-month windows",
              "primary_C_minus_D_K": "C−D, start-month windows",
              "within_month_C_minus_D_K": "C−D, within-month windows"}
    # 点估计与区间分列、并加一列 Metric：不再用 "—" 占位把两个量塞进同一列
    # （审稿人把旧版读成"列错位/重复行/缺失值"，见 12-response-to-reviewer）。
    lines = ["| Overpass | Window assignment | Metric | Estimate (K) | 95% interval (K) | Paired windows |",
             "|---|---|---|---:|---:|---:|"]
    for mode, overpass in [("day", "Day"), ("night", "Night")]:
        g = city_year[city_year["mode"].eq(mode)]
        for metric in ["primary_B_minus_C_K", "within_month_B_minus_C_K", "primary_C_minus_D_K", "within_month_C_minus_D_K"]:
            r = summary[(summary["mode"].eq(mode)) & (summary["metric"].eq(metric))].iloc[0]
            quantity = "B−C" if "B_minus_C" in metric else "C−D"
            windows = int(r.primary_paired_periods if "primary_" in metric else r.within_month_paired_periods)
            assignment = "Start month" if "primary_" in metric else "Wholly within month"
            lines.append(f"| {overpass} | {assignment} | {quantity} | {r.mean_K:.3f} | [{r.ci95_lo_K:.3f}, {r.ci95_hi_K:.3f}] | {windows:,} |")
    lines.append("\nNote. The alternative excludes any MOD11A2 window whose 8-day span crosses a calendar-month or year boundary. Primary and alternative estimates use the same city-years and calendar months within each overpass (at least 6 common months per city-year); the common panel covers 802 daytime city-years from 37 cities and 645 nighttime city-years from 33 cities, and paired-window counts are summed over those city-years. Each city is resampled once after averaging its repeated city-years, and the interval is a percentile 95% city-cluster bootstrap interval.")
    table_dir = BASE / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    (table_dir / "table7_00_window_assignment.md").write_text("# Table 7. Calendar-month assignment sensitivity\n\n" + "\n".join(lines) + "\n", encoding="utf-8")
    zh = ["| 过境 | 窗口归属 | 指标 | 点估计（K） | 95% 区间（K） | 配对窗口数 |",
          "|---|---|---|---:|---:|---:|"]
    zh_mode = {"day": "白天", "night": "夜间"}
    zh_assign = {"primary_": "起始月", "within_month_": "完全落在月内"}
    for mode in ("day", "night"):
        for metric in ["primary_B_minus_C_K", "within_month_B_minus_C_K", "primary_C_minus_D_K", "within_month_C_minus_D_K"]:
            r = summary[(summary["mode"].eq(mode)) & (summary["metric"].eq(metric))].iloc[0]
            quantity = "B−C" if "B_minus_C" in metric else "C−D"
            assignment = zh_assign["primary_" if metric.startswith("primary_") else "within_month_"]
            windows = int(r.primary_paired_periods if metric.startswith("primary_") else r.within_month_paired_periods)
            zh.append(f"| {zh_mode[mode]} | {assignment} | {quantity} | {r.mean_K:.3f} | "
                      f"[{r.ci95_lo_K:.3f}, {r.ci95_hi_K:.3f}] | {windows:,} |")
    print("\n=== 中文表 7（窗口归月敏感性） ===")
    print("\n".join(zh))
    print("ZH-NOTE-7. 替代口径排除任何跨越日历月或年份边界的 MOD11A2 窗口。两种口径在每个过境时段内使用同一批城市—年与日历月份（每城市—年至少 6 个共同月份）；共同面板包含白天 802 个城市—年（37 城）与夜间 645 个城市—年（33 城），配对窗口数为这些城市—年上的城市—合成期观测总数。每座城市先对重复城市—年求均值再重抽样一次，区间为 95% 城市聚类自助法百分位区间。")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
