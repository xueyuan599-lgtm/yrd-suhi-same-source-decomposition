"""Post-process the prespecified same-source daily/8-day decomposition.

This script does not re-estimate raster values. It derives aggregation and
fixed-season sensitivities from the full city-month panels, restricting all
methods to identical city-years/months where specified.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, theilslopes, kendalltau

BASE = Path(__file__).resolve().parents[2]
DATA = BASE / "data" / "daily"
OUT = BASE / "outputs" / "daily"
ID = ["basis", "mode", "year", "city_id", "pixel_rule", "min_daily_days", "min_8day_periods"]
METHODS = ["A_paired_daily_K", "B_independent_dates_K", "C_daily_rebuilt_8day_K", "D_native_MOD11A2_K"]
PRIMARY = {"basis": "dynamic", "pixel_rule": "20px", "min_daily_days": 6, "min_8day_periods": 2}


def sen_summary(yearly, value):
    q = yearly[["year", value]].dropna().sort_values("year")
    if len(q) < 8:
        return {"n_years": len(q), "Sen_K_per_decade": np.nan, "CI95_lo": np.nan,
                "CI95_hi": np.nan, "Kendall_tau": np.nan, "Kendall_p": np.nan}
    fit = theilslopes(q[value].to_numpy(), q.year.to_numpy(), alpha=.95)
    kt = kendalltau(q.year, q[value])
    return {"n_years": len(q), "Sen_K_per_decade": fit.slope * 10,
            "CI95_lo": fit.low_slope * 10, "CI95_hi": fit.high_slope * 10,
            "Kendall_tau": kt.statistic, "Kendall_p": kt.pvalue}


def main():
    monthly = pd.read_csv(OUT / "same_source_decomposition_monthly.csv")
    annual = pd.read_csv(OUT / "same_source_decomposition_city_year.csv")
    primary_m = monthly.query("basis == @PRIMARY['basis'] and pixel_rule == @PRIMARY['pixel_rule'] and min_daily_days == @PRIMARY['min_daily_days'] and min_8day_periods == @PRIMARY['min_8day_periods']").copy()
    primary_a = annual.query("basis == @PRIMARY['basis'] and pixel_rule == @PRIMARY['pixel_rule'] and min_daily_days == @PRIMARY['min_daily_days'] and min_8day_periods == @PRIMARY['min_8day_periods']").copy()

    # Compare equal calendar-month and pooled-date/period estimates on the same
    # months that support all four methods. The date CSV contains the matching
    # daily urban/rural summaries; the 8-day period counts are carried monthly.
    raw_rows = []
    for mode in ["night", "day"]:
        daily = pd.read_csv(DATA / f"daily_records_{mode}.csv")
        daily["year"] = daily.year.astype(int)
        daily["month"] = daily.month.astype(int)
        dm = primary_m[primary_m["mode"].eq(mode)]
        dkeep = daily[daily.u_n.ge(20) & daily.r_n.ge(20)].copy()
        dkeep["A_daily"] = dkeep.u_mean - dkeep.r_mean
        for (year, cid), g in dm.groupby(["year", "city_id"], sort=True):
            common = g.dropna(subset=METHODS)
            if len(common) < 6:
                continue
            months = set(common.month.astype(int))
            d = dkeep[(dkeep.year.eq(year)) & (dkeep.city_id.eq(cid)) & dkeep.month.isin(months)]
            a_raw = d.A_daily.mean() if len(d) else np.nan
            # B uses independent valid-date sets by zone, but the same months.
            db = daily[(daily.year.eq(year)) & (daily.city_id.eq(cid)) & daily.month.isin(months)]
            u = db[db.u_n.ge(20)]
            r = db[db.r_n.ge(20)]
            b_raw = (u.u_mean.mean() - r.r_mean.mean()) if len(u) and len(r) else np.nan
            raw_rows.append({"mode": mode, "year": int(year), "city_id": int(cid),
                             "common_months": len(common),
                             "A_month_equal_K": common.A_paired_daily_K.mean(),
                             "B_month_equal_K": common.B_independent_dates_K.mean(),
                             "C_month_equal_K": common.C_daily_rebuilt_8day_K.mean(),
                             "D_month_equal_K": common.D_native_MOD11A2_K.mean(),
                             "A_pooled_days_K": a_raw, "B_pooled_dates_K": b_raw,
                             "C_pooled_periods_K": np.average(common.C_daily_rebuilt_8day_K, weights=common.C_n_periods),
                             "D_pooled_periods_K": np.average(common.D_native_MOD11A2_K, weights=common.D_n_periods),
                             "delta_A_equal_minus_raw_K": common.A_paired_daily_K.mean() - a_raw,
                             "delta_B_equal_minus_raw_K": common.B_independent_dates_K.mean() - b_raw,
                             "delta_C_equal_minus_raw_K": common.C_daily_rebuilt_8day_K.mean() - np.average(common.C_daily_rebuilt_8day_K, weights=common.C_n_periods),
                             "delta_D_equal_minus_raw_K": common.D_native_MOD11A2_K.mean() - np.average(common.D_native_MOD11A2_K, weights=common.D_n_periods)})
    raw = pd.DataFrame(raw_rows)
    raw.to_csv(OUT / "same_source_aggregation_sensitivity_city_year.csv", index=False)
    raw_year = raw.groupby(["mode", "year"], as_index=False).agg(
        city_years=("city_id", "nunique"), mean_months=("common_months", "mean"),
        A_equal=("A_month_equal_K", "mean"), A_raw=("A_pooled_days_K", "mean"),
        B_equal=("B_month_equal_K", "mean"), B_raw=("B_pooled_dates_K", "mean"),
        C_equal=("C_month_equal_K", "mean"), C_raw=("C_pooled_periods_K", "mean"),
        D_equal=("D_month_equal_K", "mean"), D_raw=("D_pooled_periods_K", "mean"),
        dA=("delta_A_equal_minus_raw_K", "mean"), dB=("delta_B_equal_minus_raw_K", "mean"),
        dC=("delta_C_equal_minus_raw_K", "mean"), dD=("delta_D_equal_minus_raw_K", "mean"))
    raw_year.to_csv(OUT / "same_source_aggregation_sensitivity_year.csv", index=False)

    # Three-month fixed-season panels; seasonal weights do not vary over years.
    season_rows = []
    for mode in ["night", "day"]:
        q = primary_m[primary_m["mode"].eq(mode)]
        for season, months in [("JJA", [6, 7, 8]), ("DJF", [12, 1, 2])]:
            z = q[q.month.isin(months)].dropna(subset=METHODS + ["C_minus_D_paired_K"])
            for (year, cid), g in z.groupby(["year", "city_id"], sort=True):
                if set(g.month.astype(int)) != set(months):
                    continue
                vals = g[METHODS].mean()
                season_rows.append({"mode": mode, "season": season, "year": int(year), "city_id": int(cid),
                                    "n_months": 3, **vals.to_dict(),
                                    "A_minus_B_K": vals[METHODS[0]] - vals[METHODS[1]],
                                    "B_minus_C_K": vals[METHODS[1]] - vals[METHODS[2]],
                                    "C_minus_D_paired_K": g.C_minus_D_paired_K.mean()})
    season = pd.DataFrame(season_rows)
    season.to_csv(OUT / "same_source_fixed_season_city_year.csv", index=False)
    season_year = season.groupby(["mode", "season", "year"], as_index=False).agg(
        cities=("city_id", "nunique"), **{c: (c, "mean") for c in METHODS + ["A_minus_B_K", "B_minus_C_K", "C_minus_D_paired_K"]})
    season_year.to_csv(OUT / "same_source_fixed_season_year.csv", index=False)

    # Gate coverage and available-month trend diagnostics. Available-month
    # trends are explicitly secondary because city/month composition varies.
    gate_rows = []
    for keys, g in annual.groupby(["basis", "mode", "pixel_rule", "min_daily_days", "min_8day_periods"]):
        basis, mode, rule, md, mc = keys
        n12 = int((g.months_available == 12).sum())
        n6 = int(g.common_months_all_methods.ge(6).sum())
        gate_rows.append({"basis": basis, "mode": mode, "pixel_rule": rule, "min_daily_days": md,
                          "min_8day_periods": mc, "city_years": len(g), "complete_12_month_city_years": n12,
                          "common_months_ge6_city_years": n6})
    pd.DataFrame(gate_rows).to_csv(OUT / "same_source_gate_coverage.csv", index=False)

    # Re-estimate each sensitivity on the exact city-month intersection with
    # the primary rule so changes are not driven by a different calendar mix.
    variants = [("daily_min3", "dynamic", "20px", 3, 2),
                ("daily_min10", "dynamic", "20px", 10, 2),
                ("composite_min1", "dynamic", "20px", 6, 1),
                ("composite_min3", "dynamic", "20px", 6, 3),
                ("any_pixel", "dynamic", "any", 6, 2),
                ("coverage40", "dynamic", "40pct", 6, 2),
                ("coverage60", "dynamic", "60pct", 6, 2),
                ("fixed_2020_mask", "fixed2020", "20px", 6, 2)]
    match_rows = []
    match_cols = METHODS + ["C_minus_D_paired_K"]
    p_month = monthly.query("basis == 'dynamic' and pixel_rule == '20px' and min_daily_days == 6 and min_8day_periods == 2")
    for mode in ["night", "day"]:
        main_m = p_month[p_month["mode"].eq(mode)].dropna(subset=match_cols)
        for label, basis, rule, md, mc in variants:
            alt_m = monthly.query("basis == @basis and mode == @mode and pixel_rule == @rule and min_daily_days == @md and min_8day_periods == @mc").dropna(subset=match_cols)
            join_cols = ["year", "city_id", "month", "mode"] + match_cols
            joined = main_m[join_cols].merge(alt_m[join_cols],
                         on=["year", "city_id", "month", "mode"], suffixes=("_main", "_alt"))
            for (year, city), g in joined.groupby(["year", "city_id"], sort=True):
                if len(g) < 6:
                    continue
                rec = {"mode": mode, "comparison": label, "year": int(year), "city_id": int(city), "n_common_months": len(g)}
                for col in match_cols:
                    rec[col.replace("_K", "_main_K")] = g[col + "_main"].mean()
                    rec[col.replace("_K", "_alt_K")] = g[col + "_alt"].mean()
                rec["delta_AB_main_K"] = rec["A_paired_daily_main_K"] - rec["B_independent_dates_main_K"]
                rec["delta_AB_alt_K"] = rec["A_paired_daily_alt_K"] - rec["B_independent_dates_alt_K"]
                rec["delta_BC_main_K"] = rec["B_independent_dates_main_K"] - rec["C_daily_rebuilt_8day_main_K"]
                rec["delta_BC_alt_K"] = rec["B_independent_dates_alt_K"] - rec["C_daily_rebuilt_8day_alt_K"]
                match_rows.append(rec)
    matched = pd.DataFrame(match_rows)
    matched.to_csv(OUT / "same_source_matched_sensitivity_city_year.csv", index=False)
    summaries = []
    for (mode, label), g in matched.groupby(["mode", "comparison"]):
        row = {"mode": mode, "comparison": label, "city_years": len(g), "cities": g.city_id.nunique(),
               "mean_common_months": g.n_common_months.mean()}
        for col in ["delta_AB", "delta_BC", "C_minus_D_paired"]:
            for suffix in ["main_K", "alt_K"]:
                name = f"{col}_{suffix}"
                if name in g:
                    row[f"mean_{name}"] = g[name].mean()
        summaries.append(row)
    pd.DataFrame(summaries).to_csv(OUT / "same_source_matched_sensitivity_summary.csv", index=False)

    # City-cluster bootstrap intervals for the pooled common-month level
    # differences. Each city contributes one time-mean, avoiding treating
    # repeated city-years as independent observations.
    boot_rows = []
    rng = np.random.default_rng(20261004)
    metrics = {"A_minus_B": "delta_A_to_B_common_months_K",
               "B_minus_C": "delta_B_to_C_common_months_K",
               "C_minus_D": "delta_C_to_D_common_months_K"}
    for mode in ["night", "day"]:
        h = primary_a[(primary_a["mode"].eq(mode)) & primary_a.common_months_all_methods.ge(6)].copy()
        for label, col in metrics.items():
            by_city = h.groupby("city_id")[col].mean().dropna()
            vals = by_city.to_numpy()
            draws = np.array([rng.choice(vals, size=len(vals), replace=True).mean() for _ in range(10000)])
            boot_rows.append({"mode": mode, "metric": label, "cities": len(vals), "mean_K": vals.mean(),
                              "median_city_mean_K": np.median(vals), "cluster_bootstrap_95_lo_K": np.quantile(draws, .025),
                              "cluster_bootstrap_95_hi_K": np.quantile(draws, .975)})
    pd.DataFrame(boot_rows).to_csv(OUT / "same_source_cluster_bootstrap.csv", index=False)

    # Correlations across city-level means and year-level means use only the
    # common-month panel. They describe rank agreement, not equality.
    corr_rows, trend_rows = [], []
    for mode in ["night", "day"]:
        h = primary_a[(primary_a["mode"].eq(mode)) & primary_a.common_months_all_methods.ge(6)].copy()
        cols = ["A_paired_daily_common_months_mean_K", "B_independent_dates_common_months_mean_K",
                "C_daily_rebuilt_8day_common_months_mean_K", "D_native_MOD11A2_common_months_mean_K"]
        for level, grouped in [("city", h.groupby("city_id")[cols].mean()),
                               ("year", h.groupby("year")[cols].mean())]:
            for left, right, label in [(cols[0], cols[1], "A_vs_B"), (cols[1], cols[2], "B_vs_C"),
                                       (cols[2], cols[3], "C_vs_D")]:
                ok = grouped[[left, right]].dropna()
                rho, p = spearmanr(ok[left], ok[right])
                corr_rows.append({"mode": mode, "level": level, "comparison": label,
                                  "n": len(ok), "spearman_rho": rho, "p_value": p})
        # Available-month annualized trends, explicitly not the primary annual
        # estimand: only city-years with >=6 common months are retained.
        yearly = h.groupby("year", as_index=False).agg(
            cities=("city_id", "nunique"), **{c: (c, "mean") for c in cols},
            delta_A_to_B=("delta_A_to_B_common_months_K", "mean"),
            delta_B_to_C=("delta_B_to_C_common_months_K", "mean"),
            delta_C_to_D=("delta_C_to_D_common_months_K", "mean"))
        for metric in cols + ["delta_A_to_B", "delta_B_to_C", "delta_C_to_D"]:
            eligible = yearly[yearly.cities.ge(5)]
            sm = sen_summary(eligible, metric)
            trend_rows.append({"mode": mode, "metric": metric, "eligible_years": sm.pop("n_years"),
                               **sm, "note": "secondary: >=6 shared estimable months; city and month composition varies"})
    pd.DataFrame(corr_rows).to_csv(OUT / "same_source_rank_correlations.csv", index=False)
    pd.DataFrame(trend_rows).to_csv(OUT / "same_source_available_month_trends.csv", index=False)
    print("Wrote aggregation, season, gate, rank-correlation and secondary trend outputs to", OUT)


if __name__ == "__main__":
    main()
