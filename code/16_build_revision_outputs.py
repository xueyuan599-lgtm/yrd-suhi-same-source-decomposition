"""Build the three revised manuscript tables and two evidence figures."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[2]
OUT = BASE / "outputs" / "daily"
TABLES = BASE / "tables"
FIGS = BASE / "figures"
TABLES.mkdir(parents=True, exist_ok=True)
FIGS.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20261004)


def boot_mean(values, draws=10000):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if not len(v):
        return np.nan, np.nan, np.nan
    b = np.asarray([RNG.choice(v, len(v), replace=True).mean() for _ in range(draws)])
    return float(v.mean()), float(np.quantile(b, .025)), float(np.quantile(b, .975))


def city_mean(g, column):
    return g.groupby("city_id")[column].mean().dropna().to_numpy()


def write_table1():
    txt = """| Component | Operational definition | Weight and minimum support | Interpretation |
|---|---|---|---|
| Target quantity | Terra overpass urban-minus-rural mean LST on a calendar day when both zones pass the observation screen | Average daily differences within each calendar month; month estimates require at least 6 qualifying dates | Clear-sky satellite surface contrast; not an all-weather temperature difference or canopy-air UHI |
| A: same-day paired | Subtract rural from urban mean on each date valid in both zones | At least 20 valid 1-km cells in each zone on a date; equal weight across qualifying dates within month | Operational paired-date reference |
| B: independent dates | Average urban and rural LST on their own valid dates, then subtract | Same 20-cell screen and at least 6 dates for each zone in the month | Isolates non-overlapping valid-date sets relative to A |
| C: daily-rebuilt 8-day | Average the screened daily LST values at each pixel within the MOD11A2 window, then calculate the zone contrast | At least 2 qualifying windows assigned to the month in which each 8-day window starts | Adds compositing and period weighting relative to B |
| D: native MOD11A2 | Calculate the contrast from the NASA Terra 8-day product | At least 2 qualifying windows assigned by start month; compared with C on paired city-window support | Product-output comparison; not a ground-truth error |
| Annual support | Average estimable monthly contrasts with equal calendar-month weight | Primary complete-year trend requires all 12 months; partial-year summaries require at least 6 common months and are labelled secondary | Missing months are not imputed |

Note. Main spatial rule is at least 20 valid 1-km cells in each urban and rural zone. Daily MOD11A1 and native MOD11A2 use the local export's QC screen `(bits 0–1) ≤ 1`; the local daily rasters for 2020 and 2023 contain no raw QC band. Eight-day windows are assigned to their start month; Table 4 tests a rule retaining only windows wholly inside one month. Day and night are analysed separately.
"""
    (TABLES / "table1_00_estimand.md").write_text("# Table 1. Estimand and estimator definitions\n\n" + txt, encoding="utf-8")


def build_main_results():
    annual = pd.read_csv(OUT / "same_source_decomposition_city_year.csv")
    paired_boot = pd.read_csv(OUT / "same_source_cluster_bootstrap.csv").set_index(["mode", "metric"])
    q = annual.query("basis == 'dynamic' and pixel_rule == '20px' and min_daily_days == 6 and min_8day_periods == 2")
    levels = ["A_paired_daily_common_months_mean_K", "B_independent_dates_common_months_mean_K",
              "C_daily_rebuilt_8day_common_months_mean_K", "D_native_MOD11A2_common_months_mean_K"]
    steps = [("delta_A_to_B_common_months_K", "A_minus_B", "A−B"),
             ("delta_B_to_C_common_months_K", "B_minus_C", "B−C"),
             ("delta_C_to_D_common_months_K", "C_minus_D", "C−D")]
    all_cols = levels + [c for c, _, _ in steps]
    estimates, support = {}, {}
    for mode in ["day", "night"]:
        g = q[(q["mode"] == mode) & q.common_months_all_methods.ge(6)].copy()
        # ⚠ 调用顺序与旧版（A,B,C,D,A−B,B−C,C−D）逐字一致：boot_mean 用的是模块级 RNG，
        # 少调或多调一次都会让后续区间的抽样流错位、改变已报告的末位数字。
        means = {c: boot_mean(city_mean(g, c)) for c in all_cols}
        # 三步差值改用归档的城市聚类自助表，使表 3 与表 4 的主行同源。
        for col, metric, _ in steps:
            b = paired_boot.loc[(mode, metric)]
            means[col] = (float(b.mean_K), float(b.cluster_bootstrap_95_lo_K), float(b.cluster_bootstrap_95_hi_K))
        estimates[mode] = means
        support[mode] = {"city_years": len(g), "cities": g.city_id.nunique(),
                         "months": g.common_months_all_methods.mean(),
                         "complete": int((q[q["mode"] == mode].months_available == 12).sum())}

    def cell(t):
        return f"{t[0]:.3f}", f"[{t[1]:.3f}, {t[2]:.3f}]"

    # 表 2：四个估计量的水平。点估计与区间**分列**，避免宽表格被 PDF 抽取时读错。
    lv = ["| Overpass | Quantity | Estimate (K) | 95% interval (K) | City-years (cities) |",
          "|---|---|---:|---:|---:|"]
    for key, label in [("day", "Day"), ("night", "Night")]:
        s = support[key]
        for col, qname in [("A_paired_daily_common_months_mean_K", "A: same-day paired"),
                           ("B_independent_dates_common_months_mean_K", "B: independent dates"),
                           ("C_daily_rebuilt_8day_common_months_mean_K", "C: daily-rebuilt 8-day"),
                           ("D_native_MOD11A2_common_months_mean_K", "D: native MOD11A2")]:
            e, ci = cell(estimates[key][col])
            lv.append(f"| {label} | {qname} | {e} | {ci} | {s['city_years']} ({s['cities']}) |")
    lv.append("\nNote. Each estimate is a city-cluster mean: city-years are averaged within each city and cities then receive equal weight. Brackets are percentile 95% intervals from 10,000 city-cluster bootstrap draws over cities; they describe between-city sampling variation and are not a full uncertainty budget for the product or the retrieval. The panel averages 10.66 daytime and 9.53 nighttime common months, and 253 daytime and 103 nighttime city-years have all 12 months. `C−D` is computed on paired city-window observations: 8,562 daytime and 6,556 nighttime city-months, with 30,999 and 21,070 paired city-window contrasts. No complete-year balanced city panel supports a trend comparison.")
    (TABLES / "table2_00_support_and_levels.md").write_text(
        "# Table 2. Support and estimator levels\n\n" + "\n".join(lv) + "\n", encoding="utf-8")
    zh2 = ["| 过境 | 估计量 | 点估计（K） | 95% 区间（K） | 城市—年（城市数） |",
           "|---|---|---:|---:|---:|"]
    zh_names = {"A_paired_daily_common_months_mean_K": "A：同日配对",
                "B_independent_dates_common_months_mean_K": "B：独立日期",
                "C_daily_rebuilt_8day_common_months_mean_K": "C：逐日重建 8 天",
                "D_native_MOD11A2_common_months_mean_K": "D：原生 MOD11A2"}
    for key, label in [("day", "白天"), ("night", "夜间")]:
        s = support[key]
        for col, _ in [("A_paired_daily_common_months_mean_K", 0),
                       ("B_independent_dates_common_months_mean_K", 0),
                       ("C_daily_rebuilt_8day_common_months_mean_K", 0),
                       ("D_native_MOD11A2_common_months_mean_K", 0)]:
            e, ci = cell(estimates[key][col])
            zh2.append(f"| {label} | {zh_names[col]} | {e} | {ci} | {s['city_years']} ({s['cities']}) |")
    print("\n=== 中文表 2（支持度与四估计量水平） ===")
    print("\n".join(zh2))
    print("ZH-NOTE-2. 每个估计量为城市聚类均值：先在城市内对城市—年求平均，再对城市等权。方括号为按城市重抽样 10,000 次的 95% 百分位区间，反映城市间抽样变异，不是产品与反演的完整不确定度。面板平均共享 10.66（白天）/9.53（夜间）个共同月份，253（白天）/103（夜间）个城市—年拥有完整 12 个月。C−D 在配对的城市—合成期观测上计算，共 8,562（白天）/6,556（夜间）个城市—月、30,999（白天）/21,070（夜间）个配对城市—合成期对比。完整年度不平衡，不足以支持趋势比较。")

    # 表 3：三步差值（同样把点估计与区间分列）。
    df = ["| Overpass | Difference | Estimate (K) | 95% interval (K) |", "|---|---|---:|---:|"]
    for key, label in [("day", "Day"), ("night", "Night")]:
        for col, _, qname in steps:
            e, ci = cell(estimates[key][col])
            df.append(f"| {label} | {qname} | {e} | {ci} |")
    df.append("\nNote. A negative A−B means that separate urban and rural date sets give a larger monthly contrast than same-day pairing. Each interval is a percentile 95% interval from 10,000 city-cluster bootstrap draws over cities; the draws are made independently for each row, so the intervals are marginal and must not be used to test the difference between two rows. No multiplicity adjustment is applied, because the rows are descriptive.")
    (TABLES / "table3_00_step_differences.md").write_text(
        "# Table 3. Step differences\n\n" + "\n".join(df) + "\n", encoding="utf-8")
    zh3 = ["| 过境 | 差值 | 点估计（K） | 95% 区间（K） |", "|---|---|---:|---:|"]
    for key, label in [("day", "白天"), ("night", "夜间")]:
        for col, _, qname in steps:
            e, ci = cell(estimates[key][col])
            zh3.append(f"| {label} | {qname} | {e} | {ci} |")
    print("\n=== 中文表 3（三步差值） ===")
    print("\n".join(zh3))
    print("ZH-NOTE-3. A−B 为负表示城乡各自日期集给出的月度温差大于同日配对。每个区间为按城市重抽样 10,000 次的 95% 百分位区间；各行独立抽样，因此区间是边缘分布，不得用于检验两行之差；本文不做多重比较校正，因为这些行属描述性比较。")
    return estimates


def build_sensitivity():
    s = pd.read_csv(OUT / "same_source_matched_sensitivity_city_year.csv")
    scenarios = [("20px primary", None), ("≥1 pixel", "any_pixel"), ("≥40% mask pixels", "coverage40"),
                 ("≥60% mask pixels", "coverage60"), ("2020 fixed mask", "fixed_2020_mask"),
                 ("10 paired days/month", "daily_min10")]
    metrics = ["delta_AB_alt_K", "delta_BC_alt_K"]
    rows = []
    plotdata = {}
    annual = pd.read_csv(OUT / "same_source_decomposition_city_year.csv")
    primary = annual.query("basis == 'dynamic' and pixel_rule == '20px' and min_daily_days == 6 and min_8day_periods == 2")
    for mode in ["day", "night"]:
        for label, scenario in scenarios:
            if scenario is None:
                g = primary[(primary["mode"] == mode) & primary.common_months_all_methods.ge(6)].copy()
                g = g.assign(delta_AB_alt_K=g.delta_A_to_B_common_months_K,
                             delta_BC_alt_K=g.delta_B_to_C_common_months_K)
                months = g.common_months_all_methods.mean()
            else:
                g = s[(s["mode"] == mode) & (s.comparison == scenario)].copy()
                months = g.n_common_months.mean() if len(g) else np.nan
            if len(g):
                row = {"mode": mode, "scenario": label, "city_years": len(g), "cities": g.city_id.nunique(), "mean_common_months": months}
                d = {}
                for metric in metrics:
                    bycity = g.groupby("city_id")[metric].mean().dropna().to_numpy()
                    est = boot_mean(bycity)
                    row[metric + "_mean"] = est[0]
                    row[metric + "_lo"] = est[1]
                    row[metric + "_hi"] = est[2]
                    d[metric] = est
                plotdata[(mode, label)] = d
                rows.append(row)
            else:
                rows.append({"mode": mode, "scenario": label, "city_years": 0, "cities": 0,
                             "mean_common_months": np.nan, "delta_AB_alt_K_mean": np.nan,
                             "delta_AB_alt_K_lo": np.nan, "delta_AB_alt_K_hi": np.nan,
                             "delta_BC_alt_K_mean": np.nan, "delta_BC_alt_K_lo": np.nan,
                             "delta_BC_alt_K_hi": np.nan})
    out = pd.DataFrame(rows)
    # 主口径行改用归档的城市聚类自助表：与表 3 逐位相同，杜绝同一量在两表出现末位差异。
    # ⚠ 必须在生成表格文本**之前**覆盖——第一版把这段放在 lines 之后，表格仍印旧值。
    boot = pd.read_csv(OUT / "same_source_cluster_bootstrap.csv").set_index(["mode", "metric"])
    for mode in ["day", "night"]:
        i = out.index[(out["mode"] == mode) & (out.scenario == "20px primary")][0]
        for metric, col in [("delta_AB_alt_K", "A_minus_B"), ("delta_BC_alt_K", "B_minus_C")]:
            b = boot.loc[(mode, col)]
            out.loc[i, metric + "_mean"] = float(b.mean_K)
            out.loc[i, metric + "_lo"] = float(b.cluster_bootstrap_95_lo_K)
            out.loc[i, metric + "_hi"] = float(b.cluster_bootstrap_95_hi_K)
    out.to_csv(OUT / "same_source_revision_sensitivity_table.csv", index=False)
    # 长表格式：点估计与区间分列、每行一个量。旧版把两个指标塞进 7 列并用长区间串，
    # 转 PDF 后单元格折行、文本抽取读不回原字符串（2026-10-04 的保真测试抓到 16/20 条），
    # 也正是审稿人把表 2 数字读错的那类原因。
    lines = ["| Overpass | Sensitivity | Metric | Estimate (K) | 95% interval (K) | City-years (cities) |",
             "|---|---|---|---:|---:|---:|"]
    for label, _ in scenarios:
        for mode, overpass in (("day", "Day"), ("night", "Night")):
            r = out[(out["mode"] == mode) & (out.scenario == label)].iloc[0]
            for metric, qname in (("delta_AB_alt_K", "A−B"), ("delta_BC_alt_K", "B−C")):
                if not r.city_years:
                    lines.append(f"| {overpass} | {label} | {qname} | not estimable | not estimable | 0 |")
                else:
                    lines.append(f"| {overpass} | {label} | {qname} | {r[metric+'_mean']:.3f} | "
                                 f"[{r[metric+'_lo']:.3f}, {r[metric+'_hi']:.3f}] | {int(r.city_years)} ({int(r.cities)}) |")
    lines.append("\nNote. Every sensitivity is matched to the primary estimate on identical city-years and calendar months, with at least 6 common months per city-year. Intervals bootstrap cities after averaging each city's repeated city-years. The 20-pixel primary row is the same estimate as the step table, taken from the archived city-cluster bootstrap output, so the two tables cannot differ in the last digit. The 40% and 60% nighttime rows are not estimable because no matched city-year retains 6 common months. A−B and B−C are direct monthly level differences; they are not trend differences.")
    (TABLES / "table4_00_matched_sensitivities.md").write_text("# Table 4. Matched gate and mask sensitivities\n\n" + "\n".join(lines) + "\n", encoding="utf-8")
    zh4 = ["| 过境 | 敏感性口径 | 指标 | 点估计（K） | 95% 区间（K） | 城市—年（城市数） |",
           "|---|---|---|---:|---:|---:|"]
    zh_scen = {"20px primary": "20 像元主口径", "≥1 pixel": "至少 1 个像元",
               "≥40% mask pixels": "掩膜像元至少 40%", "≥60% mask pixels": "掩膜像元至少 60%",
               "2020 fixed mask": "固定 2020 掩膜", "10 paired days/month": "每月至少 10 个配对日"}
    for label, _ in scenarios:
        for mode, overpass in (("day", "白天"), ("night", "夜间")):
            r = out[(out["mode"] == mode) & (out.scenario == label)].iloc[0]
            for metric, qname in (("delta_AB_alt_K", "A−B"), ("delta_BC_alt_K", "B−C")):
                if not r.city_years:
                    zh4.append(f"| {overpass} | {zh_scen[label]} | {qname} | 不可估计 | 不可估计 | 0 |")
                else:
                    zh4.append(f"| {overpass} | {zh_scen[label]} | {qname} | {r[metric+'_mean']:.3f} | "
                               f"[{r[metric+'_lo']:.3f}, {r[metric+'_hi']:.3f}] | {int(r.city_years)} ({int(r.cities)}) |")
    print("\n=== 中文表 4（匹配敏感性，长表） ===")
    print("\n".join(zh4))
    return out, plotdata


def figures(estimates, sensitivity):
    # Step-wise level changes by day/night and sensitivity rule.
    labels = ["20px primary", "≥1 pixel", "≥40% mask pixels", "≥60% mask pixels", "2020 fixed mask", "10 paired days/month"]
    colors = {"delta_AB_alt_K": "#2878B5", "delta_BC_alt_K": "#D17A22"}
    markers = {"delta_AB_alt_K": "o", "delta_BC_alt_K": "s"}
    fig, axs = plt.subplots(1, 2, figsize=(12.4, 5.7), sharex=True)
    for ax, mode, title in zip(axs, ["day", "night"], ["Daytime", "Nighttime"]):
        ax.axvline(0, color="#555555", lw=1, ls="--")
        for i, label in enumerate(labels):
            datum = sensitivity.get((mode, label))
            if datum is None and label == "20px primary":
                # The primary estimate uses the same result as main Table 2.
                em = estimates[mode]
                datum = {"delta_AB_alt_K": em["delta_A_to_B_common_months_K"],
                         "delta_BC_alt_K": em["delta_B_to_C_common_months_K"]}
            if datum is None:
                continue
            for metric, offset in [("delta_AB_alt_K", -.12), ("delta_BC_alt_K", .12)]:
                v = datum.get(metric)
                if v is None or not np.isfinite(v[0]):
                    continue
                ax.errorbar(v[0], i+offset, xerr=[[v[0]-v[1]], [v[2]-v[0]]], color=colors[metric],
                            marker=markers[metric], capsize=2.5, lw=1.6, ms=5,
                            label={"delta_AB_alt_K": "A − B: date pairing", "delta_BC_alt_K": "B − C: 8-day reconstruction"}[metric] if i==0 else None)
        ax.set_yticks(range(len(labels)), labels)
        ax.invert_yaxis(); ax.set_title(title); ax.grid(axis="x", alpha=.2)
        ax.set_xlabel("Monthly level difference (K)")
    fig.legend(*axs[0].get_legend_handles_labels(), frameon=False, loc="lower center",
               bbox_to_anchor=(.5, -.01), ncol=2)
    fig.suptitle("The level decomposition changes with the valid-pixel rule", y=1.02, fontsize=14)
    fig.tight_layout(rect=(0, .08, 1, 1))
    fig.savefig(FIGS / "fig2_00_same_source_decomposition.png", dpi=240, bbox_inches="tight")
    fig.savefig(FIGS / "fig2_00_same_source_decomposition.pdf", bbox_inches="tight")
    plt.close(fig)

    # Annual count of supported complete 12-month city-years.
    annual = pd.read_csv(OUT / "same_source_decomposition_city_year.csv")
    q = annual.query("basis == 'dynamic' and pixel_rule == '20px' and min_daily_days == 6 and min_8day_periods == 2").copy()
    q["complete"] = q.months_available.eq(12)
    supp = q.groupby(["mode", "year"]).agg(complete=("complete", "sum"), city_years=("city_id", "nunique")).reset_index()
    fig, ax = plt.subplots(figsize=(10.5, 4.7))
    for mode, label, color, marker in [("day", "Daytime", "#D17A22", "o"), ("night", "Nighttime", "#2878B5", "s")]:
        g = supp[supp["mode"] == mode]
        ax.plot(g.year, g.complete, label=label, color=color, marker=marker, ms=4, lw=1.8)
    ax.set_ylabel("City-years with all 12 months")
    ax.set_xlabel("Year")
    ax.set_title("Strict annual support is sparse and varies over time")
    ax.set_xticks(range(2001, 2026, 2)); ax.grid(alpha=.2); ax.legend(frameon=False)
    ax.text(.01, -.23, "Counts are out of 24–37 available city-years per year; only 253 daytime and 103 nighttime city-years qualify across 2001–2025.", transform=ax.transAxes, fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGS / "fig3_00_complete_month_support.png", dpi=240, bbox_inches="tight")
    fig.savefig(FIGS / "fig3_00_complete_month_support.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    write_table1()
    estimates = build_main_results()
    sensitivity, plotdata = build_sensitivity()
    figures(estimates, plotdata)
    print("Built Table 1–3 and Figure 2–3 for the revised same-source analysis.")


if __name__ == "__main__":
    main()
