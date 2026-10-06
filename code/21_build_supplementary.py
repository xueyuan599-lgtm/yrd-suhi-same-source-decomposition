"""Assemble the Supplementary Material from the pre-revision archive.

Run after the figures and tables are frozen. Nothing is re-estimated here: the table
bodies are copied byte-for-byte from `tables/_archive_pre_major_revision/*.md` (only the
`**Table N.**` heading is rewritten to the supplement's own numbering), and the two
figures are copied from `figures/_archive_pre_major_revision/`. That keeps the rule
"SI numbers must equal their own CSVs" true by construction.

`supplementary/` is the package that would be uploaded with the manuscript; the copy at the
workspace root is the readable version for the authors.
"""

import shutil
from pathlib import Path

import pandas as pd

WS = Path(__file__).resolve().parents[2]          # = analysis/
ROOT = WS.parent                                  # = outputs/journal_rs_20260928/
TAB_ARCHIVE = WS / "tables" / "_archive_pre_major_revision"
FIG_ARCHIVE = WS / "figures" / "_archive_pre_major_revision"
SUP = ROOT / "supplementary"

# (table file, new number, new caption, source dir: "archive" or "tables")
# 顺序即输出顺序：S1–S4 是检验表，S5 是产品清单，最后列出（2026-10-04 修正：
# 原先 S5 排在 S1 之前，且未写进 Contents，读者按编号找不到入口）。
TABLES = [
    ("table4_coverage_rules.md", "S1",
     "Observation-rule sensitivity of the available-month trend diagnostic. Rows are "
     "coverage rules for the daily and 8-day zonal-mean series; the estimand here is the "
     "city-level trend of the zonal-mean contrast, not the four-estimator decomposition.",
     "archive"),
    ("table5_seasonal_composition.md", "S2",
     "Does the weighting effect survive equalising the seasonal composition? Each row "
     "recomputes both weighting conventions under one aggregation rule, so the two Sen "
     "columns differ only in the weights.",
     "archive"),
    ("table6_between_model.md", "S3",
     "Urban and rural latitude slopes on identical city-years, with a controlled slope "
     "contrast. These are associations, not causal effects.",
     "archive"),
    ("table7_rural_sensitivity.md", "S4",
     "Sensitivity of the zonal-mean urban–rural level to the rural-reference definition on "
     "the same cities and city-years.",
     "archive"),
    ("tableS5_00_data_products.md", "S5",
     "Input products and derived layers, with the version, native resolution and period used "
     "for each. All layers were accessed through the Google Earth Engine Data Catalog.",
     "tables"),
]

FIGURES = [
    ("figS1_morphology_negative_control",
     "Exploratory model-based diagnostic retained from an earlier design stage. (a) "
     "Cross-validated R² for grouped (leave-city-out) and random splits; the random split "
     "is leakage-inflated and is not a performance estimate. (b) Mean |SHAP| values for the "
     "candidate predictors. These results belong to a predictive component that is not part "
     "of the reported estimand and are included for transparency only; no main-text number "
     "depends on them."),
    ("figS2_leakage_gap",
     "Grouped versus random-split gap for the same predictive component, showing how much "
     "of the apparent skill disappears once city identity is held out. As for Figure S1, "
     "this diagnostic is not used to support any main-text claim."),
]

HEAD = """# Supplementary Material

**Manuscript:** Same-Day Pairing and Temporal Aggregation of Terra Urban–Rural Surface
Temperature Contrasts in the Yangtze River Delta

**Contents.** This supplement reports three groups of material. Tables S1–S2 give the
sensitivity of the available-month trend diagnostic to the coverage rule and to seasonal
composition. Tables S3–S4 give the sensitivity of the zonal-mean urban–rural level to the
rural-reference definition and to latitude. Table S5 lists the input products and derived
layers with the version, native resolution and period used for each. Table S6 reports the
available-month trend diagnostic itself, which the main text deliberately keeps out of the
results section: slope, interval, p value and eligible-year count for every estimator and
step. Figures S1–S2 are the exploratory model-based diagnostics kept from an earlier design
stage.

Table S6 exists because the prespecified complete-12-month trend is not estimable in this
record (Section 3.3 of the manuscript). Its rows are exploratory: the eligible-year count is
shown so that a reader can see how much support each slope rests on, and none of its slopes
supports a claim in the main text.

The trend and level checks use the zonal-mean or available-month estimand, which is not the
four-estimator same-source decomposition of Sections 2.3 and 3.1 of the manuscript. They are
reported here for completeness; the main text cites them by number only and does not restate
their values, because a level or trend computed on a zonal-mean series and a step computed
between two date sets are not interchangeable quantities.

## Figures
"""


def main() -> None:
    SUP.mkdir(exist_ok=True)
    for stem, caption in FIGURES:
        for ext in ("png", "pdf"):
            src = FIG_ARCHIVE / f"{stem}.{ext}"
            if src.exists():
                shutil.copy2(src, SUP / src.name)
    body = [HEAD]
    for stem, caption in FIGURES:
        # stem 形如 "figS1_..."：取 "S1" 作编号（不要用 stem[4:6]，那会带上下划线）。
        tag = stem[3:5].upper()
        body.append(f"![Figure {tag}]({stem}.png)\n\n"
                    f"**Figure {tag}.** {caption}\n")
    body.append("## Tables\n")
    for fname, num, caption, where in TABLES:
        src = (TAB_ARCHIVE if where == "archive" else WS / "tables") / fname
        text = src.read_text(encoding="utf-8")
        lines = text.splitlines()
        # 只改题注行（`**Table N.** ...`），表体与表注逐字保留。
        assert lines[0].startswith("**Table "), lines[0]
        lines[0] = f"**Table {num}.** {caption}"
        body.append("\n".join(lines) + "\n")
    # 表 S6：可用月份趋势诊断。正文 §3.3 只保留「预定义趋势不可估计」，数值原样搬到这里。
    tr = pd.read_csv(WS / "outputs" / "daily" / "same_source_available_month_trends.csv")
    nice = {"A_paired_daily_common_months_mean_K": "A: same-day paired",
            "B_independent_dates_common_months_mean_K": "B: independent dates",
            "C_daily_rebuilt_8day_common_months_mean_K": "C: daily-rebuilt 8-day",
            "D_native_MOD11A2_common_months_mean_K": "D: native MOD11A2",
            "delta_A_to_B": "A−B", "delta_B_to_C": "B−C", "delta_C_to_D": "C−D"}
    s6 = ["**Table S6.** Available-month trend diagnostic (exploratory).", "",
          "| Overpass | Quantity | Eligible years | Sen slope (K per decade) | 95% interval (K per decade) | Kendall τ | p |",
          "|---|---|---:|---:|---:|---:|---:|"]
    for mode, label in (("day", "Day"), ("night", "Night")):
        sub = tr[tr["mode"] == mode].set_index("metric")
        for key, qname in nice.items():
            if key not in sub.index:
                continue
            r = sub.loc[key]
            s6.append(f"| {label} | {qname} | {int(r.eligible_years)} | {r.Sen_K_per_decade:+.3f} | "
                      f"[{r.CI95_lo:+.3f}, {r.CI95_hi:+.3f}] | {r.Kendall_tau:+.3f} | {r.Kendall_p:.4f} |")
    s6 += ["", "Note. Slopes are Theil–Sen with Hamed–Rao modified Mann–Kendall p values, fitted to "
               "the city-equal annual series of available-month summaries; city and month composition "
               "changes across years, which is why the eligible-year count is shown. These slopes are "
               "not the prespecified complete-month trend estimand, for which this record has no "
               "balanced multi-year city panel.", ""]
    body.append("\n".join(s6))
    out = "\n".join(body)
    (SUP / "supplementary-material.md").write_text(out, encoding="utf-8")
    (ROOT / "supplementary-material.md").write_text(out, encoding="utf-8")
    n_rows = sum(1 for ln in out.splitlines() if ln.strip().startswith("|"))
    print(f"已写 supplementary/supplementary-material.md 与根目录副本")
    print(f"  图 {len(FIGURES)} · 表 {len(TABLES)} · 表格行 {n_rows}")
    for f in sorted(SUP.iterdir()):
        print(f"  {f.name}  {f.stat().st_size} B")


if __name__ == "__main__":
    main()
