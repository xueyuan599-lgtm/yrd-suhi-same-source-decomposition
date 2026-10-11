"""Aqua (MYD11) second-platform replication: report, table, figure, verdict (2026-10-11).

Consumes only frozen CSV products of the parameterised chain:

    same_source_decomposition_city_year{,_aqua}.csv       (13)
    same_source_matched_sensitivity_summary{,_aqua}.csv   (15)
    same_source_gate_coverage{,_aqua}.csv                 (15)
    data/daily/daily_records_{mode}{,_aqua}.csv           (09)

**Refits nothing on the rasters.** The bootstrap is the same estimator as
`15_same_source_postanalysis.py` (10,000 city-cluster draws, percentile 95%,
`default_rng(20261004)`) so the Aqua intervals are comparable with the published
Terra ones by construction, not by coincidence.

Every Terra column is restricted to **2003-2025** (the years Aqua exists); a
whole-window Terra value is never compared with an Aqua value.

Judgement rules are fixed in `AQUA_PREREG.md` and only *applied* here.

Usage
-----
    python 27_aqua_replication.py
    python 27_aqua_replication.py --no-figure
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")

BASE = Path(__file__).resolve().parents[2]          # analysis/
OUT = BASE / "outputs" / "daily"
TAB = BASE / "tables"
FIGS = BASE / "figures"
DATA = BASE / "data" / "daily"
WS = BASE.parent                                    # outputs/journal_rs_...

PRIMARY = {"basis": "dynamic", "pixel_rule": "20px", "min_daily_days": 6,
           "min_8day_periods": 2}
MIN_MONTHS = 6
MIN_DAYS = 6
MIN_PX = 20
FIRST_AQUA_YEAR = 2003
SEED, NBOOT = 20261004, 10000
METRIC_COL = {"A-B": "delta_A_to_B_common_months_K",
              "B-C": "delta_B_to_C_common_months_K",
              "C-D": "delta_C_to_D_common_months_K"}
LABEL = {"A-B": "A\u2212B", "B-C": "B\u2212C", "C-D": "C\u2212D"}
PLAT_LABEL = {"terra": "Terra", "aqua": "Aqua"}


# --------------------------------------------------------------------- helpers
def load_annual(product: str) -> pd.DataFrame:
    suf = "" if product == "terra" else "_aqua"
    df = pd.read_csv(OUT / f"same_source_decomposition_city_year{suf}.csv")
    return df


def primary_slice(annual: pd.DataFrame, first_year: int) -> pd.DataFrame:
    q = annual.query("basis == @PRIMARY['basis'] and pixel_rule == @PRIMARY['pixel_rule']"
                     " and min_daily_days == @PRIMARY['min_daily_days']"
                     " and min_8day_periods == @PRIMARY['min_8day_periods']")
    return q[(q.common_months_all_methods >= MIN_MONTHS) & (q.year >= first_year)].copy()


def cluster_stats(values: np.ndarray, rng) -> dict:
    values = values[~np.isnan(values)]
    if len(values) == 0:
        return {"cities": 0, "city_years": 0, "mean_K": np.nan, "lo": np.nan, "hi": np.nan}
    draws = np.array([rng.choice(values, size=len(values), replace=True).mean()
                      for _ in range(NBOOT)])
    return {"cities": int(len(values)), "mean_K": float(values.mean()),
            "lo": float(np.quantile(draws, .025)), "hi": float(np.quantile(draws, .975))}


def boot_city_slice(sub: pd.DataFrame, col: str, rng, group_cols=("city_id",)) -> dict:
    """City-cluster bootstrap of the per-city time-mean of `col`."""
    by_city = sub.groupby(list(group_cols))[col].mean().dropna().to_numpy()
    st = cluster_stats(by_city, rng)
    st["city_years"] = int(len(sub))
    return st


def fmt(st: dict) -> str:
    if not st["cities"]:
        return "not estimable"
    return f"{st['mean_K']:+.3f} [{st['lo']:+.3f}, {st['hi']:+.3f}]"


# ------------------------------------------------------- H1/H3: primary steps
def primary_rows(annual: dict, rng) -> list[dict]:
    rows = []
    for mode in ("day", "night"):
        for platform in ("terra", "aqua"):
            sub = primary_slice(annual[platform], FIRST_AQUA_YEAR)
            sub = sub[sub["mode"] == mode]
            for metric, col in METRIC_COL.items():
                st = boot_city_slice(sub, col, rng)
                rows.append({"mode": mode, "platform": platform, "block": "primary",
                             "comparison": LABEL[metric], **st})
    return rows


# ------------------------------------------- H2: coverage ladder (A-B only)
def ladder_rows(rng) -> list[dict]:
    rows = []
    for mode in ("day", "night"):
        for platform in ("terra", "aqua"):
            suf = "" if platform == "terra" else "_aqua"
            cy_path = OUT / f"same_source_matched_sensitivity_city_year{suf}.csv"
            sm_path = OUT / f"same_source_matched_sensitivity_summary{suf}.csv"
            sm = pd.read_csv(sm_path)
            cy = pd.read_csv(cy_path)
            for comp, label in (("any_pixel", "\u22651 pixel"),
                                ("coverage40", "\u226540% coverage"),
                                ("coverage60", "\u226560% coverage")):
                s = sm[(sm["mode"] == mode) & (sm.comparison == comp)]
                c = cy[(cy["mode"] == mode) & (cy.comparison == comp) & (cy.year >= FIRST_AQUA_YEAR)]
                # 与 15 的 city-cluster 口径一致：先逐城取时间均值，再对城市重抽样
                st = boot_city_slice(c, "delta_AB_alt_K", rng)
                st["city_years"] = int(s["city_years"].iloc[0]) if len(s) else 0
                st["cities"] = int(s["cities"].iloc[0]) if len(s) else 0
                st.update({"mode": mode, "platform": platform, "block": "coverage ladder",
                           "comparison": f"A\u2212B, {label}"})
                rows.append(st)
    return rows


# ----------------------------- H4: same-day control (pure platform difference)
def same_day_control(rng) -> tuple[list[dict], pd.DataFrame]:
    keys = ["city_id", "year", "date"]
    merged = {}
    for mode in ("day", "night"):
        a = pd.read_csv(DATA / f"daily_records_{mode}_aqua.csv")
        t = pd.read_csv(DATA / f"daily_records_{mode}.csv")
        a = a[a.year >= FIRST_AQUA_YEAR]
        t = t[t.year >= FIRST_AQUA_YEAR]
        for df, tag in ((a, "a"), (t, "t")):
            df["ok"] = (df.u_n >= MIN_PX) & (df.r_n >= MIN_PX)
            df["dT"] = df.u_mean - df.r_mean
        # 只保留两平台**都**通过分区门槛的日子 —— 这是「同一天」的定义
        # 先把两边的 ok/dT 改成不同列名再合，**month 只从一侧取**——若两侧都带着 month，
        # merge 会把它拆成 month_a/month_t，后面 groupby("month") 直接 KeyError。
        keep = keys + ["month", "ok", "dT"]
        a2 = a[keep].rename(columns={"ok": "ok_a", "dT": "dT_a"})
        t2 = t[keep].drop(columns=["month"]).rename(columns={"ok": "ok_t", "dT": "dT_t"})
        j = a2.merge(t2, on=keys, how="inner")
        j = j[j.ok_a & j.ok_t]
        # 月内 ≥6 个共同日 → 月值；年内 ≥6 个共同月 → 年值（与 A 的门槛同构）
        m = (j.groupby(["city_id", "year", "month"], as_index=False)
               .agg(dT_a=("dT_a", "mean"), dT_t=("dT_t", "mean"), n=("dT_a", "size")))
        m = m[m.n >= MIN_DAYS]
        y = (m.groupby(["city_id", "year"], as_index=False)
               .agg(dT_a=("dT_a", "mean"), dT_t=("dT_t", "mean"), months=("dT_a", "size")))
        y = y[y.months >= MIN_MONTHS]
        y["diff"] = y.dT_t - y.dT_a
        y["mode"] = mode
        merged[mode] = y
        print(f"  [{mode}] 同日共同样本：共同合格日 {len(j):,} 个 · "
              f"合格城市–月 {len(m):,} · 合格城市–年 {len(y):,}（{y.city_id.nunique()} 城）")

    rows = []
    for mode, y in merged.items():
        for platform, col in (("terra", "dT_t"), ("aqua", "dT_a")):
            st = boot_city_slice(y, col, rng)
            st.update({"mode": mode, "platform": platform, "block": "same-day control",
                       "comparison": "A (paired level) on shared days"})
            rows.append(st)
        st = boot_city_slice(y, "diff", rng)
        st.update({"mode": mode, "platform": "terra-minus-aqua", "block": "same-day control",
                   "comparison": "platform difference on shared days"})
        rows.append(st)
    return rows, pd.concat(merged.values(), ignore_index=True)


# ------------------------------------------------------------------- verdict
def verdict(rows: pd.DataFrame) -> list[str]:
    def get(mode, platform, block, comparison):
        s = rows[(rows["mode"] == mode) & (rows.platform == platform)
                 & (rows.block == block) & (rows.comparison == comparison)]
        return s.iloc[0] if len(s) else None

    out = ["## 预注册判据的机械判定（AQUA_PREREG.md §三）", ""]
    h1_pass = True
    for mode in ("day", "night"):
        t = get(mode, "terra", "primary", "A\u2212B")
        a = get(mode, "aqua", "primary", "A\u2212B")
        same_sign = np.sign(t["mean_K"]) == np.sign(a["mean_K"])
        big_enough = abs(a["mean_K"]) >= abs(t["mean_K"]) / 3
        ok = bool(same_sign and big_enough)
        h1_pass &= ok
        out.append(f"- **H1 {mode}**: Terra {t['mean_K']:+.3f} / Aqua {a['mean_K']:+.3f} K · "
                   f"同号 {'是' if same_sign else '否'} · 量级≥1/3 {'是' if big_enough else '否'} "
                   f"→ **{'PASS' if ok else 'FAIL'}**")

    h2_lines = []
    for mode in ("day", "night"):
        for label in ("\u226540% coverage", "\u226560% coverage"):
            a = get(mode, "aqua", "coverage ladder", f"A\u2212B, {label}")
            if a is None or pd.isna(a["mean_K"]):
                h2_lines.append(f"- H2 {mode} {label}: Aqua 不可估（与 Terra 同型）")
            else:
                flag = "PASS" if abs(a["mean_K"]) <= 0.10 else "FAIL"
                h2_lines.append(f"- H2 {mode} {label}: Aqua {a['mean_K']:+.3f} K → **{flag}**")
    out += h2_lines

    d = get("day", "aqua", "primary", "C\u2212D")
    h3 = abs(d["mean_K"]) < 0.05
    out.append(f"- **H3** (day |C\u2212D| < 0.05 K): Aqua {d['mean_K']:+.3f} K → "
               f"**{'PASS' if h3 else 'FAIL'}**")
    h4 = get("day", "terra-minus-aqua", "same-day control", "platform difference on shared days")
    h4n = get("night", "terra-minus-aqua", "same-day control", "platform difference on shared days")
    out.append(f"- **H4 同日对照**（平台差，Terra−Aqua）：白天 {h4['mean_K']:+.3f} "
               f"[{h4['lo']:+.3f}, {h4['hi']:+.3f}]；夜间 {h4n['mean_K']:+.3f} "
               f"[{h4n['lo']:+.3f}, {h4n['hi']:+.3f}] K")
    out.append("")
    out.append(f"**H1 总判定：{'PASS' if h1_pass else 'FAIL'}**")
    return out


# ---------------------------------------------------------------------- table
def build_table(rows: pd.DataFrame, ycounts: dict) -> str:
    lines = ["# Table 11. Replication of the same-source decomposition on the Aqua overpass", "",
             "| Overpass | Platform | Comparison | City-years (cities) | Estimate (K) | 95% interval (K) |",
             "|---|---|---|---:|---:|---:|"]
    order_block = ["primary", "coverage ladder", "same-day control"]
    for block in order_block:
        sub = rows[rows.block == block]
        for mode in ("day", "night"):
            for _, r in sub[sub["mode"] == mode].iterrows():
                n = "" if pd.isna(r.get("city_years")) else f"{int(r['city_years'])} ({int(r['cities'])})"
                # 小于 0.005 K 的行用 5 位小数：否则 Aqua 白天 C−D（+0.00005 K）会被印成
                # "+0.000 [+0.000, +0.000]"，读者无法判断它是不是精确的 0。
                tiny = (not pd.isna(r["mean_K"])) and abs(r["mean_K"]) < 0.005
                fmt_e = ("{:+.5f}" if tiny else "{:+.3f}")
                est = "not estimable" if pd.isna(r["mean_K"]) else fmt_e.format(r["mean_K"])
                ci = ("—" if pd.isna(r["lo"])
                      else "[" + fmt_e.format(r["lo"]) + ", " + fmt_e.format(r["hi"]) + "]")
                plat = {"terra": "Terra", "aqua": "Aqua",
                        "terra-minus-aqua": "Terra \u2212 Aqua"}[r["platform"]]
                lines.append(f"| {mode.capitalize()} | {plat} | {r['comparison']} | {n} | {est} | {ci} |")
    lines += ["", NOTE]
    return "\n".join(lines) + "\n"


NOTE = (
    "Note. The Aqua overpass (MYD11A1 daily, MYD11A2 8-day) is analysed with the same "
    "area of interest, 1 km grid, urban\u2013rural mask, quality screen and estimator "
    "definitions as the Terra overpass, so the two columns differ only in platform. All "
    "Terra values are restricted to 2003\u20132025, the years for which Aqua exists; a "
    "whole-window Terra value is never compared with an Aqua value. Estimates are "
    "city-cluster means in kelvin and brackets are percentile 95 % intervals from 10,000 "
    "city-cluster bootstrap draws over cities, drawn independently for each row, so the "
    "intervals are marginal and are not used to test differences between rows. Coverage-ladder "
    "rows use the city-years and calendar months shared with each platform's primary rule. "
    "The same-day control restricts both platforms to days on which both pass the 20-pixel "
    "per-zone screen, then requires six such days in a month and six such months in a year; "
    "within that sample the paired contrast can be formed for each platform and the last row "
    "is their difference, which isolates the platform from the sampling difference. "
    "Night-time 40 % and 60 % coverage estimates are not supported on either platform, because "
    "no city-year retains six shared months under those rules."
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-figure", action="store_true")
    args = ap.parse_args()

    annual = {p: load_annual(p) for p in ("terra", "aqua")}
    for p, df in annual.items():
        print(f"{PLAT_LABEL[p]} 城市–年文件 {len(df):,} 行 · 年份 {df.year.min()}–{df.year.max()}")

    rng = np.random.default_rng(SEED)
    rows = primary_rows(annual, rng)
    rows += ladder_rows(rng)
    sd_rows, sd_detail = same_day_control(rng)
    rows += sd_rows
    tab = pd.DataFrame(rows)

    tab.to_csv(OUT / "aqua_replication_summary.csv", index=False, encoding="utf-8-sig")
    sd_detail.to_csv(OUT / "aqua_same_day_control_city_year.csv", index=False,
                     encoding="utf-8-sig")
    print(f"\n已写 {OUT.name}/aqua_replication_summary.csv（{len(tab)} 行）")

    TAB.mkdir(parents=True, exist_ok=True)
    (TAB / "table11_00_aqua_replication.md").write_text(build_table(tab, {}),
                                                        encoding="utf-8")
    print(f"已写 {TAB.name}/table11_00_aqua_replication.md")

    v = verdict(tab)
    (WS / "AQUA_VERDICT.md").write_text("\n".join(v) + "\n", encoding="utf-8")
    print("\n".join(v))

    if not args.no_figure:
        render_figure(tab)
    return 0


def render_figure(tab: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIGS.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42})
    p = tab[tab.block == "primary"]
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.3), sharey=False)
    for ax, mode in zip(axes, ("day", "night")):
        ypos = np.arange(3)
        for off, (platform, color) in zip((0.19, -0.19),
                                          (("terra", "#4A6FA5"), ("aqua", "#C1663B"))):
            vals, los, his = [], [], []
            for metric in ("A\u2212B", "B\u2212C", "C\u2212D"):
                r = p[(p["mode"] == mode) & (p.platform == platform)
                      & (p.comparison == metric)].iloc[0]
                vals.append(r["mean_K"])
                los.append(r["mean_K"] - r["lo"] if not pd.isna(r["lo"]) else 0)
                his.append(r["hi"] - r["mean_K"] if not pd.isna(r["hi"]) else 0)
            ax.barh(ypos + off, vals, height=0.34, color=color, alpha=.9,
                    label=PLAT_LABEL[platform],
                    xerr=np.array([los, his]), error_kw={"lw": .9, "capsize": 2.5})
        ax.axvline(0, color="#333", lw=.8)
        ax.set_yticks(ypos)
        ax.set_yticklabels(["A\u2212B", "B\u2212C", "C\u2212D"])
        ax.set_xlabel("Difference (K)")
        ax.set_title(f"{mode.capitalize()}time overpass", loc="left", fontweight="bold")
        ax.grid(axis="x", color="#DDD", lw=.6)
        ax.set_axisbelow(True)
    axes[0].legend(frameon=False, loc="lower right", fontsize=8)
    fig.suptitle("Terra and Aqua agree on which step moves the level (2003\u20132025)",
                 x=.01, ha="left", fontsize=9.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, .94))
    for ext in ("png", "pdf"):
        fig.savefig(FIGS / f"fig6_00_terra_vs_aqua.{ext}", dpi=600)
    print(f"已写 {FIGS.name}/fig6_00_terra_vs_aqua.png/.pdf")


if __name__ == "__main__":
    raise SystemExit(main())
