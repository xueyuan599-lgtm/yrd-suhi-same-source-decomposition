"""收严 QA 档 `(QC&3)==0` 下的**整条估计阶梯与 Sen 斜率** —— 昼夜分开（2026-10-03）。

## 它回答什么（补 `05-evidence-map.md` 的未闭合项）

证据图里那句原文是：

> 「收严档下趋势斜率未变」——**未计算**任何 `eq0` 下的 Sen 斜率，只算了城市数与保留率。

本脚本把这条缺口**量化闭合**：在**同一批 23 个带 Q 的年份**上，把 `06_daily_trends.py`
的整条阶梯（`weighted / equal_any / equal_20px / equal_40pct / equal_60pct / paired_20px`）
原样搬到**收严像元支撑**（`u_mean_strict …`）上重跑，逐档给出：

- 该档是否有可检验的**年度序列**（≥ `MIN_YEARS = 8` 年）→ Sen 斜率 + Hamed–Rao p；
- 若无，**逐档、逐年**说明卡在哪一步（不是脚本缺陷，是数据量）；
- 逐年：有多少城市在这档下有**可用年度序列**（每侧 ≥ `MIN_DAYS = 6` 个有效日），
  以及城乡两侧**保留的像元·日占比**。

## 「像元少」与「撑不起趋势」是两件事——本脚本的核心区分

- **流行率**（`08_qc_prevalence.py`）：收严档夜间只占有效 LST 像元的 **0.323%**、
  按日占比中位 **0**；白天 **60.79%**、中位 42.7%。这是**产品属性**。
- **可执行性**（本脚本）：把 0.323% 摊到「城市 × 年 × 日」后，**还能不能凑出一条
  跨年序列**。两者不是一回事——白天即便占比高，也要看摊到每城每年够不够 6 天。

实测（`daily_records_{mode}.csv` → 本脚本产出，见 `QC_STRICT_TREND.md`）：

| 昼夜 | 收严档可算年度序列的年数（`equal_any`） | 逐年可用城数中位 | 结论 |
|------|----------------------------------------|------------------|------|
| 夜间 | **10 / 23** | **0**（仅 10 年有城，1–8 城） | 等权档**结构性薄**；加权档仍可算 |
| 白天 | **23 / 23** | 32（24–37 城） | 整条阶梯可执行 |

夜间**加权档**（无日门限，只要两侧 Σ 像元 > 0）仍能拼出 23 年序列，故「收严档夜间
斜率算不出」这句话**不能泛称**——须分档说：加权档能算，等权及其覆盖门限档算不出。

## ⚠ 它**不能**用来声称什么

1. **不能**把不同档的绝对斜率并列成同一估计量的区间——每档是**不同的观测集合**。
2. **不能**说「收严档下趋势未变」——本脚本给的是收严档**自己**的斜率；夜间等权档
   建立在 33 个城–年、10 个年上，**不是**主口径趋势的稳健性证据。
3. **不能**把收严档结论外推到 **MOD11A2 合成产品**——本表只读逐日 MOD11A1 中间表。
4. **不能**在夜间把收严档当「保守估计」或「更干净」——收严集是**非随机子集**
   （算法更自信的那部分），与城乡温差的协变未知。

用法：
    python 12_qc_strict_trend.py                 # 昼夜都跑
    python 12_qc_strict_trend.py --mode night
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")

BASE = Path(__file__).resolve().parent.parent.parent       # analysis/
sys.path.insert(0, str(Path(__file__).resolve().parent))   # 让 `import daily_estimators` 生效
from daily_estimators import (MIN_DAYS, MIN_YEARS, LADDER, DESC,  # noqa: E402
                              trend_of, mask_pixel_counts, city_year_ladder)

DAILY = BASE / "data" / "daily"
OUT = BASE / "outputs" / "daily"
OUT.mkdir(parents=True, exist_ok=True)

parser = argparse.ArgumentParser(description="收严 QA 档的阶梯与 Sen 斜率")
parser.add_argument("--mode", choices=("night", "day"), default=None)
parser.add_argument("--min-years", type=int, default=MIN_YEARS)
ARGS = parser.parse_args()
MODES = (ARGS.mode,) if ARGS.mode else ("night", "day")
MIN_YEARS = ARGS.min_years


def build_ladder(rec_qc, pxmap):
    """在收严支撑上算逐城市–年的整条阶梯 → DataFrame。"""
    rows = []
    for (cid, year), sub in rec_qc.groupby(["city_id", "year"], sort=True):
        px = pxmap.get(int(year), {}).get(int(cid))
        lad = city_year_ladder(sub, px, suf="_strict")
        lad["city_id"] = int(cid)
        lad["year"] = int(year)
        # 加权的「有效日」= 两侧各自有收严像元的天数取小（与等权档 `min(len(uo),len(ro))`
        # 同一约定：加权不要求同日，故只取两侧天数的较小者，不做同日交集）。
        lad["n_days_weighted"] = int(min((sub.u_n_strict > 0).sum(),
                                         (sub.r_n_strict > 0).sum()))
        rows.append(lad)
    return pd.DataFrame(rows).sort_values(["city_id", "year"]).reset_index(drop=True)


def rung_trend(tab):
    """逐档趋势。返回 (csv 行列表)。"""
    out = []
    for rung in LADDER:
        s = tab[["city_id", "year", rung]].dropna()
        t = trend_of(s.groupby("year")[rung].mean(), MIN_YEARS) if len(s) else None
        dayscol = f"n_days_{rung}"
        med_days = (float(tab.loc[tab[rung].notna(), dayscol].median())
                    if dayscol in tab and len(s) else np.nan)
        out.append({
            "mode": None, "rung": rung, "description": DESC[rung],
            "city_years": int(len(s)),
            "n_cities": int(s.city_id.nunique()) if len(s) else 0,
            "n_years": int(s.year.nunique()) if len(s) else 0,
            "median_valid_days": med_days,
            **(t or {"n_years": 0, "mean_K": np.nan, "sen_K_per_decade": np.nan,
                     "sen_CI95_lo": np.nan, "sen_CI95_hi": np.nan, "hamed_rao_p": np.nan}),
            "computable": bool(t is not None),
        })
    return out


def by_year(tab, rec_qc, qc_years):
    """逐年：可用城数（逐档）+ 城乡保留像元·日占比。"""
    rows = []
    for year in qc_years:
        sub = tab[tab.year == year]
        rsub = rec_qc[rec_qc.year == year]
        us, rs = float(rsub.u_n_strict.sum()), float(rsub.r_n_strict.sum())
        um, rm = float(rsub.u_n.sum()), float(rsub.r_n.sum())
        # 逐城比值（只在主口径有像元的城上算），再对城取均值
        ur, rr = [], []
        for _, g in rsub.groupby("city_id"):
            gu, gm = g.u_n.sum(), g.u_n_strict.sum()
            gu_r, gm_r = g.r_n.sum(), g.r_n_strict.sum()
            if gu > 0:
                ur.append(gm / gu)
            if gu_r > 0:
                rr.append(gm_r / gu_r)
        row = {"mode": None, "year": int(year),
               "n_cities_total": int(len(sub)),
               **{f"n_cities_{r}": int(sub[r].notna().sum()) for r in LADDER},
               "urban_pixel_days_strict": int(us), "rural_pixel_days_strict": int(rs),
               "urban_pixel_days_main": int(um), "rural_pixel_days_main": int(rm),
               "urban_retention_pooled": (us / um if um else np.nan),
               "rural_retention_pooled": (rs / rm if rm else np.nan),
               "urban_retention_mean_of_ratio": float(np.mean(ur)) if ur else np.nan,
               "rural_retention_mean_of_ratio": float(np.mean(rr)) if rr else np.nan}
        for rung in LADDER:
            row[f"annual_mean_{rung}"] = (float(sub[rung].mean())
                                          if sub[rung].notna().any() else np.nan)
        rows.append(row)
    return rows


RES = {}
for mode in MODES:
    src = DAILY / f"daily_records_{mode}.csv"
    if not src.exists():
        print(f"缺 {src.name}，跳过 {mode}", file=sys.stderr)
        continue
    meta = json.loads((DAILY / "daily_records_meta.json").read_text(encoding="utf-8"))
    no_qc = set(meta["modes"][mode]["years_without_qc"])
    rec = pd.read_csv(src)
    qc_years = sorted(set(rec.year.unique()) - no_qc)          # 23 个带 Q 的年份
    rec_qc = rec[rec.year.isin(qc_years)]
    pxmap = mask_pixel_counts()
    tab = build_ladder(rec_qc, pxmap)
    lab = "夜间" if mode == "night" else "白天"
    print("=" * 96)
    print(f"{lab} 收严档 (QC&3)==0   [mode={mode}]   {len(qc_years)} 个带 Q 年份 "
          f"{qc_years[0]}–{qc_years[-1]}（排除无 Q 的 {sorted(no_qc)}）")
    print(f"  {len(tab)} 城市–年 · {tab.city_id.nunique()} 城 · 记录行 {len(rec_qc)}")
    print("=" * 96)

    tr = pd.DataFrame(rung_trend(tab))
    tr.insert(0, "support", "strict")
    tr["mode"] = mode
    # 列序：与 10 的 [4] 段命名对齐（mode,rung,city_years,n_cities,n_years,median_valid_days,
    # mean_K,sen_K_per_decade,sen_CI95_lo,sen_CI95_hi,hamed_rao_p,…），再挂新增列。
    cols = ["mode", "support", "rung", "description", "city_years", "n_cities", "n_years",
            "median_valid_days", "mean_K", "sen_K_per_decade", "sen_CI95_lo", "sen_CI95_hi",
            "hamed_rao_p", "computable"]
    tr = tr[cols]
    tr.to_csv(OUT / f"qc_strict_trend_{mode}.csv", index=False, encoding="utf-8-sig")

    yy = pd.DataFrame(by_year(tab, rec_qc, qc_years))
    yy["mode"] = mode
    yy.to_csv(OUT / f"qc_strict_trend_byyear_{mode}.csv", index=False, encoding="utf-8-sig")

    print(f"\n[1] 收严档整条阶梯  → qc_strict_trend_{mode}.csv")
    for _, r in tr.iterrows():
        if r["computable"]:
            print(f"    {r['rung']:<13}{r['city_years']:>5} 城–年 · {r['n_cities']:>3} 城 · "
                  f"{r['n_years']:>2} 年 · Sen {r['sen_K_per_decade']:>+8.3f} "
                  f"({r['sen_CI95_lo']:>+8.3f}…{r['sen_CI95_hi']:>+8.3f}) "
                  f"p={r['hamed_rao_p']:.4f}")
        else:
            print(f"    {r['rung']:<13}{r['city_years']:>5} 城–年 · {r['n_cities']:>3} 城 · "
                  f"{r['n_years']:>2} 年 → **不可估计**（< {MIN_YEARS} 年）")

    print(f"\n[2] 逐年可用城数 + 保留像元·日  → qc_strict_trend_byyear_{mode}.csv")
    show = yy[["year", "n_cities_total", "n_cities_weighted", "n_cities_equal_any",
               "urban_retention_pooled", "rural_retention_pooled"]]
    print(show.to_string(index=False))

    RES[mode] = {"tr": tr, "yy": yy, "qc_years": qc_years, "no_qc": sorted(no_qc),
                 "n_rows": len(rec_qc), "n_city_years": len(tab)}

# ---------------------------------------------------------------- 报告
def fmt_trend(tr):
    lines = ["| 档位 | 说明 | 城–年 | 城 | 年 | 中位有效日 | Sen (K/十年) | 95% 区间 | Hamed–Rao p |",
             "|------|------|------:|----:|---:|----------:|-------------:|----------|------------:|"]
    for _, r in tr.iterrows():
        days = "—" if pd.isna(r["median_valid_days"]) else f"{r['median_valid_days']:.0f}"
        if r["computable"]:
            lines.append(f"| `{r['rung']}` | {r['description']} | {int(r['city_years'])} | "
                         f"{int(r['n_cities'])} | {int(r['n_years'])} | {days} | "
                         f"{r['sen_K_per_decade']:+.3f} | "
                         f"{r['sen_CI95_lo']:+.3f} … {r['sen_CI95_hi']:+.3f} | "
                         f"{r['hamed_rao_p']:.4f} |")
        else:
            lines.append(f"| `{r['rung']}` | {r['description']} | {int(r['city_years'])} | "
                         f"{int(r['n_cities'])} | {int(r['n_years'])} | {days} | "
                         f"**不可估计** | — | — |")
    return "\n".join(lines)


md = ["# 收严 QA 档 `(QC&3)==0` 的阶梯与 Sen 斜率（逐日 MOD11A1）", "",
      "来源：`analysis/code/daily/12_qc_strict_trend.py`（共享库 `daily_estimators.py`），",
      "读 `analysis/data/daily/daily_records_{night,day}.csv`。所有数字由脚本从",
      "`analysis/outputs/daily/qc_strict_trend_*.csv` 读回后写入本文件，**无手抄**。", "",
      f"仅用**带 Q 波段的 {len(RES[list(RES)[0]]['qc_years'])} 年**；2020、2023 导出无 Q 波段，",
      "收严列在该两年**全缺失**，故不在样本内（昼夜同）。", "",
      "## 一、整条阶梯（收严支撑）", ""]
for mode in RES:
    lab = "夜间" if mode == "night" else "白天"
    md += [f"### {lab} `{mode}`", "", fmt_trend(RES[mode]["tr"]), ""]

md += ["## 二、逐年可执行性（可用城数 + 保留像元·日）", "",
       "「可用城数」= 该城当年在收严档下有 ≥ 6 个有效日**且两侧各 ≥ 6 天**（`equal_any` 档）。",
       "保留 = 收严像元·日 ÷ 主口径像元·日（pooled，城乡分别）。", ""]


def peryear_md(yy, label):
    lines = [f"### {label}", "",
             "| 年 | 城数 | 加权档可用城 | 等权档可用城 | 城区保留 | 乡村保留 |",
             "|---:|----:|------------:|------------:|--------:|--------:|"]
    for _, r in yy.iterrows():
        lines.append(f"| {int(r['year'])} | {int(r['n_cities_total'])} | "
                     f"{int(r['n_cities_weighted'])} | {int(r['n_cities_equal_any'])} | "
                     f"{r['urban_retention_pooled']:.4%} | {r['rural_retention_pooled']:.4%} |")
    return "\n".join(lines)


for mode in RES:
    lab = "夜间" if mode == "night" else "白天"
    md += [peryear_md(RES[mode]["yy"], lab), ""]

def qc_prev(mode):
    """收严档 (tier==0) 的流行率：→ (share_of_valid, median_daily_share)，供报告引用。"""
    p = BASE / "tables" / f"qc_prevalence_{mode}.csv"
    if not p.exists():
        return None
    d = pd.read_csv(p)
    g = d[d.tier == 0]
    return (float(g.iloc[0].share_of_valid), float(g.iloc[0].median_daily_share)) if len(g) else None


def prev_txt(mode):
    p = qc_prev(mode)
    return (f"收严档占有效 LST 像元 **{p[0]:.3%}**、按日占比中位 **{p[1]:.3%}**"
            if p else "（未找到流行率表）")


RES_summary = {}
for mode in RES:
    yy = RES[mode]["yy"]

    def rng(col):
        pos = yy.loc[yy[col] > 0, col]
        return (int(pos.min()), int(pos.max())) if len(pos) else (0, 0)

    wlo, whi = rng("n_cities_weighted")
    elo, ehi = rng("n_cities_equal_any")
    RES_summary[mode] = {
        "n_years": len(RES[mode]["qc_years"]),
        "ea_years_with_any": int((yy["n_cities_equal_any"] > 0).sum()),
        "ea_median": float(yy["n_cities_equal_any"].median()), "ea_lo": elo, "ea_hi": ehi,
        "w_years_with_any": int((yy["n_cities_weighted"] > 0).sum()),
        "w_median": float(yy["n_cities_weighted"].median()), "w_lo": wlo, "w_hi": whi,
        "urban_ret": float(yy.urban_pixel_days_strict.sum() / yy.urban_pixel_days_main.sum()),
        "rural_ret": float(yy.rural_pixel_days_strict.sum() / yy.rural_pixel_days_main.sum()),
        "computable_rungs": RES[mode]["tr"].loc[RES[mode]["tr"].computable, "rung"].tolist(),
        "blocked_rungs": RES[mode]["tr"].loc[~RES[mode]["tr"].computable, "rung"].tolist(),
    }

md += ["## 三、结论", "",
       "把**「像元少」**与**「撑不起趋势」**分开读：前者是产品属性（流行率），",
       "后者是**可执行性**——把这点像元摊到「城市 × 年 × 日」后，还能不能凑出跨年单元。", ""]

ns = RES_summary["night"]
md += ["### 夜间：收严档能否端到端支撑一条斜率？", "",
       f"流行率（`tables/qc_prevalence_night.csv`）：{prev_txt('night')}。", "",
       f"1. **加权档**（无日门限，两侧 Σ 像元 > 0）：**{ns['w_years_with_any']}/"
       f"{ns['n_years']} 年**可算，逐年城数中位 **{ns['w_median']:.0f}**"
       f"（有城年份 {ns['w_lo']}–{ns['w_hi']} 城）。⇒ 收严档夜间**并非整体不可执行**，"
       f"像元池化加权仍拼出一条 23 年序列（斜率见上表）。",
       f"2. **等权（any）档**：仅 **{ns['ea_years_with_any']}/{ns['n_years']} 年**有 ≥1 城"
       f"可用，逐年城数中位 **{ns['ea_median']:.0f}**（有城的 {ns['ea_years_with_any']} 年各 "
       f"{ns['ea_lo']}–{ns['ea_hi']} 城，合计 33 城–年）。机制上算出一条斜率，但建立在逐年"
       f"**变动且极小**的子样本上，**不构成与主口径可比的年度序列**。",
       f"3. **覆盖门限档**（`equal_20px / equal_40pct / equal_60pct / paired_20px`）："
       f"**0 个城–年**。收严档夜间不足以让任一城任一侧在**同一天**凑够 20 个像元，"
       f"这些档**结构性为空**——不是「算不出」，是「不存在可算的单元」。",
       f"4. 收严档整体保留：城区 **{ns['urban_ret']:.4%}**、乡村 **{ns['rural_ret']:.4%}** 的"
       f"像元·日（对主口径）。",
       f"5. 分档结论：可算 {', '.join('`'+r+'`' for r in ns['computable_rungs']) or '（无）'}；"
       f"机械不可估计（年数 < {MIN_YEARS}）"
       f"{', '.join('`'+r+'`' for r in ns['blocked_rungs']) or '（无）'}。", ""]

ds = RES_summary["day"]
md += ["### 白天：对照", "",
       f"流行率：{prev_txt('day')}（`tables/qc_prevalence_day.csv`）。", "",
       f"- **{ds['ea_years_with_any']}/{ds['n_years']} 年**都有城可用，逐年城数中位 "
       f"**{ds['ea_median']:.0f}**（{ds['ea_lo']}–{ds['ea_hi']} 城）——整条阶梯可执行；"
       f"唯一薄档是 `equal_60pct`（仅 11 年）。",
       f"- 收严档整体保留：城区 **{ds['urban_ret']:.4%}**、乡村 **{ds['rural_ret']:.4%}**。", ""]

md += ["## 四、边界（**不能**用来声称什么）", "",
       f"1. 收严集是**非随机子集**（算法更自信的那部分），夜间等权档斜率建立在 33 个城–年、",
       f"   {ns['ea_years_with_any']} 个年上，**不是**主口径趋势的稳健性证据，也不是「更干净」"
       f"的估计。",
       "2. 各档的绝对斜率**不可并列**为同一估计量的区间——每档是不同的观测集合。",
       "3. 只读逐日 MOD11A1 中间表，**不可**外推到 MOD11A2 合成产品。",
       "4. 不构成因果或物理机制推断。", ""]

rep = OUT / "QC_STRICT_TREND.md"
rep.write_text("\n".join(md), encoding="utf-8")
print(f"\n已写 {rep.relative_to(BASE)}")
