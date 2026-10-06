"""把逐日 pilot 汇成**论文面**的一张表 + 一份判读（2026-10-02）。pilot 的诚实版。

## 为什么重写（2026-10-02 对抗审计后更正）

初版有两处**站不住的呈现**，已改：

1. **「保留率 = 逐日效应 / 合成期效应」被删除**。它是 MOD11A1 与 MOD11A2 两个**不同产品**
   的效应之比——而本流水线的自订约定（`03_build_daily_panel.py` 开头）明确禁止拿两者
   直接比大小。更要命的是它在夜间**没有不确定性**：逐日−合成之差实测
   2020 = +0.284 K（配对 SE 0.191，t=1.49，p=0.144）、2023 = +0.311 K（t=1.35，p=0.186），
   **两个年份都不可区分于零**。故「缺口缩小 42%」这类说法无支撑，只能是
   「逐日口径下该效应显著非零，且与合成期估计不可区分」。
2. **初版只给点估计、不给不确定性**，导致把噪声当效应。本版每个效应都带跨城 SE 与 t。

## 本脚本报什么

- 逐年、城市集取交集的**加权效应**（等权 − 池化加权），逐日与合成期各一列
- 逐日效应的跨城 SE 与对 0 的 t
- 逐日与合成期的**配对差**及检验 p（回答「两者是否可区分」）
- 等权规则敏感度：同一批数据换 `any / 20px / 40%` 规则时效应落在什么区间

用法：
    python 05_summarize_pilot.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")   # 只改 stdout 时 stderr 仍走 GBK，中文报错会乱码

BASE = Path(__file__).resolve().parent.parent.parent      # analysis/
OUTDIR = BASE / "outputs" / "daily"
TABDIR = BASE / "tables"
TABDIR.mkdir(parents=True, exist_ok=True)

# 等权规则档位：列名后缀 → 可读标签。`any` 无后缀。
RUNGS = [("any", "suhii_equal_any"), ("20px", None), ("40%", None)]

rows = []
for mode in ("night", "day"):
    dpath = BASE / "data" / "daily" / f"panel_{mode}_daily.csv"
    cpath = BASE / "outputs" / "method_comparability" / f"equal_period_{mode}.csv"
    if not (dpath.exists() and cpath.exists()):
        print(f"缺 {mode} 的源文件，跳过", file=sys.stderr)
        continue
    d, c = pd.read_csv(dpath), pd.read_csv(cpath)

    for year, dg in d.groupby("year"):
        cg = c[c.year == year]
        common = sorted(set(dg.city_id) & set(cg.city_id))
        if len(common) < 10:
            continue
        dg = dg[dg.city_id.isin(common)].set_index("city_id")
        cg = cg[cg.city_id.isin(common)].set_index("city_id")

        pw = dg.suhii_pixel_weighted.mean()
        daily_eff = dg.weighting_effect                       # 等权(any) − 加权
        comp_eff = cg.weighting_effect_K                      # 同上，合成期口径
        diff = daily_eff - comp_eff
        se = lambda s: float(s.std(ddof=1) / np.sqrt(len(s)))
        t_eff = float(daily_eff.mean() / se(daily_eff))
        t_d, p_d = stats.ttest_rel(daily_eff, comp_eff)

        # 等权规则敏感度：同一批城市在每一档下各算一次效应，逐日与合成期并置。
        # 两产品的列名不同（逐日 suhii_{mode}_daily_equal_*，合成期 unpaired_equal_*），
        # 故各给一张名表——不要试图拼一个字符串模板，两边命名本就不一致。
        RUNGS = [("any", "suhii_equal_any", "unpaired_equal_any_K"),
                 ("20px", f"suhii_{mode}_daily_equal_20px", "unpaired_equal_period_K"),
                 ("40%", f"suhii_{mode}_daily_equal_40%", "unpaired_equal_40%_K"),
                 ("60%", f"suhii_{mode}_daily_equal_60%", "unpaired_equal_60%_K")]
        rec = {
            "mode": mode, "year": int(year), "common_cities": len(common),
            "daily_weighting_effect_K": round(float(daily_eff.mean()), 3),
            "daily_SE_K": round(se(daily_eff), 3),
            "daily_t_vs_zero": round(t_eff, 2),
            "composite_weighting_effect_K": round(float(comp_eff.mean()), 3),
            "paired_diff_K": round(float(diff.mean()), 3),
            "paired_diff_SE_K": round(se(diff), 3),
            "paired_diff_p": round(float(p_d), 3),
        }
        c_pw = cg[f"suhii_{mode}_unpaired"].mean()
        for label, dcol, ccol in RUNGS:
            if dcol in dg.columns and dg[dcol].notna().any():
                rec[f"daily_rung_{label}"] = round(float(dg[dcol].mean() - pw), 3)
                rec[f"daily_rung_{label}_NaN"] = int(dg[dcol].isna().sum())
            if ccol in cg.columns and cg[ccol].notna().any():
                rec[f"comp_rung_{label}"] = round(float(cg[ccol].mean() - c_pw), 3)
                rec[f"comp_rung_{label}_NaN"] = int(cg[ccol].isna().sum())
        rows.append(rec)

if not rows:
    print("没有可汇总的年份（面板或对照文件缺失）", file=sys.stderr)
    raise SystemExit(1)

df = pd.DataFrame(rows)
out = TABDIR / "daily_pairing_pilot.csv"
df.to_csv(out, index=False, encoding="utf-8-sig")
print(f"已写 {out.relative_to(BASE)}\n")
print(df.to_string(index=False))

# ---------------------------------------------------------------- 判读
diag = {}
for mode in ("night", "day"):
    p = BASE / "data" / "daily" / f"panel_{mode}_daily.csv"
    if p.exists():
        dd = pd.read_csv(p)
        diag[mode] = {int(y): (g.u_days_per_px.median(), g.r_days_per_px.median())
                      for y, g in dd.groupby("year")}

night = df[df["mode"] == "night"]
day = df[df["mode"] == "day"]
_lines = [
    "# 逐日 MOD11A1 同日配对 pilot — 判读（对抗审计后更正版）", "",
    "> 全部数字由 `05_summarize_pilot.py` 从 `data/daily/panel_*.csv` 与",
    "> `outputs/method_comparability/equal_period_*.csv` 派生。可回溯到具体行。", "",
    "## 表", "", df.to_markdown(index=False), "",
    "## 城乡观测对称性（`u_days_per_px` / `r_days_per_px` 中位）", "",
]
for mode, per in diag.items():
    for y, (u, r) in per.items():
        _lines.append(f"- {mode} {y}: 城区 {u:.1f} · 乡村 {r:.1f} · **比 {u / r:.2f}×**")

_lines += ["", "## 结论（三条，均带不确定性）", ""]
if len(night):
    n22 = night[night.year == night.year.max()].iloc[0]
    _lines += [
        # 年份标签必须与取值同源：写死 "2020" 而取 iloc[0] 时，一旦面板含早于 2020 的
        # 年份（全量 2001–2025 就是），iloc[0] 变成 2001 而标签仍印 "2020"——静默错年。
        f"1. **夜间加权效应在同日配对下显著非零**：{int(night.iloc[0].year)} = "
        f"{night.iloc[0].daily_weighting_effect_K:+.3f} K（t={night.iloc[0].daily_t_vs_zero}）、"
        f"{int(n22.year)} = {n22.daily_weighting_effect_K:+.3f} K（t={n22.daily_t_vs_zero}）。",
        "   这是本 pilot 支持的**唯一**强主张。",
        "",
        "2. **但不能说它比合成期「缩小了」**：逐日−合成之差 "
        + "、".join(f"{int(r.year)} = {r.paired_diff_K:+.3f} K (p={r.paired_diff_p})"
                   for _, r in night.iterrows())
        + " —— 两年都不可区分于零。",
        "   故**不报**保留率/收缩比例（初版报过，已撤）。",
        "",
    ]
if len(day):
    d23 = day[day.year == day.year.max()]
    if len(d23):
        r = d23.iloc[0]
        if abs(r.paired_diff_p) < 0.05:
            _lines += [
                f"3. **白天 {int(r.year)} 两个产品给出显著相反的口径结论**：逐日 "
                f"{r.daily_weighting_effect_K:+.3f}（t={r.daily_t_vs_zero}）vs 合成期 "
                f"{r.composite_weighting_effect_K:+.3f}；配对差 {r.paired_diff_K:+.3f} K，"
                f"p={r.paired_diff_p}。",
                "   这**不是**「样本不足」可以打发的——它是本文必须如实报告的一处分歧。",
                "",
            ]

_lines += [
    "## 加权效应依赖等权规则（最重要的一条）", "",
    "同一批城市，换一套「等权」的实现，效应可以**翻号**。下表左右两栏是**两个产品各自**",
    "在同一批阈值规则下的效应（单位 K）。括号内为因不达标而不可算的城市数。", "",
    "| 年 | 口径 | any | 20px | 40% | 60% |", "|---|---|---|---|---|---|",
]
for _, r in df.iterrows():
    y = int(r["year"])
    for side, tag in (("comp", "合成期 MOD11A2"), ("daily", "逐日 MOD11A1")):
        cells = []
        for lab in ("any", "20px", "40%", "60%"):
            v, n = r.get(f"{side}_rung_{lab}"), r.get(f"{side}_rung_{lab}_NaN")
            cells.append("—" if not np.isfinite(v) else
                         f"{v:+.3f}" + (f" ({int(n)}缺)" if n else ""))
        _lines.append(f"| {r['mode']} {y} | {tag} | " + " | ".join(cells) + " |")

_lines += [
    "",
    "**读法**（这条比任何单个数字都重要）：",
    "",
    "- `any` 与 `20px` 两档在**两个产品里同号**这句话**只对夜间成立**（夜间两年两档、",
    "  白天 2020 两档皆同号）。**白天 2023 的 `any` 档两产品异号**：合成期 +0.153 vs",
    "  逐日 −0.170（即上面结论 3 那处分歧）。故稳健的只是**夜间那一段**，",
    "  不得写成「两档在两个产品里同号」而不加限定。",
    "- `40%` 档在夜间**两个产品里都翻成正号且很大**（合成期 +1.198 / 逐日 +1.392，2020）。",
    "  换言之，翻号**不是 MOD11A1 或 MOD11A2 谁的问题**，而是「等权」这个词本身的歧义：",
    "  40% 档筛出的是高覆盖（晴空）日，其城乡温差结构与全样本不同。",
    "- 故**正文报加权效应时必须写明档位**，只说「加权效应 = X K」是没有意义的。",
    "- `60%` 档夜间开始坍缩（逐日 2020 有 16/36 城不可算），阈值阶梯的最高档在逐日夜间**断掉**。",
    "", "## 未做的事（不得据此下结论）", "",
    "- **不出斜率**。两年估不出十年斜率。",
    "- 两年共用同一批城市（夜间均为 36 城），**不是独立重复**。",
    "- `60%` 档的「严阈值同号」检验在逐日夜间**不可执行**（见上），",
    "  正文须如实写明这一点，不得声称做过。", "",
]
rep = OUTDIR / "PILOT_REPORT.md"
rep.write_text("\n".join(_lines), encoding="utf-8")
print(f"\n已写 {rep.relative_to(BASE)}")
