"""逐日同日配对面板的**趋势检验**（2026-10-02）。回答最初那个斜率问题。

## 它回答什么

`54_method_comparability.py` 在**合成期**口径下给出四条序列，夜间：

| 序列 | Sen 斜率（K/十年） |
|------|------------------|
| `unpaired_pixel_weighted`（加权） | **−0.394** |
| `unpaired_equal_any`（等权） | **+0.204** |
| `unpaired_equal_20px`（等权+阈值） | +0.170 |
| `paired_equal_20px`（配对+等权） | +0.159 |

即「把加权换成等权，斜率翻号」。本脚本在**逐日**口径下重算同一阶梯。

> ⚠ **两个口径的斜率不可直接比大小**：数据源不同（MOD11A1 vs MOD11A2）、样本不同。
> 可比的是**同一口径内部**的序列之间（加权 vs 各等权档），也就是「翻号与否」这个
> 定性判断，而不是 −0.394 与某个逐日数字的差值。

## 为什么报**每一档**的翻号，而不只报 any 档

水平上已经发现加权效应**依赖等权规则**（见 `05_summarize_pilot.py` 的两产品阶梯表：
40% 档在两个产品里都翻成正号）。若斜率只在 `any` 档翻号、换成 40% 就不翻，
那「斜率翻号」这个全篇最强调的对照就是**规则选择的产物**而非数据的性质。
故本脚本对每一档都做一次翻号判定，并**并置两个产品**。

## ⚠ 实测：斜率翻号**本身也是分档的**（2026-10-02，合成期夜间，25 年）

给 54 补上 40%/60% 两档后立刻发现——**翻号只在松档成立**：

| 档 | 夜间斜率（K/十年） | 与加权同号？ |
|----|------------------|-------------|
| `weighted`（基准） | −0.394 | — |
| `equal_any` | **+0.204** | 翻号 |
| `equal_20px` | **+0.170** | 翻号 |
| `equal_40%` | **−1.145** | **不翻（比加权更负）** |
| `equal_60%` | **−2.225** | **不翻（比加权更负）** |
| `paired_20px` | +0.159 | 翻号 |

平衡面板（24 城）给出同一结论：加权 −0.310，any +0.181⇄、20px +0.195⇄、
40% −1.010、60% −2.088。手算 Sen 中位独立复核过 40% 那格（−1.145），非实现问题。

**这与水平上的发现是同一件事的两面**：40% 档的水平也翻成正号且极大
（2001 年 +6.10 K、2012 年 +9.26 K，对比加权 +3.47 / +5.67）。原因不是
「等权更保守」，而是**高覆盖要求换掉了样本**——它筛出的是城乡两侧
像元都高覆盖的晴空期，这类期次的城乡温差结构与全样本不同。

> **对正文的直接影响**：那句「把加权换成等权，斜率就翻号」**不能泛称**。
> 它只在 `≥1 像元` 与 `≥20 像元` 两档成立。若正文按 18 的结论泛化陈述，
> 审稿人只要问一句「换 40% 阈值呢」就能推翻。须改写为
> 「在 ≥1 与 ≥20 像元两种等权规则下斜率翻号；更严的覆盖要求把样本换成
> 另一类期次，斜率同号且更负」。本脚本产出的 `daily_trends_flip*.csv`
> 就是这句话的出处（单模式带模式后缀，见文末用法）。

## 平衡面板（本脚本与 54 的一处**有意不同**）

掩膜城市数逐年增长（2001 年 25 城 → 2023 年 37 城）。城集逐年变化会**自身**制造趋势。
54 用不平衡面板；本脚本**同时**报平衡面板（全部年份都有的城市交集）作为对照。
两者若给出相反的翻号结论，说明构成效应主导，正文须写明用哪一个、为什么。

用法：
    python 06_daily_trends.py                 # 昼夜都跑 → `daily_trends.csv`（合并表）
    python 06_daily_trends.py --mode night    # → `daily_trends_night.csv`（**不覆盖合并表**）
    python 06_daily_trends.py --mode day      # → `daily_trends_day.csv`

输出命名（2026-10-03 加）：单模式运行写带模式后缀的文件，只有不带 `--mode`
才写无后缀的合并表。加这条是因为实测先跑 night 再跑 day 时，后者把夜间那批
斜率整表覆盖，使已写入状态记录的夜间数字失去产出它的那一行。
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import theilslopes

sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")   # 只改 stdout 时 stderr 仍走 GBK，中文报错会乱码

BASE = Path(__file__).resolve().parent.parent.parent      # analysis/
sys.path.insert(0, str(BASE / ".deps"))
import pymannkendall as mk                                 # 与 18 / 54 同一实现

DATA = BASE / "data"
OUT = BASE / "outputs" / "daily"
OUT.mkdir(parents=True, exist_ok=True)

parser = argparse.ArgumentParser(description="逐日同日配对趋势检验")
parser.add_argument("--mode", choices=("night", "day"), default=None)
parser.add_argument("--min-years", type=int, default=8)
ARGS = parser.parse_args()

MODES = (ARGS.mode,) if ARGS.mode else ("night", "day")
MIN_YEARS = ARGS.min_years
# 输出后缀（2026-10-03 加）。**单模式运行必须带模式名**：否则先跑 `--mode night`
# 再跑 `--mode day`，后者把前者整表覆盖，夜间那批斜率就失去「产出它的那一行」。
# 2026-10-03 实跑正是如此——夜间的 14 行被日间覆盖后才补的这条。
# 不带 `--mode` 时写无后缀的合并表（含两个模式），那才是无后缀名的含义。
MOSUF = f"_{ARGS.mode}" if ARGS.mode else ""


def daily_map(mode):
    """逐日面板的 (档位, 列名)。`weighted` 是基准，其余与之比大小。"""
    return [
        ("weighted", "suhii_pixel_weighted"),
        ("equal_any", "suhii_equal_any"),
        ("equal_20px", f"suhii_{mode}_daily_equal_20px"),
        ("equal_40pct", f"suhii_{mode}_daily_equal_40%"),
        ("equal_60pct", f"suhii_{mode}_daily_equal_60%"),
        ("paired_20px", f"suhii_{mode}_daily_paired_20px"),
        ("weighting_effect", "weighting_effect"),
    ]


def comp_map(mode):
    """合成期面板的同一批档位。命名与逐日**不同**，故各给一张表，不做字符串拼接。"""
    return [
        ("weighted", f"suhii_{mode}_unpaired"),
        ("equal_any", "unpaired_equal_any_K"),
        ("equal_20px", "unpaired_equal_period_K"),
        ("equal_40pct", "unpaired_equal_40%_K"),
        ("equal_60pct", "unpaired_equal_60%_K"),
        ("paired_20px", f"suhii_{mode}_paired"),
        ("weighting_effect", "weighting_effect_K"),
    ]


DESC = {"weighted": "像元全池化（加权）", "equal_any": "等权（任一侧有像元）",
        "equal_20px": "等权 + 20 像元下限", "equal_40pct": "等权 + 40% 下限",
        "equal_60pct": "等权 + 60% 下限", "paired_20px": "同日配对 + 等权",
        "weighting_effect": "等权(any) − 加权"}


def trend_of(annual: pd.Series):
    """Theil–Sen + Hamed–Rao，与 54 逐字同参。年数不足返回 None。"""
    a = annual.dropna().sort_index()
    if len(a) < MIN_YEARS:
        return None
    sen = theilslopes(a.values, a.index.values, alpha=0.95)
    hr = mk.hamed_rao_modification_test(a.values)
    return {"years": int(len(a)),
            "grand_mean_K": float(a.mean()),
            "sen_slope_K_per_decade": float(sen.slope * 10),
            "sen_CI95_lo": float(sen.low_slope * 10),
            "sen_CI95_hi": float(sen.high_slope * 10),
            "hamed_rao_p": float(hr.p)}


rows, flip_rows = [], []
for mode in MODES:
    src = DATA / "daily" / f"panel_{mode}_daily.csv"
    if not src.exists():
        print(f"缺 {src.relative_to(BASE)}，跳过 {mode}")
        continue
    df = pd.read_csv(src)
    label_mode = "夜间" if mode == "night" else "白天"
    print("=" * 84)
    print(f"{label_mode}逐日趋势   [mode={mode}]   "
          f"{len(df)} 城市–年 · 年 {int(df.year.min())}–{int(df.year.max())} · "
          f"各年 {df.groupby('year').size().min()}–{df.groupby('year').size().max()} 城")
    print("=" * 84)

    per_year = df.groupby("year").city_id.apply(set)
    common = set.intersection(*per_year.tolist()) if len(per_year) else set()
    print(f"  平衡面板交集城市 {len(common)} 个"
          + ("  ⚠ 过少，平衡结果仅供参考" if len(common) < 10 else ""))

    for panel_name, sub in (("unbalanced", df), ("balanced", df[df.city_id.isin(common)])):
        seen = {}
        for tag, col in daily_map(mode):
            if col not in sub.columns:
                print(f"  ! 缺列 {col}", file=sys.stderr)
                continue
            t = trend_of(sub.groupby("year")[col].mean())
            if t is None:
                continue
            seen[tag] = t
            rows.append({"mode": mode, "panel": panel_name, "product": "daily",
                         "rung": tag, "description": DESC[tag], "column": col,
                         "city_years": int(sub[col].notna().sum()), **t})
        for tag, t in seen.items():
            print(f"  [{panel_name:<10}] {tag:<17}{t['sen_slope_K_per_decade']:>+8.3f} K/十年"
                  f"  (95% {t['sen_CI95_lo']:>+7.3f} … {t['sen_CI95_hi']:>+7.3f})"
                  f"  p={t['hamed_rao_p']:.3f}  n={t['years']}")
    print()

if rows:
    out = pd.DataFrame(rows)
    dst = OUT / f"daily_trends{MOSUF}.csv"
    out.to_csv(dst, index=False, encoding="utf-8-sig")
    print(f"已写 {dst.relative_to(BASE)}  {len(out)} 行")
else:
    print("（逐日面板年数不足，本次无逐年序列；翻号判定仍会对合成期执行）", file=sys.stderr)


# ------------------------------------------------- 翻号判定（两个产品并置）
def balanced_set(frame):
    """该框架**自身**全部年份都出现的城市。逐日与合成的城集不同，故各算各的——
    用逐日的交集去裁合成期会把合成期本来完整的年份无谓地削掉。"""
    per_year = frame.groupby("year").city_id.apply(set)
    return set.intersection(*per_year.tolist()) if len(per_year) else set()


def flip_table(mapping, frame):
    """对每一档判「相对 weighted 是否翻号」。→ ([(档, 斜率, 是否翻)], {档: 斜率})。"""
    got = {}
    for tag, col in mapping:
        if col in frame.columns:
            t = trend_of(frame.groupby("year")[col].mean())
            if t:
                got[tag] = t["sen_slope_K_per_decade"]
    if "weighted" not in got:
        return None, got
    w = got["weighted"]
    # `weighting_effect` 是「等权(any) − 加权」的**差值序列**，本身不是一条等权规则，
    # 与加权基准比符号没有「换口径是否翻号」的含义（实测恒为 True，见 daily_trends_flip.csv）。
    # 它与 `weighted` 一样只是参照，不参与判定。
    return [(tag, s, bool(np.sign(s) != np.sign(w)))
            for tag, s in got.items() if tag not in ("weighted", "weighting_effect")], got


print("\n" + "=" * 84)
print("【翻号判定：各等权档 vs 池化加权（斜率，K/十年）】")
print("=" * 84)
for mode in MODES:
    dsrc = DATA / "daily" / f"panel_{mode}_daily.csv"
    csrc = BASE / "outputs" / "method_comparability" / f"equal_period_{mode}.csv"
    sources = [("逐日 MOD11A1", daily_map(mode), dsrc),
               ("合成 MOD11A2", comp_map(mode), csrc)]
    if not any(p.exists() for _, _, p in sources):
        continue
    print(f"\n  {mode}：")
    for product, mapping, path in sources:
        if not path.exists():
            print(f"    {product}: 缺 {path.name}，跳过")
            continue
        frame = pd.read_csv(path)
        common = balanced_set(frame)
        for panel_name in ("unbalanced", "balanced"):
            fx = frame[frame.city_id.isin(common)] if panel_name == "balanced" else frame
            flips, got = flip_table(mapping, fx)
            if flips is None:
                print(f"    [{panel_name:<10}] {product}: 年数不足或缺加权序列，跳过")
                continue
            flip_rows.append({"mode": mode, "panel": panel_name, "product": product,
                              "weighted_slope": got["weighted"],
                              "n_balanced_cities": len(common) if panel_name == "balanced" else np.nan,
                              **{f"{t}_slope": s for t, s, _ in flips},
                              **{f"{t}_flipped": f for t, _, f in flips}})
            marks = "  ".join(f"{t}={s:+.3f}{'⇄' if f else ' '}" for t, s, f in flips)
            bal = f" · 平衡城 {len(common)}" if panel_name == "balanced" else ""
            print(f"    [{panel_name:<10}] {product}: 加权={got['weighted']:+.3f}{bal} | {marks}")
    print("    （⇄ = 与加权反号。逐日与合成**不可直接比大小**，只比各自的『翻不翻』。）")

if flip_rows:
    fp = OUT / f"daily_trends_flip{MOSUF}.csv"
    pd.DataFrame(flip_rows).to_csv(fp, index=False, encoding="utf-8-sig")
    print(f"\n已写 {fp.relative_to(BASE)}  {len(flip_rows)} 行")

if not rows and not flip_rows:
    print("没有可检验的序列（两个产品的面板都缺失或年数不足）", file=sys.stderr)
    raise SystemExit(1)

print("""
  判读：
  1. 全篇最强的那句是「加权→等权，斜率翻号」。本表逐档检验它**是否只在某一档成立**。
     若只有 `equal_any` 翻号、`equal_20px` 起就不翻，则该对照是**规则选择的产物**，
     正文必须改写为「在 ≥1 像元等权这一档下翻号」，不能泛称「换等权就翻号」。
  2. 平衡与不平衡面板若给出相反的翻号结论，说明城集构成效应主导，**不可只报一个**。
  3. `equal_60pct` 在**逐日夜间基本不可用**：`n_paired_days_60%` 中位仅 **6** 天，
     低于 MIN_DAYS=6 的城–年有 395/803。即逐日夜间的 60% 档只在 572 个城–年上存在，
     其斜率建立在**逐年变化的子样本**上，与其它档不同质，**不可与 20px/40% 并列报告**。
     （对照：逐日夜间 40% 档中位 21 天、20px 档中位 113 天，两者都可用；
     白天 60% 档中位 83 天，完全可用。这个坍缩是**夜间独有**的。）
     ⚠ 本条的**原数字是错的**：原先写「中位 4 天、51/73、19 天、118 天、86 天」，
     全部来自 2020+2023 两年 pilot（`panel_night_before.csv`，73 行）。
     2026-10-03 独立复核时发现**正文已沿用这些 pilot 数字**，遂按 803 行全量重算并
     同步改正正文。改动脚本本段时必须回 `data/daily/panel_{mode}_daily.csv` 重算，
     不要凭本注释转述。
  4. 早期年份缺日（2001 缺 17 天等）可能非随机。若不平衡面板的斜率明显被早期年份拉动，
     须在正文写明并给出剔除早期年份的敏感性。
""")
