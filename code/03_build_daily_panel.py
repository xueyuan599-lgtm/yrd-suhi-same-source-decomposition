"""逐日 MOD11A1 的**同日配对** SUHII 面板（2026-10-02）。pilot 核心问题的回答器。

## 它回答什么

合成期口径下，`54_method_comparability.py` 的四条序列是（798 个夜间城市–年）：

| 序列 | 含义 | 水平 | Sen 斜率 |
|------|------|------|---------|
| `unpaired_pixel_weighted` | 像元–期次全池化（"加权"） | 4.263 K | **−0.394** |
| `unpaired_equal_any` | 每**合成期**等权 | 2.115 K | **+0.204** |
| `unpaired_equal_20px` | 等权 + 20 像元下限 | 2.384 K | +0.170 |
| `paired_equal_20px` | 同合成期配对 + 等权 | 1.980 K | +0.159 |

**加权效应**（**第二条 − 第一条** = 等权 − 加权，与 54 的 `weighting_effect_K` 同向）
均值 **−2.148 K**，逐年在 −0.5 ~ −3.6 K 之间。把它从加权改成等权，斜率就翻号——
这是全篇最强的描述性对照。

> ⚠ 2026-10-02 更正：此行原写「第一条 − 第二条」，与代码和 54 的定义**反号**。
> 代码是对的（54 与 03 都是 `等权 − 加权`），错的是这行说明。

## ⚠ 加权效应不是单个数——它依赖等权规则（2026-10-02 实测）

「等权」有多种实现。用**同一批数据**换一套阈值规则，效应可以翻号：

| 等权规则 | 合成期夜间 2020 | 逐日夜间 2020 |
|---------|---------------|--------------|
| `any`（任一侧 ≥1 像元） | −0.683 | −0.399 |
| `20px`（≥20 像元） | −0.266 | −0.233 |
| `40%`（≥40% 掩膜像元） | **+1.198** | **+1.392** |
| `60%` | +2.327（1/36 城缺） | +0.083（16/36 城缺） |

**两个产品都翻号，量级也相当**。故这不是 MOD11A2 合成结构的产物，而是
「用什么规则定义等权」的产物——40% 档筛出的是高覆盖（晴空）日，其城乡温差结构
与全样本不同。**正文报告加权效应时必须写明用的是哪一档**，不可只说「加权效应 = X」。

> ⚠ 2026-10-02 二次更正：本表初版写 40% = **+1.008**、60% = −0.250（32/36 缺），
> 两个数都**回溯不到任何产出**。当时合成期侧根本没有 40%/60% 档——54 硬编码了
> `MIN_PX = 20`。已给 54 补上 `unpaired_equal_{20px,40%,60%}_K` 三列（纯追加，
> 旧列未动），上表全部取自 `outputs/method_comparability/equal_period_night.csv`。

本脚本在**逐日**数据上重算同一阶梯。逐日没有合成期，"期"换成"天"：
若逐日下这个 ~2 K 的加权缺口**仍然存在**，说明它不是 MOD11A2 合成结构的产物，
而是「按有效像元多寡加权」这件事本身的产物——等权口径因此更可辩护。
若缺口消失，则合成结构才是主因，正文口径须改写。

## 实测结果（2026-10-02，2020+2023，城市集取交集）

| 昼夜 | 年 | 城 | 合成期缺口 | 逐日缺口 | 逐日 t | 逐日−合成 | 配对 p |
|------|----|----|-----------|---------|--------|-----------|--------|
| 夜间 | 2020 | 36 | −0.683 | **−0.399** | −2.38 | +0.284 | 0.144 |
| 夜间 | 2023 | 36 | −1.324 | **−1.013** | −4.46 | +0.311 | 0.186 |
| 白天 | 2020 | 36 | −0.144 | −0.196 | −2.57 | −0.053 | 0.454 |
| 白天 | 2023 | 37 | +0.153 | −0.170 | −2.41 | −0.323 | <0.001 |

**能说的**：夜间加权缺口在**同日配对下显著非零**（t = −2.38 / −4.46），
故它不是「合成期结构造出来的假象」——这是本 pilot 支持的**唯一**强主张。

**不能说的**：不能说它比合成期「缩小了」。逐日−合成之差两年都是 +0.28 / +0.31 K，
**不可区分于零**（p = 0.144 / 0.186）。

> ⚠ 2026-10-02 更正：初版此处写「保留 0.58 / 0.77」并据此称「同号同量级」。
> 那个比值是 MOD11A1 与 MOD11A2 **两个不同产品**的效应之比，本文件第 58 行的自订约定
> 明文禁止这种直接比大小；且它忽略了不确定性，把噪声当成了收缩。**该说法已撤**。

白天两侧量级仅 ~0.17 K，而 2023 年**两产品显著相反**（配对差 −0.323 K，p<0.001）——
这不是「两年不足以判定」可以打发的，是正文须如实报告的一处分歧。

## 城乡观测不对称：同日配对把它压小了，但没压没（2026-10-02 实测）

⚠ **必须在同一批城–年上比**。全年的合成期中位比是 2.71×，但那是 2001–2025
全部城–年的中位；拿它去比逐日（只有 2020/2023）的 1.77× 就是本文件第 58 行
禁止的跨产品比大小。**只在 2020/2023 交集上比**才是同口径：

| 口径（2020/2023 交集的 72 个城–年） | 城区 | 乡村 | 比 |
|------|------|------|----|
| 合成期 MOD11A2（`panel_night_paired.csv`，20px 档） | 20.8 | 9.3 | **2.24×** |
| 逐日 MOD11A1（本脚本产出） | 64.9 | 35.7 | **1.82×** |

白天两侧本来就对称（合成期 0.97×、逐日 0.96×），无可压。

**同日配对把夜间的超额不对称从 1.24 压到 0.82，即压掉约 34%**——这正是它该起的作用。
但加权缺口**并没有随之消失**（仍 −0.399 / −1.013 K，显著非零）。
故：合成结构**是**不对称的来源之一（约三分之一），但不是全部；
剩余部分来自 MOD11A1 本身在城乡两侧的晴空/云掩膜差异。
**单靠同日配对不足以消除加权效应**——这是本 pilot 对正文最有用的一句。

> ⚠ 2026-10-02 更正：本节初版写「从 2.71× 压到 1.77×，解释了约 1/3」，
> 前半句是**跨样本**比较（全年的 2.71× vs 两年的 1.77×），不成立；
> 同口径下的正确数字是 2.24× → 1.82×。结论方向不变，量值已改。

（论文面表格与判读见 `05_summarize_pilot.py` 产出的 `tables/daily_pairing_pilot.csv`
与 `outputs/daily/PILOT_REPORT.md`。）

## 两条不可违反的约定

1. **城乡掩膜只读 `masks_matched.npz`**（12 的唯一产物），本脚本**不重新判定城乡**。
2. **四个口径从同一份逐日序列算出**（苹果对苹果）。不得拿 18 的合成期数字
   与逐日数字直接比大小——两者是不同数据源，只在**同一口径内**做纵向比较。

## 与合成的对照怎么做才算数

对照只在**逐年、且城市集取交集**后做。2020 掩膜 36 城、2023 掩膜 37 城，
而合成期面板 2023 年只有 36 城——不取交集会把构成差异当成效应。

用法：
    python 03_build_daily_panel.py                 # 夜间（默认）
    python 03_build_daily_panel.py --mode day
"""

import argparse
import re
import sys
from itertools import groupby
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

sys.stdout.reconfigure(encoding="utf-8")

parser = argparse.ArgumentParser(description="逐日同日配对 SUHII 面板")
parser.add_argument("--mode", choices=("night", "day"), default="night")
# 不给 --years 就跑**盘上已有的全部年份**。铺满 25 年时不必敲 25 个参数，
# 也不会因为漏写某年而静默少算（少算的年份在趋势里表现为「缺年」，
# 容易被误读成「那年没数据」）。要跑子集时显式给。
parser.add_argument("--years", type=int, nargs="+", default=None,
                    help="默认=盘上已有的全部年份")
ARGS = parser.parse_args()
MODE = ARGS.mode
TAG = MODE
LABEL = "夜间" if MODE == "night" else "白天"
COL_PAIRED = f"suhii_{MODE}_daily_paired"

FILL = -9999
SCALE_K = 0.02
# 阈值阶梯与 18 逐字同构：绝对像元下限 + 按各自像元数的比例下限
SCHEMES = [("20px", 20, 0.0), ("40%", 10, 0.40), ("60%", 10, 0.60)]
MIN_DAYS = 6

BASE = Path(__file__).resolve().parent.parent.parent      # analysis/
DATA = BASE / "data"
EXPORT = DATA / "gee_export"
DAILY = DATA / "daily"
OUT = DAILY / f"panel_{TAG}_daily.csv"
DAILY.mkdir(parents=True, exist_ok=True)

if ARGS.years is None:
    # 用正则而非 split("_")：`uhi_nightdaily_2001_h1.tif` 用 split 也能取到，
    # 但只要文件名多一段（如将来的后缀变体）下标就漂了。
    YEARS = sorted({int(m.group(1)) for p in EXPORT.glob(f"uhi_{TAG}daily_*.tif")
                    if (m := re.match(rf"uhi_{TAG}daily_(\d{{4}})_", p.name))})
    if not YEARS:
        raise SystemExit(f"盘上没有任何 uhi_{TAG}daily_*.tif —— 先跑 02 导出并 15 下载")
    print(f"  --years 未指定 → 取盘上已有的 {len(YEARS)} 个年份："
          f"{YEARS[0]}–{YEARS[-1]}")
else:
    YEARS = ARGS.years

masks = np.load(DATA / "masks_matched.npz")
H, W = int(masks["H"]), int(masks["W"])

print("=" * 78)
print(f"{LABEL}逐日同日配对 SUHII   [mode={MODE} · {len(YEARS)} 年]")
print("=" * 78)
print(f"  掩膜 {DATA.name}/masks_matched.npz（12 的唯一产物，本脚本不重判城乡）")


def load_year(year):
    """读该年的**LST**波段清单 → ([(文件, 波段号)], [日期])，均按日期序。

    只驻留一个波段（885×855×4B），366 波段全读进内存会到 GB 级。

    ⚠ 2026-10-02 起导出文件**同时含 LST 与原始 QC 两条波段**
    （band 名 `D{日期}` 与 `Q{日期}` 交替）。本脚本是**主口径**，只取 `D*`。
    若不过滤，QC 波段会被当成 LST 逐日值读进来——不报错，结果全错。
    收严档的敏感性由独立脚本读 `Q*` 波段，不走这里。
    """
    files = sorted(EXPORT.glob(f"uhi_{TAG}daily_{year}_*.tif"))
    if not files:
        return None, None
    entries, dates = [], []
    for f in files:
        with rasterio.open(f) as s:
            for bi, d in enumerate(s.descriptions, start=1):
                if d and d.startswith("D"):
                    entries.append((f, bi))
                    dates.append(d[1:])
    if len(set(dates)) != len(dates):
        dup = sorted({d for d in dates if dates.count(d) > 1})
        raise RuntimeError(f"{year}: 波段日期有重复 {dup} —— 分半段拼接错位/跨半段重叠")
    order = sorted(range(len(dates)), key=lambda i: dates[i])
    return [entries[i] for i in order], [dates[i] for i in order]


rows = []
print(f"\n{'年':<7}{'日数':<7}{'城':<6}{'配对日中位':<12}{'逐日加权':<11}{'逐日等权':<11}{'配对等权':<11}")
print("-" * 68)

for year in YEARS:
    entries, dates = load_year(year)
    if not dates:
        print(f"{year:<7}缺文件（先跑 02_export_daily_seq.py 并下载）")
        continue

    cities = sorted(int(k.split("|")[1]) for k in masks.files
                    if k.startswith(f"u{year}|") and f"r{year}|{k.split('|')[1]}" in masks.files)
    nc = len(cities)
    acc = {c: {"per_day": [], "u_sum": 0.0, "u_n": 0, "r_sum": 0.0, "r_n": 0} for c in cities}
    u_px = {c: len(masks[f"u{year}|{c}"]) for c in cities}
    r_px = {c: len(masks[f"r{year}|{c}"]) for c in cities}

    # 「像元 → 城市序号」查找表（0 = 不属于任何城）。用它 + `bincount` 取代逐城
    # fancy-indexing：原写法每波段做 2×nc 次索引（nc=36 → 72 次），每次都把 757k 的
    # float 数组拷一份，25 年 × 366 天累积到小时级；bincount 一次遍历同时算出所有城
    # 的和与计数。**口径完全不变**：仍是非填充（`(Q&3)<=1` 已由 02 掩膜）像元的均值，
    # 掩膜仍取 `masks_matched.npz`（12 的唯一产物）。2026-10-03 以 2020/2023 旧面板
    # 逐值回归验证（`--years 2020 2023` 对比 `panel_*_daily.csv`）。
    ci_u = np.zeros(H * W, dtype=np.int32)
    ci_r = np.zeros(H * W, dtype=np.int32)
    for j, c in enumerate(cities, start=1):
        ci_u[masks[f"u{year}|{c}"]] = j
        ci_r[masks[f"r{year}|{c}"]] = j

    # 按文件分组打开：同一文件开一次读完它的全部波段（50 MB 压缩 tif 开 183 次会明显变慢）
    for f, grp in groupby(entries, key=lambda e: e[0]):
        with rasterio.open(f) as src:
            for _, bi in grp:
                a = src.read(bi).astype(np.float32)
                a[a == FILL] = np.nan
                a *= SCALE_K
                flat = a.ravel()
                ok = np.isfinite(flat)
                fv = flat[ok]
                iu, ir = ci_u[ok], ci_r[ok]
                cnt_u = np.bincount(iu, minlength=nc + 1)
                sum_u = np.bincount(iu, weights=fv, minlength=nc + 1)
                cnt_r = np.bincount(ir, minlength=nc + 1)
                sum_r = np.bincount(ir, weights=fv, minlength=nc + 1)
                for j, c in enumerate(cities, start=1):
                    un, rn = int(cnt_u[j]), int(cnt_r[j])
                    d = acc[c]
                    d["u_sum"] += float(sum_u[j]); d["u_n"] += un
                    d["r_sum"] += float(sum_r[j]); d["r_n"] += rn
                    # 逐日原样存，阈值事后再套——与 18 同样的取舍（重读栅格贵，筛选便宜）
                    d["per_day"].append((float(sum_u[j] / un) if un else np.nan, un,
                                         float(sum_r[j] / rn) if rn else np.nan, rn))

    y_pw, y_eq, y_paired, n_pairs = [], [], [], []
    for c in cities:
        d = acc[c]
        if d["u_n"] == 0 or d["r_n"] == 0:
            continue
        # ① 逐像元加权（全池化）：等价于把每个 (像元,日) 都算一票
        pw = d["u_sum"] / d["u_n"] - d["r_sum"] / d["r_n"]
        # ② 每日等权（任一侧有像元即算）
        um_all = [um for um, un, _, _ in d["per_day"] if un > 0]
        rm_all = [rm for _, _, rm, rn in d["per_day"] if rn > 0]
        if len(um_all) < MIN_DAYS or len(rm_all) < MIN_DAYS:
            continue
        eq_any = float(np.mean(um_all)) - float(np.mean(rm_all))

        rec = {"city_id": c, "year": year, "mode": MODE, "n_days": len(d["per_day"]),
               "suhii_pixel_weighted": pw, "suhii_equal_any": eq_any,
               "weighting_effect": eq_any - pw,
               "u_px": u_px[c], "r_px": r_px[c],
               "u_days_per_px": d["u_n"] / u_px[c], "r_days_per_px": d["r_n"] / r_px[c]}
        for slab, minpx, minfrac in SCHEMES:
            use = [(um, rm) for um, un, rm, rn in d["per_day"]
                   if np.isfinite(um) and np.isfinite(rm)
                   and un >= max(minpx, minfrac * u_px[c])
                   and rn >= max(minpx, minfrac * r_px[c])]
            rec[f"n_paired_days_{slab}"] = len(use)
            rec[f"suhii_{MODE}_daily_paired_{slab}"] = (
                float(np.mean([um - rm for um, rm in use])) if len(use) >= MIN_DAYS else np.nan)
            # 等权 + 阈值（两侧各自达标，不要求同日）——对应 54 的 unpaired_equal_20px
            uo = [um for um, un, _, _ in d["per_day"] if un >= max(minpx, minfrac * u_px[c])]
            ro = [rm for _, _, rm, rn in d["per_day"] if rn >= max(minpx, minfrac * r_px[c])]
            rec[f"suhii_{MODE}_daily_equal_{slab}"] = (
                float(np.mean(uo)) - float(np.mean(ro)) if len(uo) >= MIN_DAYS and len(ro) >= MIN_DAYS
                else np.nan)
        rows.append(rec)
        if np.isfinite(rec[f"suhii_{MODE}_daily_paired_20px"]):
            y_pw.append(pw); y_eq.append(eq_any)
            y_paired.append(rec[f"suhii_{MODE}_daily_paired_20px"])
            n_pairs.append(rec["n_paired_days_20px"])

    if y_pw:
        print(f"{year:<7}{len(dates):<7}{len(acc):<6}{int(np.median(n_pairs)):<12}"
              f"{np.mean(y_pw):<11.3f}{np.mean(y_eq):<11.3f}{np.mean(y_paired):<11.3f}")

df = pd.DataFrame(rows)
df.to_csv(OUT, index=False, encoding="utf-8-sig")
print(f"\n已写 {OUT.relative_to(BASE)}  {len(df)} 个城市–年")

# ---------------------------------------------------------------- 判读
print(f"\n{'=' * 78}")
print("【逐日四口径（20px 口径）】")
print("=" * 78)
print(f"  {'年':<7}{'加权':>10}{'等权':>10}{'等权20px':>11}{'配对等权':>11}{'加权效应':>11}")
print("  " + "-" * 60)
for year, g in df.groupby("year"):
    pw = g.suhii_pixel_weighted.mean()
    eq = g.suhii_equal_any.mean()
    e20 = g[f"suhii_{MODE}_daily_equal_20px"].mean()
    pr = g[f"suhii_{MODE}_daily_paired_20px"].mean()
    print(f"  {year:<7}{pw:>10.3f}{eq:>10.3f}{e20:>11.3f}{pr:>11.3f}{eq - pw:>+11.3f}")

# ---------------------------------------------------------------- 与合成期对照
print(f"\n{'=' * 78}")
print("【与合成期口径对照（54_method_comparability 同源）】")
print("=" * 78)
comp_path = BASE / "outputs" / "method_comparability" / f"equal_period_{TAG}.csv"
if not comp_path.exists():
    print(f"  缺 {comp_path.relative_to(BASE)}，跳过对照")
else:
    comp = pd.read_csv(comp_path)
    print(f"  {'年':<7}{'城':<5}{'合成·加权':>11}{'合成·等权':>11}{'合成效应':>11}"
          f"{'逐日·加权':>11}{'逐日·等权':>11}{'逐日效应':>11}")
    print("  " + "-" * 80)
    cmp_rows = []
    for year, g in df.groupby("year"):
        cg = comp[comp.year == year]
        common = sorted(set(g.city_id) & set(cg.city_id))
        if not common:
            continue
        a = g[g.city_id.isin(common)]
        b = cg[cg.city_id.isin(common)]
        c_pw = b.suhii_night_unpaired.mean() if TAG == "night" else b.suhii_day_unpaired.mean()
        c_eq = b.unpaired_equal_any_K.mean()
        d_pw, d_eq = a.suhii_pixel_weighted.mean(), a.suhii_equal_any.mean()
        print(f"  {year:<7}{len(common):<5}{c_pw:>11.3f}{c_eq:>11.3f}{c_eq - c_pw:>+11.3f}"
              f"{d_pw:>11.3f}{d_eq:>11.3f}{d_eq - d_pw:>+11.3f}")
        cmp_rows.append({"mode": MODE, "year": year, "common_cities": len(common),
                         "composite_pixel_weighted_K": c_pw, "composite_equal_any_K": c_eq,
                         "composite_weighting_effect_K": c_eq - c_pw,
                         "daily_pixel_weighted_K": d_pw, "daily_equal_any_K": d_eq,
                         "daily_weighting_effect_K": d_eq - d_pw})
    OUTDIR = BASE / "outputs" / "daily"
    OUTDIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(cmp_rows).to_csv(OUTDIR / f"daily_vs_composite_{TAG}.csv",
                                  index=False, encoding="utf-8-sig")
    print(f"\n  对照表已写 outputs/daily/daily_vs_composite_{TAG}.csv")

print("""
  判读：
  1. 看**加权效应**两列。合成期 2020 = −0.683 K、2023 = −1.324 K。
     逐日若同号且显著非零，说明该缺口不是 MOD11A2 合成结构的产物。
     ⚠ 但**不要**算「逐日/合成」的保留率：那是两个产品的效应之比，
     本脚本第 58 行明文禁止；且该差的不确定性实测大到不可区分于零
     （见 `05_summarize_pilot.py` 的 paired_diff 列）。
  2. **本脚本不下斜率结论**。要斜率跑 `06_daily_trends.py`。
     加权缺口确实在逐年收缩，但收缩幅度**按年波动很大**，不是单调的：
     实测 2001–2009 均值 −2.70 K、2017–2025 均值 −1.65 K，
     而单年极值出现在 2005（−3.61）与 2015（−0.53）。
     口径：`outputs/method_comparability/equal_period_{night,day}.csv` 的
     `weighting_effect_K` 按年均值（与 `54_method_comparability.py` 同源）。
     那是**合成期**的「等权 − 像素加权」，**不是**本脚本产出的逐日列——
     两者不可混引；同一份面板里 `paired − unpaired` 是**配对效应**，
     也不是这一列（2026-10-02 实测两者数值相近而易被误当成同一个量）。
  3. `u_days_per_px` 与 `r_days_per_px` 是城乡观测对称性的体检。
     合成期口径实测（`data/panel_night_paired.csv`，scheme=20px，中位）：
     **全市集 2001–2025：城区 26.3 : 乡村 9.7 = 2.71×**；
     **2020 单年：城区 20.1 : 乡村 9.8 = 2.06×**；
     白天为**全市集 39.0 : 39.4 = 0.99×**、**2020 年 35.1 : 36.4 = 0.96×**。
     这四个数与正文 §3.2 **逐字一致**，改动任一侧都要同步另一侧。
     ⚠ 与**逐日**比值比大小前，必须**先取两者共同覆盖的年份与城集**，
     否则是跨样本比大小（第 58 行禁止）——2020+2023 同集实测为
     合成 **2.24×** → 逐日 **1.82×**，即压掉 34% 的超额不对称。
     逐日比值若显著缩小，说明逐日重采样本身改变了配对的性质，须在正文写明。
""")
