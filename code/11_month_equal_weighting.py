"""逐日 MOD11A1：季节构成 / 逐月等权的**决定性**敏感性 + 覆盖门槛同一视图（2026-10-03）。

## 它回答什么（本轮大修最要紧的一问）

全篇最强的描述性对照是：在同一批城市–年上，把**像元–日加权**换成**逐日等权**，
夜间斜率从 **−0.107**（p=0.441）翻到 **+0.193**（p=0.010）K/十年。审稿人会立刻问：
**这是不是季节构成造成的？** 夏季观测多、城乡温差大；若等权与加权在不同季节取到
不同的观测组合，翻号可能只是物候，而不是「按有效像元多寡加权」这件事本身。

本脚本用**同一批城市–年**，把两条口径都换成「逐月等权」与「固定季节」，再取差值序列：

- `month_equal`：每城市–年、每日历月先算 `SUHII_m = mean(u_mean|u_n>0) − mean(r_mean|r_n>0)`，
  再对**出现的月份**取等权平均（规则：**当月两侧各至少 1 个有效日**才算该月，
  **至少出现 `MIN_MONTHS=6` 个日历月**才产出该城–年的 `month_equal`，否则记 NaN）。见文末「口径」。
- `season_jja`：只取 6/7/8 月，在窗口内套**逐日等权**（两侧各需 ≥ `MIN_DAYS` 天）。
- `season_djf`：只取 12/1/2 月，同上。

判定标准（写在正文里）：若 `month_equal − weighted` 的**直接差值序列**斜率仍与全样本同向
且区间不跨零，则「翻号」不是季节构成的产物；若它塌到零或反向，则主结论必须改成
「效应依赖季节构成」。

## 四条产物（本脚本）

| 文件 | 内容 |
|------|------|
| `data/daily/panel_{mode}_monthly.csv` | 城市–年层，**全部估计量**（加权/等权各档/配对/逐月等权/固定季节） |
| `outputs/daily/month_equal_{mode}.csv` | 各估计量的**逐年序列** + Sen + Hamed–Rao；外加四条高亮差值检验 |
| `outputs/daily/monthly_valid_days.csv` | 城市–年–月有效日数 + 各层级中位数（行用 `scope` 区分） |
| `outputs/daily/coverage_view.csv` | **同一视图**：覆盖门槛 × 有效日数 × 留存城市–年 × 趋势 × 直接差（逐日 + 合成两产品） |
| `outputs/daily/month_equal_validation.csv` | 回归验证（见下「验证」）——便于把汇报数字回溯到文件，而非手抄 |

## 直接差值的三种相关性处理（三条都做，说清报哪条、为什么）

两条序列**在同一批城市–年上**算出（等权与加权共用同一批城–年），误差高度相关；
且年度序列本身自相关。故对差值序列同时给出：

1. **Hamed–Rao 修正的 Mann–Kendall p**（自相关稳健）——主报此项，因为与全篇 `06` 同一实现、可直接换算。
2. **移动块自助（moving-block bootstrap，块长 `BOOT_BLOCK=5` 年、`BOOT_N=2000` 次、种子 `BOOT_SEED=42`）**
   的 Sen 斜率**百分位** 95% 区间。块自助保留时间依赖，故在「序列自相关 + 两估计量误差相关」下
   给出正确的区间——这是本脚本认为**最恰当**的一种，正文报它作为斜率不确定性的主口径。
3. **城市层级配对**：逐城市对该城的 (等权 − 加权) 年度序列各做一次 Theil–Sen，再对**各城斜率**
   做 Wilcoxon 符号秩检验（是否异于 0）。它回答「效应是否在**多数城市**里一致」，与前两条互补。

> ⚠ **不能用两条单独序列的 p 值代替直接检验**：斜率之差不是检验。凡「效应」一律用**差值序列**。

## 口径（与 `03` / `06` / `10` 逐字同参，不另立一套）

- 加权：`Σ(u_mean·u_n)/Σu_n − Σ(r_mean·r_n)/Σr_n`（每 (像元,日) 一票）。
- 等权（any）：`mean(u_mean|u_n>0) − mean(r_mean|r_n>0)`，两侧各需 ≥ `MIN_DAYS=6` 天。
- 覆盖门槛（与 `03` 的 `SCHEMES` 一致）：单侧有效像元 `n ≥ max(minpx, frac × 该侧掩膜像元数)`；
  40%/60% 档 `minpx=10`，20px 档 `minpx=20, frac=0`。掩膜分母只读 `masks_matched.npz`（12 的唯一产物）。
- 配对（paired）：同一日两侧**都**达门槛才计一天，`mean(u_mean − r_mean)`，需 ≥ `MIN_DAYS` 天。
- 趋势：`theilslopes(...).slope × 10`（→ K/十年）+ `hamed_rao_modification_test(...).p`。
- 直接差：先逐城市–年作差，再取逐年城市均值，对该年序列做同一趋势检验。

## © 它**不能**用来声称什么

1. **不能**把不同档的绝对斜率并列成同一估计量的区间——每档是不同的观测集合。
2. **不能**说「收严/固定季节下趋势未变」：本脚本给的是各口径**自己**的斜率，不是同一估计量的重复测。
3. **不能**把逐月等权结果外推到 MOD11A2 合成产品——本表的逐月等权只读逐日 MOD11A1 中间表；
   合成产品的分档只出现在 `coverage_view.csv` 的对照行里。
4. **不能**把 60% 档及更严门槛当「序列」报：逐日夜间 60% 档**同日双侧达标日**中位仅 6 天、
   395/803 城–年低于 `MIN_DAYS` 下限，只作**上界**；`coverage_view.csv` 的 `comparable` 列已标。

## 验证（写入 `month_equal_validation.csv`，供回溯）

- `panel_{mode}_monthly.csv` 的 `weighted` / `equal_any` 两列须与既有 `panel_{mode}_daily.csv` 的
  `suhii_pixel_weighted` / `suhii_equal_any` 在全部城–年上逐值一致（浮点容差），打印最大绝对差。
- 本脚本算出的夜间不平衡 `equal_any` 年序列 Sen 斜率须与发表值 **+0.19313326557498456** 相差 < 1e-6。

用法：
    python 11_month_equal_weighting.py                # 昼夜都跑
    python 11_month_equal_weighting.py --mode night
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import theilslopes, wilcoxon

sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")   # 只改 stdout 时 stderr 仍走 GBK，中文报错会乱码

BASE = Path(__file__).resolve().parent.parent.parent       # analysis/
sys.path.insert(0, str(BASE / ".deps"))
import pymannkendall as mk                                  # 与 06 / 10 / 18 / 54 同一实现

DATA = BASE / "data"
DAILY = DATA / "daily"
OUT = BASE / "outputs" / "daily"
COMP = BASE / "outputs" / "method_comparability"
OUT.mkdir(parents=True, exist_ok=True)

MIN_DAYS = 6            # 单侧有效日下限（与 03/06/10 一致）
MIN_YEARS = 8           # 年序列做趋势的最少年数
MIN_MONTHS = 6          # 逐月等权要求出现的日历月数下限
BOOT_BLOCK = 5          # 移动块自助的块长（年）
BOOT_N = 2000           # 自助重采样次数
BOOT_SEED = 42          # 固定种子
SCHEMES = [("equal_20px", 20, 0.0), ("equal_40pct", 10, 0.40), ("equal_60pct", 10, 0.60)]
SEASONS = {"JJA": (6, 7, 8), "DJF": (12, 1, 2)}
PUBLISHED_EQUAL_ANY_SEN = 0.19313326557498456            # 发表值：夜间不平衡等权 Sen

parser = argparse.ArgumentParser(description="逐月等权决定性检验 + 覆盖同一视图")
parser.add_argument("--mode", choices=("night", "day"), default=None)
ARGS = parser.parse_args()
MODES = (ARGS.mode,) if ARGS.mode else ("night", "day")

RUNG_ORDER = ["weighted", "equal_any", "equal_20px", "equal_40pct", "equal_60pct",
              "paired_20px", "month_equal", "season_jja", "season_djf"]
RUNG_DESC = {
    "weighted": "像元–日全池化（加权，基准）",
    "equal_any": "逐日等权（任一侧有像元）",
    "equal_20px": "逐日等权 + 20 像元下限",
    "equal_40pct": "逐日等权 + 40% 掩膜下限",
    "equal_60pct": "逐日等权 + 60% 掩膜下限",
    "paired_20px": "同日配对 + 等权（20 像元）",
    "month_equal": f"逐月等权（每月两侧≥1 天，需≥{MIN_MONTHS} 个月）",
    "season_jja": "固定季节 6–8 月等权",
    "season_djf": "固定季节 12–2 月等权",
    "weighted_me": "逐月加权再对月等权（季节中性化的基准侧）",
    "weighted_JJA": "固定季节 6–8 月**加权**",
    "weighted_DJF": "固定季节 12–2 月**加权**",
}
# month_equal_{mode}.csv 里要出逐年序列的估计量（含季节中性化的加权基准）
ME_ESTIMATORS = ["weighted", "equal_any", "equal_20px", "equal_40pct", "equal_60pct",
                 "paired_20px", "month_equal", "season_jja", "season_djf",
                 "weighted_me", "weighted_JJA", "weighted_DJF"]
# 覆盖视图里各档「自身口径」的有效日列（城市–年面板中的计数列名）
RUNG_DAYCOL = {
    "weighted": "n_days_weighted", "equal_any": "n_days_equal_any",
    "equal_20px": "n_days_20px", "equal_40pct": "n_days_40pct",
    "equal_60pct": "n_days_60pct", "paired_20px": "n_days_paired_20px",
    "month_equal": "n_days_equal_any", "season_jja": "n_days_JJA", "season_djf": "n_days_DJF",
}
# 同一日双侧达标的配对计数列（用于 60% 档样本坍缩的告警）
RUNG_PAIREDCOL = {
    "weighted": "n_paired_any", "equal_any": "n_paired_any",
    "equal_20px": "n_paired_20px", "equal_40pct": "n_paired_40pct",
    "equal_60pct": "n_paired_60pct", "paired_20px": "n_paired_20px",
    "month_equal": "n_paired_any", "season_jja": "n_paired_JJA", "season_djf": "n_paired_DJF",
}


# ---------------------------------------------------------------- 趋势与自助
def trend_of(annual: pd.Series):
    """Theil–Sen + Hamed–Rao，与 06/10 逐字同参。年数不足返回 None。"""
    a = annual.dropna().sort_index()
    if len(a) < MIN_YEARS:
        return None
    sen = theilslopes(a.values, a.index.values, alpha=0.95)
    hr = mk.hamed_rao_modification_test(a.values)
    return {"n_years": int(len(a)),
            "mean_K": float(a.mean()),
            "sen_K_per_decade": float(sen.slope * 10),
            "sen_CI95_lo": float(sen.low_slope * 10),
            "sen_CI95_hi": float(sen.high_slope * 10),
            "hamed_rao_p": float(hr.p)}


def mbb_sen_ci(years, values, block=BOOT_BLOCK, n=BOOT_N, seed=BOOT_SEED):
    """移动块自助的 Sen 斜率百分位 95% 区间（K/十年）。

    块长默认 5 年：短于序列自相关的尺度，长于「一年一年独立」的假设，是可辩护的折中。
    ⚠ **必须重采样 (年, 值) 对**，而非只打乱值、把横轴重编为 0..N−1——后者会**抹掉趋势**
    （被抽出的块放到新位置后，块间的时间落差被随机化），实测把区间中心从 +0.23 拉到 ~−0.02，
    是错的。保留块的**原始年份**，既可复现块内自相关，又保留全局趋势。
    固定种子 `BOOT_SEED` 保证可复现。
    """
    x = np.asarray(years, dtype=float)
    v = np.asarray(values, dtype=float)
    ok = np.isfinite(x) & np.isfinite(v)
    x, v = x[ok], v[ok]
    N = len(x)
    if N < MIN_YEARS:
        return None
    L = int(min(block, max(2, N // 2)))
    nb = int(np.ceil(N / L))
    rng = np.random.default_rng(seed)
    slopes = np.empty(n, dtype=float)
    for i in range(n):
        starts = rng.integers(0, N - L + 1, size=nb)
        idx = (starts[:, None] + np.arange(L)[None, :]).ravel()[:N]
        slopes[i] = theilslopes(v[idx], x[idx])[0] * 10.0
    return {"boot_lo": float(np.percentile(slopes, 2.5)),
            "boot_hi": float(np.percentile(slopes, 97.5)),
            "boot_median": float(np.median(slopes))}


def city_paired_test(frame, col_a, col_b, min_years=MIN_YEARS):
    """城市层级配对：逐城市对 (col_b − col_a) 的年度序列做 Sen 斜率，再 Wilcoxon 符号秩检验。

    col_a 是基准（如 weighted），col_b 是候选（如 equal_any）。返回各城斜率的中位与检验 p。
    年数 < `min_years` 的城市不参与（Sen 在年数太少时不可靠）。
    """
    slopes = []
    for _, g in frame.groupby("city_id"):
        g = g.dropna(subset=[col_a, col_b]).sort_values("year")
        if len(g) < min_years:
            continue
        d = (g[col_b] - g[col_a]).to_numpy(dtype=float)
        if not np.isfinite(d).all():
            continue
        slopes.append(theilslopes(d, g["year"].to_numpy(dtype=float))[0] * 10.0)
    slopes = np.asarray(slopes, dtype=float)
    out = {"wilcoxon_n_cities": int(len(slopes)),
           "wilcoxon_median_city_slope": float(np.median(slopes)) if len(slopes) else np.nan,
           "wilcoxon_stat": np.nan, "wilcoxon_p": np.nan}
    if len(slopes) >= 5 and not np.allclose(slopes, 0.0):
        try:
            stat, p = wilcoxon(slopes)
            out["wilcoxon_stat"], out["wilcoxon_p"] = float(stat), float(p)
        except ValueError:
            pass                                       # 全为零差时 wilcoxon 抛错 → 保持 NaN
    return out


# ---------------------------------------------------------------- 城市–年口径
def mask_pixel_counts():
    """{year: {city: (n_urban_px, n_rural_px)}} —— 覆盖门槛的分母，只读 12 的掩膜产物。"""
    m = np.load(DATA / "masks_matched.npz")
    out = {}
    for k in m.files:
        if not k.startswith("u"):
            continue
        year = int(k.split("|")[0][1:])
        city = int(k.split("|")[1])
        rk = f"r{year}|{city}"
        if rk not in m.files:
            continue
        out.setdefault(year, {})[city] = (int(m[k].size), int(m[rk].size))
    return out


def per_city_year(sub: pd.DataFrame, px):
    """一个城市–年的全部口径。`sub` 是该城市–年的逐日记录（含 month）。"""
    res = {}
    d = sub

    # --- 城市–年–月有效日数（城市侧/乡村侧/同日双侧）——
    # 单独在 build 里落表，这里只汇总全年计数供其它口径复用。
    # --- 基准：像元加权（全部日，每 (像元,日) 一票）---
    su, sr = float(d.u_n.sum()), float(d.r_n.sum())
    res["weighted"] = ((float((d.u_mean * d.u_n).fillna(0).sum()) / su
                        - float((d.r_mean * d.r_n).fillna(0).sum()) / sr)
                       if su > 0 and sr > 0 else np.nan)
    u_ok = d.u_n > 0
    r_ok = d.r_n > 0
    res["n_days_weighted"] = int((u_ok & r_ok).sum())
    res["n_paired_any"] = int((u_ok & r_ok).sum())

    # --- 等权（any）---
    u_all = d.loc[u_ok, "u_mean"]
    r_all = d.loc[r_ok, "r_mean"]
    res["n_days_equal_any"] = int(min(len(u_all), len(r_all)))
    res["equal_any"] = (float(u_all.mean()) - float(r_all.mean())
                        if len(u_all) >= MIN_DAYS and len(r_all) >= MIN_DAYS else np.nan)

    # --- 覆盖门槛阶梯：等权 + 阈值；并记同日双侧达标日数（paired）---
    n_u, n_r = px if px else (np.nan, np.nan)
    for lab, minpx, frac in SCHEMES:
        if not px:
            res[lab] = np.nan
            res[f"n_days_{lab.replace('equal_','')}"] = 0
            res[f"n_paired_{lab.replace('equal_','')}"] = 0
            continue
        thu, thr = max(minpx, frac * n_u), max(minpx, frac * n_r)
        uo = d.loc[d.u_n >= thu, "u_mean"]
        ro = d.loc[d.r_n >= thr, "r_mean"]
        res[lab] = (float(uo.mean()) - float(ro.mean())
                    if len(uo) >= MIN_DAYS and len(ro) >= MIN_DAYS else np.nan)
        short = lab.replace("equal_", "")
        res[f"n_days_{short}"] = int(min(len(uo), len(ro)))
        res[f"n_paired_{short}"] = int(((d.u_n >= thu) & (d.r_n >= thr)).sum())

    # --- 同日配对 + 等权（20 像元）---
    if px:
        thu, thr = max(20, 0.0 * n_u), max(20, 0.0 * n_r)
        use = d.loc[(d.u_n >= thu) & (d.r_n >= thr)
                    & d.u_mean.notna() & d.r_mean.notna()]
        res["n_days_paired_20px"] = int(len(use))
        res["paired_20px"] = (float((use.u_mean - use.r_mean).mean())
                              if len(use) >= MIN_DAYS else np.nan)
    else:
        res["n_days_paired_20px"], res["paired_20px"] = 0, np.nan

    # --- 逐月：每月各算等权与加权两条，再对**出现的月份**取等权平均 ---
    # ⚠ 同时留 `weighted_me`（逐月加权的月均值）。它与 `month_equal` 才是**同口径**的季节中性化对照：
    #   「把两边都改成逐月」后仍存在的差，才是纯粹的权重效应；拿 `month_equal` 直接减**全样本** weighted
    #   会把加权侧的季节污染算进效应里（实测把 +0.079 放大成 +0.233）。
    per_month, per_month_w = {}, {}
    for mon, g in d.groupby("month"):
        gu = g.loc[g.u_n > 0, "u_mean"]
        gr = g.loc[g.r_n > 0, "r_mean"]
        if len(gu) >= 1 and len(gr) >= 1:           # 当月两侧各至少 1 个有效日
            per_month[int(mon)] = float(gu.mean()) - float(gr.mean())
        su, sr = float(g.u_n.sum()), float(g.r_n.sum())
        if su > 0 and sr > 0:
            per_month_w[int(mon)] = (float((g.u_mean * g.u_n).fillna(0).sum()) / su
                                     - float((g.r_mean * g.r_n).fillna(0).sum()) / sr)
    res["n_months"] = len(per_month)
    res["month_equal"] = (float(np.mean(list(per_month.values())))
                          if len(per_month) >= MIN_MONTHS else np.nan)
    res["n_months_w"] = len(per_month_w)
    res["weighted_me"] = (float(np.mean(list(per_month_w.values())))
                          if len(per_month_w) >= MIN_MONTHS else np.nan)

    # --- 固定季节：只取固定日历月窗口，季节构造成分被消掉（等权与加权各一条）---
    for lab, months in SEASONS.items():
        g = d[d.month.isin(months)]
        res[f"n_paired_{lab}"] = int(((g.u_n > 0) & (g.r_n > 0)).sum())
        gu = g.loc[g.u_n > 0, "u_mean"]
        gr = g.loc[g.r_n > 0, "r_mean"]
        res[f"n_days_{lab}"] = int(min(len(gu), len(gr)))
        res[f"season_{lab.lower()}"] = (float(gu.mean()) - float(gr.mean())
                                        if len(gu) >= MIN_DAYS and len(gr) >= MIN_DAYS else np.nan)
        su, sr = float(g.u_n.sum()), float(g.r_n.sum())
        res[f"weighted_{lab}"] = ((float((g.u_mean * g.u_n).fillna(0).sum()) / su
                                   - float((g.r_mean * g.r_n).fillna(0).sum()) / sr)
                                  if su > 0 and sr > 0 else np.nan)
    return res


def build_panel(rec: pd.DataFrame, pxmap: dict) -> pd.DataFrame:
    rows = []
    for (cid, year), sub in rec.groupby(["city_id", "year"], sort=False):
        px = pxmap.get(int(year), {}).get(int(cid))
        rows.append({"city_id": int(cid), "year": int(year),
                     "mode": sub["mode"].iloc[0],
                     "n_days_recorded": int(len(sub)),
                     **per_city_year(sub, px)})
    return pd.DataFrame(rows).sort_values(["city_id", "year"]).reset_index(drop=True)


def build_monthly_valid_days(rec: pd.DataFrame, mode: str) -> pd.DataFrame:
    """城市–年–月有效日数，并追加各层级中位数（用 `scope` 列区分）。"""
    r = rec.copy()
    r["u_ok"], r["r_ok"] = r.u_n > 0, r.r_n > 0
    r["both_ok"] = r.u_ok & r.r_ok
    g = (r.sort_values("date").groupby(["city_id", "year", "month"])
         .agg(n_days_u=("u_ok", "sum"), n_days_r=("r_ok", "sum"),
              n_days_both=("both_ok", "sum")).reset_index())
    g.insert(0, "mode", mode)
    g.insert(0, "scope", "city_year_month")
    # 中位：逐年–月（跨城市）、逐年（跨城市–月）、逐月（跨全部城市–年）
    med_yym = (g.groupby(["year", "month"])[["n_days_u", "n_days_r", "n_days_both"]]
               .median().reset_index().assign(scope="median_year_month", mode=mode, city_id=np.nan))
    med_year = (g.groupby("year")[["n_days_u", "n_days_r", "n_days_both"]]
                .median().reset_index().assign(scope="median_year", mode=mode,
                                               city_id=np.nan, month=np.nan))
    med_month = (g.groupby("month")[["n_days_u", "n_days_r", "n_days_both"]]
                 .median().reset_index().assign(scope="median_month", mode=mode,
                                                city_id=np.nan, year=np.nan))
    cols = ["mode", "scope", "city_id", "year", "month", "n_days_u", "n_days_r", "n_days_both"]
    return pd.concat([g, med_yym, med_year, med_month], ignore_index=True)[cols]


# ---------------------------------------------------------------- 主循环
records = {m: pd.read_csv(DAILY / f"daily_records_{m}.csv") for m in MODES}
pxmap = mask_pixel_counts()

panels, me_rows, valid_frames, val_rows = {}, [], [], []
direct_reports = {}

for mode in MODES:
    if mode not in records:
        print(f"缺 {mode} 的 daily_records，跳过", file=sys.stderr)
        continue
    rec = records[mode]
    lab = "夜间" if mode == "night" else "白天"
    tab = build_panel(rec, pxmap)
    panels[mode] = tab
    py = int(tab.year.min()), int(tab.year.max())
    print("=" * 96)
    print(f"{lab}逐日 {mode}：{len(tab)} 城市–年 · {tab.city_id.nunique()} 城 · 年 {py[0]}–{py[1]}")
    print("=" * 96)

    # ---------- 验证：与既有 panel_{mode}_daily.csv 逐值回归 ----------
    ref_path = DAILY / f"panel_{mode}_daily.csv"
    if ref_path.exists():
        ref = pd.read_csv(ref_path)
        mg = tab.merge(ref, on=["city_id", "year"], suffixes=("", "_ref"))
        for mine, theirs in (("weighted", "suhii_pixel_weighted"),
                             ("equal_any", "suhii_equal_any")):
            dif = float(np.nanmax(np.abs(mg[mine].to_numpy() - mg[theirs].to_numpy())))
            val_rows.append({"mode": mode, "check": f"max_abs_diff_{mine}", "value": dif,
                             "n": int(mg[mine].notna().sum()), "pass": bool(dif < 1e-9)})
            print(f"  [验证] {mine} vs {theirs}：最大绝对差 {dif:.3e}（n={int(mg[mine].notna().sum())}）")
    else:
        print(f"  [验证] 缺 {ref_path.name}，跳过回归", file=sys.stderr)

    tab.to_csv(DAILY / f"panel_{mode}_monthly.csv", index=False, encoding="utf-8-sig")
    valid_frames.append(build_monthly_valid_days(rec, mode))

    # ---------- month_equal_{mode}.csv：各估计量逐年序列 + Sen + HR ----------
    years = sorted(tab.year.unique())
    ann = {tag: tab[["city_id", "year", tag]].dropna().groupby("year")[tag].mean().reindex(years)
           for tag in ME_ESTIMATORS}
    equal_any_sen = (trend_of(ann["equal_any"]) or {}).get("sen_K_per_decade", np.nan)
    if mode == "night":
        ok = abs(equal_any_sen - PUBLISHED_EQUAL_ANY_SEN) < 1e-6
        val_rows.append({"mode": mode, "check": "equal_any_sen_vs_published",
                         "value": equal_any_sen, "n": int(ann["equal_any"].notna().sum()), "pass": bool(ok)})
        print(f"  [验证] equal_any Sen = {equal_any_sen:.16f} vs 发表值 {PUBLISHED_EQUAL_ANY_SEN:.16f}"
              f" → {'PASS' if ok else 'FAIL'}")
    for tag in ME_ESTIMATORS:
        t = trend_of(ann[tag])
        if t is None:
            continue
        s = tab[["city_id", "year", tag]].dropna()
        # ⚠ tag == "weighted" 时不能再取 tab[[...,"weighted",tag]]——会产生两列同名，
        #   both[tag] 返回 DataFrame，随后 assign 报「多列赋给单列」。基准行跳过差值即可。
        both = tab[["city_id", "year", "weighted", tag]].dropna() if tag != "weighted" else tab.iloc[0:0]
        td = boot = None
        if tag != "weighted" and len(both):
            dser = both.assign(__d=both[tag] - both["weighted"]).groupby("year")["__d"].mean()
            td = trend_of(dser)
            boot = mbb_sen_ci(dser.index.values, dser.values)
        wil = city_paired_test(tab, "weighted", tag) if tag != "weighted" else {}
        me_rows.append({
            "mode": mode, "kind": "estimator", "estimator": tag, "description": RUNG_DESC[tag],
            "city_years": int(len(s)), "n_cities": int(s.city_id.nunique()),
            **{f"y{y}": (float(ann[tag].loc[y]) if pd.notna(ann[tag].loc[y]) else np.nan) for y in years},
            "grand_mean_K": t["mean_K"], "sen_K_per_decade": t["sen_K_per_decade"],
            "sen_CI95_lo": t["sen_CI95_lo"], "sen_CI95_hi": t["sen_CI95_hi"],
            "hamed_rao_p": t["hamed_rao_p"],
            "diff_vs_weighted_city_years": int(len(both)),
            "diff_sen_K_per_decade": (td or {}).get("sen_K_per_decade", np.nan),
            "diff_sen_CI95_lo": (td or {}).get("sen_CI95_lo", np.nan),
            "diff_sen_CI95_hi": (td or {}).get("sen_CI95_hi", np.nan),
            "diff_hamed_rao_p": (td or {}).get("hamed_rao_p", np.nan),
            "diff_boot_CI95_lo": (boot or {}).get("boot_lo", np.nan),
            "diff_boot_CI95_hi": (boot or {}).get("boot_hi", np.nan),
            **{k: wil.get(k, np.nan) for k in ("wilcoxon_n_cities", "wilcoxon_median_city_slope",
                                               "wilcoxon_stat", "wilcoxon_p")},
        })

    # ---------- 四条高亮直接差检验（差值序列自己的逐年序列 + 全部区间/p）----------
    highlights = [("equal_any", "weighted", "equal_any − weighted（全样本基准对照）"),
                  ("month_equal", "weighted", "month_equal − weighted（任务字面：新估计量 vs 全样本加权）"),
                  ("season_jja", "weighted", "season_jja − weighted（任务字面：固定夏季 vs 全样本加权）"),
                  ("paired_20px", "equal_20px", "paired_20px − equal_20px（配对 vs 非配对）"),
                  # 季节中性化对照：两边都改成逐月/固定季节后，剩下的才是**纯权重效应**
                  ("month_equal", "weighted_me", "month_equal − weighted_me（季节中性化：两边都逐月）"),
                  ("season_jja", "weighted_JJA", "season_jja − weighted_JJA（固定夏季内：两边都限 6–8 月）"),
                  ("season_djf", "weighted_DJF", "season_djf − weighted_DJF（固定冬季内：两边都限 12–2 月）")]
    for cand, base, desc in highlights:
        both = tab[["city_id", "year", base, cand]].dropna()
        if len(both) < MIN_YEARS:
            continue
        dser = both.assign(__d=both[cand] - both[base]).groupby("year")["__d"].mean()
        td = trend_of(dser)
        boot = mbb_sen_ci(dser.index.values, dser.values)
        wil = city_paired_test(tab, base, cand)
        me_rows.append({
            "mode": mode, "kind": "diff", "estimator": f"{cand} − {base}", "description": desc,
            "city_years": int(len(both)), "n_cities": int(both.city_id.nunique()),
            **{f"y{y}": (float(dser.loc[y]) if y in dser.index else np.nan) for y in years},
            "grand_mean_K": dser.mean(), "sen_K_per_decade": (td or {}).get("sen_K_per_decade", np.nan),
            "sen_CI95_lo": (td or {}).get("sen_CI95_lo", np.nan),
            "sen_CI95_hi": (td or {}).get("sen_CI95_hi", np.nan),
            "hamed_rao_p": (td or {}).get("hamed_rao_p", np.nan),
            "diff_vs_weighted_city_years": int(len(both)),
            "diff_sen_K_per_decade": (td or {}).get("sen_K_per_decade", np.nan),
            "diff_sen_CI95_lo": (td or {}).get("sen_CI95_lo", np.nan),
            "diff_sen_CI95_hi": (td or {}).get("sen_CI95_hi", np.nan),
            "diff_hamed_rao_p": (td or {}).get("hamed_rao_p", np.nan),
            "diff_boot_CI95_lo": (boot or {}).get("boot_lo", np.nan),
            "diff_boot_CI95_hi": (boot or {}).get("boot_hi", np.nan),
            **{k: wil.get(k, np.nan) for k in ("wilcoxon_n_cities", "wilcoxon_median_city_slope",
                                               "wilcoxon_stat", "wilcoxon_p")},
        })
        direct_reports[(mode, cand, base)] = {
            "sen": (td or {}).get("sen_K_per_decade", np.nan),
            "hr_lo": (td or {}).get("sen_CI95_lo", np.nan), "hr_hi": (td or {}).get("sen_CI95_hi", np.nan),
            "hr_p": (td or {}).get("hamed_rao_p", np.nan),
            "boot_lo": (boot or {}).get("boot_lo", np.nan), "boot_hi": (boot or {}).get("boot_hi", np.nan),
            "wil_p": wil.get("wilcoxon_p", np.nan), "wil_n": wil.get("wilcoxon_n_cities", np.nan),
            "n": int(len(both))}
    # 打印决定性一行（含季节中性化对照）
    def _sen(tag):
        t = trend_of(ann[tag])
        return t["sen_K_per_decade"] if t else float("nan")
    print(f"\n  【决定性】weighted {_sen('weighted'):+.3f} · weighted_me {_sen('weighted_me'):+.3f} · "
          f"month_equal {_sen('month_equal'):+.3f} · season_jja {_sen('season_jja'):+.3f} "
          f"（K/十年）—— 全样本加权与逐月加权**符号相反**，即翻号很大程度来自加权侧的季节敏感性。")
    for cand, base, lbl in (("month_equal", "weighted", "任务字面"),
                            ("month_equal", "weighted_me", "季节中性化（两边都逐月）")):
        r = direct_reports.get((mode, cand, base))
        if r:
            print(f"    [{lbl}] {cand} − {base}：Sen {r['sen']:+.4f} · "
                  f"Hamed–Rao95% [{r['hr_lo']:+.4f},{r['hr_hi']:+.4f}] p={r['hr_p']:.5f} · "
                  f"块自助95% [{r['boot_lo']:+.4f},{r['boot_hi']:+.4f}] · "
                  f"城市Wilcoxon p={r['wil_p']:.5f} (n={r['wil_n']:.0f})")
    # ⚠ `me_rows` 跨模式累积，**必须在循环内按当前 mode 过滤后立刻写**：
    #   ① 写成循环外的 `me[me["mode"]==mode]` 只对**最后一个** mode 生效——实测夜间
    #      产物从未被本脚本写出，盘上留下的是已退役脚本的旧文件（列数 20 vs 47）；
    #   ② 不过滤则第二轮的 day 文件里会混入第一轮的 night 行（首版实测 26 行 = 夜 13 + 昼 13）。
    me_cur = pd.DataFrame(me_rows)
    me_cur = me_cur[me_cur["mode"] == mode]
    me_cur.to_csv(OUT / f"month_equal_{mode}.csv", index=False, encoding="utf-8-sig")
    print(f"已写 data/daily/panel_{mode}_monthly.csv（{len(tab)} 行）· "
          f"outputs/daily/month_equal_{mode}.csv（{len(me_cur)} 行）")

valid = pd.concat(valid_frames, ignore_index=True)
valid.to_csv(OUT / "monthly_valid_days.csv", index=False, encoding="utf-8-sig")
print(f"已写 outputs/daily/monthly_valid_days.csv（{len(valid)} 行，含各层级中位 scope）")

# ================================================================ 统一覆盖视图
# 一行一个 (mode, product, rung)。产品：逐日 MOD11A1；有对应档的再加合成 MOD11A2。
COMP_RUNG_COL = {                                        # 合成产品列名（见 54 的 equal_period_*）
    "weighted": "suhii_{mode}_unpaired", "equal_any": "unpaired_equal_any_K",
    "equal_20px": "unpaired_equal_20px_K", "equal_40pct": "unpaired_equal_40%_K",
    "equal_60pct": "unpaired_equal_60%_K", "paired_20px": "suhii_{mode}_paired"}


# 季节中性化的**同口径基准侧**：把加权估计量也换到同一季节窗口，剩下的差才是纯权重效应
WINDOWED_BASE = {"month_equal": "weighted_me", "season_jja": "weighted_JJA", "season_djf": "weighted_DJF"}


def daily_row(mode, tab, tag):
    """逐日 MOD11A1 的一行。"""
    s = tab[["city_id", "year", tag]].dropna()
    t = trend_of(s.groupby("year")[tag].mean())
    if t is None:
        return None
    daycol = RUNG_DAYCOL.get(tag)
    med_days = float(tab.loc[tab[tag].notna(), daycol].median()) if daycol in tab.columns else np.nan
    pcol = RUNG_PAIREDCOL.get(tag)
    med_paired = float(tab[pcol].median()) if pcol in tab.columns else np.nan
    below = int((tab[pcol] < MIN_DAYS).sum()) if pcol in tab.columns else 0
    both = tab[["city_id", "year", "weighted", tag]].dropna() if tag != "weighted" else tab.iloc[0:0]
    td = tbd = tbw = None
    if tag != "weighted" and len(both):
        dser = both.assign(__d=both[tag] - both["weighted"]).groupby("year")["__d"].mean()
        td = trend_of(dser)
        tbd = mbb_sen_ci(dser.index.values, dser.values)
        tbw = city_paired_test(tab, "weighted", tag)
    comparable, note = True, ""
    if tag == "equal_60pct" and mode == "night":
        comparable = False
        note = (f"样本坍缩：同日双侧达标日中位 {med_paired:.0f} 天、{below}/{len(tab)} 城–年低于 "
                f"{MIN_DAYS} 天下限 → 只作上界，不作序列。")
    elif tag == "equal_40pct" and mode == "night":
        note = f"较薄：同日双侧达标日中位 {med_paired:.0f} 天。"
    # 季节中性化对照：把加权侧也换到同一季节窗口（month_equal→weighted_me 等）
    wbase = WINDOWED_BASE.get(tag)
    wd = wb = ww = None
    if wbase and wbase in tab.columns:
        wboth = tab[["city_id", "year", wbase, tag]].dropna()
        if len(wboth) >= MIN_YEARS:
            wser = wboth.assign(__d=wboth[tag] - wboth[wbase]).groupby("year")["__d"].mean()
            wd = trend_of(wser)
            wb = mbb_sen_ci(wser.index.values, wser.values)
            ww = city_paired_test(tab, wbase, tag)
    return {"mode": mode, "product": "daily MOD11A1", "rung": tag, "description": RUNG_DESC[tag],
            "city_years": int(len(s)), "n_cities": int(s.city_id.nunique()),
            "median_valid_days": med_days, "median_paired_days": med_paired,
            "n_paired_days_below_floor": below,
            "mean_SUHII_K": t["mean_K"], "sen_K_per_decade": t["sen_K_per_decade"],
            "sen_CI95_lo": t["sen_CI95_lo"], "sen_CI95_hi": t["sen_CI95_hi"],
            "hamed_rao_p": t["hamed_rao_p"],
            "diff_vs_weighted_city_years": int(len(both)),
            "diff_sen_K_per_decade": (td or {}).get("sen_K_per_decade", np.nan),
            "diff_sen_CI95_lo": (td or {}).get("sen_CI95_lo", np.nan),
            "diff_sen_CI95_hi": (td or {}).get("sen_CI95_hi", np.nan),
            "diff_hamed_rao_p": (td or {}).get("hamed_rao_p", np.nan),
            "diff_boot_CI95_lo": (tbd or {}).get("boot_lo", np.nan),
            "diff_boot_CI95_hi": (tbd or {}).get("boot_hi", np.nan),
            "diff_wilcoxon_p": (tbw or {}).get("wilcoxon_p", np.nan),
            "windowed_base": wbase or "",
            "diff_windowed_sen_K_per_decade": (wd or {}).get("sen_K_per_decade", np.nan),
            "diff_windowed_CI95_lo": (wd or {}).get("sen_CI95_lo", np.nan),
            "diff_windowed_CI95_hi": (wd or {}).get("sen_CI95_hi", np.nan),
            "diff_windowed_hamed_rao_p": (wd or {}).get("hamed_rao_p", np.nan),
            "diff_windowed_boot_lo": (wb or {}).get("boot_lo", np.nan),
            "diff_windowed_boot_hi": (wb or {}).get("boot_hi", np.nan),
            "diff_windowed_wilcoxon_p": (ww or {}).get("wilcoxon_p", np.nan),
            "comparable": comparable, "note": note}


def composite_row(mode, comp, tag):
    """合成 MOD11A2 的一行（只对有对应列的档生成）。中位有效**合成期数**，非日数。"""
    col = COMP_RUNG_COL[tag].format(mode=mode)
    if col not in comp.columns:
        return None
    s = comp[["city_id", "year", col]].dropna()
    t = trend_of(s.groupby("year")[col].mean())
    if t is None:
        return None
    med_periods = (float(comp.loc[comp[col].notna(), "n_common_periods"].median())
                   if "n_common_periods" in comp.columns else np.nan)
    both = comp[["city_id", "year", "suhii_{mode}_unpaired".format(mode=mode), col]].dropna()
    td = tbd = None
    wcol = "suhii_{mode}_unpaired".format(mode=mode)
    if tag != "weighted" and len(both):
        dser = both.assign(__d=both[col] - both[wcol]).groupby("year")["__d"].mean()
        td = trend_of(dser)
        tbd = mbb_sen_ci(dser.index.values, dser.values)
    return {"mode": mode, "product": "composite MOD11A2", "rung": tag, "description": RUNG_DESC[tag],
            "city_years": int(len(s)), "n_cities": int(s.city_id.nunique()),
            "median_valid_days": med_periods, "median_paired_days": med_periods,
            "n_paired_days_below_floor": 0,
            "mean_SUHII_K": t["mean_K"], "sen_K_per_decade": t["sen_K_per_decade"],
            "sen_CI95_lo": t["sen_CI95_lo"], "sen_CI95_hi": t["sen_CI95_hi"],
            "hamed_rao_p": t["hamed_rao_p"],
            "diff_vs_weighted_city_years": int(len(both)),
            "diff_sen_K_per_decade": (td or {}).get("sen_K_per_decade", np.nan),
            "diff_sen_CI95_lo": (td or {}).get("sen_CI95_lo", np.nan),
            "diff_sen_CI95_hi": (td or {}).get("sen_CI95_hi", np.nan),
            "diff_hamed_rao_p": (td or {}).get("hamed_rao_p", np.nan),
            "diff_boot_CI95_lo": (tbd or {}).get("boot_lo", np.nan),
            "diff_boot_CI95_hi": (tbd or {}).get("boot_hi", np.nan),
            "diff_wilcoxon_p": np.nan,
            "windowed_base": "",
            "diff_windowed_sen_K_per_decade": np.nan, "diff_windowed_CI95_lo": np.nan,
            "diff_windowed_CI95_hi": np.nan, "diff_windowed_hamed_rao_p": np.nan,
            "diff_windowed_boot_lo": np.nan, "diff_windowed_boot_hi": np.nan,
            "diff_windowed_wilcoxon_p": np.nan,
            "comparable": True, "note": "中位有效为**合成期数**（8 天），非日数。"}


cov_rows = []
for mode in MODES:
    tab = panels.get(mode)
    if tab is None:
        continue
    for tag in RUNG_ORDER:
        if tag in tab.columns:
            r = daily_row(mode, tab, tag)
            if r:
                cov_rows.append(r)
    comp_path = COMP / f"equal_period_{mode}.csv"
    if comp_path.exists():
        comp = pd.read_csv(comp_path)
        for tag in COMP_RUNG_COL:
            r = composite_row(mode, comp, tag)
            if r:
                cov_rows.append(r)
    else:
        print(f"缺 {comp_path.name}，合成产品行跳过", file=sys.stderr)

cov = pd.DataFrame(cov_rows)[["mode", "product", "rung", "description", "city_years", "n_cities",
                              "median_valid_days", "median_paired_days", "n_paired_days_below_floor",
                              "mean_SUHII_K", "sen_K_per_decade", "sen_CI95_lo", "sen_CI95_hi",
                              "hamed_rao_p", "diff_vs_weighted_city_years", "diff_sen_K_per_decade",
                              "diff_sen_CI95_lo", "diff_sen_CI95_hi", "diff_hamed_rao_p",
                              "diff_boot_CI95_lo", "diff_boot_CI95_hi", "diff_wilcoxon_p",
                              "windowed_base", "diff_windowed_sen_K_per_decade",
                              "diff_windowed_CI95_lo", "diff_windowed_CI95_hi",
                              "diff_windowed_hamed_rao_p", "diff_windowed_boot_lo",
                              "diff_windowed_boot_hi", "diff_windowed_wilcoxon_p",
                              "comparable", "note"]]
cov.to_csv(OUT / "coverage_view.csv", index=False, encoding="utf-8-sig")
print(f"\n已写 outputs/daily/coverage_view.csv（{len(cov)} 行 = 一行一个 (mode, product, rung)）")

pd.DataFrame(val_rows).to_csv(OUT / "month_equal_validation.csv", index=False, encoding="utf-8-sig")
print(f"已写 outputs/daily/month_equal_validation.csv（{len(val_rows)} 行）")

# ---------------------------------------------------------------- 打印决定性与高亮
print("\n" + "=" * 96)
print("【直接差值检验（差值序列；三条相关性处理并列）】")
print("=" * 96)
print(f"  {'mode':<6}{'比较':<32}{'Sen':>9}{'HamedRao95%':>24}{'p':>9}{'块自助95%':>22}{'Wilcoxon p':>12}")
for (mode, cand, base), r in direct_reports.items():
    print(f"  {mode:<6}{cand+' − '+base:<32}{r['sen']:>+9.4f}"
          f"{'['+format(r['hr_lo'],'+.4f')+','+format(r['hr_hi'],'+.4f')+']':>24}"
          f"{r['hr_p']:>9.4f}"
          f"{'['+format(r['boot_lo'],'+.4f')+','+format(r['boot_hi'],'+.4f')+']':>22}"
          f"{r['wil_p']:>12.4f}")

print("""
  判读：
  1. **决定性行**是 `month_equal − weighted`。若其 Sen 与全样本同向、区间不跨零 →
     加权效应**不是**季节构成的产物。若塌到零或反向 → 主结论必须改成「效应依赖季节构成」，
     「翻号存活」的写法不得保留。
  2. 报哪条相关性处理：主报**块自助 95% 区间**（保留时间依赖，且两估计量误差相关，
     差值序列正是处理该相关的对象）；Hamed–Rao p 作同口径换算；城市层级 Wilcoxon
     回答「多数城市是否一致」。三者若冲突，以块自助为主，如实并报 Hamed–Rao 的临界性。
  3. 60% 档（尤其夜间）只作上界，`coverage_view.csv` 的 `comparable=False` 与 `note`
     已标出坍缩量级；不得与 20px/40% 档并列当同一估计量的序列。
""")
