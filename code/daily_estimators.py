"""逐日 SUHII 估计口径的**共享库**（2026-10-03）。`12_qc_strict_trend.py` 的唯一算术来源。

## 它是什么

从既有脚本逐字抽出的估计量定义，供下游复用，使**同一阶梯**能跑在不同的像元支撑上：

| 成分 | 逐字来源 |
|------|---------|
| `SCHEMES`（20px / 40% / 60% 阈值段）、`MIN_DAYS = 6` | `03_build_daily_panel.py` |
| 加权 / 等权 / 配对 的公式 | `03_build_daily_panel.py` 主循环 |
| `LADDER` 的 rung 名与顺序 | `06_daily_trends.py` 的 `daily_map` |
| `trend_of`（Theil–Sen ×10 + Hamed–Rao）、`mask_pixel_counts` | `10_month_equal_and_coverage.py` |

## 为什么单独一个模块

`12` 要在**两套像元支撑**（主口径 `lte1` 与收严 `eq0`）上跑**同一阶梯**。若把公式在 `12`
里再写一遍，收严档与主口径就可能悄悄漂移，且日后改阈值要改两处。故把阶梯抽成本模块，
`12` 只调 `city_year_ladder(d, px, suf)` 并传不同的支撑后缀（`""` = 主口径、`"_strict"` = 收严）。

## 它**不**解决的事

- `03` / `06` / `10` **仍各自内联自己的副本**（本轮硬约束不允许改既有脚本）。本模块是
  **新增的单一来源**，尚未回灌那三个脚本——存在多份副本漂移的风险。今后改动任一公式，
  须四处一起改，或把三个脚本改为 `import daily_estimators`。
- 不做 QA 判定、不读栅格；只对**已聚合成逐日记录**（`daily_records_{mode}.csv`）的表做算术。
"""

from pathlib import Path
import sys

import numpy as np
from scipy.stats import theilslopes

BASE = Path(__file__).resolve().parent.parent.parent       # analysis/
sys.path.insert(0, str(BASE / ".deps"))
import pymannkendall as mk                                  # 与 06 / 10 同一实现

DATA = BASE / "data"

MIN_DAYS = 6            # 等权要求每侧 ≥ 6 个有效日；与 03 / 06 / 07 / 09 / 10 同值
MIN_YEARS = 8           # 趋势最少年数；与 06 / 10 同值

# 与 03_build_daily_panel.py 的 SCHEMES 逐字相同：档位名 → (绝对像元下限, 比例下限)
SCHEMES = [("20px", 20, 0.0), ("40%", 10, 0.40), ("60%", 10, 0.60)]
# 档位名 → 06_daily_trends.py 的 rung 名
RUNG = {"20px": "equal_20px", "40%": "equal_40pct", "60%": "equal_60pct"}
# 与 06 的 daily_map 逐字相同的 rung 名与顺序（配对只做 20px 一档，与 03 一致）
LADDER = ["weighted", "equal_any", "equal_20px", "equal_40pct", "equal_60pct", "paired_20px"]
DESC = {"weighted": "像元全池化（加权）", "equal_any": "等权（任一侧有像元）",
        "equal_20px": "等权 + 20 像元下限", "equal_40pct": "等权 + 40% 下限",
        "equal_60pct": "等权 + 60% 下限", "paired_20px": "同日配对 + 等权（20 像元）"}


def trend_of(annual, min_years=MIN_YEARS):
    """Theil–Sen（×10 成每十年）+ Hamed–Rao，与 06 / 10 逐字同参。年数不足返回 None。"""
    a = annual.dropna().sort_index()
    if len(a) < min_years:
        return None
    sen = theilslopes(a.values, a.index.values, alpha=0.95)
    hr = mk.hamed_rao_modification_test(a.values)
    return {"n_years": int(len(a)),
            "mean_K": float(a.mean()),
            "sen_K_per_decade": float(sen.slope * 10),
            "sen_CI95_lo": float(sen.low_slope * 10),
            "sen_CI95_hi": float(sen.high_slope * 10),
            "hamed_rao_p": float(hr.p)}


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


def city_year_ladder(d, px, suf=""):
    """一个城市–年在**指定像元支撑**上的整条阶梯。`d` = 该城–年的逐日记录。

    `suf` 选支撑：`""` 用 `u_mean/u_n/r_mean/r_n`（主口径 lte1），`"_strict"` 用
    `u_mean_strict/u_n_strict/r_mean_strict/r_n_strict`（收严 eq0）。`px` 是该城–年的
    `(城区掩膜像元数, 乡村掩膜像元数)`，即 40%/60% 比例门槛的分母；缺则相关档落空。

    返回 {rung: 值, f"n_days_{rung}": 有效日数}。口径与 03 逐字一致：
    加权 = Σ(mean·n)/Σn 两侧相减；等权 = 两侧各取达阈日的均值再相减（不要求同日）；
    配对 = 只取**同日**两侧都达阈的行，先作差再取均值。
    """
    um, un = f"u_mean{suf}", f"u_n{suf}"
    rm, rn = f"r_mean{suf}", f"r_n{suf}"
    res = {}
    su, sr = float(d[un].sum()), float(d[rn].sum())
    res["weighted"] = (float((d[um] * d[un]).fillna(0).sum()) / su
                       - float((d[rm] * d[rn]).fillna(0).sum()) / sr
                       if su > 0 and sr > 0 else np.nan)

    ua, ra = d.loc[d[un] > 0, um], d.loc[d[rn] > 0, rm]
    res["equal_any"] = (float(ua.mean()) - float(ra.mean())
                        if len(ua) >= MIN_DAYS and len(ra) >= MIN_DAYS else np.nan)
    res["n_days_equal_any"] = int(min(len(ua), len(ra)))

    nu, nr = px if (px and np.isfinite(px[0]) and np.isfinite(px[1])) else (np.nan, np.nan)
    for slab, minpx, frac in SCHEMES:
        tag = RUNG[slab]
        if not np.isfinite(nu):
            res[tag] = np.nan
            res[f"n_days_{tag}"] = 0
            if tag == "equal_20px":
                res["paired_20px"] = np.nan
                res["n_days_paired_20px"] = 0
            continue
        thu, thr = max(minpx, frac * nu), max(minpx, frac * nr)
        uo, ro = d.loc[d[un] >= thu, um], d.loc[d[rn] >= thr, rm]
        res[tag] = (float(uo.mean()) - float(ro.mean())
                    if len(uo) >= MIN_DAYS and len(ro) >= MIN_DAYS else np.nan)
        res[f"n_days_{tag}"] = int(min(len(uo), len(ro)))
        if tag == "equal_20px":
            pair = d.loc[(d[un] >= thu) & (d[rn] >= thr), [um, rm]].dropna()
            res["paired_20px"] = (float((pair[um] - pair[rm]).mean())
                                  if len(pair) >= MIN_DAYS else np.nan)
            res["n_days_paired_20px"] = int(len(pair))
    return res
