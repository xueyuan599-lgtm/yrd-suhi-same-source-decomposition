"""逐日 MOD11A1 的**紧凑中间表** `daily_records_{mode}.csv`（2026-10-03）。三个下游分析共用它。

## 它回答什么

`03_build_daily_panel.py` 只落**城市–年**一层的四口径面板（803 行/模式）。凡是要换
**日期子集、季节、月份、覆盖规则或 QA 档**的敏感性，都得回到 **(城市, 日)** 这一层——
再读一遍 25 年栅格太贵，故先把逐日聚合**固化成一张表**，下游只读它。

本脚本把每个 `(city_id, year, date)` 上、城乡两侧在**两套 QA 判定**下的像元数与
像元均值落表：

| 列组 | 判定 | 来源 |
|------|------|------|
| `u_mean, u_n, r_mean, r_n` | 主口径 `(QC&3) <= 1`（优+次优） | `D*` 波段（导出时已服务端掩膜） |
| `u_mean_strict, u_n_strict, r_mean_strict, r_n_strict` | 收严 `(QC&3) == 0`（只要优） | `Q*` 原始 QC 波段 |

**行 = 至少一侧有 0 个以上有效像元**（主口径）。收严列与主口径列同一行，
故「同一天、同一城」两档可直接并排读。

## 三条产物属性（不是假设，是实测/输入的事实）

1. **`D` 已按 `(QC&3)<=1` 服务端掩膜**，故 `D == -9999` 恰是主口径剔除处。
   收严集是主口径集的**子集**（`(QC&3)==0 ⊂ (QC&3)<=1`），因此 `(Q&3)==0` 的像元上
   `D` 必非 `-9999`。
2. **2020 与 2023 只有 `D` 波段、`Q` 波段为 0**（同批导出的另两年）。这两年的收严列
   **必须全为缺失**：`u_mean_strict/r_mean_strict = NaN`、`u_n_strict/r_n_strict = 0`。
   本脚本**按波段名探测**这一点，不硬编码年份；缺 `Q` 的年份打印并写入
   `daily_records_meta.json`。
3. 城乡掩膜**只读** `masks_matched.npz`（12 的唯一产物），本脚本不重判城乡。

## 口径与 `03` 完全一致（这是本表唯一的意义）

`u_mean` 与 `u_n` 由 `03` 同一套 `bincount` 聚合法逐日算出：

    pixel_weighted = Σ(u_mean·u_n)/Σ u_n − Σ(r_mean·r_n)/Σ r_n
    equal_any      = mean(u_mean | u_n>0) − mean(r_mean | r_n>0)，两侧各需 ≥ 6 天

故本脚本自带**回归验证**：从 `daily_records_{mode}.csv` 重建 `panel_{mode}_daily.csv`
的 `suhii_pixel_weighted` 与 `suhii_equal_any` 两列，逐城–年比对，打印最大绝对差。

## 实测结果（2026-10-03，全量 25 年 × 昼夜，803 城市–年/模式，无缺）

| 模式 | 行数 | 城市 | 年 | max\|Δ pixel_weighted\| | max\|Δ equal_any\| |
|------|------|------|----|--------------------------|---------------------|
| 夜间 | 199 522 | 37 | 2001–2025 | **8.88e-16** | **1.14e-13** |
| 白天 | 209 572 | 37 | 2001–2025 | **4.44e-16** | **1.14e-13** |

两列**全部 803 个城市–年**逐值吻合到浮点累加位级（远优于 1e-9 的验收线），
即本中间表与 `03` 的口径**逐位同源**。`equal_any` 略大于 `pixel_weighted` 是
`mean(·)` 与 `Σ/Σ` 两次不同的浮点运算路径所致，非口径差。

**缺 Q 波段的年份**（收严列全缺失，已写入 `daily_records_meta.json`）：**2020、2023**（昼夜同）。
其余 23 年带 Q：白天收严档有实质数据（如 2019 年城区像元·日 152 156），
夜间极少（2019 年仅 285），与 `tables/qc_prevalence_night.csv` 的约 0.3% 一致。

## 它**不能**用来声称什么

- **不能**把收严档数字当作「主口径偏了多少」：收严集是**非随机子集**（算法更自信的
  那部分），与城乡温差可能相关。本表只支持「换档会不会变」，不支持「主口径的偏差量」。
- **不能**在夜间把收严档当**稳健性证据**：夜间 `(Q&3)==0` 覆盖率极低
  （见 `tables/qc_prevalence_night.csv`），有效像元·日极少，属**数据量本身不足**。
- **不能**从本表直接下趋势结论：趋势须由下游在同口径、共同城市–年份上按 `06`
  的 Theil–Sen + Hamed–Rao 约定重算。

## 内存与读法

导出是 `interleave=pixel`：**逐波段 `read` 约等于解压整幅文件**（2026-10-03 同批实测
0.518 s/波段 vs 1.75 s/整文件 332 波段，约 100×）。故本脚本与 `07`/`08` 同一取舍——
**一次读入一个文件（半年的波段），处理完即 `del`**，绝不把**全年（两文件）**或跨年
波段同时驻留。峰值约 0.5 GB/文件（332×756 675×2 B）。

用法：
    python 09_build_daily_records.py                  # 昼夜都跑
    python 09_build_daily_records.py --mode night
    python 09_build_daily_records.py --mode day --years 2020 2023
"""

import argparse
import json
import sys
from itertools import groupby
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")   # 只改 stdout 时 stderr 仍走 GBK，中文报错会乱码

FILL = -9999
SCALE_K = 0.02
MIN_DAYS = 6                     # 与 03 / 07 同值：等权要求每侧 ≥ 6 个有效日
COLS = ["city_id", "year", "mode", "date", "month",
        "u_mean", "u_n", "r_mean", "r_n",
        "u_mean_strict", "u_n_strict", "r_mean_strict", "r_n_strict"]

parser = argparse.ArgumentParser(description="逐日紧凑中间表 daily_records")
parser.add_argument("--mode", choices=("night", "day", "both"), default="both")
parser.add_argument("--years", type=int, nargs="+", default=None,
                    help="默认=盘上已有的全部年份（按 uhi_{mode}daily_*.tif 探测）")
ARGS = parser.parse_args()

MODES = ("night", "day") if ARGS.mode == "both" else (ARGS.mode,)

BASE = Path(__file__).resolve().parent.parent.parent      # analysis/
DATA = BASE / "data"
EXPORT = DATA / "gee_export"
DAILY = DATA / "daily"
DAILY.mkdir(parents=True, exist_ok=True)

masks = np.load(DATA / "masks_matched.npz")
H, W = int(masks["H"]), int(masks["W"])
NPX = H * W


def cities_of(year):
    """该年**城乡掩膜都在**的城市，升序。与 03 第 206–207 行逐字同构。"""
    out = []
    for k in masks.files:
        if k.startswith(f"u{year}|"):
            c = k.split("|")[1]
            if f"r{year}|{c}" in masks.files:
                out.append(int(c))
    return sorted(out)


def plan_year(mode, year):
    """→ [(文件, D 波段号, Q 波段号或 None, 日期串)]，按日期序。

    Q 波段**按名字探测**（`Q{date}` 是否在**同一文件**的描述里）：2020/2023 的导出
    只有 D，故它们的 `bq` 为 `None`，收严列自然落空——不靠硬编码年份。
    日期跨文件重复即报错（分半段拼接错位会静默把两天并成一天）。
    """
    plan, seen = [], set()
    for f in sorted(EXPORT.glob(f"uhi_{mode}daily_{year}_*.tif")):
        with rasterio.open(f) as s:
            desc = list(s.descriptions)
            idx = {d: i + 1 for i, d in enumerate(desc) if d}
            for d in desc:
                if d and d.startswith("D"):
                    date = d[1:]
                    if date in seen:
                        raise RuntimeError(f"{mode} {year}: 日期 {date} 在两个文件里重复")
                    seen.add(date)
                    plan.append((f, idx[d], idx.get(f"Q{date}"), date))
    plan.sort(key=lambda t: t[3])
    return plan


def build(mode, years):
    """把一个模式全部年份聚合成 `daily_records_{mode}.csv` 的 DataFrame。"""
    rows, year_info, no_qc = [], [], []
    print("=" * 84)
    print(f"逐日紧凑中间表   [mode={mode} · {len(years)} 年：{years[0]}–{years[-1]}]")
    print("=" * 84)

    for year in years:
        plan = plan_year(mode, year)
        if not plan:
            print(f"  {year:<6}缺文件，跳过")
            continue
        if not any(bq is not None for _, _, bq, _ in plan):
            no_qc.append(year)

        cities = cities_of(year)
        nc = len(cities)
        cities_arr = np.array(cities, dtype=np.int64)
        # 「像元 → 城市序号」查找表（0 = 不属于任何城）。与 03/07 同法：一次 `bincount`
        # 同时算全部城，取代逐城 fancy-indexing（后者每波段 2×nc 次整幅拷贝）。
        ci_u = np.zeros(NPX, dtype=np.int32)
        ci_r = np.zeros(NPX, dtype=np.int32)
        for j, c in enumerate(cities, start=1):
            ci_u[masks[f"u{year}|{c}"]] = j
            ci_r[masks[f"r{year}|{c}"]] = j

        n_dates = 0
        # 按文件分组：一次开、整幅读，处理完即 del（峰值 ~0.5 GB/文件，绝不跨年驻留）。
        for f, grp in groupby(plan, key=lambda t: t[0]):
            grp = list(grp)
            with rasterio.open(f) as src:
                arr = src.read()
            for _, bd, bq, date in grp:
                a = arr[bd - 1].astype(np.float32)
                a[a == FILL] = np.nan
                a *= SCALE_K
                flat = a.ravel()
                ok = np.isfinite(flat)
                fv = flat[ok]
                iu, ir = ci_u[ok], ci_r[ok]
                cu = np.bincount(iu, minlength=nc + 1)
                su = np.bincount(iu, weights=fv, minlength=nc + 1)
                cr = np.bincount(ir, minlength=nc + 1)
                sr = np.bincount(ir, weights=fv, minlength=nc + 1)
                if bq is not None:
                    qc = arr[bq - 1].ravel()
                    oks = ok & ((qc & 3) == 0)
                    fvs = flat[oks]
                    ius, irs = ci_u[oks], ci_r[oks]
                    cus = np.bincount(ius, minlength=nc + 1)
                    sus = np.bincount(ius, weights=fvs, minlength=nc + 1)
                    crs = np.bincount(irs, minlength=nc + 1)
                    srs = np.bincount(irs, weights=fvs, minlength=nc + 1)
                else:
                    z = np.zeros(nc + 1, dtype=np.float64)
                    cus = sus = crs = srs = z
                js = np.nonzero((cu[1:] > 0) | (cr[1:] > 0))[0] + 1
                if len(js) == 0:
                    continue
                n_dates += 1
                dint = int(date)
                mon = (dint // 100) % 100
                cid = cities_arr[js - 1]
                um = np.where(cu[js] > 0, su[js] / np.maximum(cu[js], 1), np.nan)
                rm = np.where(cr[js] > 0, sr[js] / np.maximum(cr[js], 1), np.nan)
                ums = np.where(cus[js] > 0, sus[js] / np.maximum(cus[js], 1), np.nan)
                rms = np.where(crs[js] > 0, srs[js] / np.maximum(crs[js], 1), np.nan)
                for k in range(len(js)):
                    j = js[k]
                    rows.append((int(cid[k]), year, mode, dint, mon,
                                 float(um[k]), int(cu[j]), float(rm[k]), int(cr[j]),
                                 float(ums[k]), int(cus[j]), float(rms[k]), int(crs[j])))
            del arr
        year_info.append({"year": year, "n_dates": n_dates, "n_cities": len(cities)})
        print(f"  {year:<6}城 {len(cities):<3}日 {n_dates:<4}"
              + ("  （无 Q 波段 → 收严列全缺失）" if year in no_qc else ""))

    df = pd.DataFrame(rows, columns=COLS)
    df = df.sort_values(["city_id", "date"]).reset_index(drop=True)
    out = DAILY / f"daily_records_{mode}.csv"
    df.to_csv(out, index=False, encoding="utf-8-sig")
    print(f"\n已写 {out.relative_to(BASE)}  {len(df)} 行 · "
          f"{df.city_id.nunique()} 城 · {int(df.year.min())}–{int(df.year.max())}")
    return df, year_info, no_qc


def reproduce(mode):
    """从中间表重建 `panel_{mode}_daily.csv` 的两列并逐城–年比对最大绝对差。"""
    rec = pd.read_csv(DAILY / f"daily_records_{mode}.csv")
    pan = pd.read_csv(DAILY / f"panel_{mode}_daily.csv")
    dif_pw, dif_eq, missing = [], [], 0
    for r in pan.itertuples():
        sub = rec[(rec.city_id == r.city_id) & (rec.year == r.year)]
        if len(sub) == 0:
            missing += 1
            continue
        # u_mean 仅在 u_n==0 时为 NaN；此时其对和的贡献应为 0，故乘积后 fillna(0)。
        pw = ((sub.u_mean * sub.u_n).fillna(0).sum() / sub.u_n.sum()
              - (sub.r_mean * sub.r_n).fillna(0).sum() / sub.r_n.sum())
        um = sub.loc[sub.u_n > 0, "u_mean"]
        rm = sub.loc[sub.r_n > 0, "r_mean"]
        if len(um) >= MIN_DAYS and len(rm) >= MIN_DAYS:
            eq = um.mean() - rm.mean()
            dif_eq.append(abs(eq - r.suhii_equal_any))
        dif_pw.append(abs(pw - r.suhii_pixel_weighted))
    d_pw = max(dif_pw) if dif_pw else np.nan
    d_eq = max(dif_eq) if dif_eq else np.nan
    print(f"  [{mode}] 城市–年 {len(pan)} · 无对应日记录的 {missing} · "
          f"max|Δpixel_weighted| = {d_pw:.3e} · max|Δequal_any| = {d_eq:.3e}")
    return {"mode": mode, "city_years": int(len(pan)), "missing": missing,
            "max_abs_diff_pixel_weighted": float(d_pw),
            "max_abs_diff_equal_any": float(d_eq)}


meta = {"modes": {}}
print("=" * 84)
print("构建")
print("=" * 84)
for mode in MODES:
    if ARGS.years is None:
        years = sorted({int(p.name.split("_")[2]) for p in EXPORT.glob(f"uhi_{mode}daily_*.tif")})
        if not years:
            print(f"盘上无 uhi_{mode}daily_*.tif，跳过", file=sys.stderr)
            continue
    else:
        years = ARGS.years
    df, yinfo, no_qc = build(mode, years)
    meta["modes"][mode] = {
        "rows": int(len(df)), "years": years,
        "years_without_qc": no_qc,
        "n_dates_per_year": yinfo,
    }

print("\n" + "=" * 84)
print("自验证：从中间表重建 panel_{mode}_daily.csv 的两列")
print("=" * 84)
for mode in MODES:
    if mode in meta["modes"]:
        meta["modes"][mode]["validation"] = reproduce(mode)

meta_out = DAILY / "daily_records_meta.json"
meta_out.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
print(f"\n已写 {meta_out.relative_to(BASE)}")

print("""
  判读：
  1. 三处下游共用本表：逐月/季节等权、覆盖门槛与有效日留存、QA 收严档——
     都只读它，不再碰栅格。改动本表列义须同步这三个下游。
  2. `u_n / r_n` 是**当日该侧有效像元数**（主口径），是覆盖门槛的分子；
     `u_mean / r_mean` 是当日该侧像元均值（K）。整年加权 = Σ(u_mean·u_n)/Σ u_n。
  3. 收严列在 **2020 / 2023 全为缺失**（那两年导出无 Q 波段），已在 meta 记录。
     夜间收严档覆盖率极低，**不得**当作稳健性证据引用。
""")
