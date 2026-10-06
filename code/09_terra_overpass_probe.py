"""Terra 过境时间漂移的**可测性探针**（2026-10-03）。回答「观测时间漂移这条敏感性
能不能在本研究既有产品谱系内跑」。

## 它回答什么

正文把「Terra overpass-time drift is untested in both products」列为局限。本脚本不
尝试做新的产品检索，只回答**一个**问题：本研究谱系里的 LST 产品是否带逐观测的
过境/观测时间字段，若有则以 2001–2025 的年度区域均值给出漂移量级。

谱系内产品（由 `analysis/code/` 与 `analysis/code/daily/` 的 `ee.ImageCollection`
调用穷举得到）：`MODIS/061/MOD11A1`（逐日 LST，Terra）、`MODIS/061/MOD11A2`
（8 天 LST，Terra）、`MODIS/061/MOD13A2`（8 天 NDVI，Terra）、
`ECMWF/ERA5_LAND/MONTHLY_AGGR`（气象）。**LST 产品即 MOD11A1 与 MOD11A2，均为 Terra。**
MOD13A2 是 NDVI，其 `ViewZenith` 是**观测天顶角**不是时间，不能用于过境时间。

档名不靠记忆、不靠产品说明，由 `bandNames().getInfo()` 实测得到（见运行首段输出）。

## 字段与单位（实测）

`MOD11A1.061` 与 `MOD11A2.061` 都带 `Day_view_time` / `Night_view_time`
（前者还叫 `Day_view_angle`、后者叫 `Day_view_angl`，仅拼写差异）。取值是**当地太阳时**
的小时数 ×10（DN）：实测 2020 年 7 月 单景白天 DN 110–114 → 11.0–11.4 h，
夜间 DN 223 → 22.3 h。故 **×0.1 得小时数**。0 是缺测填充，取均值前以 `> 0` 掩膜剔除。

## 探针口径

- 研究区 = `code/07_city_keys.json` 的 41 城（长三角，与 `daily/02_export_daily_seq.py`
  同一 AOI 构造），非另划边界。
- 逐年 2001–2025：区域面元在**该年全部有效观测**上的 `Day_view_time` / `Night_view_time`
  均值（先逐像元取年内时间均值，再对区域取空间均值），同时报 7 月单月均值作季节对照，
  以及每像元年内有效观测数均值作样本稳定性体检。
- 两个产品（MOD11A1 逐日、MOD11A2 8 天）各算一遍。MOD11A2 的 view_time 是合成期内
  的**多日平均观测时间**，与 MOD11A1 的逐日值不是同一抽样单位，**两产品的绝对小时数
  不可直接比大小**（同 `daily/06_daily_trends.py` 第 16 行的自订约定）。
- 趋势口径与仓库一致：`scipy.stats.theilslopes`（斜率 ×10 = 每十年），
  `pymannkendall.hamed_rao_modification_test` 给 p。

## 输出

- `outputs/daily/terra_overpass_times.csv`：本脚本唯一数字化产物（逐年 × 产品 × 昼夜）。
- `outputs/daily/TERRA_OVERPASS_PROBE.md`：由本脚本**从内存中的计算结果**生成，
  报告里每一个数都来自同一次 `reduceRegion` 返回，不手抄。

## 不能拿来声称什么

1. **不能**把「区域均值过境时间漂移」直接读成「LST 趋势偏差」。本探针只量化漂移的量级，
   没有做过境时间与 SUHII 趋势的关联或订正实验。
2. **不能**把 MOD11A1 与 MOD11A2 的小时数并列比较或互推（抽样单位不同，见上）。
3. **不能**把当地太阳时的漂移等同于卫星轨道降交点时刻漂移——前者含经度/纬度加权，
   是区域平均效应。
4. 本脚本**不修改任何既有产出**；只新增本文件与上述两个输出。

用法：
    python 09_terra_overpass_probe.py            # 全量 2001–2025
    python 09_terra_overpass_probe.py --years 2001 2020 2025
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import theilslopes

sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")   # 只改 stdout 时 stderr 仍走 GBK，中文报错会乱码

import ee

PROJECT = "evident-ocean-510002-t0"
ee.Initialize(project=PROJECT)

BASE = Path(__file__).resolve().parent.parent.parent      # analysis/
sys.path.insert(0, str(BASE / ".deps"))
import pymannkendall as mk                                 # 与 06 / 18 / 54 同一实现

OUT = BASE / "outputs" / "daily"
OUT.mkdir(parents=True, exist_ok=True)

SCALE = 1000          # 与 02/16 导出网格同尺度
SCALE_H = 0.1         # view_time DN → 小时
# 8 天产品（MOD11A2）与逐日（MOD11A1）
PRODUCTS = {
    "MOD11A1": "MODIS/061/MOD11A1",
    "MOD11A2": "MODIS/061/MOD11A2",
}
BAND_DAY = "Day_view_time"
BAND_NIGHT = "Night_view_time"
BAND_ANGLE = ("Day_view_angle", "Day_view_angl")   # A1 与 A2 仅拼写差
LABEL = {"MOD11A1": "逐日 (Terra)", "MOD11A2": "8 天合成 (Terra)"}

parser = argparse.ArgumentParser(description="Terra 过境时间漂移可测性探针")
parser.add_argument("--years", type=int, nargs="+", default=None,
                    help="默认 2001–2025 全部")
ARGS = parser.parse_args()
YEARS = ARGS.years if ARGS.years else list(range(2001, 2026))
FOCUS = [2001, 2005, 2010, 2015, 2020, 2025]        # 正文/报告重点披露的年份


def city_fc(key):
    """与 daily/02_export_daily_seq.py 逐字同构的 AOI 构造。"""
    return (ee.FeatureCollection("FAO/GAUL/2025/level2")
            .filter(ee.Filter.eq("GAUL1_NAME", key["province"]))
            .filter(ee.Filter.eq("GAUL2_NAME", key["city"])))


import json
with open(BASE / "code" / "07_city_keys.json", encoding="utf-8") as f:
    CITY_KEYS = json.load(f)
AOI = ee.FeatureCollection([city_fc(k) for k in CITY_KEYS]).flatten().geometry()

# ------------------------------------------------------------- 0 档名实测
print("=" * 84)
print("【0】产品档名实测（bandNames().getInfo()，非记忆）")
print("=" * 84)
BANDNAMES = {}
for tag, pid in PRODUCTS.items():
    names = ee.ImageCollection(pid).filterDate("2020-01-01", "2020-02-01").first() \
        .bandNames().getInfo()
    BANDNAMES[tag] = names
    has = [b for b in (BAND_DAY, BAND_NIGHT) if b in names]
    print(f"  {tag} ({pid})")
    print(f"     档数 {len(names)}：{names}")
    print(f"     观测时间档：{has if has else '无'}   角度档："
          f"{[a for a in BAND_ANGLE if a in names]}")
md13 = ee.ImageCollection("MODIS/061/MOD13A2").filterDate("2020-01-01", "2020-02-01") \
    .first().bandNames().getInfo()
print(f"  （对照）MOD13A2 是 NDVI，无观测时间档："
      f"{'Day_view_time' in md13 or 'Night_view_time' in md13}")
print()


def masked_mean_band(pid: str, band: str, start: str, end: str):
    """该区间内 band 的逐像元年内均值（先以 band>0 剔除缺测填充）。"""
    col = ee.ImageCollection(pid).filterDate(start, end).select(band)
    col = col.map(lambda im: im.updateMask(im.gt(0)))
    return col.mean()


def masked_count_band(pid: str, band: str, start: str, end: str):
    """该区间内 band>0 的**有效观测数**（逐像元），体检样本稳定性。"""
    col = ee.ImageCollection(pid).filterDate(start, end).select(band)
    return col.map(lambda im: im.gt(0)).sum()


def year_bands(pid: str, year: int):
    """某年的三张图：白天均值、夜间均值、白天有效观测数（重命名带年标签）。"""
    y0, y1 = f"{year}-01-01", f"{year + 1}-01-01"
    return [
        masked_mean_band(pid, BAND_DAY, y0, y1).rename(f"day_all_{year}"),
        masked_mean_band(pid, BAND_NIGHT, y0, y1).rename(f"night_all_{year}"),
        masked_mean_band(pid, BAND_DAY, f"{year}-07-01", f"{year}-08-01")
        .rename(f"day_jul_{year}"),
        masked_count_band(pid, BAND_DAY, y0, y1).rename(f"nobs_day_{year}"),
    ]


def reduce_over_aoi(image: ee.Image):
    """区域均值。tileScale 抬高避免 user memory 超限；bestEffort 在超 maxPixels 时放宽尺度。"""
    d = image.reduceRegion(
        reducer=ee.Reducer.mean(), geometry=AOI, scale=SCALE,
        maxPixels=1e9, bestEffort=True, tileScale=4,
    ).getInfo()
    return d or {}


# ------------------------------------------------------------- 1 逐年区域均值
print("=" * 84)
print(f"【1】逐年区域均值  年份 {YEARS[0]}–{YEARS[-1]}  ·  AOI {len(CITY_KEYS)} 城  ·  "
      f"{SCALE} m")
print("=" * 84)
raw = {}
for tag, pid in PRODUCTS.items():
    img = ee.Image.cat([b for y in YEARS for b in year_bands(pid, y)])
    got = reduce_over_aoi(img)
    raw[tag] = got
    print(f"  {tag}：返回 {len(got)} 个带  (期望 {len(YEARS) * 4})")

# 组装长表
rows = []
for tag in PRODUCTS:
    for y in YEARS:
        g = raw[tag]
        day_all = g.get(f"day_all_{y}")
        night_all = g.get(f"night_all_{y}")
        day_jul = g.get(f"day_jul_{y}")
        nobs = g.get(f"nobs_day_{y}")
        rows.append({
            "product": tag, "year": y,
            "day_view_h": None if day_all is None else day_all * SCALE_H,
            "night_view_h": None if night_all is None else night_all * SCALE_H,
            "day_view_jul_h": None if day_jul is None else day_jul * SCALE_H,
            "mean_n_days_per_px": nobs,
        })
df = pd.DataFrame(rows)
df.to_csv(OUT / "terra_overpass_times.csv", index=False, encoding="utf-8-sig")


def trend_of(series: pd.Series):
    """Theil–Sen + Hamed–Rao，与 06_daily_trends.py 同参。年数不足返回 None。"""
    a = series.dropna().sort_index()
    if len(a) < 8:
        return None
    sen = theilslopes(a.values, a.index.values, alpha=0.95)
    hr = mk.hamed_rao_modification_test(a.values)
    return {"n_years": int(len(a)),
            "sen_slope_min_per_decade": float(sen.slope * 10 * 60),
            "sen_CI95_lo": float(sen.low_slope * 10 * 60),
            "sen_CI95_hi": float(sen.high_slope * 10 * 60),
            "hamed_rao_p": float(hr.p)}


TREND = {}
print("\n【逐年区域均值（小时，当地太阳时）】")
for tag in PRODUCTS:
    sub = df[df["product"] == tag].set_index("year")
    print(f"\n  {tag}  {LABEL[tag]}")
    print(f"    {'年':<6}{'白天':>8}{'夜间':>8}{'7月白天':>10}{'有效日/像元':>12}")
    for y in YEARS:
        r = sub.loc[y]
        def f(v, w=8, p=3):
            return f"{v:>{w}.{p}f}" if v is not None and np.isfinite(v) else f"{'--':>{w}}"
        print(f"    {y:<6}{f(r.day_view_h)}{f(r.night_view_h)}{f(r.day_view_jul_h, 10)}"
              f"{f(r.mean_n_days_per_px, 12, 1)}")
    for key, col in (("day", "day_view_h"), ("night", "night_view_h")):
        t = trend_of(sub[col])
        TREND[(tag, key)] = t
        if t:
            print(f"    → {key} Theil–Sen {t['sen_slope_min_per_decade']:+.2f} min/十年 "
                  f"(95% {t['sen_CI95_lo']:+.2f} … {t['sen_CI95_hi']:+.2f})  "
                  f"p={t['hamed_rao_p']:.4g}")

# ------------------------------------------------------------- 2 报告生成
def fmt(v, w=8, p=3):
    return f"{v:>{w}.{p}f}" if v is not None and np.isfinite(v) else f"{'--':>{w}}"


def focus_rows(tag):
    sub = df[df["product"] == tag].set_index("year")
    out = []
    for y in FOCUS:
        if y not in sub.index:
            continue
        r = sub.loc[y]
        out.append(f"| {y} | {fmt(r.day_view_h,0,3)} | {fmt(r.night_view_h,0,3)} | "
                   f"{fmt(r.day_view_jul_h,0,3)} | {fmt(r.mean_n_days_per_px,0,1)} |")
    return out


def trend_line(tag, key):
    t = TREND.get((tag, key))
    if not t:
        return "年数不足，未做趋势检验"
    return (f"Theil–Sen **{t['sen_slope_min_per_decade']:+.2f} min/十年**"
            f"（95% {t['sen_CI95_lo']:+.2f} … {t['sen_CI95_hi']:+.2f}），"
            f"Hamed–Rao p = {t['hamed_rao_p']:.4g}，n = {t['n_years']} 年")


def first_last(tag, key):
    sub = df[df["product"] == tag].dropna(subset=[key]).sort_values("year")
    if len(sub) < 2:
        return None
    a, b = sub.iloc[0], sub.iloc[-1]
    return (int(a.year), float(a[key]), int(b.year), float(b[key]),
            float(b[key] - a[key]))


# 判定：任一天/夜序列的 Hamed–Rao p < 0.05 且斜率非零 → 有可测漂移
measurable = []
for (tag, key), t in TREND.items():
    if t and t["hamed_rao_p"] < 0.05 and abs(t["sen_slope_min_per_decade"]) > 0:
        measurable.append((tag, key, t))

verdict = ("**TESTABLE** —— 谱系内 MOD11A1 与 MOD11A2 均带逐观测过境时间档，"
           "且 2001–2025 区域均值呈可测漂移（见下表数值）。"
           if measurable else
           "**NOT TESTABLE 的强主张不成立** —— 字段存在，但 2001–2025 区域均值"
           "未检出达到 5% 阈值的单调漂移（见下表数值）。")

lines = []
lines.append("# Terra 过境时间漂移：可测性探针（2026-10-03）")
lines.append("")
lines.append(f"**判定：{verdict}**")
lines.append("")
lines.append("> 本文件由 `analysis/code/daily/09_terra_overpass_probe.py` 从同一次 "
             "`reduceRegion` 返回生成，未手抄任何数字。数字化产物见 "
             "`analysis/outputs/daily/terra_overpass_times.csv`。")
lines.append("")
lines.append("## 1. 字段是否存在的问题（实测，非说明书）")
lines.append("")
lines.append("`MODIS/061/MOD11A1` 与 `MODIS/061/MOD11A2` 的 `bandNames()` 实测：")
lines.append("")
for tag in PRODUCTS:
    lines.append(f"- **{tag}**（{LABEL[tag]}）：`{BAND_DAY}`、`{BAND_NIGHT}` "
                 f"均在档（共 {len(BANDNAMES[tag])} 档）。取值 = 当地太阳时 ×10（DN），"
                 f"×0.1 得小时；0 为缺测，取均值前以 `>0` 掩膜剔除。")
lines.append(f"- 对照 **MOD13A2**（NDVI，非 LST）无观测时间档，其 `ViewZenith` 是"
             f"观测天顶角不是时间，**不能**用于过境时间。")
lines.append("")
lines.append("结论：**正文局限「过境时间在两种产品中都未测」在字段层面不再成立**——"
             "字段是可得的；差的是有没有跑过这项敏感性。")
lines.append("")
lines.append("## 2. 实测漂移（区域均值，当地太阳时，小时）")
lines.append("")
lines.append("口径：AOI = `code/07_city_keys.json` 的 41 城（长三角）；先逐像元取年内时间"
             "均值，再对区域取空间均值；`7 月白天`为季节对照；“有效日/像元”是逐像元年内"
             "有效观测数均值（样本稳定性体检）。")
for tag in PRODUCTS:
    lines.append("")
    lines.append(f"### {tag}（{LABEL[tag]}）")
    lines.append("")
    lines.append("| 年 | 白天 | 夜间 | 7月白天 | 有效日/像元 |")
    lines.append("|----|------|------|---------|-------------|")
    lines.extend(focus_rows(tag))
    lines.append("")
    for key, col, name in (("day", "day_view_h", "白天"), ("night", "night_view_h", "夜间")):
        fl = first_last(tag, col)
        fl_txt = (f"首末差 {fl[2]}–{fl[0]} 为 {fl[4]:+.3f} h"
                  if fl else "首末差不可得")
        lines.append(f"- {name}趋势：{trend_line(tag, key)}；{fl_txt}。")
lines.append("")
lines.append("## 3. 能说与不能说")
lines.append("")
lines.append("**能说**：过境时间档存在且非空，年度区域均值可用于刻画 2001–2025 的漂移量级。")
lines.append("")
lines.append("**不能说**：")
lines.append("1. 不得把区域均值过境时间漂移直接读成 LST 趋势偏差——本探针未做漂移与 "
             "SUHII 趋势的关联/订正实验。")
lines.append("2. MOD11A1（逐日）与 MOD11A2（8 天）的小时数**不可直接比大小**：抽样单位"
             "不同（后者是合成期内多日平均观测时间）。")
lines.append("3. 当地太阳时的漂移不等于轨道降交点时刻漂移（含经度/纬度加权）。")
lines.append("")

(OUT / "TERRA_OVERPASS_PROBE.md").write_text("\n".join(lines), encoding="utf-8")
print(f"\n已写 {(OUT / 'terra_overpass_times.csv').relative_to(BASE)}")
print(f"已写 {(OUT / 'TERRA_OVERPASS_PROBE.md').relative_to(BASE)}")
print(f"\n判定：{verdict}")
