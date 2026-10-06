"""补导：**逐日** MOD11A1 LST 序列（2026-10-02）。逐日 pilot 的数据源。

与 `16_export_nightseq.py` / `16b_export_dayseq.py` 逐字同构，**只有三处不同**：

| 项 | 16/16b（合成期） | 本脚本（逐日） |
|----|-----------------|---------------|
| 集合 | `MODIS/061/MOD11A2` → 46 期/年 | `MODIS/061/MOD11A1` → 365/366 天/年 |
| 年份 | 2001–2025 | **仅 2020 / 2023**（pilot） |
| 输出名 | `uhi_{night,day}seq_{年}.tif` | `uhi_{mode}daily_{年}.tif` |

region / crs / scale / QC 口径 / FILL / DN×0.02 **全部与 16 逐字相同** ——
网格必须重合，否则本地叠加静默错位。

## 为什么不走服务器端 zonal reduce

原计划是「服务器端归约 → 只导城市–日结果表（约 6 万行）」，但那条路需要把
`12_build_panel.py` 的城乡掩膜先建成 EE asset，而 asset 上传在本机**被阻断**
（无 GCS 桶、CLI 只收 gs://，详见 `01_build_mask_asset.py` 抬头）。

改走栅格路线的理由：**掩膜零漂移**——本地直接读 `masks_matched.npz`，与 18/19 同一份
掩膜、同一条代码路径，不引入第二份实现。

**体量已实测**（2026-10-02，非估算）：2020 = 210.9 MB、2023 = 243.6 MB，两年合计
**454.5 MB**（day 332.7 + night 122.8）。按年均 227 MB 外推，铺满 2001–2025 约
**5.7 GB**——而非先前的 28 GB 估算（该估算高了约 5 倍，据此得出的「此路不可行」结论
**作废**）。即栅格路线本身可以支撑全量；asset 问题从「必须解决」降为「可选优化」。

## 为什么按「半年两段」切

366 个波段塞进一个 `Export.image.toDrive` 是 16 脚本的 8 倍，任务重、超时风险高。
切成上/下半年各约 183 波段，本地再沿波段轴拼接。拼接处的日期由波段描述
（`D{YYYYMMDD}`）决定，不依赖顺序假设。

## QC 口径（**分析决策，已于 10-02 定档并写入状态记录**）

`bitwiseAnd(3).lte(1)`（优+次优），昼夜同口径，与 08/16/16b 逐字一致。
**不收严**：白天虽有「优」档（QC 位 00）而夜间没有，但收严会让逐日新序列与
已发表面板不可比；要收须两套一起换，不在本 pilot 内。

用法：
    python 02_export_daily_seq.py --mode night --dry-run
    python 02_export_daily_seq.py --mode night
    python 02_export_daily_seq.py --mode day
"""

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

import ee

sys.stdout.reconfigure(encoding="utf-8")

PROJECT = "evident-ocean-510002-t0"
ee.Initialize(project=PROJECT)

parser = argparse.ArgumentParser(description="逐日 MOD11A1 LST 序列导出")
parser.add_argument("--mode", choices=("night", "day"), default="night")
parser.add_argument("--dry-run", action="store_true", help="只预检，不提交")
parser.add_argument("--years", type=int, nargs="+", default=[2020, 2023])
# 带上原始 QC 波段：任务数不变，但本地可复算**任何** QC 档位（含收严）。
# 代价是每景多一条 int16 波段 → 数据量与耗时可能近翻倍，故先用 A/B 探针实测。
parser.add_argument("--with-qc", action="store_true",
                    help="每日多导一条原始 QC 波段（band 名 Q{YYYYMMDD}）")
parser.add_argument("--prefix", default=None,
                    help="覆盖输出文件名前缀（探针用，避免混进主流程的 glob）")
ARGS = parser.parse_args()

MODE = ARGS.mode
TAG = MODE
# 昼夜差异只在这三处，与 16/16b 的约定一致
QC_BAND = "QC_Night" if MODE == "night" else "QC_Day"
LST_BAND = "LST_Night_1km" if MODE == "night" else "LST_Day_1km"
LABEL = "夜间" if MODE == "night" else "白天"

SCALE = 1000
FOLDER = "GEE_uhi_export"
TRACK = f"02_export_tasks_{TAG}.json"   # 探针版本在下方按 SUF 改名
FILL = -9999
# 半年一段：波段数与任务重量都减半。
# ⚠ `filterDate` 的**结束日是排他的**。原先写 `("01-01","06-30")` / `("07-01","12-31")`
# 会静默丢掉 6/30 与 12/31 两天（实测 2020 只取到 364 天、2023 只取到 363 天，
# 一年正好差 2 天）。结束日一律写**下一天**，跨年段写次年 1 月 1 日。
HALVES = ("h1", "h2")


def chunk_range(year: int, half: str):
    """返回该半年段的 (起, 止-排他)。"""
    if half == "h1":
        return f"{year}-01-01", f"{year}-07-01"
    return f"{year}-07-01", f"{year + 1}-01-01"

BASE = Path(__file__).resolve().parent.parent.parent      # analysis/
with open(BASE / "code" / "07_city_keys.json", encoding="utf-8") as f:
    CITY_KEYS = json.load(f)

# 探针（--prefix）走独立清单与独立任务表，不污染主流程的 glob 与监视器
SUF = "_probe" if ARGS.prefix else ""
if SUF:
    TRACK = f"02_export_tasks_{TAG}{SUF}.json"
MANIFEST = BASE / "data" / "daily" / f"daily_dates_{TAG}{SUF}.json"
MANIFEST.parent.mkdir(parents=True, exist_ok=True)


def city_fc(key):
    return (
        ee.FeatureCollection("FAO/GAUL/2025/level2")
        .filter(ee.Filter.eq("GAUL1_NAME", key["province"]))
        .filter(ee.Filter.eq("GAUL2_NAME", key["city"]))
    )


aoi = ee.FeatureCollection([city_fc(k) for k in CITY_KEYS]).flatten().geometry()


def ymd(ms: int) -> str:
    """毫秒时间戳 → YYYYMMDD（UTC）。与 16/16b 逐字相同。"""
    return dt.datetime.fromtimestamp(ms / 1000, dt.timezone.utc).strftime("%Y%m%d")


def daily_seq(year: int, start: str, end: str):
    """该区间的逐日影像与日期（start / end 是**完整**的 YYYY-MM-DD，end 排他）。

    ⚠ **必须 sort 后再取**：ImageCollection 默认顺序不保证是时间序，
    不排序会让波段与日期错位，而且**不会报错**——配对结果全错但不自知（16 脚本同款告诫）。
    """
    col = (ee.ImageCollection("MODIS/061/MOD11A1")
           .filterDate(start, end)
           .sort("system:time_start"))

    times = col.aggregate_array("system:time_start").getInfo()
    dates = [ymd(t) for t in times]
    if len(set(dates)) != len(dates):
        raise RuntimeError(f"{year} {start}~{end} 日期有重复：{len(dates)} 天 / "
                           f"{len(set(dates))} 个不同日期")

    lst = col.toList(col.size())

    def one(i):
        img = ee.Image(lst.get(i))
        valid = img.select(QC_BAND).bitwiseAnd(3).lte(1)   # 与 08/16/16b 逐字相同
        bands = [(img.select(LST_BAND)
                  .updateMask(valid)
                  .unmask(FILL)
                  .toInt16()
                  .rename(f"D{dates[i]}"))]
        if ARGS.with_qc:
            # 原始 QC（未加掩膜）另存一景，供本地复算收严档。
            # 与 LST 同为 int16（EE 多波段导出要求同类型），QC 取值 0–255 无损。
            bands.append(img.select(QC_BAND).unmask(FILL).toInt16().rename(f"Q{dates[i]}"))
        return bands[0] if len(bands) == 1 else ee.Image.cat(bands)

    return ee.Image.cat([one(i) for i in range(len(dates))]), dates


# ------------------------------------------------------------- 预检
print("=" * 78)
print(f"{LABEL}逐日 MOD11A1 序列导出预检   [mode={MODE}]")
print("=" * 78)
print(f"  QC {QC_BAND} bitwiseAnd(3).lte(1)  ·  存原始 DN(int16)，本地 ×0.02 得 K")
print(f"  年份 {len(ARGS.years)} 个：{ARGS.years[0]}–{ARGS.years[-1]}"
      f"  ·  分 {len(HALVES)} 段/年 → {len(ARGS.years) * len(HALVES)} 个任务")
print(f"  网格 {SCALE} m 同 08/16"
      + (f"  ·  **附带原始 QC 波段**（每景波段数 ×2）" if ARGS.with_qc else "")
      + (f"  ·  前缀 {ARGS.prefix}" if ARGS.prefix else "") + "\n")

built, alldates = {}, {}
for y in ARGS.years:
    for half in HALVES:
        start, end = chunk_range(y, half)
        img, dates = daily_seq(y, start, end)
        built[(y, half)] = img
        alldates[(y, half)] = dates
        flag = "   ⚠ 天数异常" if not (170 <= len(dates) <= 185) else ""
        print(f"  {y} {half}  {len(dates):>3} 天   {dates[0]} … {dates[-1]}{flag}")
    n_year = sum(len(v) for (yy, _), v in alldates.items() if yy == y)
    expect = 366 if y % 4 == 0 and (y % 100 != 0 or y % 400 == 0) else 365
    print(f"  → {y} 全年 {n_year} 天（平/闰年应有 {expect}）"
          + ("" if n_year == expect else "   ⚠ 缺日"))

# 跨昼夜对照：同一年必须落在同一批日期上（日期逻辑若与另一模式有分毫差异，配对会全错）
other = BASE / "data" / "daily" / f"daily_dates_{'day' if MODE == 'night' else 'night'}.json"
if other.exists():
    prev = json.loads(other.read_text(encoding="utf-8"))
    # ⚠ JSON 的键**只能是字符串**，而 alldates 的键是 (年, 半段) 元组。
    # 曾经直接写 `k in prev`，元组永远不等于字符串键 → bad 恒为空 →
    # 每次都打印「✓ 全部一致」，检查形同虚设（2026-10-02 审计发现并修）。
    # 必须先规范化成同一套键再比。
    prev_norm = {tuple(k.split("|")): v for k, v in prev.items()}
    prev_norm = {(int(y), h): v for (y, h), v in prev_norm.items()}
    common = sorted(set(alldates) & set(prev_norm))
    bad = [k for k in common if prev_norm[k] != alldates[k]]
    print(f"\n与另一模式逐日对照：{'✓ 全部一致' if not bad else f'✗ {bad} 不一致'}"
          f"（比对 {len(common)} 个半年段）")
    if not common:
        print("  ⚠ 无重叠半年段可比——对照**未实际发生**，不得当作已校验")
    if bad:
        raise SystemExit(1)
else:
    print(f"\n（另一模式清单尚未生成：{other.name}，跳过跨昼夜对照）")

# 清单**合并**写入：分批提交时后一批不得冲掉前一批的日期与任务 ID
def merged(path: Path, new: dict) -> dict:
    old = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    return {**old, **new}


all_dates = merged(MANIFEST, {f"{y}|{h}": v for (y, h), v in alldates.items()})
MANIFEST.write_text(json.dumps(all_dates, ensure_ascii=False), encoding="utf-8")
print(f"\n日期清单已写 {MANIFEST.relative_to(BASE)}（累计 {len(all_dates)} 个半年段）")

if ARGS.dry_run:
    print("\n[--dry-run] 未提交。")
    raise SystemExit(0)

# ------------------------------------------------------------- 提交
print(f"\n{'=' * 78}\n提交\n{'=' * 78}")
tasks = {}
for (y, half), img in built.items():
    # --prefix 只给基名，仍按 年_半段 区分，否则两半会互相覆盖
    name = (f"{ARGS.prefix}_{y}_{half}" if ARGS.prefix
            else f"uhi_{TAG}daily_{y}_{half}")
    task = ee.batch.Export.image.toDrive(
        image=img,
        description=name,
        folder=FOLDER,
        fileNamePrefix=name,
        region=aoi,
        scale=SCALE,
        crs="EPSG:4326",
        maxPixels=1e10,
        fileFormat="GeoTIFF",
    )
    task.start()
    tasks[name] = task.id
    print(f"  {name:<26}{task.id}")

TRACK_PATH = BASE / "code" / "daily" / TRACK
all_tasks = merged(TRACK_PATH, tasks)
TRACK_PATH.write_text(json.dumps(all_tasks, ensure_ascii=False, indent=2), encoding="utf-8")
print(f"\n已提交 {len(tasks)} 个任务 → Drive 文件夹 {FOLDER}")
print(f"任务 ID 存 code/daily/{TRACK}（累计 {len(all_tasks)} 个，含往批）")
