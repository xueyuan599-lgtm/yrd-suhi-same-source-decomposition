"""平台映射的**唯一源**（2026-10-11）。

Terra 与 Aqua 的 LST 谱系在集合名、输出名前缀、清单文件名上平行但**不共用文件**——
本模块把这些差异收在一处，导出脚本与本地分析脚本一律 `from products import get`，
不在各自文件里各写一份字典（否则改一处漏一处，且是静默的）。

## 为什么用不同的输出前缀

`03` / `09` 用 `uhi_{tag}daily_*.tif` 做 glob 圈定 Terra 的输入。Aqua 若沿用同一前缀，
两个平台的文件会**同时命中同一个 glob**，被无差别读进同一条流水线——
不报错，但结果是把两个平台混成一份。故 Aqua 一律 `mya_` 前缀 + `_aqua` 输出后缀。

## 已实测的事实（2026-10-11，直接查 EE）

| 集合 | 波段名 | 起始 |
|---|---|---|
| `MODIS/061/MOD11A1` | LST_Day_1km / LST_Night_1km / QC_Day / QC_Night | 2000-02-24 |
| `MODIS/061/MYD11A1` | **与 MOD11A1 逐字相同** | **2002-07-04** |
| `MODIS/061/MYD11A2` | 与 MOD11A2 逐字相同（含 8 天合成期） | **2002-07-04** |

Aqua 起始晚于 Terra 一年半 → **完整年自 2003 起**；`years` 字段即由此而来。
2002-07～2002-12 的半段在完整年口径下不可用，故不导。
"""

from __future__ import annotations

import sys
from pathlib import Path

PRODUCTS: dict[str, dict] = {
    "terra": {
        "collection_daily": "MODIS/061/MOD11A1",
        "collection_8day": "MODIS/061/MOD11A2",
        "stem_daily": "uhi_{tag}daily",
        "stem_8day": "uhi_{tag}seq",
        "suffix": "",                       # 既有产物名一概不动
        "manifest": "daily_dates_{tag}.json",
        "track_daily": "02_export_tasks_{tag}.json",
        "label": "Terra",
        "label_daily": "MOD11A1",
        "label_8day": "MOD11A2",
        "years": (2001, 2025),
    },
    "aqua": {
        "collection_daily": "MODIS/061/MYD11A1",
        "collection_8day": "MODIS/061/MYD11A2",
        "stem_daily": "mya_{tag}daily",
        "stem_8day": "mya_{tag}seq",
        "suffix": "_aqua",
        "manifest": "daily_dates_aqua_{tag}.json",
        "track_daily": "02_export_tasks_aqua_{tag}.json",
        "label": "Aqua",
        "label_daily": "MYD11A1",
        "label_8day": "MYD11A2",
        "years": (2003, 2025),
    },
}

CHOICES = tuple(PRODUCTS)


def get(product: str) -> dict:
    """取某平台的配置；未知平台立即报错（不静默退回 terra）。"""
    if product not in PRODUCTS:
        raise SystemExit(f"未知 product={product!r}；可选 {CHOICES}")
    return PRODUCTS[product]


def add_sys_path() -> None:
    """让 `code/` 与 `code/daily/` 下的脚本都能 `import products`。"""
    here = Path(__file__).resolve().parent
    for p in (str(here), str(here / "daily")):
        if p not in sys.path:
            sys.path.insert(0, p)
