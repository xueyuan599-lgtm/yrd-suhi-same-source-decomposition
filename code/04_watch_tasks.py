"""盯住逐日导出的 8 个 EE 任务，全部离开 READY/RUNNING 即退出。

只读 `02_export_tasks_{night,day}.json`，**不重新提交**。
退出码：0 = 全部 COMPLETED；1 = 有 FAILED/CANCELLED；2 = 超时仍在跑。

用法：
    python 04_watch_tasks.py            # 每 45 s 轮询一次
    python 04_watch_tasks.py --once     # 只查一次就退出
"""

import json
import sys
import time
from pathlib import Path

import ee

sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")   # 只改 stdout 时 stderr 仍走 GBK，中文报错会乱码
ee.Initialize(project="evident-ocean-510002-t0")

BASE = Path(__file__).resolve().parent.parent.parent
ONCE = "--once" in sys.argv
POLL_S, TIMEOUT_S = 45, 5400          # 最多盯 90 分钟

tasks, missing_modes = {}, []
for mode in ("night", "day"):
    p = Path(__file__).resolve().parent / f"02_export_tasks_{mode}.json"
    if p.exists():
        tasks.update(json.loads(p.read_text(encoding="utf-8")))
    else:
        missing_modes.append(mode)
if not tasks:
    print("没有任务清单（先跑 02_export_daily_seq.py）", file=sys.stderr)
    raise SystemExit(1)
if missing_modes:
    # 只跑了一半的清单是可达到的中间态（02 按 --mode 分次提交、分次写清单）。
    # 静默忽略的后果是「只监视子集，却宣布全部结束并返回 0」——与
    # `10_export_cityid.py` 记录的那次事故同型，故这里必须吵。
    print(f"⚠ 缺 {missing_modes} 的任务清单：本次只监视 {len(tasks)} 个任务，"
          f"跑完**不代表**昼夜两批都结束", file=sys.stderr, flush=True)

TERMINAL = {"COMPLETED", "FAILED", "CANCELLED", "UNKNOWN"}
t0 = time.time()
last = {}

while True:
    states = {s["id"]: s for s in ee.data.getTaskStatus(list(tasks.values()))}
    pending = 0
    for name, tid in tasks.items():
        st = states.get(tid, {})
        state = st.get("state", "MISSING")
        # 只在状态变化时打印，避免每轮刷屏（stdout 是事件流）
        if last.get(name) != state:
            err = st.get("error_message", "")
            print(f"[{time.strftime('%H:%M:%S')}] {name} → {state}"
                  + (f"  ERR={str(err)[:120]}" if err else ""), flush=True)
            last[name] = state
        if state not in TERMINAL:
            pending += 1

    if pending == 0:
        bad = [n for n, t in tasks.items() if last.get(n) != "COMPLETED"]
        scope = f"{len(tasks) - len(bad)} 完成 / {len(bad)} 未完成"
        if missing_modes:
            print(f"已监视的 {len(tasks)} 个任务全部结束：{scope}，"
                  f"但缺 {missing_modes} 的清单，**不能宣布全部结束**", flush=True)
        else:
            print(f"全部结束：{scope}" + (f" {bad}" if bad else ""), flush=True)
        raise SystemExit(0 if (not bad and not missing_modes) else 1)

    if ONCE:
        print(f"仍有 {pending} 个在跑", flush=True)
        raise SystemExit(2)
    if time.time() - t0 > TIMEOUT_S:
        print(f"盯了 {TIMEOUT_S // 60} 分钟仍有 {pending} 个在跑，退出", flush=True)
        raise SystemExit(2)
    time.sleep(POLL_S)
