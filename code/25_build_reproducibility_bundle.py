"""Assemble the submission-side reproducibility package (response to reviewer point 5).

This is **not** an archive of the manuscript: it is the code + frozen inputs that regenerate
every number, table and figure reported in the paper, which the journal asks for. Files are
copied byte-for-byte and a SHA-256 manifest is written, so a reader can verify what they hold.

Output: `submission/reproducibility/`
    README.md            what each script produces and the rerun order
    environment.txt      `pip freeze` of the interpreter that actually ran the analysis
    code/                analysis/code/daily/*.py as run
    data/                the frozen CSVs that back Tables 2–10 and the daily records
    MANIFEST.sha256      SHA-256 of every bundled file
"""

import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

WS = Path(__file__).resolve().parents[2]
ROOT = WS.parent
BUNDLE = ROOT / "submission" / "reproducibility"
CODE_DIR = WS / "code" / "daily"
ARCH = WS / "tables" / "_archive_pre_major_revision"

DATA_FILES = [
    (WS / "outputs" / "daily" / "same_source_cluster_bootstrap.csv", "data"),
    (WS / "outputs" / "daily" / "same_source_decomposition_city_year.csv", "data"),
    (WS / "outputs" / "daily" / "same_source_revision_sensitivity_table.csv", "data"),
    (WS / "outputs" / "daily" / "same_source_same_month_window_sensitivity.csv", "data"),
    (WS / "outputs" / "daily" / "same_source_available_month_trends.csv", "data"),
    (WS / "outputs" / "daily" / "support_zone_size.csv", "data"),
    (WS / "outputs" / "daily" / "support_single_sided.csv", "data"),
    (WS / "outputs" / "daily" / "support_missing_units.csv", "data"),
    (ARCH / "qc_tier_sensitivity_day.csv", "data"),
    (ARCH / "qc_tier_sensitivity_night.csv", "data"),
    (WS / "data" / "daily" / "daily_records_day.csv", "data"),
    (WS / "data" / "daily" / "daily_records_night.csv", "data"),
    (WS / "data" / "city_covariates.csv", "data"),
    (WS / "data" / "daily" / "zones_manifest.json", "data"),
]

README = """# Reproducibility package

Everything in the manuscript is regenerated from the frozen inputs in `data/` by the numbered
scripts in `code/`. Run them from the `analysis/` directory in this order:

    16_build_revision_outputs.py     Table 2, Table 3, Table 4   (and the legacy Figure 2/3 pass)
    17_same_month_window_sensitivity.py   Table 7
    18_render_publication_figures.py Figure 2, Figure 3
    19_render_city_maps.py           Figure 4
    20_render_mechanism_example.py   Figure 5
    21_build_supplementary.py       Supplementary Tables S1-S6 and Figures S1-S2
    22_support_asymmetry.py         single-sided-day summary; Table S5 (products)
    23_build_support_tables.py      Table 6, Table 8, Table 10
    24_support_diagnostics.py       Table 5, Table 9, missing-unit and zone-size diagnostics

The two manuscript-gate scripts live one level up in `my-paper/`: `build_draft.py` assembles
`draft.md` from `draft-parts/`, `captions.md` and `analysis/tables/`; `build_formatted.py`
produces the formatted manuscript; `my-paper-zh/assemble_zh.py` checks that every numeric token
in the Chinese text matches the English source; `build_docx.py` renders both DOCX files.

`environment.txt` is the `pip freeze` of the interpreter that produced the reported numbers.
`MANIFEST.sha256` lists the SHA-256 digest of every file in this package, so a reader can verify
that their copy is unmodified.

Note on the two daily-record files: they are the frozen per-city-day urban and rural summaries
(valid-pixel counts and zone means) from which every estimator, table and figure is derived.
They are included because they, rather than the raw 5.7 GB raster archive, are what the reported
numbers depend on.
"""


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    (BUNDLE / "code").mkdir(parents=True, exist_ok=True)
    (BUNDLE / "data").mkdir(parents=True, exist_ok=True)
    bundled = []
    for src in sorted(CODE_DIR.glob("*.py")):
        dst = BUNDLE / "code" / src.name
        shutil.copy2(src, dst)
        bundled.append(dst)
    for src, sub in DATA_FILES:
        if not src.exists():
            raise SystemExit(f"缺少输入文件：{src}")
        dst = BUNDLE / sub / src.name
        shutil.copy2(src, dst)
        bundled.append(dst)
    (BUNDLE / "README.md").write_text(README, encoding="utf-8")
    freeze = subprocess.run([sys.executable, "-m", "pip", "freeze"],
                            capture_output=True, text=True).stdout
    (BUNDLE / "environment.txt").write_text(
        f"# python {sys.version.split()[0]} · {sys.executable}\n" + freeze, encoding="utf-8")
    manifest = [f"{sha256(p)}  {p.relative_to(BUNDLE).as_posix()}"
                for p in sorted(BUNDLE.rglob("*")) if p.is_file()]
    (BUNDLE / "MANIFEST.sha256").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    total = sum(p.stat().st_size for p in BUNDLE.rglob("*") if p.is_file())
    print(f"已写 {BUNDLE.relative_to(ROOT)}：{len(bundled)} 个文件 + README/环境/清单")
    print(f"  合计 {total/1e6:.1f} MB · 清单 {len(manifest)} 行")


if __name__ == "__main__":
    main()
