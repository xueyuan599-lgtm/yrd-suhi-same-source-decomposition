# Reproducibility package - `Same-Day Pairing and Temporal Aggregation of Terra Urban-Rural Surface Temperature Contrasts in the Yangtze River Delta`

**Author:** Yuan Xu, School of Data Science, Zhejiang University of Finance and Economics, Hangzhou 310018, China; ORCID 0009-0003-9590-102X

**Contact:** xuyuan2@zufe.edu.cn

This repository holds the analysis code and the frozen derived tables that regenerate every number, table and figure in the manuscript. code/ is the numbered analysis scripts as run; data/ holds the frozen inputs they read (per-city-day urban and rural summaries from MOD11A1, the city-year decomposition, the sensitivity outputs, and the city covariates); nvironment.txt is the `pip freeze` of the interpreter that produced the reported numbers; MANIFEST.sha256 lists the SHA-256 digest of every file so a reader can verify their copy.

**License:** Creative Commons Attribution 4.0 International (CC BY 4.0) - see LICENSE. Please cite the manuscript when reusing this material.

---

# Reproducibility package

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

