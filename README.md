# CytoAtlas: PI Report

Interactive report for the Pan-Disease Single-Cell Cytokine Activity Atlas.

**Live report:** [https://psychemistz.github.io/2cytoatlas-report/](https://psychemistz.github.io/2cytoatlas-report/)

## Contents

| File | Description |
|------|-------------|
| [`index.html`](index.html) | Interactive HTML report with Plotly charts (GitHub Pages entry point) |
| [`methodology.html`](methodology.html) | LinCytoSig methodology supplement |
| [`REPORT.md`](REPORT.md) | Markdown version of the full report |
| [`figures/`](figures/) | 15 static figures (PNG, 300 DPI) + summary statistics |
| [`scripts/`](scripts/) | Python scripts that generated the report and figures |

## Figures

| # | Figure | Description |
|---|--------|-------------|
| 1 | Dataset overview | Cell counts, signature matrices, validation layers |
| 2 | Correlation summary | Spearman rho distributions across atlases |
| 3 | Good/bad correlations | Top/bottom 15 targets per atlas (CytoSig) |
| 4 | Bio targets heatmap | Biologically important targets across atlases |
| 5 | Representative scatter | Good/bad scatter plots (CIMA) |
| 6 | Cross-atlas consistency | Target consistency profiles |
| 7 | Validation levels | Aggregation level effect on correlation |
| 8 | Method comparison | 10-way: CytoSig, LinCytoSig x8, SecAct |
| 9 | LinCytoSig vs CytoSig | Matched target scatter |
| 10 | LinCytoSig advantage | Cell-type-specific advantage analysis |
| 11 | SecAct novel signatures | Novel high-correlation targets |
| 12 | Bulk validation | GTEx/TCGA bulk RNA-seq validation |
| 13 | LinCytoSig specificity | Wins vs losses by cell type |
| 14 | Celltype scatter | Cell-type-level scatter examples |
| 15 | Summary table | Complete validation statistics |

## Documentation

| Document | Description |
|----------|-------------|
| [`METHODS.md`](METHODS.md) | Per-figure and per-table provenance: data source, statistical method, script reference, interpretation |
| [`DATA_DICTIONARY.md`](DATA_DICTIONARY.md) | Column definitions for every CSV and JSON file consumed by the report scripts |
| [`PROVENANCE.md`](PROVENANCE.md) | Pipeline diagram and data flow from raw H5AD to final report |
| [`methodology.html`](methodology.html) | LinCytoSig methodology supplement |

## Generation Scripts

The scripts in `scripts/` reproduce all figures and the interactive HTML from the CytoAtlas analysis results. They require access to the HPC data environment and are included here for methodological transparency.

- `generate_report_figures.py` — generates 15 static matplotlib figures (PNG + PDF)
- `generate_interactive_report.py` — generates the interactive HTML with embedded Plotly.js charts
