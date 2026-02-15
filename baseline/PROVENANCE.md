# Data Provenance

How every number in the PI report traces back to raw data.

---

## Pipeline Overview

```
Raw H5AD files (4 datasets, ~29M cells)
│
├─ scripts/01_cima_activity.py ──────────────► results/cima/
├─ scripts/02_inflam_activity.py ────────────► results/inflammation/
├─ scripts/03_scatlas_analysis.py ───────────► results/scatlas/
├─ scripts/18_parse10m_activity.py ──────────► results/parse10m/
├─ scripts/19_tahoe_activity.py ─────────────► results/tahoe/
└─ scripts/20_spatial_activity.py ───────────► results/spatial/
     │
     ▼
Validation pipelines
├─ scripts/10_atlas_validation_pipeline.py ──► results/cross_sample_validation/
├─ scripts/11_donor_level_pipeline.py ───────► results/cross_sample_validation/correlations/
├─ scripts/15_bulk_validation.py ────────────► results/cross_sample_validation/correlations/all_correlations.csv
└─ scripts/16_resampled_validation.py ───────► results/atlas_validation/
     │
     ▼
Preprocessing for visualization
├─ scripts/14_preprocess_bulk_validation.py ─► visualization/data/validation/donor_scatter/
│                                              visualization/data/validation/celltype_scatter/
├─ scripts/17_preprocess_validation_summary.py ► visualization/data/validation/level_comparison.json
└─ scripts/compute_lincytosig_gene_filter.py ► results/lincytosig_gene_filter/
                                               visualization/data/validation/method_comparison*.json
                                               visualization/data/validation/best_lincytosig_selection.json
     │
     ▼
Report generation (this repo)
├─ scripts/generate_report_figures.py ───────► figures/fig1–fig15 (PNG + PDF)
│                                              figures/summary_statistics.csv
└─ scripts/generate_interactive_report.py ───► index.html (self-contained with embedded Plotly.js)
```

All paths are relative to the upstream `2cytoatlas` repository root unless otherwise noted.

---

## Per-Figure Data Lineage

Each figure traces through three stages: upstream pipeline script(s), intermediate data file(s), and the report generation function.

| Figure | Upstream Script(s) | Intermediate File(s) | Report Function |
|--------|-------------------|----------------------|-----------------|
| 1 | None (hardcoded counts) | — | `fig1_dataset_overview()` |
| 2 | 11, 15 | `correlations/{atlas}_correlations.csv`, `correlations/all_correlations.csv` | `fig2_correlation_summary()` |
| 3 | 11, 15 | `correlations/{atlas}_correlations.csv`, `correlations/all_correlations.csv` | `fig3_good_bad_correlations()` |
| 4 | 11, 15 | `correlations/{atlas}_correlations.csv`, `correlations/all_correlations.csv` | `fig4_bio_targets_heatmap()` |
| 5 | 14 | `validation/donor_scatter/cima_cytosig.json` | `fig5_representative_scatter()` |
| 6 | 11, 15, `compute_lincytosig_gene_filter.py` | `correlations/{atlas}_correlations.csv`, `lincytosig_gene_filter/*_lincyto_filt_correlations.csv`, `validation/best_lincytosig_selection.json` | `fig6_cross_atlas_consistency()` |
| 7 | 11, 15 | `correlations/{atlas}_correlations.csv`, `correlations/all_correlations.csv` | `fig7_validation_levels()` |
| 8 | `compute_lincytosig_gene_filter.py` | `validation/method_comparison_8way_all.json` | `fig8_method_comparison()` |
| 9 | `compute_lincytosig_gene_filter.py` | `validation/method_comparison.json` | `fig9_lincytosig_vs_cytosig()` |
| 10 | `compute_lincytosig_gene_filter.py` | `validation/method_comparison.json` | `fig10_lincytosig_advantage()` |
| 11 | `compute_lincytosig_gene_filter.py` | `validation/method_comparison.json` | `fig11_celltype_delta_rho()` |
| 12 | 15 | `correlations/all_correlations.csv` | `fig12_bulk_validation()` |
| 13 | `compute_lincytosig_gene_filter.py` | `validation/method_comparison.json` | `fig13_lincytosig_specificity()` |
| 14 | 14 | `validation/celltype_scatter/{atlas}_{level}_cytosig.json` | `fig14_celltype_scatter_examples()` |
| 15 | 11, 15 | `correlations/{atlas}_correlations.csv`, `correlations/all_correlations.csv` | `fig15_summary_table()` |

Script numbers refer to `scripts/{num}_*.py` in the main `2cytoatlas` repo.

---

## Upstream Script Reference

| Script | Input | Output | Purpose |
|--------|-------|--------|---------|
| `01_cima_activity.py` | CIMA H5AD (6.5M cells) | `results/cima/` activity H5ADs + CSVs | Ridge regression activity inference for CIMA |
| `02_inflam_activity.py` | Inflammation H5AD (6.3M cells, 3 cohorts) | `results/inflammation/` | Activity inference for Inflammation Atlas |
| `03_scatlas_analysis.py` | scAtlas H5ADs (6.4M cells) | `results/scatlas/` | Activity inference for scAtlas normal + cancer |
| `10_atlas_validation_pipeline.py` | Activity H5ADs from 01-03 | `results/cross_sample_validation/` | Multi-level validation framework |
| `11_donor_level_pipeline.py` | Activity H5ADs from 01-03 | `correlations/{atlas}_correlations.csv` | Donor-level Spearman correlations per target |
| `14_preprocess_bulk_validation.py` | Correlation CSVs, activity H5ADs | `validation/donor_scatter/*.json`, `validation/celltype_scatter/*.json` | Split scatter data into per-atlas JSON files |
| `15_bulk_validation.py` | GTEx + TCGA expression data | `correlations/all_correlations.csv` | Bulk RNA-seq activity inference + correlations |
| `16_resampled_validation.py` | Activity H5ADs | `results/atlas_validation/` resampled H5ADs | Bootstrap resampled pseudobulk for CIs |
| `17_preprocess_validation_summary.py` | Correlation CSVs, method comparison JSONs | `validation/level_comparison.json` | Aggregation-level comparison preprocessing |
| `compute_lincytosig_gene_filter.py` | LinCytoSig signatures, correlation CSVs | `lincytosig_gene_filter/*_lincyto_filt_correlations.csv`, `validation/method_comparison*.json`, `validation/best_lincytosig_selection.json` | Gene-filtered LinCytoSig + 10-way method comparison |

---

## How to Verify a Specific Number

**Example: "CIMA CytoSig median ρ = 0.114" (REPORT.md Section 4.1)**

1. **Find the intermediate file:**
   ```
   results/cross_sample_validation/correlations/cima_correlations.csv
   ```

2. **Filter to the relevant subset:**
   ```python
   import pandas as pd
   df = pd.read_csv('cima_correlations.csv')
   sub = df[(df['atlas'] == 'cima') & (df['level'] == 'donor_only') & (df['signature'] == 'cytosig')]
   print(sub['spearman_rho'].median())  # → 0.114
   ```

3. **Trace back to the pipeline script:**
   - `cima_correlations.csv` was produced by `scripts/11_donor_level_pipeline.py`
   - Which consumed activity scores from `scripts/01_cima_activity.py`
   - Which ran ridge regression on the CIMA H5AD: `/data/Jiang_Lab/Data/Seongyong/CIMA/Cell_Atlas/CIMA_RNA_6484974cells_36326genes_compressed.h5ad`

4. **Verify in the report script:**
   - `fig15_summary_table()` in `generate_report_figures.py` (line 1306) computes the same median
   - `prepare_summary_table()` in `generate_interactive_report.py` (line 161) computes it for the HTML

**Example: "LinCytoSig Best (comb+filt) wins 63 out of 136 matched targets in CIMA"**

1. **Find the intermediate file:**
   ```
   visualization/data/validation/method_comparison.json
   ```

2. **Compute from source:**
   ```python
   import json, numpy as np
   mc = json.load(open('method_comparison.json'))
   cat_key = 'cima_celltype'  # CIMA category
   wins = 0
   for lin_target, match_info in mc['matched_targets'].items():
       cyto = match_info.get('cytosig')
       if cyto is None: continue
       lin_rho = mc['lincytosig']['rhos'].get(lin_target, {}).get(cat_key)
       cyto_rho = mc['cytosig']['rhos'].get(cyto, {}).get(cat_key)
       if lin_rho is not None and cyto_rho is not None:
           if lin_rho > cyto_rho + 0.05:
               wins += 1
   print(wins)  # → 63
   ```

3. **Trace to the report:** `fig9_lincytosig_vs_cytosig()` in `generate_report_figures.py` (line 779) computes these same win counts.
