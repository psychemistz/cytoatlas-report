# Methods and Figure Provenance

Per-figure and per-table provenance for the PI report. Every visualization, table, and inline statistic is traced to its data source, statistical method, and generation code.

See also: [PROVENANCE.md](PROVENANCE.md) for pipeline overview, [DATA_DICTIONARY.md](DATA_DICTIONARY.md) for column definitions.

---

## Figures

### Figure 1: Dataset Overview

| Field | Content |
|-------|---------|
| **Shows** | Three-panel summary: (A) cell counts per dataset, (B) signature matrix sizes, (C) validation layer counts |
| **Data source** | Hardcoded values from dataset metadata. Cell counts: CIMA 6.5M, Inflammation 6.3M, scAtlas 6.4M, parse_10M 9.7M. Signature counts: CytoSig 43, LinCytoSig 178, SecAct 1,170. Layer counts: number of atlas x level combinations per aggregation type |
| **Method** | No computation; descriptive bar charts of known quantities |
| **Script** | `generate_report_figures.py` → `fig1_dataset_overview()` (lines 180–224) |
| **Interpretation** | Orients the reader to the scale of the project (29M cells, 1,213 signatures, 4 validation layers). Cell counts come from `.n_obs` of the source H5AD files |

---

### Figure 2: Correlation Summary Boxplot

| Field | Content |
|-------|---------|
| **Shows** | Boxplots of Spearman rho distributions at donor-level pseudobulk, organized in a 2x2 grid by signature type (CytoSig, LinCytoSig Best, LinCytoSig, SecAct), with 6 atlases per panel |
| **Data source** | `results/cross_sample_validation/correlations/{atlas}_correlations.csv` (CIMA, Inflammation, scAtlas) and `correlations/all_correlations.csv` (GTEx, TCGA). LinCytoSig Best data from `results/lincytosig_gene_filter/*_lincyto_filt_correlations.csv` via `best_lincytosig_selection.json` |
| **Method** | Filter to `level == 'donor_only'` (or `donor_organ` for scAtlas). Group by atlas. Box shows Q1–Q3, whiskers extend to 1.5×IQR, outliers hidden (`showfliers=False`). Each data point is one target's Spearman rho |
| **Script** | `generate_report_figures.py` → `fig2_correlation_summary()` (lines 230–281) |
| **Interpretation** | Compare signature performance across atlases. Higher boxes = better validation. SecAct typically has the highest median; Inflammation Atlas shows more spread due to disease-driven variance |

---

### Figure 3: Good/Bad Correlations (CytoSig)

| Field | Content |
|-------|---------|
| **Shows** | Horizontal bar charts of top 15 (best) and bottom 15 (worst) CytoSig targets by Spearman rho for 3 atlases: CIMA (donor), Inflammation Main (donor), scAtlas Normal (donor×organ). Bars colored by cytokine family |
| **Data source** | `correlations/{atlas}_correlations.csv` filtered to `signature == 'cytosig'` at the appropriate donor level |
| **Method** | Sort targets by `spearman_rho` descending. Take `.head(15)` for good, `.tail(15)` for bad. Color by membership in 7 cytokine families (Interferon, TGF-beta, Interleukin, TNF, Growth Factor, Chemokine, Colony-Stimulating) |
| **Script** | `generate_report_figures.py` → `fig3_good_bad_correlations()` (lines 287–345) |
| **Interpretation** | IL1B, TNFA, VEGFA consistently top-ranked. CD40L, TRAIL consistently bottom-ranked (membrane-bound or post-transcriptionally regulated targets) |

---

### Figure 4: Biologically Important Targets Heatmap

| Field | Content |
|-------|---------|
| **Shows** | Heatmap of Spearman rho for biologically important targets across 6 atlases, shown in 3 panels: CytoSig (all 43 targets), LinCytoSig Best (43 targets), SecAct (subset of matched targets). Color scale: RdBu_r, range -0.5 to 0.8 |
| **Data source** | Merged correlation CSVs (same as Fig 2), filtered to donor-level per atlas |
| **Method** | For each target × atlas combination: look up `spearman_rho` at donor level. Build matrix (targets × atlases). Display as `imshow` with RdBu_r colormap. Annotate cells with rho values where non-null |
| **Script** | `generate_report_figures.py` → `fig4_bio_targets_heatmap()` (lines 351–422) |
| **Interpretation** | Reveals atlas-specific patterns. TGFB1 strong in scAtlas but variable in CIMA. CXCL12 shows high rho in scAtlas (organ-level chemokine gradient). Missing cells indicate the target was not measured in that atlas |

---

### Figure 5: Representative Scatter Plots (CIMA)

| Field | Content |
|-------|---------|
| **Shows** | 8 scatter plots (2×4 grid): 4 good correlations (top row) and 4 poor correlations (bottom row) for CIMA donor-level CytoSig. Each point = one donor. X-axis: mean target gene expression. Y-axis: predicted activity z-score |
| **Data source** | `visualization/data/validation/donor_scatter/cima_cytosig.json`. Targets selected from `cima_correlations.csv` by sorting `spearman_rho` |
| **Method** | Good targets: top 4 by rho (rho > 0.3). Bad targets: bottom 4 by rho. Each scatter shows donor-level expression vs activity with a linear regression line. Colored green (good), red (bad), or gray (neutral) based on rho threshold |
| **Script** | `generate_report_figures.py` → `fig5_representative_scatter()` (lines 428–502) |
| **Interpretation** | Good correlations (e.g., IL1B: rho=0.67) show clear positive trends. Poor correlations (e.g., CD40L: rho=-0.48) show inverse or no relationship, explained by membrane-bound biology |

---

### Figure 6: Cross-Atlas Consistency

| Field | Content |
|-------|---------|
| **Shows** | Line chart tracking 14 key cytokine targets across 6 atlases. Solid lines = CytoSig, dashed lines = LinCytoSig Best. X-axis: atlas (CIMA → Inflammation → scAtlas Normal → scAtlas Cancer → GTEx → TCGA). Y-axis: Spearman rho. Lines colored by cytokine family |
| **Data source** | Merged correlation CSVs (CytoSig) and gene-filtered LinCytoSig Best CSVs, both at donor level per atlas |
| **Method** | For each of 14 targets (IFNG, IL1B, TNFA, TGFB1, IL6, IL10, IL17A, IL4, BMP2, EGF, HGF, VEGFA, CXCL12, GMCSF): extract donor-level rho per atlas. Plot as connected lines across atlases |
| **Script** | `generate_report_figures.py` → `fig6_cross_atlas_consistency()` (lines 508–580) |
| **Interpretation** | Consistent positive slopes (IL1B, TNFA) indicate generalizable signatures. Divergent patterns (IL4, IL17A) reflect disease-specific biology. GTEx/TCGA generally higher than single-cell atlases |

---

### Figure 7: Validation Levels Comparison

| Field | Content |
|-------|---------|
| **Shows** | 2×2 grid of boxplots comparing Spearman rho at different aggregation levels (donor → donor×L1 → donor×L2 → ...) for 4 single-cell atlases: CIMA, Inflammation, scAtlas Normal, scAtlas Cancer. Blue = CytoSig, brown = LinCytoSig Best, side by side |
| **Data source** | Merged correlation CSVs filtered to each atlas and each aggregation level. LinCytoSig Best loaded via gene-filtered CSVs |
| **Method** | For each atlas × level × signature: extract rho array. Plot as paired boxplots (CytoSig blue, LinCytoSig Best brown) with median labels. CIMA has 5 levels (donor, L1–L4); Inflammation has 3 (donor, L1, L2); scAtlas has 3 (donor×organ, celltype1, celltype2) |
| **Script** | `generate_report_figures.py` → `fig7_validation_levels()` (lines 586–687) |
| **Interpretation** | Finer cell-type stratification consistently decreases median rho because per-group sample sizes shrink. Donor-level gives highest correlations. The drop quantifies the bias-variance tradeoff of pseudobulk aggregation |

---

### Figure 8: 10-Way Method Comparison

| Field | Content |
|-------|---------|
| **Shows** | 2×2 grid of boxplots comparing 10 methods (CytoSig, 8 LinCytoSig variants, SecAct) across 4 combined atlases. Each box = distribution of matched-cytokine Spearman rhos for that method in that atlas |
| **Data source** | `visualization/data/validation/method_comparison_8way_all.json` (pre-computed by `compute_lincytosig_gene_filter.py`). Contains arrays of rhos for ~20 matched cytokines per method per atlas |
| **Method** | Load pre-computed rho arrays. Display as boxplots with outliers shown. Median labels annotated above each box. Methods ordered: CytoSig → LinCytoSig variants → SecAct |
| **Script** | `generate_report_figures.py` → `fig8_method_comparison()` (lines 693–773) |
| **Interpretation** | SecAct consistently highest median. CytoSig outperforms most LinCytoSig variants at donor level. Gene filtering improves LinCytoSig in 3 of 4 atlases. See REPORT.md Section 5.1 for detailed analysis |

---

### Figure 9: LinCytoSig vs CytoSig Matched Scatter

| Field | Content |
|-------|---------|
| **Shows** | Scatter plots (one per atlas category) with CytoSig rho on x-axis and LinCytoSig rho on y-axis. Each point = one matched target. Diagonal line = equal performance. Points colored: amber = LinCytoSig wins (by >0.05), blue = CytoSig wins, gray = tie |
| **Data source** | `visualization/data/validation/method_comparison.json`. Uses `matched_targets` mapping to pair LinCytoSig targets with CytoSig counterparts, then `lincytosig.rhos` and `cytosig.rhos` for correlation values per category |
| **Method** | For each matched pair: look up rho in both methods for each category. Plot as scatter. Win/loss threshold: absolute difference > 0.05. Top 3 winners and losers labeled with arrows |
| **Script** | `generate_report_figures.py` → `fig9_lincytosig_vs_cytosig()` (lines 779–890) |
| **Interpretation** | Points above diagonal favor LinCytoSig. Win counts shown per panel. CytoSig wins overall, but specific cell types (Basophil, NK, DC) favor LinCytoSig |

---

### Figure 10: LinCytoSig Advantage by Cell Type

| Field | Content |
|-------|---------|
| **Shows** | Two-panel chart. (A) Mean Δρ (LinCytoSig − CytoSig) per cell type, sorted descending. (B) Win rate per cell type (% of comparisons where LinCytoSig or CytoSig wins). Amber = LinCytoSig better, blue = CytoSig better |
| **Data source** | `visualization/data/validation/method_comparison.json`. Extracts cell type from LinCytoSig target names (format `CellType__Cytokine`), computes per-celltype mean difference and win/loss counts across all atlas categories |
| **Method** | For each LinCytoSig target: extract cell type prefix (before `__`). Compute `lin_rho - cyto_rho` for all matched pairs across all categories. Aggregate by cell type: mean difference, win count (diff > 0.05), loss count (diff < -0.05) |
| **Script** | `generate_report_figures.py` → `fig10_lincytosig_advantage()` (lines 896–994) |
| **Interpretation** | Basophil, NK, DC benefit most from cell-type specificity (+0.18 to +0.21 Δρ). Lymphatic Endothelial, Adipocyte, Osteocyte suffer most (too few experiments in LinCytoSig database) |

---

### Figure 11: Cell-Type Δρ (LinCytoSig − CytoSig)

| Field | Content |
|-------|---------|
| **Shows** | Horizontal bar chart of mean Δρ (LinCytoSig − CytoSig) per cell type with SEM error bars, aggregated across all 4 atlases at donor × celltype level. Win/loss counts annotated per bar |
| **Data source** | `visualization/data/validation/method_comparison.json` (same as Fig 10 but different aggregation and presentation) |
| **Method** | Same cell-type extraction as Fig 10. Compute per-cell-type: mean Δρ, SEM (std/sqrt(n)), win count, loss count. Sort by mean Δρ descending. Global summary: overall mean Δρ, count of cell types favoring each method |
| **Script** | `generate_report_figures.py` → `fig11_celltype_delta_rho()` (lines 1000–1083) |
| **Interpretation** | Complements Figure 10 with error bars and sample sizes. Shows statistical confidence of per-cell-type advantages. Summary box gives the global picture |

---

### Figure 12: Bulk RNA-seq Validation (GTEx & TCGA)

| Field | Content |
|-------|---------|
| **Shows** | 2×3 grid: rows = GTEx, TCGA; columns = CytoSig, LinCytoSig, SecAct. For small target sets (≤60): ranked bar chart of Spearman rhos. For large sets (>60, i.e., SecAct): histogram with median vertical line |
| **Data source** | `correlations/all_correlations.csv` filtered to `atlas in ('gtex', 'tcga')`, `level == 'donor_only'`, per signature type |
| **Method** | Filter, sort by rho descending. Color bars by rho threshold (green > 0.2, red < -0.1, gray otherwise). For histograms: 40 bins, median line overlaid |
| **Script** | `generate_report_figures.py` → `fig12_bulk_validation()` (lines 1089–1138) |
| **Interpretation** | Bulk RNA-seq provides the cleanest validation signal (92-100% significance). SecAct achieves highest median rho in both GTEx (0.395) and TCGA (0.415) |

---

### Figure 13: LinCytoSig Specificity Deep Dive

| Field | Content |
|-------|---------|
| **Shows** | Two-panel horizontal bar chart: (A) top 20 LinCytoSig advantages over CytoSig, (B) top 20 CytoSig advantages over LinCytoSig. Both methods shown as overlapping bars for each matched target |
| **Data source** | `visualization/data/validation/method_comparison.json`. Iterates all matched targets across all categories, computes `diff = lin_rho - cyto_rho` |
| **Method** | Collect all matched-pair rhos across all categories. Rank by difference. Take top 20 (LinCytoSig wins) and bottom 20 (CytoSig wins). Display as overlapping horizontal bars with Δ labels |
| **Script** | `generate_report_figures.py` → `fig13_lincytosig_specificity()` (lines 1144–1228) |
| **Interpretation** | Identifies the most extreme cases of method advantage. LinCytoSig wins tend to be immune cell types with distinctive biology. CytoSig wins tend to be cell types with insufficient LinCytoSig training data |

---

### Figure 14: Cell-Type-Level Scatter Examples

| Field | Content |
|-------|---------|
| **Shows** | 4×3 grid of scatter plots for 4 key targets (IFNG, TGFB1, IL1B, IL6) across 3 atlases (CIMA L1, Inflammation Main L1, scAtlas Normal Celltype1). Each point = cell-type × donor pseudobulk. Points optionally colored by cell type |
| **Data source** | `visualization/data/validation/celltype_scatter/{atlas}_{level}_cytosig.json` |
| **Method** | Load scatter JSON for each atlas × level. Extract target data (points, rho, optional celltype indices). Plot scatter with regression line. Title shows target name, atlas, and rho value |
| **Script** | `generate_report_figures.py` → `fig14_celltype_scatter_examples()` (lines 1234–1299) |
| **Interpretation** | Cell-type-level scatter has more data points than donor-level but more variance. Strong targets (IFNG, IL1B) maintain positive trends. Weaker targets show cell-type-dependent behavior |

---

### Figure 15: Summary Statistics Table

| Field | Content |
|-------|---------|
| **Shows** | Publication-ready table with validation statistics per atlas × signature combination. Columns: Atlas, Signature, N Targets, Median ρ, Mean ρ, Std ρ, Min ρ, Max ρ, % Significant, % Positive |
| **Data source** | Merged correlation CSVs at donor level per atlas (same as Fig 2) |
| **Method** | For each atlas × signature: compute `median()`, `mean()`, `std()`, `min()`, `max()` of `spearman_rho`. Count targets with `spearman_pval < 0.05` for significance rate. Count `rho > 0` for positive rate. Render as matplotlib table with colored rows by signature type |
| **Script** | `generate_report_figures.py` → `fig15_summary_table()` (lines 1305–1383). Also exported as `figures/summary_statistics.csv` |
| **Interpretation** | The definitive summary table. Key numbers (e.g., "CIMA CytoSig median ρ = 0.114") are directly readable and traceable to the correlation CSVs |

---

## Tables in REPORT.md

### Table: Section 1.2 — Processing Scale

| Field | Content |
|-------|---------|
| **Shows** | Dataset name, cell count, sample count, processing time, GPU used |
| **Data source** | Hardcoded from `.n_obs` of source H5AD files and SLURM job logs |
| **Method** | No computation; descriptive metadata |

### Table: Section 2.1 — Datasets and Scale

| Field | Content |
|-------|---------|
| **Shows** | Detailed dataset catalog: cells, donors/samples, cell types, source |
| **Data source** | H5AD file metadata (`.n_obs`, `.obs['donor_id'].nunique()`, cell type annotation counts) |
| **Method** | No computation; descriptive metadata from data exploration |

### Table: Section 2.3 — Signature Matrices

| Field | Content |
|-------|---------|
| **Shows** | CytoSig (43 targets), LinCytoSig (178), SecAct (1,170): construction method and source |
| **Data source** | `secactpy.load_cytosig()` and `secactpy.load_secact()` for target counts. LinCytoSig count from `method_comparison.json` target list |
| **Method** | No computation; descriptive table of signature matrix properties |

### Table: Section 4.1 — Key Numbers (CytoSig)

| Field | Content |
|-------|---------|
| **Shows** | Per-atlas CytoSig donor-level: Median ρ, % Significant, % Positive |
| **Data source** | `correlations/{atlas}_correlations.csv` and `correlations/all_correlations.csv` |
| **Method** | Filter to `signature == 'cytosig'`, `level` per `LEVEL_MAP`. Compute `spearman_rho.median()`, `(spearman_pval < 0.05).mean() * 100`, `(spearman_rho > 0).mean() * 100` |
| **Script** | Same computation as `fig15_summary_table()` and `prepare_summary_table()` |

### Table: Section 4.4 — Aggregation-Level Statistics

| Field | Content |
|-------|---------|
| **Shows** | Per-atlas CytoSig median ρ at each aggregation level (donor, +L1, +L2, +L3, +L4) |
| **Data source** | Same correlation CSVs, grouped by `level` column |
| **Method** | Filter to `signature == 'cytosig'` and each aggregation level. Compute `spearman_rho.median()` per group |
| **Script** | Same logic as `fig7_validation_levels()` and `prepare_levels_data()` |

### Table: Section 4.5 — Comprehensive Validation (3 tables)

| Field | Content |
|-------|---------|
| **Shows** | Three tables (CytoSig, LinCytoSig, SecAct): per-atlas N Targets, Median ρ, Mean ρ, % Significant, % Positive |
| **Data source** | `correlations/{atlas}_correlations.csv` and `correlations/all_correlations.csv` |
| **Method** | Same as Section 4.1 but expanded to all three signature types |
| **Script** | Same computation as `fig15_summary_table()` |

### Table: Section 5.1 — 10-Way Method Comparison

| Field | Content |
|-------|---------|
| **Shows** | Per-atlas median Spearman ρ for 10 methods (CytoSig, 8 LinCytoSig variants, SecAct) across 4 atlases |
| **Data source** | `visualization/data/validation/method_comparison_8way_all.json` |
| **Method** | Median of each method's rho array per atlas. Same data as Figure 8 boxplots |
| **Script** | `prepare_method_comparison_boxplot()` in `generate_interactive_report.py` |

### Table: Section 5.2 — LinCytoSig vs CytoSig Win/Loss Counts

| Field | Content |
|-------|---------|
| **Shows** | Per-atlas win/loss/tie counts for LinCytoSig vs CytoSig (threshold: ±0.05) |
| **Data source** | `visualization/data/validation/method_comparison.json` |
| **Method** | Same as Figure 9. For each matched pair per category: count `lin_rho > cyto_rho + 0.05` as win, `cyto_rho > lin_rho + 0.05` as loss, else tie |
| **Script** | Same computation as `fig9_lincytosig_vs_cytosig()` |

---

## Key Statistics Trace

Inline numbers from REPORT.md prose, traced to their data source.

| Statistic | Section | Source File | Filter / Computation |
|-----------|---------|-------------|---------------------|
| "~29 million cells" | Exec Summary | H5AD `.n_obs` sums | 6.5M + 6.3M + 6.4M + 9.7M |
| "1,213 signatures" | Exec Summary | `secactpy` | 43 (CytoSig) + 1,170 (SecAct) |
| "178 cell-type-specific" | Exec Summary | LinCytoSig signature count | 45 cell types × 1-13 cytokines |
| "ρ=0.6-0.9 for IL1B, TNFA, VEGFA, TGFB" | Exec Summary | `correlations/*.csv` | Per-target rho at donor level across atlases |
| "CIMA median ρ = 0.114" | 4.1 | `cima_correlations.csv` | `atlas=='cima', level=='donor_only', signature=='cytosig'` → `.spearman_rho.median()` |
| "72.1% significant" | 4.1 | `cima_correlations.csv` | `(spearman_pval < 0.05).mean() * 100` |
| "IL1B ρ = 0.67 in CIMA" | 4.2 | `cima_correlations.csv` | `target=='IL1B', level=='donor_only', signature=='cytosig'` → `spearman_rho` |
| "CD40L ρ = -0.48 in CIMA" | 4.2 | `cima_correlations.csv` | `target=='CD40L', level=='donor_only', signature=='cytosig'` |
| "CIMA 0.114 → 0.005 from donor to L4" | 4.4 | `cima_correlations.csv` | `atlas=='cima', signature=='cytosig'` at levels `donor_only` and `donor_l4` |
| "Basophil +0.21 mean Δρ" | 5.2 | `method_comparison.json` | Mean of `lin_rho - cyto_rho` for all Basophil__ targets across categories |
| "SecAct median ρ = 0.395 in GTEx" | 4.5 | `all_correlations.csv` | `atlas=='gtex', level=='donor_only', signature=='secact'` → `.spearman_rho.median()` |
| "97.1% positive in TCGA" | 5.4 | `all_correlations.csv` | `atlas=='tcga', signature=='secact'` → `(spearman_rho > 0).mean() * 100` |
| "262 REST endpoints" | 1.1 | FastAPI OpenAPI spec | Count from `/docs` endpoint |
| "11.4K LOC" | 1.1 | `cloc` or `wc -l` on `cytoatlas-api/frontend/src/` | TypeScript + TSX source lines |

---

## Reproducibility

To regenerate all figures from the intermediate data files:

```bash
cd /data/parks34/projects/2cytoatlas
source ~/bin/myconda && conda activate secactpy

# Static figures (PNG + PDF)
python scripts/05_figures.py --all
# or equivalently:
python report/generate_report_figures.py

# Interactive HTML report
python report/generate_interactive_report.py
```

To regenerate the intermediate data files from raw H5AD (requires GPU node, ~8 hours):

```bash
sbatch scripts/slurm/run_all.sh --main
```

See `PROVENANCE.md` for the complete pipeline diagram.
