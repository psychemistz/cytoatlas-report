# Data Dictionary

Column and field definitions for every data file consumed by the report generation scripts.

---

## 1. Correlation CSVs

**Files:** `results/cross_sample_validation/correlations/{atlas}_correlations.csv`, `correlations/all_correlations.csv`

Per-atlas donor-level and cell-type-level Spearman correlations between predicted activity and target gene expression.

| Column | Type | Description |
|--------|------|-------------|
| `target` | str | Signature target name (e.g., `IFNG`, `IL1B`, `Macrophage__IFNG` for LinCytoSig) |
| `gene` | str | HGNC gene symbol used for expression measurement (e.g., `IFNG`, `INHBA` for Activin A) |
| `celltype` | str | Cell type used for aggregation. `all` for donor-only level; specific type name otherwise |
| `spearman_rho` | float | Spearman rank correlation between mean expression and predicted activity across donors/samples |
| `spearman_pval` | float | Two-sided p-value for the Spearman correlation |
| `n_samples` | int | Number of donors or pseudobulk samples in the correlation |
| `mean_expr` | float | Mean of the target gene expression across samples |
| `std_expr` | float | Standard deviation of target gene expression |
| `mean_activity` | float | Mean of predicted activity z-scores across samples |
| `std_activity` | float | Standard deviation of predicted activity z-scores |
| `lincytosig_celltype` | str or null | For LinCytoSig: the cell type used in the signature (e.g., `Macrophage`). Null for CytoSig/SecAct |
| `matched_atlas_celltypes` | str or null | Atlas cell types that mapped to the LinCytoSig cell type. Null for CytoSig/SecAct |
| `atlas` | str | Atlas identifier: `cima`, `inflammation_main`, `inflammation_val`, `inflammation_ext`, `scatlas_normal`, `scatlas_cancer`, `gtex`, `tcga` |
| `level` | str | Aggregation level: `donor_only`, `donor_l1`, `donor_l2`, `donor_l3`, `donor_l4`, `donor_organ`, `donor_organ_celltype1`, `donor_organ_celltype2` |
| `signature` | str | Signature matrix used: `cytosig`, `lincytosig`, `secact` |

**Notes:**
- `all_correlations.csv` contains GTEx and TCGA bulk RNA-seq results; per-atlas CSVs contain single-cell results.
- The report scripts merge `inflammation_main`, `inflammation_val`, and `inflammation_ext` into a single `inflammation` atlas.

---

## 2. LinCytoSig Gene-Filtered Correlation CSVs

**Files:** `results/lincytosig_gene_filter/{atlas}_lincyto_filt_correlations.csv`

Same schema as Section 1, but computed using gene-filtered LinCytoSig signatures (restricted to CytoSig's ~4,881 gene space). The `target` column contains `CellType__Cytokine` format names (e.g., `Macrophage__IFNG`). The `signature` column is `lincytosig` in the raw files; the report scripts remap it to `lincyto_best` after selecting the best cell-type variant per cytokine.

---

## 3. Method Comparison JSON

**File:** `visualization/data/validation/method_comparison.json`

Matched-target comparison of CytoSig vs LinCytoSig at the cell-type level.

```
{
  "categories": [
    {"key": "cima_celltype", "label": "CIMA (celltype)"},
    {"key": "inflammation_celltype", "label": "Inflammation Atlas (celltype)"},
    ...
  ],
  "cytosig": {
    "rhos": {"IFNG": {"cima_celltype": 0.45, ...}, ...},
    "targets": ["Activin A", "BDNF", ...]
  },
  "lincytosig": {
    "rhos": {"Macrophage__IFNG": {"cima_celltype": 0.52, ...}, ...},
    "targets": ["B_Cell__IL13", ...]
  },
  "secact": {
    "rhos": {"IFNG": {"cima_celltype": 0.40, ...}, ...},
    "targets": [...]
  },
  "matched_targets": {
    "Macrophage__IFNG": {"cytosig": "IFNG", "secact": null},
    ...
  },
  "n_matched": 136
}
```

| Field | Description |
|-------|-------------|
| `categories` | List of atlas × level combinations. Each has `key` (used as lookup) and `label` (display name) |
| `cytosig.rhos` | Dict mapping target name to dict of category_key → Spearman rho |
| `lincytosig.rhos` | Same structure, keyed by `CellType__Cytokine` names |
| `secact.rhos` | Same structure for SecAct targets |
| `matched_targets` | Dict mapping each LinCytoSig target to its CytoSig/SecAct counterparts (or null) |
| `n_matched` | Number of matched target pairs |

**Used by:** Figures 9, 10, 11, 13.

---

## 4. 10-Way Method Comparison JSON

**File:** `visualization/data/validation/method_comparison_8way_all.json`

Pre-computed Spearman rho arrays for 10 methods across 4 combined atlases.

```
{
  "CIMA": {
    "cytokines": ["IFNG", "IL1B", ...],
    "cytosig": [0.45, 0.67, ...],
    "lincyto_orig": [0.32, ...],
    "lincyto_filt": [0.38, ...],
    "lincyto_best_orig": [...],
    "lincyto_best_filt": [...],
    "lincyto_best_gtex": [...],
    "lincyto_best_tcga": [...],
    "lincyto_best_gtex_filt": [...],
    "lincyto_best_tcga_filt": [...],
    "secact": [0.50, ...]
  },
  "Inflammation Atlas": {...},
  "scAtlas Normal": {...},
  "scAtlas Cancer": {...}
}
```

| Field | Description |
|-------|-------------|
| `cytokines` | Ordered list of matched cytokine names (same across all methods) |
| `cytosig` | Array of Spearman rhos for CytoSig, one per cytokine |
| `lincyto_orig` | LinCytoSig (no gene filter) rhos |
| `lincyto_filt` | LinCytoSig (gene-filtered to CytoSig gene space) rhos |
| `lincyto_best_orig` | Best cell-type variant per cytokine (combined GTEx+TCGA selection, all genes) |
| `lincyto_best_filt` | Best cell-type variant (combined selection, gene-filtered) |
| `lincyto_best_gtex` | Best variant selected by GTEx correlation only |
| `lincyto_best_tcga` | Best variant selected by TCGA correlation only |
| `lincyto_best_gtex_filt` | GTEx-selected, gene-filtered |
| `lincyto_best_tcga_filt` | TCGA-selected, gene-filtered |
| `secact` | SecAct rhos for matched targets |

**Used by:** Figure 8.

---

## 5. Best LinCytoSig Selection JSON

**File:** `visualization/data/validation/best_lincytosig_selection.json`

Maps each of the 43 CytoSig cytokines to the best-performing LinCytoSig cell-type variant.

```
{
  "best_orig": {"IFNG": "Macrophage__IFNG", ...},
  "best_filt": {"IFNG": "NK_Cell__IFNG", ...},
  "best_orig_gtex": {"IFNG": "Macrophage__IFNG", ...},
  "best_orig_tcga": {"IFNG": "T_Cell__IFNG", ...}
}
```

| Field | Description |
|-------|-------------|
| `best_orig` | Best cell-type variant per cytokine selected via combined GTEx+TCGA correlation, all genes |
| `best_filt` | Best variant via combined bulk, restricted to CytoSig gene space |
| `best_orig_gtex` | Best variant selected by GTEx correlation alone |
| `best_orig_tcga` | Best variant selected by TCGA correlation alone |

Values are `CellType__Cytokine` strings (e.g., `Macrophage__IFNG`).

**Used by:** `load_lincyto_best_data()` in both report scripts. The `best_filt` mapping is used as `lincyto_best` throughout the report.

---

## 6. Donor Scatter JSONs

**Files:** `visualization/data/validation/donor_scatter/{atlas}_{signature}.json`

Per-target scatter plot data (mean expression vs predicted activity) at donor level.

```
{
  "IFNG": {
    "gene": "IFNG",
    "rho": 0.452,
    "pval": 1.2e-05,
    "n": 421,
    "points": [[expr_1, activity_1], [expr_2, activity_2], ...]
  },
  ...
}
```

| Field | Description |
|-------|-------------|
| `gene` | HGNC gene symbol |
| `rho` | Spearman correlation coefficient |
| `pval` | Spearman p-value |
| `n` | Number of donors/samples |
| `points` | Array of `[mean_expression, predicted_activity]` pairs, one per donor |

**Used by:** Figures 5 (CIMA donor scatter).

---

## 7. Cell-Type Scatter JSONs

**Files:** `visualization/data/validation/celltype_scatter/{atlas}_{level}_{signature}.json`

Same schema as donor scatter (Section 6), but each point represents a cell-type x donor pseudobulk aggregate. May include an optional `celltypes` array and a third element per point indexing into it.

**Used by:** Figure 14 (cell-type scatter examples).

---

## 8. Level Comparison JSON

**File:** `visualization/data/validation/level_comparison.json`

Three-way matched comparison (CytoSig vs LinCytoSig vs SecAct) at each cell-type aggregation level.

```
{
  "CIMA": [
    {
      "level": "L1",
      "n_matched_2way": 136,
      "n_matched_3way": 22,
      "vs_cyto_win": 63,
      "vs_cyto_loss": 64,
      "vs_cyto_mean": -0.012,
      "cyto_median": 0.062,
      "lin_median": 0.055,
      "vs_sec_win": 45,
      "vs_sec_loss": 78,
      "vs_sec_mean": -0.08,
      "sec_median": 0.10,
      "cyto_rhos": [...],
      "lin_rhos": [...],
      "sec_rhos": [...]
    },
    ...
  ],
  ...
}
```

| Field | Description |
|-------|-------------|
| `level` | Aggregation level label (L1, L2, L3, L4, Celltype1, Celltype2) |
| `n_matched_2way` | Number of targets matched between CytoSig and LinCytoSig |
| `n_matched_3way` | Number of targets matched across all three methods |
| `vs_cyto_win` | Count of targets where LinCytoSig rho > CytoSig rho |
| `vs_cyto_loss` | Count where CytoSig rho > LinCytoSig rho |
| `vs_cyto_mean` | Mean difference (LinCytoSig rho - CytoSig rho) |
| `cyto_median` | Median Spearman rho for CytoSig at this level |
| `lin_median` | Median Spearman rho for LinCytoSig at this level |
| `vs_sec_*` | Same comparison vs SecAct (may be null if <10 matched targets) |
| `sec_median` | Median Spearman rho for SecAct at this level |
| `cyto_rhos`, `lin_rhos`, `sec_rhos` | Full arrays of per-target Spearman rhos |

**Used by:** Interactive report Section 5.2 (`prepare_level_comparison_data()`).

---

## 9. Output Files (This Repo)

### `figures/summary_statistics.csv`

Generated by `fig15_summary_table()`. One row per atlas x signature combination.

| Column | Description |
|--------|-------------|
| `Atlas` | Atlas display name (e.g., `Cima`, `Scatlas Normal`) |
| `Signature` | Signature display name (`CytoSig`, `LinCytoSig Best (comb+filt)`, `LinCytoSig`, `SecAct`) |
| `N Targets` | Number of targets with non-null Spearman rho |
| `Median ρ` | Median Spearman rho |
| `Mean ρ` | Mean Spearman rho |
| `Std ρ` | Standard deviation of Spearman rho |
| `Min ρ` | Minimum Spearman rho |
| `Max ρ` | Maximum Spearman rho |
| `% Significant` | Percentage of targets with p < 0.05 |
| `% Positive` | Percentage of targets with rho > 0 |

### `figures/lincytosig_variants_44cytokines.csv`

Per-cytokine Spearman rhos across all 10 methods and 4 atlases.

| Column | Description |
|--------|-------------|
| `Cytokine` | CytoSig target name |
| `CT(comb)` | Best cell-type variant selected via combined GTEx+TCGA |
| `CT(GTEx)` | Best cell-type variant selected via GTEx |
| `CT(TCGA)` | Best cell-type variant selected via TCGA |
| `{Atlas}_CytoSig` | CytoSig Spearman rho for this cytokine in this atlas |
| `{Atlas}_Best(comb,orig)` | Best combined, all genes |
| `{Atlas}_Best(comb,filt)` | Best combined, gene-filtered |
| `{Atlas}_Best(GTEx,orig)` | GTEx-selected, all genes |
| `{Atlas}_Best(TCGA,orig)` | TCGA-selected, all genes |
| `{Atlas}_Best(GTEx,filt)` | GTEx-selected, gene-filtered |
| `{Atlas}_Best(TCGA,filt)` | TCGA-selected, gene-filtered |
| `Mean_*` | Mean across the 4 atlases for each method |

Atlas prefixes: `CIMA_`, `Inflam_`, `scN_` (scAtlas Normal), `scC_` (scAtlas Cancer).
