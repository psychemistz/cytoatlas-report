---
title: "Dataset Analytics: Cross-Sample Correlation Design"
date: "2026-02-13"
status: draft
purpose: "Document each dataset's structure, identify cleaning decisions, and design a reproducible pseudobulk correlation pipeline"
---

# Dataset Analytics: Cross-Sample Correlation Design

## Goal

Compute clean **cross-sample (cross-donor) Spearman correlations** between:
- **X**: Pseudobulk expression of target gene (one value per sample unit)
- **Y**: Predicted activity z-score from ridge regression (one value per sample unit)

Every data-cleaning decision flows from this goal: each sample unit must be independent, identifiable, and comparable.

---

## 1. Dataset Inventory

### Summary Table

| Dataset | Type | Cells | Donors | Samples | Sample unit for correlation | Key issue |
|---------|------|-------|--------|---------|----------------------------|-----------|
| CIMA | scRNA-seq | 6,484,974 | 421 | 421 | donor (1:1) | Clean — one sample per donor |
| Inflammation Main | scRNA-seq | 4,918,140 | ? | 817 | sampleID | No donorID column; donor–sample relationship unknown |
| Inflammation Val | scRNA-seq | 849,922 | ? | 144 | sampleID | Same issue; cell type labels are predicted (Level2pred) |
| Inflammation Ext | scRNA-seq | 572,872 | ? | 86 | sampleID | Same issue; different gene space (37K vs 23K genes) |
| scAtlas Normal | scRNA-seq | 2,293,951 | 317 | 706 | donor × organ | 115 donors have multiple samples from different organs |
| scAtlas Cancer | scRNA-seq | 4,146,975 | 717 | 1,062 | donor × tissue | 82 donors have multiple samples (tumor + adjacent + metastasis) |
| GTEx | Bulk RNA-seq | — | 946 | 19,788 | sample (tissue biopsy) | Multi-tissue per donor (1–39 samples/donor); NOT pseudobulk |
| TCGA | Bulk RNA-seq | — | 10,274 | 11,069 | sample (tumor biopsy) | Mostly 1:1 donor:sample; some donors have normal+tumor pairs |

---

## 2. Per-Dataset Analysis

### 2.1 CIMA (Cell Atlas)

**Source:** Q. Zhang et al., *Nature*, 2024
**File:** `CIMA_RNA_6484974cells_36326genes_compressed.h5ad`

| Property | Value |
|----------|-------|
| Cells | 6,484,974 |
| Genes | 36,326 |
| Donors | 421 |
| Samples | 421 (1:1 with donors) |
| Cell type levels | L1 (6), L2 (27), L3 (38), L4 (72) |
| Population | Healthy donors |
| Metadata | age, sex, BMI, blood biochemistry, metabolomics |

**Obs columns:** `sample`, `cell_type_l1`, `cell_type_l2`, `cell_type_l3`, `cell_type_l4`, `age`, `sex`, `BMI`, ...

**Note:** Raw data contains NO `LowQuality_cells` or `Doublets` — already QC'd upstream. No additional cell exclusion needed.

#### Cell Type Hierarchy

**L1 (6 types) — all present in all 421 donors:**

| L1 type | Total cells | Donors | Cells/donor (mean) |
|---------|------------|--------|-------------------|
| CD4_T | 2,320,103 | 421 | 5,511 |
| CD8_T | 1,363,182 | 421 | 3,238 |
| unconventional_T | 977,267 | 421 | 2,321 |
| Myeloid | 754,770 | 421 | 1,793 |
| B | 566,569 | 420 | 1,349 |
| ILC | 503,083 | 421 | 1,195 |

**L2 (27 types) — top types well-represented, bottom types sparse:**

| L2 type | Cells/donor (mean) | Donors | Status |
|---------|-------------------|--------|--------|
| CD4_naive | 2,720 | 421 | Excellent |
| Mono | 1,535 | 421 | Excellent |
| CD8_naive | 1,389 | 421 | Excellent |
| CD4_helper | 1,328 | 421 | Excellent |
| CD8_CTL | 1,002 | 421 | Excellent |
| ... (12 more types with >100 cells/donor) | | | |
| pDC | 79 | 420 | Adequate |
| CD56_bright_NK | 60 | 418 | Adequate |
| Total_Plasma | 52 | 420 | Marginal |
| MK | 40 | 420 | Marginal |
| Proliferative_T | 39 | 421 | Marginal |
| Immature_T | 24 | 402 | Sparse |
| Proliferative_NK | 14 | 400 | Sparse |
| HSPC | 10 | 411 | Sparse |
| ILC2 | 7 | 406 | Sparse |

#### Pseudobulk Strategy

**Approach:** Compute all 5 levels. Correlation quality will naturally demonstrate where signal degrades.

| Level | Groupby | Pseudobulk groups | Min cells/group | Groups <30 cells |
|-------|---------|-------------------|-----------------|-----------------|
| Donor-only | `sample` | 421 | 3,181 | 0 (0%) |
| Donor × L1 | `sample` × `cell_type_l1` | 2,525 | 45 | 0 (0%) |
| Donor × L2 | `sample` × `cell_type_l2` | 10,251 | 10 | 1,288 (12.6%) |
| Donor × L3 | `sample` × `cell_type_l3` | 14,304 | 10 | 2,156 (15.1%) |
| Donor × L4 | `sample` × `cell_type_l4` | 24,741 | 10 | 5,158 (20.8%) |

**Correlation computation per level:**
- **Donor-only:** One rho per target. Spearman across 421 donors. Fully independent.
- **Donor × celltype (L1-L4):** For each cell type separately, Spearman across donors that have ≥min_cells of that type. Each correlation uses independent donors. One rho per (target × cell type).

**Note on existing processed data:** One donor (CIMA_H365) lacks B cells above threshold at L1, giving 2,525 instead of 2,526 groups. This is correct behavior.

#### Existing Processed Data Verification

**Location:** `/data/parks34/projects/2cytoatlas/results/cross_sample_validation/cima/`

| Level | Pseudobulk | CytoSig (43) | SecAct (1,170) | Shapes match? |
|-------|-----------|-------------|---------------|--------------|
| Donor-only | 421 × 36,326 | 421 × 43 | 421 × 1,170 | Yes |
| Donor × L1 | 2,525 × 36,326 | 2,525 × 43 | 2,525 × 1,170 | Yes |
| Donor × L2 | 10,251 × 36,326 | 10,251 × 43 | 10,251 × 1,170 | Yes |
| Donor × L3 | 14,304 × 36,326 | 14,304 × 43 | 14,304 × 1,170 | Yes |
| Donor × L4 | 24,741 × 36,326 | 24,741 × 43 | 24,741 × 1,170 | Yes |

Resampled (bootstrap) versions also exist for L1-L4.

**Min cells at generation:** None (`01_cima_activity.py` includes all groups). **Min cells at validation:** 10 (`11_donor_level_pipeline.py` filters at correlation time). The "Min cells/group" column above shows the smallest observed group in the data, not a threshold. L2-L4 include some small pseudobulk profiles (10-29 cells) whose noise will be visible in correlation distributions.

#### Decision: CONFIRMED (with caveat on thresholds)

- Use all 5 levels (donor-only + L1-L4) with all 3 signature types
- For correlation: within each cell type, correlate across donors (independent)
- L3/L4 noise will be apparent in the results, demonstrating the resolution limit naturally
- **Threshold note:** The original activity script (`01_cima_activity.py`) has **no min_cells filter** — all sample×celltype groups are included regardless of size. The downstream validation pipeline (`11_donor_level_pipeline.py`, `validation/`) applies min_cells=10 at the correlation stage. This is acceptable: the pseudobulk H5ADs preserve all groups, and filtering happens at analysis time. See Section 3.8 for cross-pipeline threshold reconciliation.

---

### 2.2 Inflammation Atlas

**Source:** Inflammatory Atlas consortium
**Files:**
- Main: `INFLAMMATION_ATLAS_main_afterQC.h5ad` (4,918,140 cells, 817 samples)
- Validation: `INFLAMMATION_ATLAS_validation_afterQC.h5ad` (849,922 cells, 144 samples)
- External: `INFLAMMATION_ATLAS_external_afterQC.h5ad` (572,872 cells, 86 samples)
- **Combined:** 6,340,934 cells, 1,047 samples

| Property | Main | Validation | External |
|----------|------|------------|----------|
| Cells | 4,918,140 | 849,922 | 572,872 |
| Samples (sampleID) | 817 | 144 | 86 |
| Genes | 22,838 | 22,838 | 37,169 |
| Cell type L1 | 17 (Level1) | 15 (Level1pred) | 15 (Level1pred) |
| Cell type L2 | 66 (Level2) | 53 (Level2pred) | 54 (Level2pred) |

**Obs columns:** `libraryID`, `sampleID`, `Level1`/`Level1pred`, `Level2`/`Level2pred`

**Gene IDs:** Ensembl IDs (ENSG*) as var_names, HGNC symbols in `var['symbol']` (22,826 unique symbols, 0 NaN). Gene matching for activity inference and correlation uses `var['symbol']`.

#### Cohort Overlap Analysis

| Metric | Main ↔ Val | Main ↔ Ext | Val ↔ Ext |
|--------|-----------|-----------|----------|
| Shared sampleIDs | 0 | 0 | 0 |
| Shared libraryIDs | 217 | 0 | 0 |
| Shared donors | 0 | 0 | 0 |

The 217 shared libraryIDs between Main and Val come from multiplexed 10x experiments (CellPlex/hashtag). The same physical library contains cells from multiple donors, but different donors were assigned to Main vs Val during the cohort split. Zero shared cells or donors — the shared libraryIDs mean shared sequencing run only (mild batch correlation, not data duplication).

#### Decision: Analyze Separately (Option C — CONFIRMED)

**Main (817 samples)** as primary for cross-donor correlation analysis:
- Expert-curated cell type annotations (Level1, Level2)
- 22,838 genes (shared with Val)
- Sufficient power: 817 samples, 30,227 pseudobulk groups at L2

**Validation (144 samples)** as internal replication:
- Predicted annotations (Level2pred) — noisier
- Same gene space as Main
- Checks whether correlation patterns hold

**External (86 samples)** as independent validation:
- Predicted annotations
- Different gene space (37,169 vs 22,838)
- Smallest cohort, truly independent

Report clearly as "samples" not "donors" until donor identity is resolved.

#### Critical Issue: No `donorID` Column

- Only `sampleID` is available
- Unknown whether multiple sampleIDs can map to the same donor (e.g., pre/post treatment, multiple timepoints)
- If some donors contributed multiple samples, treating each sampleID as independent inflates correlation N and may introduce within-donor correlations
- **Mitigation:** Treat each sampleID as an independent unit; report N as "samples" not "donors"

#### Cell Exclusion Issue (REQUIRES REGENERATION)

The Level1 annotation in the raw data contains non-cell categories that **must be excluded** before pseudobulk aggregation:

**Main dataset Level1 types (17 total):**

| Type | Status | Notes |
|------|--------|-------|
| B, DC, ILC, Mono, Plasma, Platelets, Progenitors, RBC | Keep | Real cell types |
| T_CD4_Naive, T_CD4_NonNaive, T_CD8_Naive, T_CD8_NonNaive, UTC | Keep | Real cell types |
| Cycling_cells, pDC | Keep | Real cell types |
| **Doublets** | **EXCLUDE** | Artifact: 780 sample-groups at L1 |
| **LowQuality_cells** | **EXCLUDE** | Artifact: 776 sample-groups at L1 |

**Level2 also contains:**
- `Doublets` — exclude
- `LowQuality_cells` — exclude
- `NK_lowRibocontent` — keep (real cell type with quality caveat, 956 correlation rows)

**Impact of NOT excluding (current state):**
- **donor_only level:** Doublet and LowQuality cells are **mixed into every donor's pseudobulk profile**, contaminating expression data. Effect is moderate (~1-5% of cells per donor) but systematic.
- **L1/L2 levels:** They form their own cell type strata. Correlations computed for "Doublets" and "LowQuality_cells" cell types are meaningless (1,922 correlation rows each in `inflammation_main_correlations.csv`).

#### Existing Processed Data Verification

**Location:** `/data/parks34/projects/2cytoatlas/results/cross_sample_validation/inflammation_main/`

| Level | Pseudobulk | CytoSig (43) | LinCytoSig (178) | SecAct (1,170) | Match? |
|-------|-----------|-------------|-----------------|---------------|--------|
| donor_only | 817 × 22,838 | 817 × 43 | 817 × 178 | 817 × 1,170 | Yes |
| donor × L1 | 9,771 × 22,838 | 9,771 × 43 | 9,771 × 178 | 9,771 × 1,170 | Yes |
| donor × L2 | 28,671 × 22,838 | 28,671 × 43 | 28,671 × 178 | 28,671 × 1,170 | Yes |

**Correlation statistics (AFTER cell exclusion + ridge_batch + min_samples=10, regenerated 2026-02-13):**

| Level | Signature | Median ρ | N targets | n_samples range |
|-------|-----------|---------|-----------|----------------|
| donor_only | CytoSig | 0.323 | 33 | 817 |
| donor_only | LinCytoSig | 0.178 | 123 | 817 |
| donor_only | SecAct | 0.173 | 805 | 817 |
| donor_l1 | CytoSig | 0.079 | 33 | — |
| donor_l1 | LinCytoSig | 0.067 | 123 | — |
| donor_l1 | SecAct | 0.075 | 805 | — |
| donor_l2 | CytoSig | 0.037 | 33 | — |
| donor_l2 | LinCytoSig | 0.048 | 123 | — |
| donor_l2 | SecAct | 0.040 | 805 | — |

**Changes from regeneration:**
- L1: 15 cell types (was 17 — Doublets, LowQuality_cells removed). 9,771 groups (was 11,327).
- L2: 63 cell types (was 66). 28,671 groups (was 30,227).
- Total correlation rows: 72,358 (reduced from ~76K — artifact cell type rows removed).
- donor_only CytoSig median ρ: 0.321 → 0.323 (slight increase from cleaner profiles).
- donor_only SecAct median ρ: 0.188 → 0.173 (shifted slightly with ridge_batch).
- L2 min n_samples now 10 (was 5) — min_samples=10 enforced consistently.
- Activity inferred with `ridge_batch` (GPU-accelerated CuPy) instead of `ridge`.

**Note:** Inflammation donor_only median ρ (0.323 CytoSig) is substantially higher than CIMA (0.114 CytoSig). This is expected: inflammation samples span healthy + 20 disease conditions, creating much wider biological variance that inflates cross-sample correlations. CIMA is healthy-only with less variance.

#### Comparison with CIMA (Gold Standard)

| Property | CIMA | Inflammation Main |
|----------|------|------------------|
| Sample unit | donor (1:1) | sampleID (donor mapping unknown) |
| Donor/Sample count | 421 | 817 |
| Gene ID format | HGNC symbols (var_names) | Ensembl IDs (var_names), symbols in var['symbol'] |
| Cell exclusion needed | No (pre-QC'd upstream) | **Yes (Doublets, LowQuality_cells)** |
| Annotation quality | Curated (L1–L4) | Curated L1/L2 (Main); predicted (Val/Ext) |
| Cell type levels | L1(6), L2(27), L3(38), L4(72) | L1(17→15 after exclusion), L2(66→63) |
| Population | Healthy only | Healthy + 20 disease conditions |

#### Pipeline Comparison: 02_inflam_activity.py vs 11_donor_level_pipeline.py

Two separate scripts generate inflammation pseudobulk with **different normalization approaches**:

| Property | `02_inflam_activity.py` | `11_donor_level_pipeline.py` |
|----------|------------------------|------------------------------|
| Output location | `results/inflammation/` | `results/cross_sample_validation/inflammation_main/` |
| Normalization | Sum raw counts → CPM → log2(x+1) | Per-cell CPM → log1p → mean across cells |
| Min cells filter | None (all groups) | 10 (at generation) |
| Cell exclusion | None | Doublets + LowQuality_cells via `exclude_celltypes` config (**FIXED**) |
| Used for | Disease differential, treatment response | **Cross-donor correlation validation** |

The cross-sample validation data (from `11_donor_level_pipeline.py`) is what we use for cross-donor correlations. Both CIMA and Inflammation were processed through the same pipeline, ensuring internal consistency.

**Normalization note:** The mean-of-log approach (averaging per-cell log1p(CPM) values) differs from the standard pseudobulk approach (summing raw counts, then normalizing). For Spearman correlation (rank-based), both approaches produce valid rankings. Since all single-cell atlases (CIMA, Inflammation, scAtlas) use the same 11_donor_level_pipeline.py, they are internally consistent.

---

### 2.3 scAtlas Normal

**Source:** Q. Shi et al., *Nature*, 2025
**File:** `igt_s9_fine_counts.h5ad`

| Property | Value |
|----------|-------|
| Cells | 2,293,951 |
| Genes | 21,812 |
| Donors | 317 |
| Samples | 706 |
| Organ systems | 10 |
| Tissues | 35 |
| Cell type levels | majorCluster (8), subCluster (102), cellType1 (468), cellType2 (1,069) |
| Population | Healthy / normal tissue |

**Obs columns:** `donorID`, `sampleID`, `system`, `tissue`, `region`, `majorCluster`, `subCluster`, `cellType1`, `cellType2`, `sex`, `age`, ...

**Multi-sample donors:**
- 202 donors: 1 sample (63.7%)
- 53 donors: 2 samples (16.7%)
- 62 donors: 3+ samples (19.6%, up to 27 samples from TSP14 — blood, liver, lung, colon, skin, thymus, etc.)

**Top tissues by donor count:**
Breast (124), Lung (97), Colon (65), Heart (52), Liver (43), SmallIntestine (30), Spleen (30), Ovary (28), Thymus (23)

#### 2.3.1 Cell Type Annotation Hierarchy

Four annotation columns with very different quality and granularity:

| Column | Unique values | Quality | Use case |
|--------|--------------|---------|----------|
| `majorCluster` | 8 | Clean, standardized | Coarse grouping (Epithelial, Stromal, Myeloid, CD8T, CD4T, Endothelial, B, ILC) |
| `subCluster` | 102 | Clean, structured naming (`CD8T02_Tem_GZMK`, `M07_Mph_FCN1`, `Epi_Breast`) | **Recommended for cell-type analysis** |
| `cellType1` | 468 | **Inconsistent** — 25 case-duplicate groups, 7+ singular/plural pairs | Legacy; do not use without normalization |
| `cellType2` | 1,069 | **Inconsistent** — 50 case-duplicate groups | Most granular; too sparse for correlation |

#### 2.3.2 Cell Type Annotation Inconsistencies

Full audit of `cellType1` problems with cell counts:

**Case duplicates** (25 groups, major examples):
- "Myeloid" (134K) vs "myeloid" (18K)
- "Fibroblasts" (157K) vs "fibroblasts" (11K)
- "Macrophage" (28K) vs "macrophage" (1.2K)
- "NK cell" (7.6K) vs "nk cell" (varies by organ)

**Singular/plural duplicates** (7+ pairs):
- "t cell" (11.7K) vs "t cells" (83.8K)
- "fibroblast" (31.3K) vs "fibroblasts" (168.5K)
- "macrophage" (29.2K) vs "macrophages" (5.4K)
- "b cell" (18.3K) vs "b cells" (5.7K)

**Decision:** Use `subCluster` (102 values, already standardized) for all cell-type-stratified analyses. The `subCluster` column uses structured naming (`PREFIX##_Subtype_MARKER`) and has zero case/plural issues. Already adopted by `03_scatlas_analysis.py`.

#### 2.3.3 Cell Type Distribution Across Organs

Organ-specificity spectrum of `cellType1` (after case normalization):

| Category | Cell types | Percentage |
|----------|-----------|------------|
| Organ-specific (1 tissue only) | 333 | 73.5% |
| Local (2–5 tissues) | 75 | 16.0% |
| Regional (6–10 tissues) | 28 | 6.0% |
| Ubiquitous (11+ tissues) | 17 | 3.6% |

**Top ubiquitous types:**

| cellType1 (normalized) | Tissues |
|------------------------|---------|
| macrophage | 20 |
| fibroblast | 15 |
| Monocyte | 15 |
| NK cell | 15 |
| T cell | 15 |
| plasma cell | 14 |
| endothelial cell | 14 |

**Key insight:** Most immune/stromal `subCluster`s (M05_Mo_CD14, S01_Fb_PI16, CD8T02_Tem_GZMK, etc.) appear across many organs. Epithelial `subCluster`s (Epi_Breast, Epi_Lung, etc.) are organ-specific by definition. Cross-organ analysis is possible for immune/stromal cells but not epithelial.

#### 2.3.4 Donor × Tissue Sparsity Analysis

**Pseudobulk quality at donor × tissue level:**
- 706 donor × tissue combinations with meaningful cells
- Min cells per group: 30, Median: 2,364, Max: 28,619
- All 706 groups are adequate for pseudobulk (no min_cells filtering needed at this level)

**Per-tissue donor counts** (35 tissues):

| Tier | Criteria | Tissues | Donor counts |
|------|----------|---------|-------------|
| A | ≥30 donors | 7 | Breast 124, Lung 97, Colon 65, Heart 52, Liver 43, SmallIntestine 30, Spleen 30 |
| B | 20–29 donors | 5 | Ovary 28, Thymus 23, Kidney 22, Blood 21, Skin 20 |
| C | 10–19 donors | 7 | Stomach 18, BoneMarrow 16, Eye 14, Pancreas 13, Esophagus 12, Brain 11, Uterus 10 |
| D | <10 donors | 16 | Prostate 8, Bladder 7, Gallbladder 6, Trachea 5, Testis 5, ... (down to 1) |

**Donor independence:**
- 202 single-tissue donors (63.7%)
- 53 donors with 2 tissues (16.7%)
- 62 donors with 3+ tissues (19.6%, up to 27 tissues)
- Multi-tissue donors create non-independence across per-tissue analyses (same issue as GTEx)

**Threshold sweep for Level 1 (by_organ):**

| min_samples | Tissues passing | Total samples | Tiers included |
|-------------|-----------------|---------------|----------------|
| ≥30 | 7 | 441 | A only |
| ≥25 | 8 | 469 | A + Ovary |
| ≥20 | 12 | 556 | A + B |
| ≥15 | 14 | 590 | A + B + BoneMarrow, Eye |
| ≥10 | 19 | 652 | A + B + C |

**Organ × subCluster sparsity** (for Level 3 feasibility):

Of 3,570 total tissue × subCluster combinations (any cells), viable combinations with ≥N donors having ≥10 cells per donor:

| Donor threshold | Viable combos |
|-----------------|---------------|
| ≥10 donors | 200 |
| ≥20 donors | 74 |
| ≥30 donors | 39 |

**Per-tissue subCluster viability** (donors with ≥10 cells, min_cells=10):

| Tissue (Tier A) | subClusters ≥10 donors | ≥20 donors | ≥30 donors |
|-----------------|----------------------|------------|------------|
| Lung | 38 | 26 | 16 |
| Spleen | 34 | 8 | 0 |
| Breast | 32 | 26 | 23 |
| Liver | 21 | 7 | 0 |
| Colon | 5 | 0 | 0 |
| SmallIntestine | 3 | 0 | 0 |
| Heart | 0 | 0 | 0 |

| Tissue (Tier B) | subClusters ≥10 donors | ≥20 donors | ≥30 donors |
|-----------------|----------------------|------------|------------|
| Skin | 10 | 0 | 0 |
| Blood | 4 | 0 | 0 |
| Ovary | 0 | 0 | 0 |
| Thymus | 0 | 0 | 0 |
| Kidney | 0 | 0 | 0 |

**Level 3 conclusion:** Only **Breast** (26 subClusters) and **Lung** (26) are viable at min_samples=20. At min_samples=10, Spleen (34) and Liver (21) become exploratory candidates. Heart, Colon, SmallIntestine, and all Tier B tissues lack subCluster depth entirely.

#### 2.3.5 Correlation Strategy Design

**The GTEx analogy:** scAtlas Normal is structurally analogous to GTEx — multiple tissues per donor, need per-tissue stratification for independence. The key difference: scAtlas has cell-type resolution that GTEx lacks.

**The CIMA analogy:** CIMA uses per-celltype correlation across donors, but all cells are PBMC (one compartment). scAtlas has 35 organs, so "macrophage" in Liver ≠ "macrophage" in Lung biologically.

**Multi-level strategy** (mirroring GTEx/TCGA three-level design):

| Level | Unit | N | Mean-centering | Independence | Status | Analogy |
|-------|------|---|----------------|--------------|--------|---------|
| 0. donor_only | donor | 317 | Global | Fully independent | Supplementary | GTEx pooled |
| 1. by_organ (PRIMARY) | per-organ donors | 20–124 per organ | Within-organ | Fully independent within organ | **Primary** | **GTEx by_tissue** |
| 2. donor_organ | donor × tissue | ~706 | Global | Partial (115 multi-tissue) | Supplementary | TCGA primary_only |
| 3. by_organ_subCluster | per-organ per-subCluster | Breast/Lung only | Within-organ | Fully independent | Exploratory | CIMA per-celltype |

**Level 0 (donor_only):** Sum all cells per donor → 317 profiles. Fully independent but mixes tissue biology. Donors with different tissue coverage are not comparable. Supplementary with explicit caveat. Pipeline support: `donor_col='donorID'` in `12_cross_sample_correlation.py` config (sampleID is per donor×organ at 706; donorID aggregates to 317 true donors). Generated output: pseudobulk 317 × 21,812 (n_cells range 48–145,977, mean 7,236), activity H5ADs for cytosig (43), lincytosig (178), secact (1,170).

**Level 1 (by_organ) — PRIMARY:** For each tissue with ≥20 donors (12 tissues, Tier A + B), aggregate all cells per donor within that tissue, run within-organ mean-centered ridge regression, correlate across donors. This is the direct analogue of GTEx by_tissue. Each donor contributes at most one point per tissue → fully independent within each tissue correlation. Report Tier A (≥30 donors, 7 tissues) as high-confidence and Tier B (20–29 donors, 5 tissues) as adequate-confidence separately.

**Level 2 (donor_organ):** All 706 donor × tissue profiles with global mean-centering. 115 multi-tissue donors create non-independence (same issue as GTEx pooled). Useful as a sensitivity check.

**Level 3 (by_organ_subCluster) — Exploratory:** Feasible only for Breast (26 subClusters at ≥20 donors) and Lung (26 subClusters). At min_samples=10, Spleen (34 subClusters) and Liver (21) become additional candidates. Heart, Colon, SmallIntestine, and all Tier B tissues lack subCluster depth. Not primary analysis, but demonstrates where cell-type resolution adds or loses signal within the best-sampled tissues.

**Threshold rationale:**
- min_samples=20 for Level 1 — at N=20, Spearman has reasonable power to detect rho ≥ 0.45 (p < 0.05); gains 5 tissues (+26% samples) over strict ≥30 cutoff while avoiding the wide confidence intervals below N=20
- min_samples=20 for Level 3 primary (Breast/Lung), min_samples=10 for Level 3 exploratory (Spleen/Liver)
- No min_cells filter needed at donor × tissue level (all 706 groups have ≥30 cells)

**Pipeline implementation:** `12_cross_sample_correlation.py` generates `donor_organ` pseudobulk with within-tissue mean-centered ridge regression. Level 0 (donor_only) now generated via `donor_col='donorID'` config (required because `sampleID` is per donor×organ, not per donor). Level 1 (by_organ) requires only adding per-tissue stratified correlation to `12_`/`13_` — no new pseudobulk regeneration needed. This mirrors the pattern in `15_bulk_validation.py` (GTEx by_tissue).

#### 2.3.6 Cross-Platform Tissue Comparison

Of the 12 scAtlas tissues passing the min_samples=20 threshold, **11 have direct GTEx counterparts** (Thymus is scAtlas-only):

| scAtlas Normal | GTEx (SMTS) | scAtlas donors | GTEx samples | Match |
|----------------|-------------|----------------|--------------|-------|
| Breast | Breast | 124 | 514 | Exact |
| Lung | Lung | 97 | 604 | Exact |
| Colon | Colon | 65 | 925 | Exact |
| Heart | Heart | 52 | 913 | Exact |
| Liver | Liver | 43 | 282 | Exact |
| SmallIntestine | Small Intestine | 30 | 226 | Semantic |
| Spleen | Spleen | 30 | 277 | Exact |
| Ovary | Ovary | 28 | 193 | Exact |
| Kidney | Kidney | 22 | 115 | Exact |
| Blood | Blood | 21 | 1,130 | Exact |
| Skin | Skin | 20 | 2,057 | Exact |
| Thymus | — | 23 | — | No GTEx match |

Per-tissue correlations can be directly compared: same tissue, same biological question, independent data, different technology (pseudobulk scRNA-seq vs bulk RNA-seq). This gives 11 independent tissue-level comparison points — a strong cross-platform validation opportunity.

---

### 2.4 scAtlas Cancer

**Source:** Q. Shi et al., *Nature*, 2025
**File:** `PanCancer_igt_s9_fine_counts.h5ad`

| Property | Value |
|----------|-------|
| Cells | 4,146,975 |
| Genes | 21,812 |
| Donors | 717 |
| Samples | 1,062 |
| Cancer types | 29 |
| Tissue types | Tumor, Adjacent, Metastasis, Blood, PreLesion, PleuralFluids |
| Cell type levels | cellType1 (162), cellType2 (153) |

**Obs columns:** `donorID`, `sampleID`, `cancerType`, `sub_cancerType`, `tissue`, `treatment`, `treatmentResponse`, `treatmentPhase`, `cellType1`, `cellType2`, ...

**Multi-sample donors (82 donors with >2 samples):**
- All 82 have the SAME cancer type across samples
- Multiple samples come from: tumor vs adjacent normal, multiple tumor regions, blood, metastasis
- Example: CRC patients with 5-8 biopsies (multiple tumor regions + 1 adjacent)
- Example: TNBC patients with blood + metastasis at pre/post/progression timepoints

**Top cancer types by donor count:**
HCC (88), PAAD (75), CRC (65), ESCA (60), LUAD (59), BRCA (49), HNSC (46), NPC (39), STAD (39), KIRC (35), ICC (30)

**Tissue breakdown for donors with ≤2 samples (635 donors):**
- Tumor only: 438
- Unlabeled: 82
- Tumor + Adjacent: 63
- Metastasis only: 53
- PreLesion only: 38
- Tumor + Blood: 22
- Others: 12

**Cleaning decisions:**

1. **Tissue type filtering:** Should we include adjacent normal, metastasis, blood, pre-lesion samples alongside tumor samples in the same correlation?
   - **Problem:** Mixing tumor and adjacent normal samples from the SAME donor creates within-donor pairs that are not independent but are biologically very different
   - **Recommendation:** Stratify by tissue type. Primary analysis on Tumor-only. Separate analyses for Adjacent, Metastasis.

2. **Cancer-type stratification:**
   - **Option A:** Pool all cancer types (N=717 donors, higher power, but heterogeneous)
   - **Option B:** Per-cancer-type correlations (HCC: N=88, PAAD: N=75, etc.)
   - **Recommendation:** Both. Pooled as primary (more power), per-cancer as sensitivity analysis. Exclude cancer types with <20 donors.

3. **Multi-region tumor samples:**
   - For donors with multiple tumor biopsies (e.g., CRC with 5-7 tumor samples): average across tumor samples per donor, or pick one
   - **Recommendation:** Average across same-tissue samples per donor → one value per (donor × tissue_type × cancer_type)

---

### 2.5 GTEx (Bulk RNA-seq)

**Source:** GTEx Consortium, v11
**File:** `GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.parquet` (4.1 GB)

| Property | Value |
|----------|-------|
| Samples | 19,788 |
| Donors | 946 |
| Genes (raw) | 74,628 ENSG IDs (all GENCODE features: protein-coding + lncRNA + pseudogenes + small RNA) |
| Genes (mapped) | 72,930 unique HGNC symbols (after dedup; ~20K protein-coding, rest non-coding/pseudo) |
| Tissues (SMTS) | 30 |
| Tissue subtypes (SMTSD) | 68 |
| Data type | TPM (Transcripts Per Million) |
| Transform | log2(TPM + 1) |

**This is true bulk RNA-seq, NOT pseudobulk.** No cell-type information.

#### Expression Value Characterization

The file contains **TPM** values computed by RNASeQC v2.4.3.

**Evidence (in log2(TPM+1) space after transform):**
- Value range: 0.0 to 19.4 (corresponding to 0–697K TPM)
- Mean: ~0.87, median: 0.0
- **51% of values are zero** — expected for tissue-specific gene expression (many genes silent in any given tissue)
- **No negative values** — clean TPM, no batch correction artifacts (unlike TCGA RSEM)

**Pipeline handling (`15_bulk_validation.py`):**
- No preprocessing needed (TPM is always ≥ 0)
- Applies `log2(TPM + 1)` directly
- Data format tracked in output H5AD metadata (`uns['data_format']`, `uns['transform']`)

#### Gene ID Mapping

**Source gene IDs:** Versioned ENSG (e.g., `ENSG00000223972.6`)
**Probemap:** GENCODE v23 (`gencode.v23.annotation.gene.probemap`, 60,498 entries → 58,581 unique gene symbols)

**3-tier mapping cascade** (in `_map_ensg_to_symbols()`):
1. Probemap lookup with full versioned ENSG ID
2. Probemap lookup with base ENSG (strip `.version`)
3. Fallback to `Description` column in parquet (contains HGNC symbols directly)

**Mapping results:**
- 25,687 mapped via probemap (versioned ENSG match)
- 26,841 mapped via probemap (base ENSG, version stripped)
- 22,100 mapped via Description column fallback (newer genes not in GENCODE v23 probemap)
- **0 unmapped** — 100% mapping rate (unlike a probemap-only approach which would lose 22K genes)
- 74,628 mapped rows → 72,930 unique symbols after averaging 231 duplicate symbols

**Why 72,930 genes and not ~20K?** The raw GTEx v11 parquet contains ALL GENCODE features (protein-coding + lncRNA + pseudogenes + small RNAs). Only ~20K are protein-coding. The non-coding genes are harmlessly carried through — they get silently excluded during the gene intersection step of ridge regression, which keeps only genes present in both expression data AND signature matrix. Effective gene counts used for inference: CytoSig 4,859, SecAct 7,451.

**Note:** The download script (`15a_download_bulk_data.sh`) downloads GENCODE v23 probemap and GTEx **v8** sample attributes, but the pipeline uses GTEx **v11** parquet and v11 sample attributes (obtained separately from GTEx Portal). The GENCODE v23 probemap remains valid for v11 gene mapping.

#### Donor Independence

| Samples/donor | Donors | Cumulative % |
|---------------|--------|-------------|
| 1–5 | 12 | 1.3% |
| 6–10 | 48 | 6.3% |
| 11–15 | 131 | 20.2% |
| 16–20 | 230 | 44.5% |
| 21–25 | 267 | 72.7% |
| 26–30 | 143 | 87.8% |
| 31–35 | 67 | 94.9% |
| 36–39 | 48 | 100% |

- **All 946 donors contribute multiple tissues** (only 1 donor has a single sample)
- Median: ~20 samples per donor
- Maximum: 39 samples (donor GTEX-15ER7)
- **Every correlation in the pooled analysis is inflated** by within-donor tissue correlation

#### Tissue Breakdown (30 tissues, SMTS level)

| Tissue | Samples | Tissue | Samples |
|--------|---------|--------|---------|
| Brain | 3,234 | Stomach | 491 |
| Skin | 2,057 | Testis | 414 |
| Esophagus | 1,578 | Pancreas | 384 |
| Blood Vessel | 1,431 | Pituitary | 313 |
| Adipose Tissue | 1,301 | Adrenal Gland | 295 |
| Blood | 1,130 | Liver | 282 |
| Colon | 925 | Prostate | 282 |
| Heart | 913 | Spleen | 277 |
| Muscle | 818 | Small Intestine | 226 |
| Thyroid | 684 | Ovary | 193 |
| Nerve | 670 | Salivary Gland | 181 |
| Lung | 604 | Vagina | 170 |
| Breast | 514 | Uterus | 153 |
| | | Kidney | 115 |
| | | Bladder | 77 |
| | | Cervix Uteri | 47 |
| | | Fallopian Tube | 29 |

All 30 tissues have ≥29 samples. At min_samples=30, Fallopian Tube is excluded → **29 tissues in stratified analysis**.

**Tissue detail (SMTSD level):** 68 subtypes (e.g., Brain splits into 13 regions: Cortex, Cerebellum, Hippocampus, etc.). The by_tissue stratification uses SMTS (coarse), not SMTSD.

#### Decision: Compute Two Levels

**Level 1 — All samples pooled (donor_only):**
- N = 19,788 samples from 946 donors
- **Severe non-independence:** every donor contributes ~20 samples across tissues
- Exists: `gtex_donor_only_*` ✓
- Useful as a sensitivity check but NOT a valid independent correlation

**Level 2 — Per-tissue (by_tissue):**
- Within each of 29 tissues (≥30 samples), correlate across donors
- Each donor contributes at most one sample per tissue → **truly independent within tissue**
- N per tissue: 47–3,234 (after dropping Fallopian Tube)
- Exists: `gtex_by_tissue_*` ✓

**Recommendation:** Level 2 (per-tissue) as primary — biologically meaningful and statistically valid. Level 1 (pooled) as supplementary with explicit caveat about non-independence.

#### Existing Processed Data Verification

**Location:** `/data/parks34/projects/2cytoatlas/results/cross_sample_validation/gtex/`

| Level | Expression | CytoSig (43) | LinCytoSig (178) | SecAct (1,170) | Match? |
|-------|-----------|-------------|-----------------|---------------|--------|
| donor_only (all) | 19,788 × 72,930 | 19,788 × 43 | 19,788 × 178 | 19,788 × 1,170 | Yes |
| by_tissue (29 tissues) | 19,759 × 72,930 | 19,759 × 43 | 19,759 × 178 | 19,759 × 1,170 | Yes |

**Gene coverage per signature:**

| Signature | Common genes | Total sig genes | Coverage |
|-----------|-------------|----------------|---------|
| CytoSig | 4,859 | 4,881 | 99.5% |
| LinCytoSig | 18,639 | 19,918 | 93.6% |
| SecAct | 7,451 | 7,919 | 94.1% |

**Note on by_tissue sample count:** 19,759 = 19,788 − 29 (Fallopian Tube samples excluded by min_samples=30 filter).

**Metadata:** Regenerated 2026-02-13 with `data_format='TPM'` and `transform='log2(TPM+1)'` in `.uns`.

#### Correlation Approach Comparison

`gtex_correlations.csv` (regenerated 2026-02-13 with `data_format` metadata):

| Approach | N | Median rho (CytoSig) | Median rho (SecAct) | Independence | Notes |
|----------|---|---------------------|--------------------|--------------|----|
| donor_only (pooled) | 19,788 | 0.211 | 0.394 | **No** | Cross-tissue variation inflates rho |
| by_tissue "all" | 19,759 | 0.105 | 0.192 | **No** | Within-tissue centered, same donors across tissues |

**Why donor_only rho is inflated:** Tissue differences dominate. Liver samples have high albumin expression AND high albumin activity; brain samples have low both. This cross-tissue variation creates strong correlations that don't reflect within-tissue donor variation.

**Why by_tissue "all" is different from donor_only:** Activity was inferred with within-tissue mean centering (each tissue group centered separately before ridge regression), removing tissue-level variation. The "all" row then correlates across all samples using these tissue-adjusted values.

**Proper approach: per-tissue correlations.**
- Donors are truly independent within each tissue (at most one sample per donor per tissue)
- Each ridge regression is mean-centered within tissue, asking: "Within lung, does the donor with higher X expression also show higher X activity?"
- Report median rho across tissues, with per-tissue rho distributions

**Self-referential consideration:** Both expression and activity come from the same matrix. However, the activity z-score is derived from ridge regression across the full signature matrix (~4,859 genes for CytoSig). The target gene contributes only 1/4,859 of the signal, so the self-referential bias is minimal.

---

### 2.6 TCGA (Bulk RNA-seq)

**Source:** TCGA PanCancer
**File:** `EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv`

| Property | Value |
|----------|-------|
| Total samples | 11,069 |
| Unique donors | 10,274 |
| Genes | 20,501 (after symbol\|entrezID parsing, dedup by averaging) |
| Cancer types | 34 (33 named + 1 "Unknown") |
| Data type | EBPlusPlus batch-adjusted RSEM normalized counts (see below) |
| Transform | clip(0) → log2(RSEM + 1) |

**True bulk RNA-seq, NOT pseudobulk.** No cell-type information.

#### Expression Value Characterization

The file contains **EBPlusPlus batch-adjusted RSEM normalized counts** — not TPM, not FPKM, not raw counts.

**Evidence:**
- Column sums vary widely (759K to 16.4M per sample); TPM would sum to exactly 1M
- Values range from -0.89 to ~717K (mean ~990, median ~56)
- **0.55% of values are negative** — impossible for raw counts or TPM; confirms EBPlusPlus batch correction artifact introduces small negatives

**Pipeline handling (`15_bulk_validation.py`):**
- Clips negative values to 0 before log transform (removes batch correction artifacts)
- Applies `log2(x + 1)` on clipped values
- Logs the count and percentage of clipped values for traceability
- Data format tracked in output H5AD metadata (`uns['data_format']`, `uns['transform']`)

**Gene ID mapping:** `symbol|entrezID` format (e.g., `TP53|7157`). Pipeline splits on `|`, takes the symbol. No probemap needed. Genes with `?` symbol are dropped. Duplicate symbols averaged via `groupby().mean()`.

#### Sample Type Breakdown

| Code | Type | Count | Notes |
|------|------|-------|-------|
| 01 | Primary Tumor | 9,706 | Main analysis target |
| 11 | Solid Tissue Normal | 737 | Matched normals — same donor as tumor |
| 06 | Metastatic | 395 | Different biology from primary |
| 03 | Primary Blood Cancer | 173 | AML — no solid tissue |
| 02 | Recurrent Tumor | 46 | Same donor as primary |
| 05 | Additional Primary | 11 | Rare |
| 07 | Additional Metastatic | 1 | Rare |

#### Donor Independence

| Samples/donor | Donors | Notes |
|---------------|--------|-------|
| 1 | 9,489 | Fully independent |
| 2 | 775 | Most are tumor + matched normal |
| 3 | 10 | Tumor + normal + recurrent/metastatic |

- 712 donors have BOTH tumor (01) and normal (11) → not independent
- Filtering to Primary Tumor (01) + Blood Cancer (03) gives 9,879 samples from 9,875 donors (~1:1)

#### Cancer Types (33 named, by_cancer level)

| Cancer type | Samples | Cancer type | Samples |
|-------------|---------|-------------|---------|
| Breast Invasive Carcinoma | 1,212 | Pheochromocytoma & Paraganglioma | 185 |
| Kidney Clear Cell Carcinoma | 603 | Pancreatic Adenocarcinoma | 183 |
| Lung Adenocarcinoma | 573 | Glioblastoma Multiforme | 172 |
| Thyroid Carcinoma | 571 | Acute Myeloid Leukemia | 167 |
| Head & Neck Squamous Cell Carcinoma | 564 | Testicular Germ Cell Tumor | 154 |
| Lung Squamous Cell Carcinoma | 549 | Thymoma | 121 |
| Prostate Adenocarcinoma | 548 | Rectum Adenocarcinoma | 103 |
| Brain Lower Grade Glioma | 527 | Kidney Chromophobe | 91 |
| Skin Cutaneous Melanoma | 470 | Mesothelioma | 87 |
| Stomach Adenocarcinoma | 447 | Uveal Melanoma | 79 |
| Bladder Urothelial Carcinoma | 426 | Adrenocortical Cancer | 77 |
| Liver Hepatocellular Carcinoma | 422 | Uterine Carcinosarcoma | 57 |
| Colon Adenocarcinoma | 331 | Diffuse Large B-Cell Lymphoma | 47 |
| Kidney Papillary Cell Carcinoma | 321 | Cholangiocarcinoma | 45 |
| Cervical & Endocervical Cancer | 309 | | |
| Ovarian Serous Cystadenocarcinoma | 307 | | |
| Sarcoma | 264 | | |
| Uterine Corpus Endometrioid Carcinoma | 204 | | |
| Esophageal Carcinoma | 193 | | |

All 33 named cancer types have ≥30 samples in donor_only (all sample types). After filtering to primary-only (Level 2), all 33 still have ≥30 samples (smallest: Cholangiocarcinoma with 36).

#### Decision: Compute All Three Levels

**Level 1 — All samples pooled (as-is):**
- N = 11,069 samples (all sample types, all cancer types)
- 785 donors contribute >1 sample → mild non-independence
- Exists: `tcga_donor_only_*` ✓

**Level 2 — Primary tumor + blood cancer only:**
- Filter to sample types 01 (Primary Tumor) + 03 (Primary Blood Cancer)
- N = 9,879 samples from 9,875 donors (~1:1, nearly perfectly independent)
- Removes matched normals, metastatic, recurrent
- Exists: `tcga_primary_only_*` ✓

**Level 3 — Per-cancer type (from Level 2):**
- Within each of 33 cancer types (≥30 primary-only samples), correlate across donors
- N per type: 36–1,092 (33 of 34 pass min_samples=30; "Unknown" excluded)
- Independent within each cancer type
- Exists: `tcga_primary_by_cancer_*` ✓

#### Existing Processed Data Verification

**Location:** `/data/parks34/projects/2cytoatlas/results/cross_sample_validation/tcga/`

| Level | Expression | CytoSig (43) | LinCytoSig (178) | SecAct (1,170) | Match? |
|-------|-----------|-------------|-----------------|---------------|--------|
| donor_only (all) | 11,069 × 20,501 | 11,069 × 43 | 11,069 × 178 | 11,069 × 1,170 | Yes |
| by_cancer (no Unknown) | 10,409 × 20,501 | 10,409 × 43 | 10,409 × 178 | 10,409 × 1,170 | Yes |
| primary_only (01+03) | 9,879 × 20,501 | 9,879 × 43 | 9,879 × 178 | 9,879 × 1,170 | Yes |
| primary_by_cancer (33 types) | 9,236 × 20,501 | 9,236 × 43 | 9,236 × 178 | 9,236 × 1,170 | Yes |

**Regenerated 2026-02-13:** All 4 levels regenerated with negative clipping (`clip(0)` before `log2(RSEM+1)`). Previous data had 1,247,092 negative values (0.55%) from EBPlusPlus batch correction artifacts propagated through log2 transform. Metadata now includes `data_format='EBPlusPlus_RSEM_normalized_counts'` and `transform='clip(0) -> log2(RSEM+1)'`.

**Correlation summary (post-clipping, 2026-02-13):**

| Level | CytoSig median ρ | LinCytoSig median ρ | SecAct median ρ | N targets (CytoSig/LinCytoSig/SecAct) |
|-------|------------------|--------------------|-----------------|------------------------------------|
| donor_only | 0.238 | 0.206 | 0.415 | 41 / 152 / 1,085 |
| by_cancer | 0.178 | 0.152 | 0.289 | 41 / 152 / 1,085 |
| primary_only | 0.225 | 0.203 | 0.406 | 41 / 152 / 1,085 |
| primary_by_cancer | 0.147 | 0.151 | 0.279 | 41 / 152 / 1,085 |

**Note on existing `by_cancer`:** Removes "Unknown" cancer type (660 samples) but does NOT filter by sample type — it still includes normals, metastatic, and recurrent within each cancer type.

**Note on `primary_by_cancer`:** Filters to primary tumors (01) + blood cancer (03) first, then stratifies by cancer type with min_samples=30. This gives 33 groups (9,236 samples) — the one excluded group is "Unknown" cancer type which has fewer than 30 primary-only samples. Sample counts per cancer type differ from the donor_only table above (e.g., Breast: 1,092 vs 1,212; Cholangiocarcinoma: 36 vs 45) because normals, metastatic, and recurrent samples are excluded.

---

## 3. Cross-Cutting Cleaning Decisions

### 3.1 Defining the Independent Sample Unit

The core question for every dataset:

| Dataset | Independent unit | N | Fully independent? |
|---------|-----------------|---|-------------------|
| CIMA | donor | 421 | Yes |
| Inflammation Main | sampleID | 817 | Unknown (need to verify donor identity) |
| scAtlas Normal | donor × organ | 706 | No (115 donors share across organs) |
| scAtlas Cancer | donor × tissue_type | ~800 | No (82 donors share across tissue types) |
| GTEx | tissue biopsy | 19,788 | No (946 donors × multiple tissues) |
| TCGA | tumor biopsy | ~10,000 | Mostly yes (after removing matched normals) |

### 3.2 Cell Exclusion (Single-Cell Datasets)

Before pseudobulk aggregation, exclude:
- `LowQuality_cells` (Inflammation) — CIMA is already QC'd upstream; no such cells exist
- `Doublets` (Inflammation, scAtlas) — CIMA is already QC'd upstream; no such cells exist
- Cells with `predicted_doublet == True` (scAtlas)

### 3.3 Minimum Cell Count per Pseudobulk Group

A pseudobulk aggregate with too few cells is noisy. Two-stage approach:
- **Generation stage** (`01_`, `02_`, `03_`): No filtering — all groups preserved in H5ADs
- **Validation stage** (`11_`, `validation/`): min_cells=10 applied at correlation time
- See Section 3.8 for full per-script threshold audit

### 3.4 Minimum Samples per Correlation

Spearman correlation requires sufficient N for reliable estimates:
- **Single-cell validation** (`validation/config.py`): min_samples=**10** — appropriate for cell-type strata where some types appear in few donors
- **scAtlas Normal by_organ** (`12_`/`13_`): min_samples=**20** — balances tissue coverage (12 tissues, +26% samples vs ≥30) against Spearman power (detects rho ≥ 0.45 at p < 0.05). Report Tier A (≥30) and Tier B (20–29) separately.
- **Bulk stratified** (`15_bulk_validation.py`, `15b_tcga_primary_filter.py`): min_samples=**30** — appropriate for per-tissue/per-cancer groups which are larger
- **Report N alongside every rho** so readers can assess reliability

### 3.8 Threshold Reconciliation Across Pipeline Scripts

**Problem:** Different pipeline stages use different thresholds. Audit of all scripts:

**min_cells (pseudobulk aggregation — minimum cells per sample×celltype group):**

| Script | min_cells | Notes |
|--------|-----------|-------|
| `01_cima_activity.py` | **None** | All groups included regardless of size |
| `02_inflam_activity.py` | **None** | All groups included |
| `03_scatlas_analysis.py` | **10** (default), **1** (actual calls at lines 1340, 1357) | Inconsistent: default vs actual usage |
| `11_donor_level_pipeline.py` | **10** | Applied at validation stage |
| `validation/config.py` | **10** (standard), **50** (resampled) | Two-tier |
| `validation/01_generate_pseudobulk.py` | **10** | |
| `validation/run_validation.py` | **10** | |

**min_samples (correlation — minimum samples per Spearman correlation):**

| Script | min_samples | Notes |
|--------|-------------|-------|
| `validation/config.py` | **10** | |
| `validation/run_validation.py` | **10** | |
| `15_bulk_validation.py` (stratified) | **30** | GTEx by-tissue, TCGA by-cancer |
| `15b_tcga_primary_filter.py` (stratified) | **30** | TCGA primary by-cancer |
| scAtlas Normal by_organ (`12_`/`13_`) | **20** | Per-tissue correlation (12 tissues, Tier A+B) |
| scAtlas Normal by_organ_subCluster | **20** (Breast/Lung), **10** (Spleen/Liver) | Level 3 exploratory |

**Assessment:** The two-stage design is actually correct: activity scripts (`01_`, `02_`, `03_`) generate pseudobulk with all groups (no filtering), and downstream validation scripts (`11_`, `validation/`) apply min_cells=10 at correlation time. This means the pseudobulk H5ADs are complete and reusable — filtering decisions are deferred to analysis.

**Remaining item to review:**
- `03_scatlas_analysis.py` passes min_cells=1 in some calls (lines 1340, 1357), overriding its default of 10. This is in the tumor vs adjacent comparison where all groups are needed regardless of size. Should be verified as intentional.

**Resolved:** min_samples differences (10 vs 30) are context-appropriate and documented in Section 3.4. The two-stage pseudobulk design (generate all, filter at analysis) is confirmed as correct architecture.

### 3.5 Gene ID Harmonization

| Dataset | Gene ID format | Mapping needed |
|---------|---------------|----------------|
| CIMA | HGNC symbols (var_names) | Direct match to signatures |
| Inflammation | Ensembl IDs (var_names), HGNC in var['symbol'] | Symbol column used for matching (22,826 unique symbols) |
| scAtlas Normal | HGNC symbols | Direct match |
| scAtlas Cancer | HGNC symbols | Direct match |
| GTEx | ENSG (versioned) | GENCODE v23 probemap → HGNC (3-tier: versioned → unversioned → Description column) |
| TCGA | `symbol\|entrezID` (e.g., `TP53\|7157`) | Split on `\|`, take symbol; no probemap needed |

CytoSig target aliases: TNFA→TNF, GMCSF→CSF2, IFNA→IFNA1, etc. Must resolve before correlation.

### 3.5b Duplicate Gene Symbol Handling

After gene ID mapping, multiple source IDs can map to the same HGNC symbol. Both bulk pipelines resolve this with `groupby(level=0).mean()`.

**Duplicate counts per dataset:**

| Dataset | Raw IDs | After `?` removal | Unique symbols | Duplicated symbols | Duplicated rows |
|---------|---------|-------------------|----------------|-------------------|----------------|
| GTEx | 74,628 ENSG | 72,930 mapped | ~72,700 | 231 | 1,538 |
| TCGA | 20,531 | 20,502 | 20,501 | 1 (SLC35E2) | 2 |
| CIMA / scAtlas | N/A | N/A | N/A | N/A | N/A (already HGNC symbols) |
| Inflammation | 22,838 ENSG | N/A | 22,826 symbols via var['symbol'] | 12 (symbol collisions) | Handled by `11_donor_level_pipeline.py` via gene_symbols parameter |

**GTEx duplicates are overwhelmingly non-coding RNAs:**
- Y_RNA (727 ENSG IDs → 1 symbol), Metazoa_SRP (167), U3 (43), U6 (33), SNORA70 (25), ...
- These are absent from CytoSig/SecAct signature matrices, so averaging has no impact on activity inference

**Resolution method:** `pandas.DataFrame.groupby(level=0).mean()`
- Averages expression values across all rows mapping to the same symbol
- **NaN handling:** `skipna=True` by default — NaN values are excluded from the mean, not propagated. If all values for a gene+sample are NaN, the result is NaN (later replaced by `np.nan_to_num(Y, nan=0.0)` during ridge regression input preparation).
- This is conservative: averaging smooths rather than inflates, and the affected genes are mostly non-coding.

**TCGA has essentially no duplicates** (1 gene: SLC35E2 with 2 entrezIDs). The averaging is trivially correct.

**Single-cell datasets** use HGNC symbols natively — no gene mapping or deduplication needed.

### 3.6 Expression Normalization

**Single-cell datasets — two pipelines with different normalization:**

*Activity scripts (`01_cima_activity.py`, `02_inflam_activity.py`):*
1. Sum raw counts per group (sample × celltype)
2. CPM normalize on summed counts: (counts / total_counts) × 1e6
3. Log2 transform: log2(CPM + 1)
4. Used for: disease differential, treatment response analyses

*Donor-level pipeline (`11_donor_level_pipeline.py`):*
1. Per-cell: CPM normalize → log1p (natural log)
2. Mean of per-cell log1p(CPM) values across cells in group
3. Used for: **cross-donor correlation validation** (the main validation analysis)
4. All single-cell atlases use this same pipeline → internally consistent

Note: Mean-of-log ≠ log-of-mean (Jensen's inequality), but for Spearman correlation (rank-based), both approaches produce valid rankings. The key requirement is consistency across atlases, which is satisfied.

**Bulk datasets (format-specific preprocessing):**

| Dataset | Source format | Preprocessing | Transform |
|---------|-------------|---------------|-----------|
| GTEx v11 | TPM (always ≥ 0) | None needed | `log2(TPM + 1)` |
| TCGA PanCancer | EBPlusPlus batch-adjusted RSEM normalized counts | Clip negatives to 0 (0.55% of values are small negatives from batch correction) | `log2(RSEM + 1)` |
| TOIL fallback | Uniformly recomputed TPM | None needed | `log2(TPM + 1)` |

**Why this is sufficient for our analysis:**
- **Spearman correlation** is rank-based — only rank order matters, not absolute scale
- **Ridge regression** mean-centers per gene across samples before inference — absolute scale cancels out
- Both formats produce comparable downstream results because activity z-scores depend on relative variation across samples, not absolute normalization

**Pipeline implementation** (`15_bulk_validation.py`): Each dataset's `data_format` is tracked through loading → expression H5AD → activity H5AD metadata (`uns['data_format']`, `uns['transform']`) for full provenance.

### 3.7 Confounders (Age, Sex, Batch, Tumor Purity)

Not a concern for expression-activity correlation validation. Both the target gene expression (X) and predicted activity (Y) are derived from the **same expression matrix**, so any factor that shifts expression also shifts the activity prediction proportionally. Biological confounders (age, sex, tumor purity) create natural variation that the model should track — this is signal, not noise. Technical confounders (batch, sequencing center) are partially handled by ridge regression's mean-centering step.

The one caveat: confounders can inflate **absolute rho magnitude** by widening the expression range (e.g., TCGA's heterogeneous tumor purity likely contributes to its higher rho vs. CIMA's homogeneous healthy donors). This affects cross-dataset comparison of rho values but does not create false positive/negative correlations for individual targets.

Confounders would become relevant for downstream analyses correlating activity with **independent** measurements (protein levels, clinical outcomes, treatment response).

---

## 4. Proposed Reprocessing Pipeline

### Phase 1: Standardized Pseudobulk Generation

For each single-cell dataset:
```
1. Load H5AD (backed mode)
2. Exclude low-quality cells and doublets (Section 3.2)
3. For each aggregation level:
   a. Define groupby keys (donor, donor×celltype, donor×organ, etc.)
   b. Sum counts for ALL groups (no min_cells filtering at this stage)
   c. CPM normalize → log1p
   d. Save as clean pseudobulk H5AD with n_cells per group in metadata
4. Run ridge regression on pseudobulk → activity z-scores
```
Note: min_cells filtering is deferred to Phase 2 (correlation), keeping H5ADs reusable.

### Phase 2: Standardized Correlation Computation

For each dataset × aggregation level:
```
1. Load pseudobulk expression + activity
2. For each target signature:
   a. Resolve target → gene name
   b. Extract expression and activity vectors
   c. Remove NaN/Inf, require ≥ min_samples
   d. Require non-zero variance in both vectors
   e. Compute Spearman rho + p-value
   f. Record N, mean, std for both vectors
3. Apply BH-FDR correction per dataset × signature group
4. Save correlation table
```

### Phase 3: Reporting

For each dataset, report:
- N (samples used), not just "number of significant correlations"
- Independence status (fully independent vs partially correlated samples)
- Aggregation level used
- Filtering applied (cells excluded, groups dropped)

---

## 5. Action Items

### Immediate (before reprocessing)

- [x] **Inflammation Atlas cohort decision:** Analyze Main/Val/Ext separately. Main as primary (curated annotations), Val for internal replication, Ext for independent validation.
- [x] **Inflammation Atlas overlap analysis:** 0 shared sampleIDs across all cohorts. 217 shared libraryIDs between Main↔Val (multiplexed sequencing, no shared cells/donors). Safe to analyze separately.
- [ ] **Inflammation Atlas donor identity:** Investigate sampleID naming convention to determine if donors can be identified. Low priority — treat sampleIDs as independent for now.
- [x] **scAtlas Normal:** Decided on primary aggregation level — by_organ (per-tissue) as primary, donor_organ as supplementary, donor_only as supplementary
- [ ] **scAtlas Cancer:** Decide on tissue-type filtering (tumor-only vs all)
- [x] **scAtlas Normal:** Use `subCluster` (102 values) instead of `cellType1` (468 values) for cell-type-stratified analyses
- [x] **scAtlas Normal:** Document `cellType1` annotation inconsistencies (25 case groups, 7+ plural pairs)
- [x] **scAtlas Normal:** Threshold decided — min_samples=20 for Level 1 by_organ (12 tissues); report Tier A (≥30) and Tier B (20–29) separately
- [x] **scAtlas Normal:** Level 3 (by_organ_subCluster) scoped — exploratory only; Breast/Lung at min_samples=20, Spleen/Liver at min_samples=10
- [x] **scAtlas Normal:** GTEx cross-platform tissue mapping — 11 of 12 tissues have direct GTEx matches (Thymus excluded)
- [x] **scAtlas Normal:** Per-tissue stratified correlation — already implemented via `donor_organ` pseudobulk. `13_cross_sample_correlation_analysis.py` detects `tissue` column and computes per-tissue Spearman automatically. Results in `scatlas_normal_correlations.csv` (7 Tier A + 5 Tier B tissues)
- [x] **scAtlas Normal:** Generate donor_only level — used `donor_col='donorID'` (not CIMA pattern since sampleID ≠ donor). Generated: pseudobulk 317 × 21,812, cytosig 317 × 43, lincytosig 317 × 178, secact 317 × 1,170
- [x] **GTEx:** Decided on stratification — per-tissue (by_tissue) as primary, pooled (donor_only) as supplementary with non-independence caveat
- [x] **TCGA:** Decided on sample filtering — primary tumor (01) + blood cancer (03); implemented in `15b_tcga_primary_filter.py`

### Inflammation Regeneration (PRIORITY)

The existing inflammation cross_sample_validation data contains Doublets and LowQuality_cells in the pseudobulk. This requires regeneration:

**Step 1: Add cell exclusion to pseudobulk generation scripts** ✅ DONE
- Modified `11_donor_level_pipeline.py`: Added `exclude_celltypes` config to inflammation atlas entries. Exclusion logic marks cells as NaN (preserving positional alignment with adata.X) so they are skipped during batch accumulation.
- Modified `12_cross_sample_correlation.py`: Same `exclude_celltypes` config. Exclusion sets `sample_col` to NaN → cells become `__INVALID__` in group keys. This script handles all levels including `donor_only`.
- Both scripts restored from `archive/scripts/early_pipeline/` to `scripts/`
- Main: `{'Level1': ['Doublets', 'LowQuality_cells']}`
- Val/Ext: `{'Level1pred': ['Doublets', 'LowQuality_cells']}` (defensive — these files are already QC'd and contain neither category)
- Does NOT exclude `NK_lowRibocontent` (real cell type with quality caveat)

**Step 2: Regenerate pseudobulk + activity for all 3 cohorts** ✅ DONE (2026-02-13)
- inflammation_main: Pseudobulk regenerated via SLURM (482,218 cells excluded = 9.8%). Activity re-inferred with `ridge_batch` on A100 GPU (22 min).
- inflammation_val: Pseudobulk unchanged (already QC'd upstream). Activity re-inferred with `ridge_batch` (16.5 min).
- inflammation_ext: Pseudobulk unchanged (already QC'd upstream). Activity re-inferred with `ridge_batch` (17.7 min).
- Initial SLURM job (11654934) failed due to missing `libcublas.so.12`; re-executed directly on GPU node after `module load CUDA/12.8.1 cuDNN`.

**Step 3: Recompute correlations** ✅ DONE (2026-02-13)
- Fixed min_samples threshold in `13_cross_sample_correlation_analysis.py`: per-celltype changed from 5 → 10 to match overall threshold
- Added TCGA primary-only levels (`primary_only`, `primary_by_cancer`) to `13_` ATLAS_CONFIGS
- Both `12_` and `13_` restored from archive to `scripts/`
- Executed: `python scripts/13_cross_sample_correlation_analysis.py --atlas inflammation_main inflammation_val inflammation_ext`
- Results: inflammation_main 72,358 rows, inflammation_val 52,090 rows, inflammation_ext 56,325 rows

**Step 4: Update resampled validation**
```bash
python scripts/16_resampled_validation.py --atlas inflammation_main
python scripts/16_resampled_validation.py --atlas inflammation_val
python scripts/16_resampled_validation.py --atlas inflammation_ext
```

**Observed changes after regeneration (inflammation_main):**
- donor_only: CytoSig median ρ 0.321 → 0.323, SecAct 0.188 → 0.173 (cleaner profiles + ridge_batch)
- L1: 15 cell types instead of 17 (Doublets and LowQuality_cells removed)
- L2: 63 cell types instead of 65 (Doublets and LowQuality_cells removed)
- Correlation CSVs: ~3,844 fewer rows (1,922 Doublets + 1,922 LowQuality_cells)
- Median rho at donor_only may shift slightly (cleaner expression profiles)

### Other Reprocessing

- [x] Standardize cell exclusion criteria across all single-cell datasets (Section 3.2) — verified: CIMA pre-QC'd (no bad cells), Inflammation has `exclude_celltypes` for Doublets/LowQuality_cells, scAtlas Normal/Cancer have `predicted_doublet` all False (doublets removed upstream). No additional exclusion needed.
- [x] Keep two-stage threshold design: no min_cells at generation, min_cells=10 at validation (Section 3.3, 3.8) — confirmed across `12_`, `validation/config.py`
- [x] Keep context-appropriate min_samples: 10 for single-cell celltype strata, 20 for scAtlas by_organ, 30 for bulk stratified (Section 3.4, 3.8) — confirmed across `12_`, `13_`, `15_`
- [x] Regenerate TCGA data with negative clipping (Section 3.6) — `15_bulk_validation.py --dataset tcga --force --backend cupy` + `15b_tcga_primary_filter.py --force --backend cupy`. All 4 levels regenerated with `clip(0)` before `log2(RSEM+1)`.
- [x] Regenerate GTEx H5ADs to add `data_format` metadata — `15_bulk_validation.py --dataset gtex --force --backend cupy`. Both donor_only and by_tissue regenerated with `data_format` field.
- [x] Generate TCGA Level 2 (primary-only) and Level 3 (per-cancer from primary) — `15b_tcga_primary_filter.py --force` (9,879 primary-only samples; 9,236 in 33 cancer-type strata)
- [x] Recompute correlations for updated datasets — `13_cross_sample_correlation_analysis.py --atlas tcga gtex`
- [ ] Rebuild figures and tables from new correlations

### Documentation

- [x] Record inflammation cohort overlap analysis and decision rationale
- [x] Document cell exclusion issue in inflammation data
- [x] Report sample independence status for each dataset — documented in Section 3.1 table (CIMA: fully independent, Inflammation: unknown donor identity, scAtlas: partial non-independence from multi-tissue donors, GTEx: severe non-independence, TCGA: mostly independent after primary filtering)
- [ ] Create a provenance chain: raw H5AD → filtered cells → pseudobulk → activity → correlation → figure (partially done: `baseline/PROVENANCE.md` exists)

---

## 6. Impact on Existing Report

If reprocessing changes the numbers, these sections of the baseline report need updating:

| Report section | What changes |
|----------------|-------------|
| 1.2 Processing Scale | Cell counts, donor/sample counts |
| 2.1 Datasets and Scale | scAtlas: 781 donors → 1,034 donors (317 normal + 717 cancer) |
| 4.1 Summary statistics | All median rho values |
| 4.2 Per-atlas boxplots | All correlation distributions |
| 4.3 Cross-atlas consistency | Consistency profiles may shift |
| 4.5 Aggregation levels | Results depend on cleaning choices |
| 5.x Method comparison | All 10-way comparison numbers |
| Figures 2-15 | Most figures need regeneration |

**Known discrepancy:** The report states scAtlas has "781 donors" — actual data shows 317 (normal) + 717 (cancer) = 1,034 donors. This needs investigation: either the data has changed since the report was generated, or the 781 figure used a different counting method.
