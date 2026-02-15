# CytoAtlas: Pan-Disease Single-Cell Cytokine Activity Atlas

**February 14, 2026**

---

## Executive Summary

CytoAtlas is a comprehensive computational resource that maps cytokine and secreted protein signaling activity across **~29 million human cells and ~31,000 bulk RNA-seq samples** from six independent datasets spanning two bulk RNA-seq resources (GTEx, TCGA) and four single-cell compendia: CIMA (6.5M healthy donor cells), Inflammation Atlas (6.3M disease cells across 3 cohorts; Main cohort used for validation), scAtlas (6.4M organ and cancer cells), and parse_10M (9.7M cytokine-perturbed cells). The system uses **linear ridge regression** against experimentally derived signature matrices to infer activity — producing fully interpretable, conditional z-scores rather than black-box predictions. Each prediction traces back to a weighted combination of known gene-to-cytokine relationships, with permutation-based z-scores providing significance estimates.

**Key results:**
- 1,213 signatures (43 CytoSig cytokines + 1,170 SecAct secreted proteins), plus 178 cell-type-specific LinCytoSig variants, validated across 6 independent datasets
- Spearman correlations between predicted activity and target gene expression reach ρ=0.6–0.8 at independent levels for well-characterized cytokines (IL1B, TNFA, VEGFA, TGFB family), exceeding ρ=0.9 in specific tissue/cell-type strata
- Cross-dataset consistency demonstrates that signatures generalize across CIMA, Inflammation Atlas Main, scAtlas, GTEx, and TCGA
- Cell-type-specific signatures (LinCytoSig) improve prediction for select immune cell types (Basophil, NK, DC: +0.18-0.21 Δρ) but generally underperform global CytoSig for non-immune cell types
- SecAct provides the broadest validated coverage with 805–1,161 targets per dataset (varying by gene overlap), achieving the highest median correlations in 5 of 6 datasets (independence-corrected median ρ=0.19–0.46)

---

## 1. System Architecture and Design Rationale

### 1.1 Architecture and Processing

**Linear interpretability over complex models.**
Ridge regression (L2-regularized linear regression) was chosen deliberately over methods like autoencoders, graph neural networks, or foundation models. The resulting activity z-scores are **conditional on the specific genes in the signature matrix**, meaning every prediction can be traced to a weighted combination of known gene responses. This is critical for biological interpretation — a scientist can ask "which genes drive the IFNG activity score in this sample?" and get a direct answer.

**Reproducibility through separation of concerns.** The system is divided into independent components, each chosen for the constraints of HPC/SLURM infrastructure:

| Component | Technology | Purpose | Rationale |
|-----------|-----------|---------|-----------|
| **Pipeline** | Python + CuPy (GPU) | Activity inference | 10–34x speedup over NumPy; batch-streams H5AD files (500K–1M cells/batch) with projection matrix held on GPU; automatic CPU fallback when GPU unavailable |
| **Storage** | DuckDB (3 databases, 68 tables) | Columnar analytics | Single-file databases require no server — essential on HPC where database servers are unavailable; each database regenerates independently without affecting others |
| **API** | FastAPI (262 endpoints) | RESTful data access | Async I/O for concurrent DuckDB queries; automatic OpenAPI documentation; Pydantic request validation |
| **Frontend** | React 19 + TypeScript | Interactive exploration (12 pages) | Migrated from 25K-line vanilla JS SPA to 11.4K lines (54% reduction) with type safety, component reuse, and lazy-loaded routing |

**Processing scale.** Ridge regression (λ=5×10⁵) is applied using `secactpy.ridge()` against each signature matrix. For single-cell data, expression is first aggregated to pseudobulk (donor or donor×celltype level), then genes are intersected with the signature matrix (CytoSig: ~4,860 genes; SecAct: ~7,450 genes). The resulting z-scored activity coefficients are compared to target gene expression via Spearman correlation across donors.

| Dataset | Cells/Samples | Processing Time | Hardware |
|---------|---------------|-----------------|----------|
| GTEx | 19,788 bulk samples | ~10min | A100 80GB |
| TCGA | 11,069 bulk samples | ~10min | A100 80GB |
| CIMA | 6.5M cells | ~2h | A100 80GB |
| Inflammation Atlas (main/val/ext) | 6.3M cells | ~2h | A100 80GB |
| scAtlas Normal | 2.3M cells | ~1h | A100 80GB |
| scAtlas Cancer | 4.1M cells | ~1h | A100 80GB |
| parse_10M | 9.7M cells | ~3h | A100 80GB |

**Total:** ~29M single cells + ~31K bulk RNA-seq samples, processed through ridge regression against 3 signature matrices (CytoSig: 43 cytokines, LinCytoSig: 178 cell-type-specific, SecAct: 1,170 secreted proteins). **Processing Time** = wall-clock time for full activity inference on a single NVIDIA A100 GPU (80 GB VRAM). For bulk datasets (GTEx/TCGA), ridge regression is applied with within-tissue/within-cancer mean centering to remove tissue-level variation. See Section 2.1 for per-dataset details.

> **System Architecture** (`fig_system_architecture.png`): Top-to-bottom layered system design diagram. Users interact via a React SPA through Nginx and FastAPI (17 routers, JWT auth). The API serves a Data Query Service (DuckDB, 3 databases, 80+ tables) and an AI Chat Service with dual LLM (Mistral-Small-24B via vLLM + Claude fallback), RAG (LanceDB + MiniLM), and 22 data tools. An offline GPU pipeline (SLURM/A100) performs batch activity inference into the storage layer.

### 1.2 Validation Strategy

CytoAtlas validates at **four aggregation levels**, each testing whether predicted activity correlates with target gene expression (Spearman ρ) across independent samples:

| Level | Description | Datasets | Report Section |
|-------|-------------|----------|----------------|
| **Donor pseudobulk** | One value per donor, averaging across cell types | CIMA, Inflammation Atlas Main, scAtlas Normal/Cancer | §4.1, §4.3 |
| **Donor × cell-type** | Stratified by cell type within each donor | CIMA, Inflammation Atlas Main, scAtlas Normal/Cancer | §4.7 |
| **Per-tissue / per-cancer** | Median-of-medians across tissues or cancer types | GTEx (29 tissues), TCGA (33 cancer types) | §4.2, §4.3 |
| **Cross-platform** | Bulk vs pseudobulk concordance per tissue/cancer | GTEx vs scAtlas Normal, TCGA vs scAtlas Cancer | §4.4 |

All statistics use **independence-corrected** values — preventing inflation from repeated measures across tissues, cancer types, or cell types. CytoSig vs SecAct comparisons use Mann-Whitney U (total) and Wilcoxon signed-rank (32 matched targets) with BH-FDR correction. See Section 3.3 for the validation philosophy and Section 4 for full results.

> **Why independence correction matters:** Pooling across tissues or cancer types yields higher correlations than within-stratum values. For example, GTEx pooled CytoSig median ρ (0.211) is 40% higher than the independence-corrected by-tissue value (0.151); SecAct shows +26% inflation (0.394 vs 0.314). All results in this report use the corrected values. For a detailed comparison of pooled vs independent levels, including inflation magnitude and finer cell-type stratification, see the [Section 4.1 statistical supplement](stats_section_4.1.html).

> **Figure 1** (`fig1_dataset_overview.png`): Data sources, activity inference pipeline, and validation analyses.

---

## 2. Dataset Catalog

### 2.1 Datasets and Scale ([detailed analytics](../reports/weekly/dataset_analytics.html))

| # | Dataset | Type | Cells/Samples | Donors | Cell Types | Source |
|---|---------|------|---------------|--------|------------|--------|
| 1 | **GTEx** | Bulk RNA-seq | 19,788 samples | 946 donors | — | GTEx v11 (30 tissues) |
| 2 | **TCGA** | Bulk RNA-seq | 11,069 samples | 10,274 donors | — | PanCancer (33 cancer types) |
| 3 | **CIMA** | scRNA-seq | 6,484,974 | 421 donors | 27 L2 / 100+ L3 | J. Yin et al., *Science*, 2026 |
| 4 | **Inflammation Atlas Main** | scRNA-seq | 4,918,140 | 817 samples\* | 66+ | Jimenez-Gracia et al., *Nature Medicine*, 2026 |
| 5 | **Inflammation Atlas Val** | scRNA-seq | 849,922 | 144 samples\* | 66+ | Validation cohort |
| 6 | **Inflammation Atlas Ext** | scRNA-seq | 572,872 | 86 samples\* | 66+ | External cohort |
| 7 | **scAtlas Normal** | scRNA-seq | 2,293,951 | 317 donors | 102 subCluster | 35 organs |
| 8 | **scAtlas Cancer** | scRNA-seq | 4,146,975 | 717 donors (601 tumor-only) | 162 cellType1 | 29 cancer types |
| 9 | **parse_10M** | scRNA-seq | 9,697,974 | 12 donors × 90 cytokines (+PBS control) | 18 PBMC types | Cytokine perturbation |

**Grand total: ~29 million single cells + ~31K bulk samples from 6 independent studies (9 datasets), 100+ cell types**

*\* Inflammation Atlas does not provide donor-level identifiers; the 817/144/86 values are sample counts. The donor–sample relationship is unknown, so correlations use sampleID as the independent unit.*

### 2.2 Disease and Condition Categories

**CIMA (421 healthy donors):** Healthy population atlas with paired blood biochemistry (19 markers: ALT, AST, glucose, lipid panel, etc.) and plasma metabolomics (1,549 features). Enables age, BMI, sex, and smoking correlations with cytokine activity.

**Inflammation Atlas (19 diseases):**
- Autoimmune: RA, SLE, Sjogren's, PSA
- IBD: Crohn's disease, Ulcerative Colitis
- Infectious: COVID-19, Sepsis, HIV, HBV
- Cancer: BRCA, CRC, HNSCC, NPC
- Other: COPD, Cirrhosis, MS, Asthma, Atopic Dermatitis

**scAtlas Normal (317 donors):** 35 organs, 12 tissues with ≥20 donors for per-organ stratification (Breast 124, Lung 97, Colon 65, Heart 52, Liver 43, etc.)

**scAtlas Cancer (717 donors, 601 tumor-only):** 29 cancer types, 11 with ≥20 tumor-only donors for per-cancer stratification (HCC 88, PAAD 58, CRC 51, ESCA 48, HNSC 39, LUAD 36, NPC 36, KIRC 31, BRCA 30, ICC 29, STAD 27)

**parse_10M perturbations:** 90 cytokines × 12 donors (+PBS control; perturbation resource for validation, not ground truth — true ground truth is unattainable for inferential activity scores)


### 2.3 Signature Matrices

| Matrix | Targets | Construction | Source |
|--------|---------|-------------|--------|
| **CytoSig** | 43 cytokines | Median log2FC across all experimental bulk RNA-seq | Jiang et al. |
| **LinCytoSig** | 178 (45 cell types × 1-13 cytokines) | Cell-type-stratified median from CytoSig database | This work |
| **SecAct** | 1,170 secreted proteins | Median global Moran's I across 1,000 Visium datasets | Ru et al., *Nature Methods*, 2026 (in press) |

---

## 3. Scientific Value Proposition

### 3.1 What Makes CytoAtlas Different from Deep Learning Approaches?

Most single-cell analysis tools use complex models (variational autoencoders, graph neural networks, transformer-based foundation models) that produce **aggregated, non-linear representations** difficult to interpret biologically. CytoAtlas takes the opposite approach:

| Property | CytoAtlas (Ridge Regression) | Typical DL Approach |
|----------|------|------|
| **Model** | Linear (z = Xβ + ε) | Non-linear (multi-layer NN) |
| **Interpretability** | Every gene's contribution is a coefficient | Feature importance approximated post-hoc |
| **Conditionality** | Activity conditional on specific gene set | Latent space mixes all features |
| **Confidence** | Permutation-based z-scores with CI | Often point estimates only |
| **Generalization** | Tested across 6 independent datasets (bulk + single-cell) | Often tested on held-out splits of same cohort |
| **Bias** | Transparent — limited by signature matrix genes | Hidden in architecture and training data |

**The key insight:** CytoAtlas is not trying to replace DNN-based tools. It provides a linear, interpretable readout: when CytoAtlas reports "IFNG activity is elevated in CD8+ T cells from RA patients," the contributing signature genes and their weights can be examined directly in those cells.

### 3.2 What Scientific Questions Does CytoAtlas Answer?

1. **Which cytokines are active in which cell types across diseases?** — IL1B/TNFA in monocytes/macrophages, IFNG in CD8+ T and NK cells, IL17A in Th17, VEGFA in endothelial/tumor cells, TGFB family in stromal cells — quantified across 19 diseases, 35 organs, and 15 cancer types.
2. **Are cytokine activities consistent across independent cohorts?** — Yes. IL1B, TNFA, VEGFA, and TGFB family show consistent positive correlations across all 6 validation datasets (Figure 6).
3. **Does cell-type-specific biology matter for cytokine inference?** — For select targets, yes: LinCytoSig improves IL6×Macrophage, VEGFA×Endothelial, and IL2×CD8T prediction in normal tissue, but global CytoSig wins overall and the advantage does not transfer to cancer (Figures 11–12).
4. **Which secreted proteins beyond cytokines show validated activity?** — SecAct (1,170 targets) achieves the highest correlations in 5 of 6 datasets (median ρ=0.19–0.46), with validated targets including INHBA/Activin A (ρ=0.91 in TCGA Colon), CXCL12 (ρ=0.94 in scAtlas Normal Fibroblast), and BMP family.
5. **Can we predict treatment response from cytokine activity?** — We are incorporating cytokine-blocking therapy outcomes from bulk RNA-seq to test whether predicted cytokine activity associates with therapy response. Additionally, Inflammation Atlas responder/non-responder labels (208 samples across 6 diseases: RA, PS, PSA, CD, UC, SLE) enable treatment response prediction using cytokine activity profiles as features.

### 3.3 Validation Philosophy

CytoAtlas validates against a direct principle: **if CytoSig predicts high IFNG activity for a sample, that sample should have high IFNG gene expression.** This expression-activity correlation is computed via Spearman rank correlation across donors/samples.

This validation only captures signatures where the target gene itself is expressed. Signatures that act primarily through downstream effectors without upregulating the ligand gene itself would not be captured by this metric.

---

## 4. Validation Results

Validation uses a simple, conservative principle: if a method predicts high activity for a target, that target's gene expression should be high. All Spearman ρ values below use **independence-corrected** statistics at each dataset's primary independent level — preventing inflation from repeated measures across tissues, cancer types, or cell types. For GTEx and TCGA, each target's representative ρ is the **median across per-tissue/per-cancer values** (median-of-medians); for single-cell atlases, each target has one ρ at the independent donor/tumor level. Full detail tables and supplementary levels are in the linked HTML supplements.

### 4.1 Overall Performance Summary ([full details](stats_section_4.1.html))

> **Figure 2** (`fig2_summary_table.png`): Complete validation statistics table.

| Dataset (Independent Level) | Method | N Targets | Median ρ | % Sig (BH q<0.05) | % Pos |
|-----------------------------|--------|-----------|----------|-------------------|-------|
| GTEx (by_tissue, 29 tissues) | CytoSig | 43 | 0.151 | 81.4% | 83.7% |
| | SecAct | 1,132 | 0.314 | 88.8% | 90.5% |
| TCGA (primary_by_cancer, 33 types) | CytoSig | 41 | 0.214 | 68.3% | 92.7% |
| | SecAct | 1,085 | 0.357 | 81.5% | 95.8% |
| CIMA (donor_only, 421 donors) | CytoSig | 43 | 0.114 | 69.8% | 58.1% |
| | SecAct | 1,161 | 0.191 | 75.0% | 76.1% |
| Inflam Main (donor_only, 817 samples) | CytoSig | 33 | 0.323 | 78.8% | 69.7% |
| | SecAct | 805 | 0.173 | 84.8% | 72.7% |
| scAtlas Normal (donor_only, 317 donors) | CytoSig | 43 | 0.173 | 74.4% | 69.8% |
| | SecAct | 1,140 | 0.455 | 83.9% | 87.8% |
| scAtlas Cancer (tumor_only, 601 donors) | CytoSig | 43 | 0.184 | 81.4% | 83.7% |
| | SecAct | 1,140 | 0.399 | 92.3% | 95.0% |

**Column values:** Median Spearman ρ between predicted activity and target gene expression. For GTEx/TCGA, values are median-of-medians (each target's ρ = median across tissues/cancers, then median across targets). % Sig = Benjamini-Hochberg FDR-corrected significance rate (q < 0.05). **Independent Level** = aggregation where samples are fully independent (see [analytics](../reports/weekly/dataset_analytics.html)).

**Summary ranking:** SecAct > CytoSig in 5 of 6 datasets (GTEx, TCGA, CIMA, scAtlas Normal, scAtlas Cancer). CytoSig leads Inflammation Main (0.323 vs 0.173 SecAct), the only disease-focused single-cell dataset. SecAct achieves the highest median ρ in scAtlas Normal (0.455) and TCGA (0.357). CIMA shows the lowest correlations overall (0.114 CytoSig, 0.191 SecAct).

- Inflammation Val (144 samples) and Ext (86 samples) replicate directionally; full results in [supplement](stats_section_4.1.html)

### 4.2 Cross-Dataset Comparison: CytoSig vs SecAct ([statistical methods](stats_section_4.1.html#cross-dataset-comparison))

> **Figure 3** (`fig3_cross_dataset_boxplot.png`): CytoSig vs SecAct Spearman ρ distributions across atlases, with Total/Matched tabs.

**Independence-corrected approach:** Each target gets one representative ρ at the independent level. For datasets with per-stratum correlations (GTEx 29 tissues, TCGA 33 cancers), the representative ρ is the **median across strata** (median-of-medians). For donor-only datasets (CIMA, Inflammation Main, scAtlas Normal/Cancer), each target already has one ρ. This yields n = 33–43 CytoSig values and n = 805–1,161 SecAct values per dataset (Inflammation Main has fewer targets due to stricter expression filtering in a smaller cohort).

**Total comparison (Mann-Whitney U):** Compares all CytoSig targets vs all SecAct targets.

**Matched comparison (Wilcoxon signed-rank):** Restricts to 32 shared targets (30 for TCGA due to missing IFNL1/IL36A; 27 for Inflammation Main due to smaller target set). After full alias resolution (TNFA→TNF, GMCSF→CSF2, CD40L→CD40LG, TRAIL→TNFSF10, TWEAK→TNFSF12, Activin A→INHBA, GCSF→CSF3, MCSF→CSF1, IFNL→IFNL1, IL36→IL36A), the CytoSig–SecAct overlap is 32 targets — 10 more than the previous 22-target list which only matched raw cytokine names.

| Dataset | Mann-Whitney U p | Wilcoxon p | N pairs | Direction |
|---------|-----------------|-----------|---------|-----------|
| GTEx | 4.76 × 10⁻⁴ | 3.54 × 10⁻⁵ | 32 | SecAct > CytoSig |
| TCGA | 2.85 × 10⁻³ | 3.24 × 10⁻⁶ | 30 | SecAct > CytoSig |
| CIMA | 3.18 × 10⁻² | 2.28 × 10⁻² | 32 | SecAct > CytoSig |
| Inflammation Main | 0.548 (ns) | 0.141 (ns) | 27 | ns (mixed) |
| scAtlas Normal | 1.04 × 10⁻⁴ | 3.54 × 10⁻⁵ | 32 | SecAct > CytoSig |
| scAtlas Cancer | 1.06 × 10⁻⁵ | 3.54 × 10⁻⁵ | 32 | SecAct > CytoSig |

SecAct significantly outperforms CytoSig in 5 of 6 datasets. The Inflammation Main exception is non-significant on both tests (p = 0.548 total, p = 0.141 matched); on total targets CytoSig has higher median ρ (0.32 vs 0.17), but on the 27 matched targets SecAct has higher median ρ (0.44 vs 0.35). The total comparison includes ~1,170 SecAct targets vs 33 CytoSig targets, so differences in the non-overlapping targets drive the result. The matched comparison on equal footing shows no significant difference in this dataset.

### 4.3 Per-Tissue and Per-Cancer Stratified Validation ([full results](stats_section_4.1.html#per-tissue-stratified))

> **Figure 4a** (`fig4a_gtex_per_tissue.png`): GTEx per-tissue stratified validation (29 tissues).
> **Figure 4b** (`fig4b_tcga_per_cancer.png`): TCGA per-cancer stratified validation (33 cancer types).

Four datasets have per-stratum breakdowns:
- **GTEx**: 29 tissues (47–3,234 donors each)
- **TCGA**: 33 cancer types (36–1,092 each)
- **scAtlas Normal**: 19 organs (10–124 donors each; 7 Tier A ≥30, 12 Tier B <30)
- **scAtlas Cancer**: 17 cancer types (10–88 tumors each; 9 Tier A ≥30, 8 Tier B <30)

**Statistical tests per stratum:** Mann-Whitney U (CytoSig vs SecAct total targets within each tissue/cancer) and Wilcoxon signed-rank (32 shared targets). BH-FDR correction applied within each dataset (29 tests for GTEx, 33 for TCGA, etc.). Full per-stratum results with FDR q-values are in the [supplement](stats_section_4.1.html#per-tissue-stratified).

**Key findings:**
- **Matched targets (strongest test):** On matched targets (23 pairs per GTEx tissue, 19–21 per TCGA cancer after expression filtering), SecAct wins direction in 29/29 GTEx tissues and 32/33 TCGA cancer types, with 23/29 and 31/33 reaching significance (q < 0.05). The sole exception is Acute Myeloid Leukemia (Δ = −0.08, q = 0.83). This consistency across 61 of 62 strata indicates the SecAct advantage holds at the individual tissue/cancer level, not only in aggregate
- **Total targets:** SecAct wins direction in 28/29 GTEx tissues (21 significant) and 30/33 TCGA cancers (15 significant). The few CytoSig-favored strata on total targets (Brain in GTEx; Kidney Chromophobe, Ovarian, Uveal Melanoma in TCGA) are all non-significant
- **Per-target signature performance:** Because this comparison uses the same 23–32 matched cytokines, the difference reflects per-target signature performance rather than target count. The largest Δ values occur in Small Intestine (+0.55), Vagina (+0.34), Salivary Gland (+0.30), Pancreatic (+0.34), Head & Neck (+0.27), and Lung Adeno (+0.25); the smallest in Heart (+0.06), Skin (+0.06), Brain (+0.06), Uveal Melanoma (+0.07), and Kidney Chromophobe (+0.10)
- scAtlas strata with Tier B sample sizes (<30 donors) are shown with reduced confidence

### 4.4 Cross-Platform Comparison: Bulk vs Pseudobulk

> **Figure 4c** (`fig4c_cross_platform_concordance.png`): Cross-platform concordance scatter (GTEx vs scAtlas Normal, TCGA vs scAtlas Cancer).

This section tests whether expression-activity relationships replicate across measurement technologies. For each matching tissue or cancer type, per-target Spearman ρ from bulk RNA-seq is compared to the same target's ρ from single-cell pseudobulk. Per-stratum Wilcoxon signed-rank tests (paired by target) with BH-FDR correction are used to test whether ρ values differ systematically between platforms.

**GTEx vs scAtlas Normal** (13 matching tissues: Blood, Breast, Colon, Esophagus, Heart, Kidney, Liver, Lung, Ovary, Skin, Small Intestine, Spleen, Uterus): Per-tissue CytoSig concordance ranges from ρ = −0.10 (Uterus) to 0.49 (Breast), with most tissues showing moderate positive concordance.

**TCGA vs scAtlas Cancer** (11 matching cancer types: BRCA, CRC, ESCA, HCC, HNSC, KIRC, LUAD, OV, PAAD, PRAD, STAD): Per-cancer concordance ranges from ρ = 0.23 (OV) to 0.83 (PAAD), generally higher than the GTEx tissue concordance values.

**Key finding:** Using all targets, SecAct shows significant bulk–pseudobulk differences in most strata (11/13 GTEx tissues, 5/11 TCGA cancers), while CytoSig shows almost none (1/13, 0/11). When restricted to the same 32 shared targets (Matched tabs), both methods show no significant platform differences (CytoSig: 0/13, 0/11; SecAct: 0/13, 1/11). Matched and unmatched SecAct targets show the same per-target platform shift (mean |Δ| = 0.298 vs 0.302, Mann–Whitney p = 0.82), but SecAct's ~1,000 paired targets per tissue provide 25× more observations than CytoSig's ~40, allowing detection of a small systematic shift (Δ ≈ 0.03) that is not detectable with CytoSig's sample size. Full per-stratum results in [statistical methods](stats_section_4.1.html#cross-platform-comparison).

### 4.5 Best and Worst Correlated Targets ([per-atlas details](stats_section_4.1.html#per-tissue-stratified))

> **Figure 5** (`fig5_good_bad_correlations_cytosig.png`): Top 15 and bottom 15 targets per atlas.

All 32 matched targets have data in ≥4 datasets. Using donor-level aggregation (matching Figure 5), two categories emerge (mean ρ across datasets):

**Consistent in Both CytoSig and SecAct (16 targets, both mean ρ > 0.25):**

| Target | CytoSig Mean ρ | Range | SecAct Mean ρ |
|--------|----------------|-------|---------------|
| IL1B | +0.56 | +0.24 to +0.72 | +0.58 |
| TNFA | +0.50 | +0.26 to +0.63 | +0.46 |
| IFNG | +0.44 | +0.30 to +0.62 | +0.35 |
| IL1A | +0.43 | +0.03 to +0.71 | +0.33 |
| IL27 | +0.43 | +0.19 to +0.56 | +0.40 |
| TGFB3 | +0.39 | +0.19 to +0.53 | +0.42 |
| IL6 | +0.38 | +0.17 to +0.53 | +0.48 |
| OSM | +0.38 | +0.06 to +0.49 | +0.57 |
| LIF | +0.37 | +0.20 to +0.62 | +0.50 |
| IL10 | +0.35 | +0.06 to +0.55 | +0.56 |
| CXCL12 | +0.34 | +0.10 to +0.57 | +0.59 |
| TGFB1 | +0.34 | +0.05 to +0.56 | +0.41 |
| BMP4 | +0.33 | −0.02 to +0.61 | +0.43 |
| BMP2 | +0.31 | +0.19 to +0.41 | +0.45 |
| Activin A | +0.29 | +0.12 to +0.54 | +0.56 |
| GMCSF | +0.26 | +0.01 to +0.46 | +0.37 |

**SecAct-Only: CytoSig Near-Zero, SecAct Positive (11 targets):**

| Target | CytoSig Mean ρ | Range | SecAct Mean ρ | Δ |
|--------|----------------|-------|---------------|------|
| LTA | −0.02 | −0.33 to +0.26 | **+0.53** | +0.55 |
| HGF | +0.06 | −0.29 to +0.40 | **+0.58** | +0.51 |
| TWEAK | −0.02 | −0.22 to +0.11 | **+0.44** | +0.47 |
| IL15 | +0.11 | −0.05 to +0.43 | **+0.57** | +0.47 |
| BMP6 | +0.04 | −0.40 to +0.26 | **+0.49** | +0.45 |
| TRAIL | −0.00 | −0.54 to +0.58 | **+0.44** | +0.44 |
| CD40L | +0.02 | −0.55 to +0.57 | **+0.46** | +0.43 |
| FGF2 | +0.04 | −0.23 to +0.29 | **+0.46** | +0.42 |
| MCSF | +0.16 | −0.26 to +0.50 | **+0.40** | +0.24 |
| IL21 | −0.02 | −0.22 to +0.09 | **+0.22** | +0.24 |
| BDNF | +0.11 | −0.07 to +0.20 | **+0.33** | +0.21 |

**Remaining 5 targets** do not fit either category: VEGFA and IFNL are CytoSig-only (CytoSig mean +0.38/+0.21, SecAct < 0.2), GDF11 and GCSF are borderline (CytoSig mean +0.23 each, SecAct +0.43/+0.35), and IL36 shows near-zero correlations in both methods.

**Key insight:** 16 of 32 matched targets show mean ρ > 0.25 for both CytoSig and SecAct. For 11 additional targets, CytoSig averages near zero while SecAct achieves mean ρ > 0.2 (e.g., TWEAK: CytoSig −0.02 vs SecAct +0.44; LTA: −0.02 vs +0.53; HGF: +0.06 vs +0.58). Note: this section uses donor-level aggregation (matching Figure 5's interactive chart), not the independence-corrected median-of-medians used in §4.1–4.3. Full per-atlas top-15/bottom-15 lists in [supplement](stats_section_4.1.html#per-tissue-stratified).

### 4.6 Cross-Atlas Consistency

> **Figure 6** (`fig6_cross_atlas_consistency.png`): Key targets tracked across 6 atlases at independent levels.

14 key CytoSig targets tracked across 6 datasets at independence-corrected levels reveal three consistency tiers:

| Tier | Targets | Pattern |
|------|---------|---------|
| **Universal** (ρ > 0.1 in all 6 datasets) | IL1B, TNFA, IFNG, IL6, BMP2, VEGFA | Robust across all cohorts and platforms (mean ρ 0.31–0.58) |
| **Mostly consistent** (4–5 of 6 positive) | IL10, TGFB1, CXCL12, GMCSF, HGF | Occasional near-zero outliers; TGFB1 slightly negative only in scAtlas Normal (−0.05) |
| **Context-dependent** (≤3 of 6 positive) | IL4, IL17A, EGF | Sign changes across cohorts; IL4 and IL17A absent from Inflammation Main |

**Key insight:** Mean |ρ| across 14 targets does not separate by platform: GTEx 0.26, TCGA 0.33, CIMA 0.26, Inflammation Main 0.44, scAtlas Normal 0.36, scAtlas Cancer 0.37. The universal tier (IL1B, TNFA, IFNG, IL6) shows ρ > 0.1 in all 6 datasets; context-dependent targets (IL4, IL17A) show sign changes and are absent from Inflammation Main.

### 4.7 Effect of Aggregation Level ([statistical methods](stats_section_4.1.html#aggregation-level))

> **Figure 7** (interactive): CytoSig vs SecAct boxplots across aggregation levels for each single-cell dataset (Total/Matched tabs, atlas dropdown). Mann-Whitney U (Total) and Wilcoxon signed-rank (Matched) tests with BH-FDR correction.

**Median Spearman ρ at coarsest and finest aggregation levels:**

| Dataset | Coarsest Level | CytoSig | SecAct | Finest Level | CytoSig | SecAct |
|---------|---------------|---------|--------|-------------|---------|--------|
| CIMA | Donor Only | 0.114 | 0.191 | Donor × L4 (73 types) | 0.005 | 0.052 |
| Inflammation Main | Donor Only | 0.323 | 0.173 | Donor × L2 (65 types) | 0.044 | 0.048 |
| scAtlas Normal | Donor × Organ | 0.150 | 0.295 | Donor × Organ × CT2 (356 types) | 0.059 | 0.129 |
| scAtlas Cancer | All Tumor Cells | 0.184 | 0.399 | Per Cancer Type × CT1 (~120 types) | 0.033 | 0.171 |

**Signal retention (finest / coarsest median ρ):**

| Dataset | CytoSig Retention | SecAct Retention | SecAct Advantage |
|---------|-------------------|------------------|-----------------|
| CIMA | 4% | 27% | 6.2× more robust |
| Inflammation Main | 14% | 28% | 2.0× |
| scAtlas Normal | 39% | 44% | 1.1× |
| scAtlas Cancer | 18% | 43% | 2.4× |

**Key insights:**
- **Beyond ~L2 annotation depth, % positive correlations approach 50% and signal/null ratio drops below 2×.** Inflammation Main retains signal at L2 (553 samples/group, signal/null = 3.3×) while scAtlas Normal shows weaker correlations at its shallowest stratification (22 samples/group)
- **In aggregate, SecAct's larger target pool shows higher retention** (27–44% vs CytoSig's 4–39%), but this reflects unequal target counts (~1,170 vs 43). On the 32 matched targets, CytoSig retains more signal in 2 of 4 datasets
- **All datasets show monotonic decline** except scAtlas Cancer CytoSig, which *increases* from All Tumor Cells (0.184) to Per Cancer Type (0.223) before dropping to CT1 (0.033)
- **CytoSig and SecAct converge in some datasets** — Inflammation Main L2: 0.044 vs 0.048 — but not universally (scAtlas Cancer CT1: 0.033 vs 0.171)

Per-level Mann-Whitney/Wilcoxon tests with BH-FDR correction for all datasets are in the [supplement](stats_section_4.1.html#aggregation-level).

---

## 5. CytoSig vs LinCytoSig vs SecAct Comparison

Can cell-type-specific signatures outperform global CytoSig for specific targets? We test 7 biologically motivated celltype–cytokine pairs (Macrophage×IL6, HUVEC×VEGFA, T Cell×IL2, Macrophage×IFNG, Macrophage×IL10, Macrophage×TNFA, Fibroblast×TGFB1) across all 6 datasets, comparing CytoSig, LinCytoSig, and SecAct.

For the full 10-way method comparison (8 LinCytoSig strategies + CytoSig + SecAct), see the [Open Issues & Analysis document](lincytosig_issues.html#full10way).

### 5.1 Donor-Level: 7 Representative Targets

> **Figure 11** (interactive): CytoSig vs LinCytoSig vs SecAct grouped bar chart for 7 representative celltype–cytokine pairs at donor level. LinCytoSig uses the biologically matched cell-type signature (e.g., Macrophage__IL6 for IL6).

**Key findings (donor-level, 7 targets):**

- **LinCytoSig wins the median in 2 of 6 datasets** (Inflammation Atlas Main, scAtlas Normal) — datasets with rich immune cell diversity.
- **Per-target:** LinCytoSig wins 15/42 comparisons vs CytoSig's 11/42 (SecAct: 16/42). IL2 and IFNG are LinCytoSig's strongest targets; Fibroblast__TGFB1 is the main drag.
- **SecAct wins the median in 4 of 6 datasets** (GTEx, TCGA, CIMA, scAtlas Cancer), driven by its broad coverage of secreted proteins.

### 5.2 Celltype-Level Evaluation

> **Figure 12** (interactive): CytoSig vs LinCytoSig vs SecAct evaluated on matched cell-type pseudobulk (e.g., IL6 on macrophage pseudobulk) in scAtlas Normal and scAtlas Cancer. Only these two datasets have sufficient cell-type diversity for all 7 targets.

**Key findings (celltype-level):**

- **LinCytoSig excels where biology matches:** IL6×Macrophage, VEGFA×Endothelial, and IL2×CD8T in scAtlas Normal — cytokines with strong cell-type specificity to the tested cell type.
- **Advantage does NOT transfer to cancer:** LinCytoSig loses for IL6 and VEGFA in scAtlas Cancer, likely because tumor-associated cells have different transcriptional programs.
- **CytoSig wins TNFA and TGFB1 consistently** across all datasets at celltype level — the global signature captures these signaling responses better even on the "right" cell type.

CIMA and Inflammation Atlas Main lack macrophage, endothelial, and fibroblast annotations — only 2–3 of 7 targets evaluable using Myeloid/Mono proxies. See [full celltype analysis](lincytosig_issues.html#celltype7).

### 5.3 Limitations: Experimental Bias in the CytoSig Database

The CytoSig database has systematic experimental bias: it preferentially acquires specific cell types (macrophage, fibroblast, cancer lines) and specific cytokines (IFNG, TGFB1, IL6). Cell types with fewer than 10 experiments per cytokine produce noisy signatures, making it difficult to systematically select the "best" cell-type-specific signature for each target.

- **45% of cytokines get different "best" cell type** depending on the selection dataset (GTEx vs TCGA), indicating that best-selection is unstable (see [Issue 1](lincytosig_issues.html#issue1)).
- **The approach cannot be validated generally** with current data — but shows promise for specific high-quality targets where the biological match between signature origin and evaluation context is strong.

**The cell-type-specific approach has merit for well-characterized targets (IL6×Macrophage, VEGFA×Endothelial, IFNG×Macrophage) but the current CytoSig database cannot support systematic evaluation.** Future work could use perturbation data (parse_10M: 90 cytokines × 12 donors × 18 cell types) or expanded cell-type databases to revisit this question. See [full issues & analysis](lincytosig_issues.html).

---

## 6. Key Takeaways for Scientific Discovery

### 6.1 What CytoAtlas Enables That Other Tools Cannot

1. **Quantitative cytokine activity per cell type per disease:** Not just "which genes are differentially expressed" but "which cytokine signaling pathways are active, and how much."

2. **Cross-disease comparison:** The same 43 CytoSig signatures measured identically across 20 diseases, 35 organs, and 15 cancer types — enabling systematic comparison.

3. **Perturbation validation:** parse_10M provides 90 cytokine perturbations × 12 donors × 18 cell types for evaluating whether predicted activity matches known perturbation conditions.


### 6.2 Limitations (Honest Assessment)

1. **Linear model limitation:** Ridge regression cannot capture non-linear interactions between cytokines. If IFNG and TNFA synergize, CytoAtlas scores them independently.

2. **Transcriptomics-only:** Post-translational regulation (protein stability, secretion, receptor binding) is invisible. Some targets (e.g., CD40L) show negative correlations in specific datasets, which may relate to post-transcriptional biology not captured by mRNA-level inference.

3. **Signature matrix bias:** CytoSig is derived from bulk RNA-seq experiments. Cell types underrepresented in the source database (rare cell types, tissue-resident cells) have weaker signatures.

4. **Validation metric limitation:** Expression-activity correlation only validates targets where the gene itself is expressed. Downstream effector validation would require perturbation data (which we have for 90 cytokines via parse_10M).

### 6.3 Future Directions

1. **scGPT cohort integration** (~35M cells) — pending
2. **cellxgene Census integration** — pending
3. **Treatment response biomarkers** — using Inflammation Atlas responder/non-responder labels

---

## Figure Index

| Figure | File | Description |
|--------|------|-------------|
| 1 | `fig1_dataset_overview.png` | Dataset scale, signature matrices, validation layers |
| 2 | `fig2_summary_table.png` | Complete validation statistics table |
| 3 | `fig3_cross_dataset_boxplot.png` | CytoSig vs SecAct ρ distributions (independence-corrected) |
| 4a | `fig4a_gtex_per_tissue.png` | GTEx per-tissue stratified validation (29 tissues) |
| 4b | `fig4b_tcga_per_cancer.png` | TCGA per-cancer stratified validation (33 cancer types) |
| 4c | `fig4c_cross_platform_concordance.png` | Cross-platform concordance scatter (GTEx vs scAtlas) |
| 5 | `fig5_good_bad_correlations_cytosig.png` | Top/bottom 15 targets per atlas (CytoSig) |
| 6 | `fig6_cross_atlas_consistency.png` | Cross-atlas target consistency profiles |
| 7 | `fig7_validation_levels.png` | Aggregation level effect on correlation |
| 8 | `fig8_bio_targets_heatmap.png` | Heatmap of biologically important targets |
| 9 | `fig9_representative_scatter_cima.png` | Representative good/bad scatter plots (CIMA) |
| 10 | `fig10_method_comparison.png` | 10-way method comparison (CytoSig, LinCytoSig ×8, SecAct) |
| 11 | `fig11_lincytosig_vs_cytosig_scatter.png` | Matched target scatter (Lin vs Cyto) |
| 12 | `fig12_lincytosig_advantage_by_celltype.png` | Cell-type-specific LinCytoSig advantage |
| 13 | `fig13_secact_novel_signatures.png` | Novel SecAct high-correlation targets |
| 14 | `fig14_bulk_validation.png` | GTEx/TCGA bulk RNA-seq validation |
| 15 | `fig15_lincytosig_specificity.png` | LinCytoSig wins vs CytoSig wins |
| 16 | `fig16_celltype_scatter_examples.png` | Cell-type-level scatter examples |

All figures saved at: `report/figures/` (in the upstream `2cytoatlas` repository)
(Both PNG at 300 DPI and PDF vector formats)

---

## Appendix: Technical Specifications

### A. Computational Infrastructure
- **GPU:** NVIDIA A100 80GB (SLURM gpu partition)
- **Memory:** 256-512GB host RAM per node
- **Storage:** ~50GB results, 590MB DuckDB (atlas), 3-5GB DuckDB (perturbation), 2-4GB DuckDB (spatial)
- **Pipeline:** 24 Python scripts, 18 pipeline subpackages (~18.7K lines)
- **API:** 262 REST endpoints across 17 routers
- **Frontend:** 12 pages, 122 source files, 11.4K LOC

### B. Statistical Methods
- **Activity inference:** Ridge regression (λ=5×10⁵, z-score normalization, permutation-based significance)
- **Correlation:** Spearman rank correlation (robust to outliers)
- **Multiple testing:** Benjamini-Hochberg FDR (q < 0.05)
- **Bootstrap:** 100-1000 resampling iterations for confidence intervals
- **Differential:** Wilcoxon rank-sum test with effect size

### C. Data Availability
- Web interface: `http://[server]:8000/static/`
- API documentation: `http://[server]:8000/docs`
- All results: `results/` (in the upstream `2cytoatlas` repository)
- DuckDB databases: `cytoatlas-api/` (in the upstream `2cytoatlas` repository)
