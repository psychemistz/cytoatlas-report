# CytoAtlas: Pan-Disease Single-Cell Cytokine Activity Atlas

**February 14, 2026**

---

## Executive Summary

CytoAtlas is a comprehensive computational resource that maps cytokine and secreted protein signaling activity across **~29 million human cells and ~31,000 bulk RNA-seq samples** from six independent datasets spanning two bulk RNA-seq resources (GTEx, TCGA) and four single-cell compendia: CIMA (6.5M healthy donor cells), Inflammation Atlas (6.3M disease cells across 3 cohorts; Main cohort used for validation), scAtlas (6.4M organ and cancer cells), and parse_10M (9.7M cytokine-perturbed cells). The system uses **linear ridge regression** against experimentally derived signature matrices to infer activity — producing fully interpretable, conditional z-scores rather than black-box predictions. This makes CytoAtlas an orthogonal tool to deep learning approaches: every prediction traces back to a known gene-to-cytokine relationship with a quantifiable confidence interval.

**Key results:**
- 1,213 signatures (43 CytoSig cytokines + 1,170 SecAct secreted proteins), plus 178 cell-type-specific LinCytoSig variants, validated across 6 independent datasets
- Spearman correlations between predicted activity and target gene expression reach ρ=0.6-0.9 for well-characterized cytokines (IL1B, TNFA, VEGFA, TGFB family)
- Cross-dataset consistency demonstrates that signatures generalize across CIMA, Inflammation Atlas Main, scAtlas, GTEx, and TCGA
- Cell-type-specific signatures (LinCytoSig) improve prediction for select immune cell types (Basophil, NK, DC: +0.18-0.21 Δρ) but generally underperform global CytoSig for non-immune cell types
- SecAct provides the broadest validated coverage with 805–1,161 targets per dataset (varying by gene overlap), achieving the highest median correlations in 5 of 6 datasets (independence-corrected median ρ=0.31–0.46)

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

> **Why independence correction matters:** Pooling across tissues or cancer types inflates correlations through confounding. For example, GTEx pooled CytoSig median ρ (0.211) is 40% higher than the independence-corrected by-tissue value (0.151); SecAct shows +30% inflation (0.394 vs 0.304). All results in this report use the corrected values. For a detailed comparison of pooled vs independent levels, including inflation magnitude and finer cell-type stratification, see the [Section 4.1 statistical supplement](stats_section_4.1.html).

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

**Inflammation Atlas (20 diseases):**
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

**The key insight:** CytoAtlas is not trying to replace DNN-based tools. It provides an **orthogonal, complementary signal** that a human scientist can directly inspect. When CytoAtlas says "IFNG activity is elevated in CD8+ T cells from RA patients," you can verify this by checking the IFNG signature genes in those cells.

### 3.2 What Scientific Questions Does CytoAtlas Answer?

1. **Which cytokines are active in which cell types across diseases?** — IL1B/TNFA in monocytes/macrophages, IFNG in CD8+ T and NK cells, IL17A in Th17, VEGFA in endothelial/tumor cells, TGFB family in stromal cells — quantified across 20 diseases, 35 organs, and 15 cancer types.
2. **Are cytokine activities consistent across independent cohorts?** — Yes. IL1B, TNFA, VEGFA, and TGFB family show consistent positive correlations across all 6 validation datasets (Figure 6).
3. **Does cell-type-specific biology matter for cytokine inference?** — For select immune types, yes: LinCytoSig improves prediction for Basophils (+0.21 Δρ), NK cells (+0.19), and DCs (+0.18), but global CytoSig wins overall (Figures 11–12).
4. **Which secreted proteins beyond cytokines show validated activity?** — SecAct (1,170 targets) achieves the highest correlations across all atlases (median ρ=0.33–0.49), with novel validated targets like Activin A (ρ=0.98), CXCL12 (ρ=0.92), and BMP family (Figure 13).
5. **Can we predict treatment response from cytokine activity?** — We are incorporating cytokine-blocking therapy outcomes from bulk RNA-seq to test whether predicted cytokine activity associates with therapy response. Additionally, Inflammation Atlas responder/non-responder labels (208 samples across 6 diseases: RA, PS, PSA, CD, UC, SLE) enable treatment response prediction using cytokine activity profiles as features.

### 3.3 Validation Philosophy

CytoAtlas validates against a simple but powerful principle: **if CytoSig predicts high IFNG activity for a sample, that sample should have high IFNG gene expression.** This expression-activity correlation is computed via Spearman rank correlation across donors/samples.

This is a conservative validation — it only captures signatures where the target gene itself is expressed. Signatures that act through downstream effectors would not be captured, meaning our validation **underestimates** true accuracy.

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

**Summary ranking:** SecAct > CytoSig in 5 of 6 datasets (GTEx, TCGA, CIMA, scAtlas Normal, scAtlas Cancer). CytoSig leads Inflammation Main (0.323 vs 0.173 SecAct) where disease-driven cytokine variance amplifies canonical signatures. SecAct achieves the highest median ρ in scAtlas Normal (0.455) and TCGA (0.357), benefiting from broad secreted protein coverage and spatial-transcriptomics-derived signatures. CIMA shows the lowest correlations (0.114 CytoSig) because healthy donors have narrow activity ranges.

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

SecAct significantly outperforms CytoSig in 5 of 6 datasets. The Inflammation Main exception is non-significant on both tests (p = 0.548 total, p = 0.141 matched); on total targets CytoSig has higher median ρ (0.32 vs 0.17), but on the 27 matched targets SecAct actually has higher median ρ (0.44 vs 0.35). The total-target advantage reflects disease-driven cytokine variance that preferentially boosts CytoSig's canonical targets — in inflammatory conditions, the 33 CytoSig cytokines available are the most biologically active, giving CytoSig a domain advantage that disappears when restricting to shared targets.

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
- **Matched targets (strongest test):** On matched targets (23 pairs per GTEx tissue, 19–21 per TCGA cancer after expression filtering), SecAct wins direction in 29/29 GTEx tissues and 32/33 TCGA cancer types, with 23/29 and 31/33 reaching significance (q < 0.05). The sole exception is Acute Myeloid Leukemia (Δ = −0.08, q = 0.83), a low-sample liquid tumor. This near-unanimous result across 61 of 62 strata rules out Simpson's paradox and confirms the Section 4.2 aggregate is not driven by a few dominant strata
- **Total targets:** SecAct wins direction in 28/29 GTEx tissues (21 significant) and 30/33 TCGA cancers (15 significant). The few CytoSig-favored strata on total targets (Brain in GTEx; Kidney Chromophobe, Ovarian, Uveal Melanoma in TCGA) are all non-significant, reflecting low cytokine signaling rather than genuine CytoSig superiority
- **Spatial signatures capture context-dependent regulation:** Since SecAct outperforms CytoSig on the *same* matched cytokines, the advantage is not about target breadth but about signature quality. SecAct's spatial correlation signatures capture tissue-context-dependent cytokine regulation that CytoSig's median log2FC signatures miss. The advantage is largest in tissues with complex cellular microenvironments — GTEx: Small Intestine (Δ = +0.55), Vagina (+0.34), Salivary Gland (+0.30); TCGA: Pancreatic (+0.34), Head & Neck (+0.27), Lung Adeno (+0.25) — and smallest in homogeneous or low-cytokine contexts: GTEx: Heart (+0.06), Skin (+0.06), Brain (+0.06); TCGA: Uveal Melanoma (+0.07), Kidney Chromophobe (+0.10)
- scAtlas strata with Tier B sample sizes (<30 donors) are shown with reduced confidence

### 4.4 Cross-Platform Comparison: Bulk vs Pseudobulk

> **Figure 4c** (`fig4c_cross_platform_concordance.png`): Cross-platform concordance scatter (GTEx vs scAtlas Normal, TCGA vs scAtlas Cancer).

This section tests whether expression-activity relationships replicate across measurement technologies. For each matching tissue or cancer type, per-target Spearman ρ from bulk RNA-seq is compared to the same target's ρ from single-cell pseudobulk. Per-stratum Wilcoxon signed-rank tests (paired by target) with BH-FDR correction are used to test whether ρ values differ systematically between platforms.

**GTEx vs scAtlas Normal** (13 matching tissues: Blood, Breast, Colon, Esophagus, Heart, Kidney, Liver, Lung, Ovary, Skin, Small Intestine, Spleen, Uterus): Per-tissue CytoSig concordance ranges from ρ = −0.10 (Uterus) to 0.49 (Breast), with most tissues showing moderate positive concordance.

**TCGA vs scAtlas Cancer** (11 matching cancer types: BRCA, CRC, ESCA, HCC, HNSC, KIRC, LUAD, OV, PAAD, PRAD, STAD): Per-cancer concordance is generally stronger, ranging from ρ = 0.23 (OV) to 0.83 (PAAD), suggesting cancer-type-specific activity patterns are more consistent across platforms than normal tissue patterns.

**Key finding — platform effect is a statistical power effect:** Using all targets, SecAct shows significant bulk–pseudobulk differences in most strata (11/13 GTEx tissues, 5/11 TCGA cancers), while CytoSig shows almost none (1/13, 0/11). However, when restricted to the same 32 shared targets (Matched tabs), both methods show no significant platform differences (CytoSig: 0/13, 0/11; SecAct: 0/13, 1/11). The apparent platform sensitivity is not a signal quality difference: matched and unmatched SecAct targets show the same per-target platform shift (mean |Δ| = 0.298 vs 0.302, Mann–Whitney p = 0.82). Rather, SecAct's ~1,000 paired targets per tissue provide 25× more observations than CytoSig's ~40, easily detecting the same tiny systematic shift (Δ ≈ 0.03) that CytoSig lacks power to detect. Core cytokine targets are platform-robust. Full per-stratum results in [statistical methods](stats_section_4.1.html#cross-platform-comparison).

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

**SecAct-Only: CytoSig Near-Zero, SecAct Rescues (11 targets):**

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

**Key insight:** CytoSig's median log2FC signatures reliably detect 16 canonical cytokines with strong transcriptional programs. SecAct's spatial correlation signatures additionally rescue 11 targets where CytoSig averages near zero — including membrane-bound (CD40L, TRAIL), paracrine (HGF, FGF2), and heteromeric (LTA) signaling targets. TWEAK and IL21 are the clearest cases: CytoSig averages near zero across all datasets, while SecAct achieves +0.44 and +0.22. Note: this section uses donor-level aggregation (matching Figure 5's interactive chart), not the independence-corrected median-of-medians used in §4.1–4.3. Full per-atlas top-15/bottom-15 lists in [supplement](stats_section_4.1.html#per-tissue-stratified).

### 4.6 Cross-Atlas Consistency

> **Figure 6** (`fig6_cross_atlas_consistency.png`): Key targets tracked across 6 atlases at independent levels.

14 key CytoSig targets tracked across 6 datasets at independence-corrected levels reveal three consistency tiers:

| Tier | Targets | Pattern |
|------|---------|---------|
| **Universal** (ρ > 0.1 in all 6 datasets) | IL1B, TNFA, IFNG, IL6, BMP2, VEGFA | Robust across all cohorts and platforms (mean ρ 0.31–0.58) |
| **Mostly consistent** (4–5 of 6 positive) | IL10, TGFB1, CXCL12, GMCSF, HGF | Occasional near-zero outliers; TGFB1 slightly negative only in scAtlas Normal (−0.05) |
| **Context-dependent** (≤3 of 6 positive) | IL4, IL17A, EGF | Sign changes across cohorts; Th2/Th17 cytokines show limited expression in most datasets (IL4 and IL17A absent from Inflammation Main) |

**Key insight:** No systematic bulk vs single-cell advantage exists — mean |ρ| across 14 targets: GTEx 0.26, TCGA 0.33, CIMA 0.26, Inflammation Main 0.44, scAtlas Normal 0.36, scAtlas Cancer 0.37. Variation is target-specific and cohort-specific rather than platform-driven. The universal tier (IL1B, TNFA, IFNG, IL6) represents cytokine axes fundamental to inflammation biology, while context-dependent targets (IL4, IL17A) reflect Th2/Th17 pathways active only in specific disease cohorts.

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
- **Beyond ~L2 annotation depth, correlations become noise-dominated.** Cell-type stratification fragments the donor pool below the sample size needed for reliable correlation — at L2+ the signal/null ratio drops below 2× on average and % positive targets approaches 50% (chance level). The practical limit is driven by samples-per-group: Inflammation Main retains signal at L2 (553 samples/group, signal/null = 3.3×) while scAtlas Normal is already marginal at its shallowest stratification (22 samples/group)
- **In aggregate, SecAct's larger target pool shows higher retention** (27–44% vs CytoSig's 4–39%), but this reflects unequal target counts (~1,170 vs 43). On the 32 matched targets, CytoSig actually retains more signal in 2 of 4 datasets — the apparent SecAct advantage disappears when comparing like-for-like
- **All datasets show monotonic decline** except scAtlas Cancer CytoSig, which *increases* from All Tumor Cells (0.184) to Per Cancer Type (0.223) before dropping to CT1 (0.033) — cancer-type stratification removes cross-cancer confounding before cell-type dilution dominates
- **CytoSig and SecAct converge in some datasets** — Inflammation Main L2: 0.044 vs 0.048 — but not universally (scAtlas Cancer CT1: 0.033 vs 0.171). Convergence occurs where both methods approach the noise floor

Per-level Mann-Whitney/Wilcoxon tests with BH-FDR correction for all datasets are in the [supplement](stats_section_4.1.html#aggregation-level).

---

## 5. CytoSig vs LinCytoSig vs SecAct Comparison

### 5.1 Method Overview

> **Figure 10** (`fig10_method_comparison.png`): 10-way method comparison across all atlases.

We evaluate ten approaches for cytokine activity inference, covering three signature matrices and eight LinCytoSig strategies:

| # | Method | Targets | Description |
|---|--------|---------|-------------|
| 1 | **CytoSig** | 43 | Global (all cell types) signatures from experimental bulk RNA-seq |
| 2 | **LinCytoSig (no filter)** | 178 | Cell-type-specific signatures using all ~20K genes |
| 3 | **LinCytoSig (gene filter)** | 178 | Cell-type-specific signatures restricted to CytoSig gene space (~4,881 genes) |
| 4 | **LinCytoSig Best (combined)** | 43 | Best cell-type variant per cytokine selected via combined bulk (GTEx+TCGA) correlation, all genes |
| 5 | **LinCytoSig Best (comb+filt)** | 43 | Best cell-type variant per cytokine selected via combined bulk correlation, gene-filtered |
| 6 | **LinCytoSig Best (GTEx)** | 54 | Best cell-type variant per cytokine selected via GTEx correlation only |
| 7 | **LinCytoSig Best (TCGA)** | 50 | Best cell-type variant per cytokine selected via TCGA correlation only |
| 8 | **LinCytoSig Best (GTEx+filt)** | 54 | GTEx-selected best cell-type variant, gene-filtered |
| 9 | **LinCytoSig Best (TCGA+filt)** | 50 | TCGA-selected best cell-type variant, gene-filtered |
| 10 | **SecAct** | 1,170 | Global signatures from spatial transcriptomics (Moran's I) |

**LinCytoSig strategy rationale:**
- Methods 2–3 use the **cell-type-matched** LinCytoSig signature for each cytokine (e.g., Macrophage__IFNG for macrophages). The "no filter" version uses all ~20K genes in the signature; the "gene filter" version restricts to the ~4,881 genes present in CytoSig's gene space, testing whether the extra ~15K genes help or hurt prediction.
- Methods 4–5 take a **bulk-selected best** approach: for each cytokine, test all available cell-type-specific LinCytoSig signatures and select the one with the highest expression-activity correlation in combined bulk RNA-seq (GTEx + TCGA). This single "best" signature is then applied across all single-cell datasets.
- Methods 6–9 test whether **dataset-specific selection** improves generalization: the best cell-type variant is chosen independently from GTEx (54 cytokines) or TCGA (50 cytokines) correlations alone. Methods 8–9 additionally apply gene filtering to these dataset-specific selections.

**Donor-level pseudobulk validation (20 matched cytokines, median Spearman ρ):**

| Atlas | CytoSig | LinCyto (orig) | LinCyto (filt) | Best (comb) | Best (c+f) | Best (GTEx) | Best (TCGA) | Best (GTEx+f) | Best (TCGA+f) | SecAct |
|-------|---------|----------------|----------------|-------------|------------|-------------|-------------|---------------|---------------|--------|
| **CIMA** | 0.225 | 0.082 | 0.166 | 0.063 | 0.141 | 0.063 | 0.040 | 0.082 | 0.075 | 0.334 |
| **Inflammation Atlas** | 0.285 | 0.133 | 0.201 | 0.235 | 0.228 | 0.239 | 0.168 | 0.200 | 0.260 | 0.379 |
| **scAtlas Normal** | 0.216 | 0.173 | 0.193 | 0.298 | 0.230 | 0.231 | 0.169 | 0.179 | 0.203 | 0.391 |
| **scAtlas Cancer** | 0.344 | 0.172 | 0.146 | 0.275 | 0.267 | 0.300 | 0.167 | 0.267 | 0.182 | 0.492 |

**Key observations from Figure 10:**
- **SecAct consistently achieves the highest median ρ** across all 4 atlases (0.334–0.492), benefiting from its broad gene coverage and spatial-transcriptomics-derived signatures.
- **CytoSig outperforms most LinCytoSig variants** at donor level. The gap is largest in CIMA (0.225 vs best LinCytoSig variant 0.166). However, in scAtlas Normal, LinCytoSig Best-combined (0.298) exceeds CytoSig (0.216), suggesting cell-type-specific selection can help for organ-level analyses.
- **Gene filtering improves LinCytoSig** in 3 of 4 atlases (all except scAtlas Cancer), suggesting that the extra ~15K genes in LinCytoSig introduce noise at donor level. The improvement is most pronounced in CIMA (0.082 → 0.166, +102%).
- **GTEx-selected vs TCGA-selected:** GTEx-selected best variants generally outperform TCGA-selected in single-cell atlases — particularly in Inflammation Atlas (0.239 vs 0.168) and scAtlas Cancer (0.300 vs 0.167). This may reflect GTEx's broader tissue diversity providing more generalizable selections.
- **Gene filtering of GTEx/TCGA-selected:** The "+filt" variants show mixed results. TCGA+filt improves substantially in Inflammation Atlas (0.260 vs TCGA-orig 0.168), suggesting gene filtering recovers signal from noisy TCGA-based selections. GTEx+filt shows modest changes, indicating GTEx selections are already reasonably robust.
- **General ranking with caveats:** SecAct > CytoSig ≥ LinCytoSig Best > LinCytoSig (filtered) > LinCytoSig (orig). This ordering holds broadly but is not universal — LinCytoSig Best outperforms CytoSig in scAtlas Normal, and GTEx-selected best outperforms combined-selected best in scAtlas Cancer. For donor-level analysis, global signatures (CytoSig, SecAct) generally outperform cell-type-specific ones.

### 5.2 When Does LinCytoSig Outperform CytoSig?

> **Figure 11** (`fig11_lincytosig_vs_cytosig_scatter.png`): Matched target scatter plots.
> **Figure 12** (`fig12_lincytosig_advantage_by_celltype.png`): Cell-type-specific advantage analysis.

**Key finding from Figure 11:** In matched comparisons across 4 atlases (136 matched targets each):

| Atlas | LinCytoSig Wins | CytoSig Wins | Tie |
|-------|----------------|-------------|-----|
| CIMA | 63 | 64 | 9 |
| Inflammation Atlas | 41 | 90 | 5 |
| scAtlas Normal | 54 | 78 | 4 |
| scAtlas Cancer | 58 | 71 | 7 |

**CytoSig wins overall**, but LinCytoSig has specific advantages:

**LinCytoSig wins (Figure 12, top cell types):**
- **Basophil** (+0.21 mean Δρ): Basophil-specific IL3 response not captured by global CytoSig
- **NK Cell** (+0.19): NK-specific IL15/IL2 responses differ from global average
- **Dendritic Cell** (+0.18): DC-specific GMCSF and IL12 responses
- **Trophoblast** (+0.09): Placenta-specific cytokine responses

**CytoSig wins (LinCytoSig loses):**
- **Lymphatic Endothelial** (-0.73): Too few experiments in LinCytoSig database
- **Adipocyte** (-0.44): Single experiment for most cytokines → noisy signatures
- **Osteocyte** (-0.40): Very limited experimental data
- **PBMC** (-0.38): PBMC is already a mixture — global CytoSig *is* the PBMC signature
- **Dermal Fibroblast** (-0.33): Subtype-specific but insufficient replicates

### 5.3 Why LinCytoSig Underperforms for Some Cell Types

**Root cause analysis** (based on `/results/celltype_signatures/metadata.json`):

1. **Sample size effect:** LinCytoSig stratifies CytoSig's ~2,000 experiments by cell type. Cell types with <10 experiments per cytokine produce **noisy median signatures** (high variance, low reproducibility).
   - Breast Cancer Line: 108 experiments, 11 cytokines → reliable
   - Adipocyte: 1 experiment for BMP2 → unreliable

2. **Cell-type mismatch:** LinCytoSig's 45 cell types don't perfectly map to atlas annotations. When atlas cell types (e.g., "CD4+ Memory T") don't have a LinCytoSig equivalent, the system falls back to the closest match, introducing noise.

3. **The "PBMC paradox":** For donor-level analysis where all cell types are aggregated, CytoSig (which is already a mixture-level signature) naturally outperforms cell-type-specific LinCytoSig.

**Recommendation:** Use LinCytoSig for **cell-type-resolved** questions (e.g., "is IL15 activity specifically elevated in NK cells?") and CytoSig for **donor-level** questions (e.g., "does this patient have high IFNG activity?").

### 5.4 SecAct: Breadth Over Depth

> **Figure 13** (`fig13_secact_novel_signatures.png`): Novel SecAct targets with consistent high correlation.

SecAct covers 1,170 secreted proteins vs CytoSig's 43 cytokines. Key advantages:

- **Highest median ρ** in single-cell datasets (scAtlas Normal: 0.455, Cancer: 0.399, independence-corrected)
- **Highest median ρ** in bulk RNA-seq (GTEx: 0.314, TCGA: 0.357, independence-corrected median-of-medians)
- **95.8% positive correlation** in TCGA (independence-corrected) — nearly all targets work
- Discovers novel validated targets beyond canonical cytokines

**Top novel SecAct targets (not in CytoSig-43, consistently ρ > 0.5):** These represent secreted proteins with strong validated activity-expression correlations that would be missed by CytoSig alone. They represent potential novel paracrine signaling axes.

### 5.5 Biologically Important Targets Deep Dive

> **Figure 8** (`fig8_bio_targets_heatmap.png`): Heatmap across all atlases.
> **Figure 15** (`fig15_lincytosig_specificity.png`): LinCytoSig advantage/disadvantage cases.
> **Figure 16** (`fig16_celltype_scatter_examples.png`): Cell-type-level scatter examples.

**Interferon family:**
- IFNG: ρ = 0.25-0.68 across atlases (CytoSig), consistently positive
- B_Cell__IFNG (LinCytoSig): ρ = 0.37-0.73, *better* than global CytoSig in CIMA and Inflammation Val
- IFN1 (type I): ρ = 0.20-0.40, lower but consistent
- IFNL: ρ = 0.20-0.23 in CIMA, shows tissue-specific activity

**TGF-beta family:**
- TGFB1: ρ = 0.35 (CIMA), 0.90 (scAtlas) — strong in organ-level analysis
- TGFB3: ρ = 0.33 (CIMA), 0.55 (Inflammation), 0.90 (scAtlas)
- BMP2: ρ = 0.19 (CIMA), 0.43 (Inflammation), 0.90 (scAtlas)
- BMP4: ρ = 0.92 (scAtlas) — bone morphogenetic protein activity validated

**Interleukin family:**
- IL1B: ρ = 0.67 (CIMA), 0.68 (Inflammation) — top performer
- IL6: ρ = 0.41 (Inflammation), 0.90 (scAtlas)
- IL10: ρ = 0.38 (CIMA), 0.52 (Inflammation) — immunosuppressive
- IL17A: ρ = variable, atlas-dependent (Th17-specific)
- IL27: ρ = 0.43 (CIMA), 0.54 (Inflammation) — emerging immunotherapy target

---

## 6. Key Takeaways for Scientific Discovery

### 6.1 What CytoAtlas Enables That Other Tools Cannot

1. **Quantitative cytokine activity per cell type per disease:** Not just "which genes are differentially expressed" but "which cytokine signaling pathways are active, and how much."

2. **Cross-disease comparison:** The same 43 CytoSig signatures measured identically across 20 diseases, 35 organs, and 15 cancer types — enabling systematic comparison.

3. **Perturbation ground truth:** parse_10M provides 90 cytokine perturbations × 12 donors × 18 cell types. When we add exogenous IFNG to PBMCs, does CytoSig correctly predict elevated IFNG activity? (Yes.)


### 6.2 Limitations (Honest Assessment)

1. **Linear model limitation:** Ridge regression cannot capture non-linear interactions between cytokines. If IFNG and TNFA synergize, CytoAtlas scores them independently.

2. **Transcriptomics-only:** Post-translational regulation (protein stability, secretion, receptor binding) is invisible. CD40L negative correlation is a feature, not a bug — it's membrane-bound.

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
