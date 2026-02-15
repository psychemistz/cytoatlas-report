# CLAUDE.md

## Multi-Session Coordination

Multiple Claude Code sessions may be working in this repo concurrently. Follow these rules:

1. **Before editing any file**, read `COORDINATION.md` and check the File Lock Table. If another session owns the file, do NOT edit it.
2. **Before starting work**, claim your files in the File Lock Table by editing `COORDINATION.md`. Update your session status to reflect what you're working on.
3. **After finishing with a file**, release it in the File Lock Table (set Owner back to `—`).
4. **Never edit a file owned by the other session.** If you need changes in a locked file, post a message in the Messages section of `COORDINATION.md` and wait.
5. **`COORDINATION.md` is shared** — both sessions can read/write it, but keep edits to your own section and the Messages area.

## File Concordance (Three-Location Sync)

This repo has files that must stay in sync with the upstream pipeline's report directory. **After any edit to a synced file, copy it to the sync target. After regenerating HTML, copy that too.**

| Source of Truth (edit here) | Sync Target | What |
|-----------------------------|-------------|------|
| `scripts/generate_interactive_report.py` | `/vf/users/parks34/projects/2cytoatlas/report/generate_interactive_report.py` | Report generator script |
| `scripts/generate_report_figures.py` | `/vf/users/parks34/projects/2cytoatlas/report/generate_report_figures.py` | Static figure generator |
| `baseline/index.html` | `/vf/users/parks34/projects/2cytoatlas/report/REPORT.html` | Generated interactive report |

**Sync workflow after editing a script:**
```bash
# 1. Regenerate HTML
python scripts/generate_interactive_report.py
# 2. Copy HTML to both locations
cp /data/parks34/projects/2cytoatlas/report/REPORT.html baseline/index.html
# 3. Copy script to upstream
cp scripts/generate_interactive_report.py /vf/users/parks34/projects/2cytoatlas/report/generate_interactive_report.py
```

**Remote:** After committing, push to remote to keep GitHub Pages current.

## Git Commit Guidelines

- Do NOT include `Co-Authored-By` lines in commit messages.

## Repository Structure

This is the reporting/documentation companion to the upstream `2cytoatlas` pipeline repository (`../2cytoatlas`).

```
baseline/                           # Baseline report outputs
  index.html                        # Interactive HTML report (~13 MB, generated)
  REPORT.md                         # Markdown baseline report
  stats_section_4.1.html            # Stats supplement: all Section 4 statistical methods
scripts/
  generate_interactive_report.py    # Main report generator (~2900 lines, produces index.html)
  generate_report_figures.py        # Static figure generator
reports/weekly/                     # Weekly dataset analytics and design documents
```

### Key Documents

| Document | Description |
|----------|-------------|
| `reports/weekly/DATASET_ANALYTICS.md` | Cross-sample correlation design: dataset inventory, filtering, thresholds, gene ID handling, normalization, confounders |
| `baseline/REPORT.md` | Baseline validation report with existing results |
| `SESSION_HANDOFF.md` | Session handoff file for context continuity across sessions |

### Upstream Pipeline Reference

Pipeline scripts live in `../2cytoatlas/scripts/`. Key scripts for bulk validation:

| Script | Description |
|--------|-------------|
| `15_bulk_validation.py` | GTEx/TCGA bulk RNA-seq activity + correlations (format-aware: TPM vs RSEM) |
| `15b_tcga_primary_filter.py` | TCGA primary tumor filtering + per-cancer activity inference |
| `validation/config.py` | Shared thresholds (min_cells=10, min_samples=10/30) |

### Key Design Decisions

- **Two-stage threshold design**: Activity scripts generate pseudobulk with ALL groups (no min_cells filtering); validation scripts apply min_cells=10 at correlation time
- **Format-specific preprocessing**: GTEx uses TPM (log2(TPM+1)), TCGA uses EBPlusPlus RSEM normalized counts (clip negatives → log2(RSEM+1))
- **TCGA gene IDs**: `symbol|entrezID` format parsed by string split (no probemap needed)
- **GTEx gene IDs**: Versioned ENSG with 3-tier probemap fallback
- **Inflammation Atlas → Inflammation Main**: Only `inflammation_main` cohort (817 donors) is used in Section 4+ validation. The old `merge_inflammation_atlases()` that combined main/val/ext is removed. Pre-computed JSON files still use "Inflammation Atlas" keys — remapped at load time via `.replace()`.
- **scAtlas Cancer levels**: Uses `tumor_only`, `tumor_by_cancer`, `tumor_by_cancer_celltype1` (NOT `donor_organ` variants)

## Interactive Report Architecture

### Report Sections (Section 4: Validation Results)

| Section | Title | Figure | Interactive Features |
|---------|-------|--------|---------------------|
| 4.1 | Overall Performance Summary | Table 1 | Sortable table |
| 4.2 | Cross-Dataset Comparison | Figure 2 | Total/Matched tabs |
| 4.3 | Per-Tissue/Per-Cancer Stratified | Figure 3 | Total/Matched tabs, GTEx/TCGA dropdown |
| 4.4 | Cross-Platform Comparison | Figure 4 | GTEx/TCGA toggle, CytoSig/SecAct + Matched tabs |
| 4.5 | Best/Worst Correlated Targets | Figure 5 | Signature + atlas dropdowns |
| 4.6 | Cross-Atlas Consistency | Figure 6 | Line chart with legend toggle |
| 4.7 | Effect of Aggregation Level | Figure 7 | Total/Matched tabs, atlas dropdown |
| 4.8 | Representative Scatter Plots | Figure 8 | Atlas + target + signature dropdowns |
| 4.9 | Biologically Important Targets Heatmap | Figure 9 | CytoSig/SecAct tabs |
| 4.10 | Per-Target Correlation Rankings | Figure 10 | Dataset + signature dropdowns |

### Statistical Testing Pattern

All CytoSig vs SecAct comparisons use two tests:
- **Total (Mann-Whitney U)**: All CytoSig targets (~43) vs all SecAct targets (~1,170). Unpaired, unequal n.
- **Matched (Wilcoxon signed-rank)**: 32 shared targets (22 direct + 10 alias-resolved via ALIAS_MAP). Paired by target.
- **Multiple testing**: BH-FDR correction via `statsmodels.stats.multitest.multipletests(pvals, method='fdr_bh')`
- **JavaScript helpers**: `sigStars(p)` for star annotations, `formatPval(p)` for raw p-values, `formatQval(q)` for BH-corrected q-values

### Key Constants (in generate_interactive_report.py)

```python
ATLAS_ORDER = ['gtex', 'tcga', 'cima', 'inflammation_main', 'scatlas_normal', 'scatlas_cancer']
ATLAS_LABELS = ['GTEx', 'TCGA', 'CIMA', 'Inflammation Main', 'scAtlas (Normal)', 'scAtlas (Cancer)']
MATCHED_TARGETS = [22 direct] + [10 alias keys]      # CytoSig names
MATCHED_TARGETS_SECACT = [22 direct] + [10 alias values]  # SecAct gene symbols
```

### Embedded JavaScript in Python f-strings

The HTML template uses Python f-strings, so JavaScript curly braces must be **doubled**: `{{` and `}}`. Example:
```python
html += f"""
allLevels.forEach(function(level) {{
  if (sigData[level]) {{
    sigData[level].rhos.forEach(function(v) {{ y.push(v); }});
  }}
}});
"""
```

## GPU Node Environment

When running GPU-dependent computation (ridge_batch, CuPy, H5AD streaming) on the GPU node:

```bash
module load CUDA/12.8.1 cuDNN
source ~/bin/myconda && conda activate secactpy
```

This is required for `libcublas.so.12` and CuPy to work. The SLURM job that failed (11654934) was missing these module loads.
