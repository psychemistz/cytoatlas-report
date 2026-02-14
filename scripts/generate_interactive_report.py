#!/usr/bin/env python3
"""
CytoAtlas PI Report — Interactive HTML Generation
===================================================
Generates REPORT.html with embedded Plotly.js interactive figures.

Changes from previous version:
  0) Remove floating TOC
  1) Update dataset references
  2) Fix parse_10M description (not ground truth)
  3) Update SecAct reference (Ru et al.)
  4) Replace DNN-based → DL-based
  5) Embed summary table as HTML
  6) Interactive Plotly boxplot (all methods × all atlases)
  7) Add TWEAK duplicate explanation
  8) Add gene mapping verification note
  9) Consistent atlas ordering (bulk, CIMA, Inflammation, scAtlas)
 10) Interactive consistency plot
 11) Add LinCytoSig/SecAct to aggregation levels
 12) Interactive heatmap with tabs
 13) Interactive bulk validation with dropdown
 14) Interactive cell-type scatter with dropdown
"""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

warnings.filterwarnings('ignore')

# ─── Paths ────────────────────────────────────────────────────────────────────
BASE = Path('/data/parks34/projects/2cytoatlas')
CORR_DIR = BASE / 'results' / 'cross_sample_validation' / 'correlations'
VIZ_DIR = BASE / 'visualization' / 'data' / 'validation'
LINCYTO_FILT_DIR = BASE / 'results' / 'lincytosig_gene_filter'
BEST_SELECTION = VIZ_DIR / 'best_lincytosig_selection.json'
SCATTER_DONOR = VIZ_DIR / 'donor_scatter'
SCATTER_CT = VIZ_DIR / 'celltype_scatter'
LEVEL_COMPARISON_PATH = VIZ_DIR / 'level_comparison.json'
REPORT_DIR = BASE / 'report'

# ─── Atlas ordering (item 9): bulk, CIMA, Inflammation, scAtlas ───────────────
ATLAS_ORDER = [
    'gtex', 'tcga', 'cima', 'inflammation_main',
    'scatlas_normal', 'scatlas_cancer',
]
ATLAS_LABELS = [
    'GTEx', 'TCGA', 'CIMA', 'Inflammation Main',
    'scAtlas (Normal)', 'scAtlas (Cancer)',
]

LEVEL_MAP = {
    'gtex': 'donor_only', 'tcga': 'donor_only',
    'cima': 'donor_only', 'inflammation_main': 'donor_only',
    'scatlas_normal': 'donor_only', 'scatlas_cancer': 'tumor_only',
}

# Independence-corrected levels (per-tissue/cancer for GTEx/TCGA)
INDEPENDENT_LEVEL_MAP = {
    'gtex': ('by_tissue', True),          # (level, needs_median_of_medians)
    'tcga': ('primary_by_cancer', True),
    'cima': ('donor_only', False),
    'inflammation_main': ('donor_only', False),
    'scatlas_normal': ('donor_only', False),
    'scatlas_cancer': ('tumor_only', False),
}

SIG_COLORS = {
    'cytosig': '#2563EB',
    'lincyto_best': '#92400E',
    'lincytosig': '#D97706',
    'secact': '#059669',
}

BIO_FAMILIES = {
    'Interferon': ['IFNG', 'IFN1', 'IFNL'],
    'TGF-beta': ['TGFB1', 'TGFB2', 'TGFB3', 'BMP2', 'BMP4', 'BMP6', 'GDF11'],
    'Interleukin': ['IL1A', 'IL1B', 'IL2', 'IL4', 'IL6', 'IL10', 'IL12', 'IL13',
                    'IL15', 'IL17A', 'IL21', 'IL22', 'IL27', 'IL33'],
    'TNF': ['TNFA', 'LTA', 'TRAIL', 'TWEAK', 'CD40L'],
    'Growth Factor': ['EGF', 'FGF2', 'HGF', 'VEGFA', 'PDGFB', 'IGF1'],
    'Chemokine': ['CXCL12'],
    'Colony-Stimulating': ['GMCSF', 'GCSF', 'MCSF'],
}
TARGET_TO_FAMILY = {}
for fam, targets in BIO_FAMILIES.items():
    for t in targets:
        TARGET_TO_FAMILY[t] = fam

FAMILY_COLORS = {
    'Interferon': '#DC2626', 'TGF-beta': '#2563EB', 'Interleukin': '#059669',
    'TNF': '#D97706', 'Growth Factor': '#8B5CF6', 'Chemokine': '#EC4899',
    'Colony-Stimulating': '#6366F1', 'Other': '#9CA3AF',
}


# ═══════════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ═══════════════════════════════════════════════════════════════════════════════

def load_all_correlations():
    dfs = []
    for f in CORR_DIR.glob('*_correlations.csv'):
        if f.name == 'all_correlations.csv' or 'resampled' in f.name or 'summary' in f.name:
            continue
        dfs.append(pd.read_csv(f, low_memory=False))
    bulk = pd.read_csv(CORR_DIR / 'all_correlations.csv', low_memory=False)
    dfs.append(bulk)
    merged = pd.concat(dfs, ignore_index=True)
    merged = merged.drop_duplicates(subset=['target', 'celltype', 'atlas', 'level', 'signature'])
    return merged


def load_method_comparison():
    with open(VIZ_DIR / 'method_comparison.json') as f:
        return json.load(f)


def load_scatter(directory, filename):
    fp = directory / filename
    if fp.exists():
        with open(fp) as f:
            return json.load(f)
    return None


def merge_inflammation_atlases(df):
    """Replace inflammation_main/val/ext with single 'inflammation' atlas."""
    mask = df['atlas'].isin(['inflammation_main', 'inflammation_val', 'inflammation_ext'])
    inflam = df[mask].copy()
    inflam['atlas'] = 'inflammation'
    return pd.concat([df[~mask], inflam], ignore_index=True)


def load_lincyto_best_data():
    """Load LinCytoSig Best (comb+filt) as signature='lincyto_best'.

    Reads gene-filtered CSVs, keeps only 43 best_filt targets,
    renames targets from CellType__Cytokine to Cytokine.
    """
    if not BEST_SELECTION.exists():
        return pd.DataFrame()
    with open(BEST_SELECTION) as f:
        best_filt = json.load(f)['best_filt']
    reverse_map = {v: k for k, v in best_filt.items()}
    best_targets = set(best_filt.values())
    dfs = []
    for fp in LINCYTO_FILT_DIR.glob('*_lincyto_filt_correlations.csv'):
        df = pd.read_csv(fp, low_memory=False)
        df = df[df['target'].isin(best_targets)].copy()
        df['lincyto_best_source'] = df['target']
        df['target'] = df['target'].map(reverse_map)
        df['signature'] = 'lincyto_best'
        dfs.append(df)
    if not dfs:
        return pd.DataFrame()
    merged = pd.concat(dfs, ignore_index=True)
    for col in ['celltype', 'lincytosig_celltype', 'matched_atlas_celltypes']:
        if col not in merged.columns:
            merged[col] = None
    return merged


# ═══════════════════════════════════════════════════════════════════════════════
# DATA PREPARATION FOR INTERACTIVE FIGURES
# ═══════════════════════════════════════════════════════════════════════════════

def prepare_summary_table(df):
    """Prepare summary statistics as list of dicts for HTML table (Section 4.1).

    Uses independence-corrected methodology:
    - GTEx/TCGA: median-of-medians (one rho per target = median across tissues/cancers)
    - Other datasets: direct donor_only/tumor_only values
    Two methods: CytoSig, SecAct.

    Output:
        List of dicts embedded as JSON in the interactive HTML.
    """
    from statsmodels.stats.multitest import multipletests

    rows = []
    for sig_type in ['cytosig', 'secact']:
        for atlas in ATLAS_ORDER:
            level, needs_mom = INDEPENDENT_LEVEL_MAP[atlas]
            label = ATLAS_LABELS[ATLAS_ORDER.index(atlas)]
            sub = df[(df['atlas'] == atlas) & (df['level'] == level) & (df['signature'] == sig_type)]
            if len(sub) == 0:
                continue

            if needs_mom:
                # Median-of-medians: exclude 'all'/'unmappable', per-target median
                sub = sub[~sub['celltype'].isin(['all', 'unmappable'])]
                n_strata = sub['celltype'].nunique()
                if 'tissue' in level:
                    primary_level = f'by_tissue ({n_strata} tissues)'
                elif 'cancer' in level:
                    primary_level = f'{level} ({n_strata} types)'
                else:
                    primary_level = f'{level} ({n_strata} groups)'
                target_medians = sub.groupby('target')['spearman_rho'].median().dropna()
                rhos = target_medians
                n_sig = 0
                if 'spearman_pval' in sub.columns:
                    target_pvals = sub.groupby('target')['spearman_pval'].median().dropna()
                    if len(target_pvals) > 0:
                        _, fdr_qvals, _, _ = multipletests(target_pvals.values, method='fdr_bh')
                        n_sig = int((fdr_qvals < 0.05).sum())
            else:
                rhos = sub['spearman_rho'].dropna()
                # Build primary_level label with sample/donor count
                n_samples_val = int(sub['n_samples'].iloc[0]) if 'n_samples' in sub.columns and len(sub) > 0 else None
                if n_samples_val:
                    primary_level = f'{level} ({n_samples_val} donors)'
                else:
                    primary_level = level
                n_sig = 0
                if 'spearman_pval' in sub.columns:
                    pvals = sub['spearman_pval'].dropna().values
                    if len(pvals) > 0:
                        _, fdr_qvals, _, _ = multipletests(pvals, method='fdr_bh')
                        n_sig = int((fdr_qvals < 0.05).sum())

            n_pos = (rhos > 0).sum()
            n_total = len(rhos)
            sig_label = sig_type.upper()
            rows.append({
                'atlas': label,
                'signature': sig_label,
                'primary_level': primary_level,
                'n_targets': int(n_total),
                'median_rho': round(float(rhos.median()), 3),
                'mean_rho': round(float(rhos.mean()), 3),
                'std_rho': round(float(rhos.std()), 3),
                'min_rho': round(float(rhos.min()), 3),
                'max_rho': round(float(rhos.max()), 3),
                'pct_sig': round(100 * n_sig / n_total, 1) if n_total > 0 else 0,
                'pct_pos': round(100 * n_pos / n_total, 1) if n_total > 0 else 0,
            })
    return rows


# 22 direct matches (same target name in CytoSig and SecAct)
MATCHED_TARGETS_DIRECT = [
    'BDNF', 'BMP2', 'BMP4', 'BMP6', 'CXCL12', 'FGF2', 'GDF11', 'HGF',
    'IFNG', 'IL10', 'IL15', 'IL1A', 'IL1B', 'IL21', 'IL27', 'IL6',
    'LIF', 'LTA', 'OSM', 'TGFB1', 'TGFB3', 'VEGFA',
]

# 10 alias matches (CytoSig common name → SecAct HGNC gene symbol)
ALIAS_MAP = {
    'Activin A': 'INHBA', 'CD40L': 'CD40LG', 'GCSF': 'CSF3', 'GMCSF': 'CSF2',
    'IFNL': 'IFNL1', 'IL36': 'IL36A', 'MCSF': 'CSF1', 'TNFA': 'TNF',
    'TRAIL': 'TNFSF10', 'TWEAK': 'TNFSF12',
}

# All 32 matched targets: CytoSig names + alias CytoSig names
MATCHED_TARGETS = MATCHED_TARGETS_DIRECT + list(ALIAS_MAP.keys())
# Corresponding SecAct names
MATCHED_TARGETS_SECACT = MATCHED_TARGETS_DIRECT + list(ALIAS_MAP.values())


def prepare_boxplot_data(df):
    """Prepare rho arrays for section 4.2 boxplot with two tabs (Total / Matched).

    Uses independence-corrected methodology:
    - GTEx/TCGA: median-of-medians (one rho per target = median across tissues/cancers)
    - CIMA/Inflammation Main/scAtlas: direct donor_only/tumor_only rho values
    Matched comparison uses 32 shared targets (22 direct + 10 alias-resolved).

    Output:
        Dict of {atlas_label: {cytosig, secact, cytosig_matched, secact_matched,
        stats_total, stats_matched}} embedded as JSON in the interactive HTML.
        Stats include BH-FDR corrected q-values across 6 atlases.
    """
    from statsmodels.stats.multitest import multipletests
    reverse_alias = {v: k for k, v in ALIAS_MAP.items()}
    result = {}
    for atlas in ATLAS_ORDER:
        level, needs_mom = INDEPENDENT_LEVEL_MAP[atlas]
        label = ATLAS_LABELS[ATLAS_ORDER.index(atlas)]
        result[label] = {}

        # --- Total tab: all targets ---
        for sig_type in ['cytosig', 'secact']:
            sub = df[(df['atlas'] == atlas) & (df['level'] == level) & (df['signature'] == sig_type)]
            if needs_mom:
                # Median-of-medians: exclude 'all'/'unmappable', one rho per target
                sub = sub[~sub['celltype'].isin(['all', 'unmappable'])]
                rhos = sub.groupby('target')['spearman_rho'].median().dropna().tolist()
            else:
                rhos = sub['spearman_rho'].dropna().tolist()
            result[label][sig_type] = [round(r, 4) for r in rhos]

        # --- Matched targets (32 shared between CytoSig and SecAct) ---
        cytosig_sub = df[(df['atlas'] == atlas) & (df['level'] == level) & (df['signature'] == 'cytosig')]
        secact_sub = df[(df['atlas'] == atlas) & (df['level'] == level) & (df['signature'] == 'secact')]
        if needs_mom:
            cytosig_sub = cytosig_sub[~cytosig_sub['celltype'].isin(['all', 'unmappable'])]
            secact_sub = secact_sub[~secact_sub['celltype'].isin(['all', 'unmappable'])]

        cytosig_matched = cytosig_sub[cytosig_sub['target'].isin(MATCHED_TARGETS)]
        secact_matched = secact_sub[secact_sub['target'].isin(MATCHED_TARGETS_SECACT)]

        if needs_mom:
            cs_matched_rhos = cytosig_matched.groupby('target')['spearman_rho'].median().dropna().tolist()
            sa_matched_rhos = secact_matched.groupby('target')['spearman_rho'].median().dropna().tolist()
        else:
            cs_matched_rhos = cytosig_matched['spearman_rho'].dropna().tolist()
            sa_matched_rhos = secact_matched['spearman_rho'].dropna().tolist()

        result[label]['cytosig_matched'] = [round(r, 4) for r in cs_matched_rhos]
        result[label]['secact_matched'] = [round(r, 4) for r in sa_matched_rhos]

        # --- Statistical tests ---
        cyto_rhos = np.array(result[label]['cytosig'])
        sec_rhos = np.array(result[label]['secact'])

        # Total: Mann-Whitney U (unpaired, different n)
        if len(cyto_rhos) >= 3 and len(sec_rhos) >= 3:
            u_stat, u_pval = scipy_stats.mannwhitneyu(cyto_rhos, sec_rhos, alternative='two-sided')
            result[label]['stats_total'] = {
                'test': 'Mann-Whitney U',
                'stat': round(float(u_stat), 2),
                'p_value': float(u_pval),
                'n_cytosig': len(cyto_rhos),
                'n_secact': len(sec_rhos),
            }
        else:
            result[label]['stats_total'] = None

        # Matched: Wilcoxon signed-rank (paired by target, 32 shared)
        # Compute per-target representative rho (median for mom, mean for direct)
        if needs_mom:
            cyto_m = cytosig_matched.groupby('target')['spearman_rho'].median().dropna()
            sec_m_raw = secact_matched.groupby('target')['spearman_rho'].median().dropna()
        else:
            cyto_m = cytosig_matched.groupby('target')['spearman_rho'].mean().dropna()
            sec_m_raw = secact_matched.groupby('target')['spearman_rho'].mean().dropna()

        # Map SecAct alias names back to CytoSig names for pairing
        sec_m = pd.Series({
            reverse_alias.get(t, t): v for t, v in sec_m_raw.items()
        })
        common = sorted(set(cyto_m.index) & set(sec_m.index))
        if len(common) >= 6:
            c_vals = np.array([float(cyto_m.loc[t]) for t in common])
            s_vals = np.array([float(sec_m.loc[t]) for t in common])
            w_stat, w_pval = scipy_stats.wilcoxon(c_vals, s_vals, alternative='two-sided')
            result[label]['stats_matched'] = {
                'test': 'Wilcoxon signed-rank',
                'stat': round(float(w_stat), 2),
                'p_value': float(w_pval),
                'n_pairs': len(common),
            }
        else:
            result[label]['stats_matched'] = None

    # --- BH-FDR correction across 6 atlases ---
    atlas_labels = [ATLAS_LABELS[ATLAS_ORDER.index(a)] for a in ATLAS_ORDER]

    # Mann-Whitney
    mw_entries = [(lbl, result[lbl]['stats_total']['p_value'])
                  for lbl in atlas_labels if result[lbl].get('stats_total')]
    if mw_entries:
        mw_lbls, mw_ps = zip(*mw_entries)
        _, mw_qs, _, _ = multipletests(list(mw_ps), method='fdr_bh')
        for lbl, q in zip(mw_lbls, mw_qs):
            result[lbl]['stats_total']['q_bh'] = round(float(q), 6)

    # Wilcoxon
    wx_entries = [(lbl, result[lbl]['stats_matched']['p_value'])
                  for lbl in atlas_labels if result[lbl].get('stats_matched')]
    if wx_entries:
        wx_lbls, wx_ps = zip(*wx_entries)
        _, wx_qs, _, _ = multipletests(list(wx_ps), method='fdr_bh')
        for lbl, q in zip(wx_lbls, wx_qs):
            result[lbl]['stats_matched']['q_bh'] = round(float(q), 6)

    return result


def prepare_stratified_data(df):
    """Per-tissue/per-cancer CytoSig vs SecAct comparison for GTEx and TCGA.

    For each tissue (GTEx by_tissue) or cancer type (TCGA primary_by_cancer),
    computes CytoSig vs SecAct rho distributions and statistical tests:
      - Total: Mann-Whitney U on all targets
      - Matched: Wilcoxon signed-rank on 32 shared targets
    BH-FDR correction applied across all strata within each dataset.
    """
    from statsmodels.stats.multitest import multipletests

    reverse_alias = {v: k for k, v in ALIAS_MAP.items()}

    datasets = {
        'GTEx': ('gtex', 'by_tissue'),
        'TCGA': ('tcga', 'primary_by_cancer'),
    }

    result = {}
    for ds_label, (atlas, level) in datasets.items():
        sub = df[(df['atlas'] == atlas) & (df['level'] == level)]
        strata = sorted([s for s in sub['celltype'].unique()
                         if s not in ('all', 'unmappable')])

        strata_data = {}
        mw_pvals = []
        wilcox_pvals = []
        strata_order = []

        for stratum in strata:
            st = sub[sub['celltype'] == stratum]

            # Total: all CytoSig vs all SecAct targets
            cytosig_rhos = st[st['signature'] == 'cytosig']['spearman_rho'].dropna()
            secact_rhos = st[st['signature'] == 'secact']['spearman_rho'].dropna()

            if len(cytosig_rhos) < 3 or len(secact_rhos) < 3:
                continue

            _, mw_p = scipy_stats.mannwhitneyu(
                cytosig_rhos.values, secact_rhos.values, alternative='two-sided')

            # Matched: 32 shared targets, paired via alias resolution
            cytosig_m_df = st[(st['signature'] == 'cytosig') &
                             (st['target'].isin(MATCHED_TARGETS))]
            secact_m_df = st[(st['signature'] == 'secact') &
                            (st['target'].isin(MATCHED_TARGETS_SECACT))]

            cyto_m = cytosig_m_df.set_index('target')['spearman_rho'].dropna()
            sec_m_raw = secact_m_df.set_index('target')['spearman_rho'].dropna()
            sec_m = pd.Series({reverse_alias.get(t, t): v
                               for t, v in sec_m_raw.items()})
            common = sorted(set(cyto_m.index) & set(sec_m.index))

            wilcox_p = None
            if len(common) >= 6:
                c_vals = np.array([float(cyto_m.loc[t]) for t in common])
                s_vals = np.array([float(sec_m.loc[t]) for t in common])
                _, wilcox_p = scipy_stats.wilcoxon(
                    c_vals, s_vals, alternative='two-sided')

            strata_order.append(stratum)
            mw_pvals.append(float(mw_p))
            wilcox_pvals.append(float(wilcox_p) if wilcox_p is not None
                                else None)

            strata_data[stratum] = {
                'n_cytosig': int(len(cytosig_rhos)),
                'n_secact': int(len(secact_rhos)),
                'cytosig': [round(r, 4) for r in cytosig_rhos.tolist()],
                'secact': [round(r, 4) for r in secact_rhos.tolist()],
                'cytosig_matched': [round(float(cyto_m.loc[t]), 4)
                                    for t in common],
                'secact_matched': [round(float(sec_m.loc[t]), 4)
                                   for t in common],
                'median_cytosig': round(float(cytosig_rhos.median()), 4),
                'median_secact': round(float(secact_rhos.median()), 4),
                'n_pairs': len(common),
                'stats_total': {
                    'test': 'Mann-Whitney U',
                    'p_raw': float(mw_p),
                    'n_cytosig': int(len(cytosig_rhos)),
                    'n_secact': int(len(secact_rhos)),
                },
                'stats_matched': {
                    'test': 'Wilcoxon signed-rank',
                    'p_raw': float(wilcox_p),
                    'n_pairs': len(common),
                } if wilcox_p is not None else None,
            }

        # BH-FDR correction for Mann-Whitney
        if mw_pvals:
            _, mw_qvals, _, _ = multipletests(mw_pvals, method='fdr_bh')
            for i, stratum in enumerate(strata_order):
                strata_data[stratum]['stats_total']['q_bh'] = round(
                    float(mw_qvals[i]), 6)

        # BH-FDR correction for Wilcoxon (only valid p-values)
        valid_idx = [i for i, p in enumerate(wilcox_pvals) if p is not None]
        if valid_idx:
            valid_pvals = [wilcox_pvals[i] for i in valid_idx]
            _, wilcox_qvals, _, _ = multipletests(valid_pvals, method='fdr_bh')
            for j, idx in enumerate(valid_idx):
                strata_data[strata_order[idx]]['stats_matched']['q_bh'] = round(
                    float(wilcox_qvals[j]), 6)

        # Sort by n_secact descending (proxy for data richness)
        sorted_strata = sorted(strata_order,
                               key=lambda s: strata_data[s]['n_secact'],
                               reverse=True)

        result[ds_label] = {
            'strata': sorted_strata,
            'data': strata_data,
        }

    return result


def prepare_method_comparison_boxplot(df):
    """Prepare rho arrays for section 5.1 boxplot — 10-way method comparison.

    Loads pre-computed comparison data (ridge regression on pseudobulk)
    for 4 combined atlases: CIMA, Inflammation Main, scAtlas Normal/Cancer.

    Ten methods:
      1. CytoSig — 43 cytokines, 4,881 curated genes (global)
      2. LinCytoSig (orig) — cell-type-matched signatures, all ~20K genes
      3. LinCytoSig (gene-filtered) — cell-type-matched, restricted to CytoSig 4,881 genes
      4. LinCytoSig (best-bulk, orig) — best per cytokine via combined GTEx+TCGA, all genes
      5. LinCytoSig (best-bulk, filtered) — best combined bulk, CytoSig genes only
      6. LinCytoSig (best GTEx) — best per cytokine via GTEx only, all genes
      7. LinCytoSig (best TCGA) — best per cytokine via TCGA only, all genes
      8. LinCytoSig (best GTEx, filtered) — GTEx-selected, CytoSig genes only
      9. LinCytoSig (best TCGA, filtered) — TCGA-selected, CytoSig genes only
     10. SecAct — 1,170 secreted proteins
    """
    eightway_path = VIZ_DIR / 'method_comparison_8way_all.json'
    if not eightway_path.exists():
        eightway_path = VIZ_DIR / 'method_comparison_6way_all.json'
    if not eightway_path.exists():
        print(f'  WARNING: no comparison JSON found — falling back to empty')
        return {}

    with open(eightway_path) as f:
        data = json.load(f)

    # Map atlas keys in JSON to display labels
    atlas_key_to_label = {
        'CIMA': 'CIMA',
        'Inflammation Atlas': 'Inflammation Main',
        'scAtlas Normal': 'scAtlas (Normal)',
        'scAtlas Cancer': 'scAtlas (Cancer)',
    }

    # Method keys match the JSON keys directly
    method_keys = ['cytosig', 'lincyto_orig', 'lincyto_filt',
                   'lincyto_best_orig', 'lincyto_best_filt',
                   'lincyto_best_gtex', 'lincyto_best_tcga',
                   'lincyto_best_gtex_filt', 'lincyto_best_tcga_filt', 'secact']

    result = {}
    for src_key, label in atlas_key_to_label.items():
        if src_key not in data:
            continue
        atlas_data = data[src_key]
        entry = {}
        for key in method_keys:
            vals = atlas_data.get(key, [])
            if vals:
                entry[key] = [round(v, 4) for v in vals]
        result[label] = entry
    return result


def prepare_consistency_data(df):
    """Prepare data for cross-atlas consistency line chart.

    Data source:
        Merged correlation CSVs for CytoSig and SecAct at donor level.
        14 key targets: IFNG, IL1B, TNFA, TGFB1, IL6, IL10, IL17A,
        IL4, BMP2, EGF, HGF, VEGFA, CXCL12, GMCSF.
    Method:
        For each target x signature: extract donor-level rho at each atlas
        in ATLAS_ORDER. CytoSig entries keyed by target name, SecAct by
        target_secact. Each entry includes rho array, cytokine family,
        family color, and sig_type.
    Output:
        Dict of {target_key: {rhos: [...], family, color, sig_type}}
        for the interactive consistency line chart (Section 4.5).
    """
    key_targets = ['IFNG', 'IL1B', 'TNFA', 'TGFB1', 'IL6', 'IL10', 'IL17A',
                   'IL4', 'BMP2', 'EGF', 'HGF', 'VEGFA', 'CXCL12', 'GMCSF']
    result = {}
    for sig_type in ['cytosig', 'secact']:
        sig_df = df[df['signature'] == sig_type]
        for target in key_targets:
            rhos = []
            for atlas in ATLAS_ORDER:
                level, needs_mom = INDEPENDENT_LEVEL_MAP[atlas]
                sub = sig_df[(sig_df['target'] == target) & (sig_df['atlas'] == atlas)
                              & (sig_df['level'] == level)]
                if needs_mom:
                    sub = sub[~sub['celltype'].isin(['all', 'unmappable'])]
                    if len(sub) > 0:
                        rhos.append(round(float(sub['spearman_rho'].median()), 4))
                    else:
                        rhos.append(None)
                elif len(sub) > 0:
                    rhos.append(round(float(sub['spearman_rho'].values[0]), 4))
                else:
                    rhos.append(None)
            family = TARGET_TO_FAMILY.get(target, 'Other')
            entry_key = target if sig_type == 'cytosig' else f'{target}_secact'
            result[entry_key] = {
                'rhos': rhos, 'family': family,
                'color': FAMILY_COLORS.get(family, '#9CA3AF'),
                'sig_type': sig_type,
            }
    return result


def prepare_heatmap_data(df):
    """Prepare heatmap matrices for CytoSig and SecAct using independence-corrected rho."""
    result = {}
    for sig_type in ['cytosig', 'secact']:
        sub = df[df['signature'] == sig_type].copy()
        sig_atlas_order = ATLAS_ORDER
        sig_atlas_labels = ATLAS_LABELS

        # Collect rows at independence-corrected level
        rows = []
        for atlas in sig_atlas_order:
            level, needs_mom = INDEPENDENT_LEVEL_MAP[atlas]
            atlas_sub = sub[(sub['atlas'] == atlas) & (sub['level'] == level)]
            if needs_mom:
                atlas_sub = atlas_sub[~atlas_sub['celltype'].isin(['all', 'unmappable'])]
                # Per-target median across tissues/cancers
                target_medians = atlas_sub.groupby('target')['spearman_rho'].median().reset_index()
                target_medians['atlas'] = atlas
                target_medians['level'] = level
                for _, row in target_medians.iterrows():
                    rows.append(row)
            else:
                for _, row in atlas_sub.iterrows():
                    rows.append(row)
        if not rows:
            result[sig_type] = {'targets': [], 'matrix': [], 'atlases': sig_atlas_labels}
            continue
        sub_df = pd.DataFrame(rows)

        if sig_type == 'cytosig':
            targets = sorted(sub_df['target'].unique())
        else:  # secact — show matched targets (using SecAct names) + top additional
            matched_set = set(MATCHED_TARGETS_SECACT)
            matched = sorted([t for t in sub_df['target'].unique() if t in matched_set])
            remaining = sub_df[~sub_df['target'].isin(matched_set)]
            if len(remaining) > 0:
                median_rhos = remaining.groupby('target')['spearman_rho'].median().sort_values(ascending=False)
                top_extra = list(median_rhos.head(25).index)
            else:
                top_extra = []
            targets = matched + sorted(top_extra)

        matrix = []
        for target in targets:
            row_vals = []
            for atlas in sig_atlas_order:
                level, _ = INDEPENDENT_LEVEL_MAP[atlas]
                match = sub_df[(sub_df['target'] == target) & (sub_df['atlas'] == atlas)]
                if len(match) > 0:
                    row_vals.append(round(float(match['spearman_rho'].values[0]), 3))
                else:
                    row_vals.append(None)
            matrix.append(row_vals)

        result[sig_type] = {'targets': targets, 'matrix': matrix, 'atlases': sig_atlas_labels}
    return result


def prepare_levels_data(df):
    """Prepare aggregation level comparison data for CytoSig and SecAct (section 4.6).

    Data source:
        Merged correlation CSVs filtered by atlas and aggregation level.
        4 single-cell atlases: CIMA (5 levels), Inflammation (3),
        scAtlas Normal (3), scAtlas Cancer (3).
    Method:
        For each atlas x level: extract CytoSig/SecAct rho arrays + matched
        (32 shared targets, paired via alias resolution).
        Statistical tests per level:
          - Total: Mann-Whitney U (CytoSig vs SecAct, all targets)
          - Matched: Wilcoxon signed-rank (32 paired shared targets)
        BH-FDR correction applied across levels within each atlas.
    Output:
        Nested dict with keys: cytosig, secact, cytosig_matched,
        secact_matched, stats (per-level test results with q_bh).
    """
    from statsmodels.stats.multitest import multipletests

    reverse_alias = {v: k for k, v in ALIAS_MAP.items()}

    configs = [
        ('cima', ['donor_only', 'donor_l1', 'donor_l2', 'donor_l3', 'donor_l4'], 'CIMA'),
        ('inflammation_main', ['donor_only', 'donor_l1', 'donor_l2'], 'Inflammation Main'),
        ('scatlas_normal', ['donor_organ', 'donor_organ_celltype1', 'donor_organ_celltype2'], 'scAtlas Normal'),
        ('scatlas_cancer', ['tumor_only', 'tumor_by_cancer', 'tumor_by_cancer_celltype1'], 'scAtlas Cancer'),
    ]
    level_labels = {
        'donor_only': 'Donor Only', 'donor_l1': 'Donor x L1', 'donor_l2': 'Donor x L2',
        'donor_l3': 'Donor x L3', 'donor_l4': 'Donor x L4',
        'donor_organ': 'Donor x Organ', 'donor_organ_celltype1': 'Donor x Organ x CT1',
        'donor_organ_celltype2': 'Donor x Organ x CT2',
        'tumor_only': 'Tumor Only', 'tumor_by_cancer': 'Tumor x Cancer',
        'tumor_by_cancer_celltype1': 'Tumor x Cancer x CT1',
    }
    result = {}
    for atlas, levels, title in configs:
        atlas_sub = df[df['atlas'] == atlas]
        cytosig_all = atlas_sub[atlas_sub['signature'] == 'cytosig']
        secact_all = atlas_sub[atlas_sub['signature'] == 'secact']

        cyto_data, sec_data = {}, {}
        cyto_m_data, sec_m_data = {}, {}
        mw_pvals, wilcox_pvals, level_order = [], [], []

        for level in levels:
            ll = level_labels.get(level, level)
            level_order.append(ll)

            # --- Total: all targets ---
            c_rhos = cytosig_all[cytosig_all['level'] == level]['spearman_rho'].dropna()
            s_rhos = secact_all[secact_all['level'] == level]['spearman_rho'].dropna()

            if len(c_rhos) > 0:
                cyto_data[ll] = {
                    'rhos': [round(r, 4) for r in c_rhos.tolist()],
                    'median': round(float(c_rhos.median()), 4),
                    'n': int(len(c_rhos)),
                }
            if len(s_rhos) > 0:
                sec_data[ll] = {
                    'rhos': [round(r, 4) for r in s_rhos.tolist()],
                    'median': round(float(s_rhos.median()), 4),
                    'n': int(len(s_rhos)),
                }

            # Mann-Whitney U test (Total)
            mw_p = None
            if len(c_rhos) >= 3 and len(s_rhos) >= 3:
                _, mw_p = scipy_stats.mannwhitneyu(
                    c_rhos.values, s_rhos.values, alternative='two-sided')
            mw_pvals.append(float(mw_p) if mw_p is not None else None)

            # --- Matched: 32 shared targets, paired via alias resolution ---
            cm_df = cytosig_all[(cytosig_all['level'] == level) &
                                (cytosig_all['target'].isin(MATCHED_TARGETS))]
            sm_df = secact_all[(secact_all['level'] == level) &
                               (secact_all['target'].isin(MATCHED_TARGETS_SECACT))]

            # Median across celltypes per target for pairing
            cyto_m = cm_df.groupby('target')['spearman_rho'].median().dropna()
            sec_m_raw = sm_df.groupby('target')['spearman_rho'].median().dropna()
            sec_m = pd.Series({reverse_alias.get(t, t): v
                               for t, v in sec_m_raw.items()})
            common = sorted(set(cyto_m.index) & set(sec_m.index))

            if len(common) > 0:
                c_vals = [float(cyto_m.loc[t]) for t in common]
                s_vals = [float(sec_m.loc[t]) for t in common]
                cyto_m_data[ll] = {
                    'rhos': [round(v, 4) for v in c_vals],
                    'median': round(float(np.median(c_vals)), 4),
                    'n': len(common),
                }
                sec_m_data[ll] = {
                    'rhos': [round(v, 4) for v in s_vals],
                    'median': round(float(np.median(s_vals)), 4),
                    'n': len(common),
                }

            # Wilcoxon signed-rank test (Matched, paired)
            wilcox_p = None
            if len(common) >= 6:
                _, wilcox_p = scipy_stats.wilcoxon(
                    np.array(c_vals), np.array(s_vals), alternative='two-sided')
            wilcox_pvals.append(float(wilcox_p) if wilcox_p is not None else None)

        # --- BH-FDR correction across levels within this atlas ---
        stats = {ll: {'total': None, 'matched': None} for ll in level_order}

        # Mann-Whitney
        valid_mw = [(i, p) for i, p in enumerate(mw_pvals) if p is not None]
        if valid_mw:
            mw_idx, mw_ps = zip(*valid_mw)
            _, mw_qs, _, _ = multipletests(list(mw_ps), method='fdr_bh')
            for j, idx in enumerate(mw_idx):
                stats[level_order[idx]]['total'] = {
                    'test': 'Mann-Whitney U',
                    'p_raw': round(mw_ps[j], 6),
                    'q_bh': round(float(mw_qs[j]), 6),
                }

        # Wilcoxon
        valid_wx = [(i, p) for i, p in enumerate(wilcox_pvals) if p is not None]
        if valid_wx:
            wx_idx, wx_ps = zip(*valid_wx)
            _, wx_qs, _, _ = multipletests(list(wx_ps), method='fdr_bh')
            for j, idx in enumerate(wx_idx):
                stats[level_order[idx]]['matched'] = {
                    'test': 'Wilcoxon signed-rank',
                    'p_raw': round(wx_ps[j], 6),
                    'q_bh': round(float(wx_qs[j]), 6),
                }

        result[title] = {
            'cytosig': cyto_data,
            'secact': sec_data,
            'cytosig_matched': cyto_m_data,
            'secact_matched': sec_m_data,
            'stats': stats,
        }
    return result


def prepare_bulk_validation_data(df):
    """Prepare validation data for all datasets — CytoSig and SecAct.
    For SecAct: show matched CytoSig targets + top additional targets to match count.
    """
    BULK_ATLAS_MAP = {
        'gtex': 'GTEx',
        'tcga': 'TCGA',
        'cima': 'CIMA',
        'inflammation_main': 'Inflammation Main',
        'scatlas_normal': 'scAtlas (Normal)',
        'scatlas_cancer': 'scAtlas (Cancer)',
    }
    result = {}
    for atlas, label in BULK_ATLAS_MAP.items():
        level, needs_mom = INDEPENDENT_LEVEL_MAP[atlas]
        result[label] = {}

        # CytoSig
        cyto_sub = df[(df['atlas'] == atlas) & (df['signature'] == 'cytosig') & (df['level'] == level)]
        if needs_mom:
            cyto_sub = cyto_sub[~cyto_sub['celltype'].isin(['all', 'unmappable'])]
            cyto_agg = cyto_sub.groupby('target')['spearman_rho'].median().dropna().reset_index()
            cyto_agg = cyto_agg.sort_values('spearman_rho', ascending=False)
        else:
            cyto_agg = cyto_sub[['target', 'spearman_rho']].dropna().sort_values('spearman_rho', ascending=False)

        if len(cyto_agg) > 0:
            cyto_targets = set(cyto_agg['target'].tolist())
            n_cyto = len(cyto_targets)
            result[label]['cytosig'] = {
                'targets': cyto_agg['target'].tolist(),
                'rhos': [round(r, 4) for r in cyto_agg['spearman_rho'].tolist()],
                'median': round(float(cyto_agg['spearman_rho'].median()), 3),
                'n': int(len(cyto_agg)),
            }
        else:
            cyto_targets = set()
            n_cyto = 43

        # SecAct — matched targets first, then top additional to match CytoSig count
        secact_sub = df[(df['atlas'] == atlas) & (df['signature'] == 'secact') & (df['level'] == level)]
        if needs_mom:
            secact_sub = secact_sub[~secact_sub['celltype'].isin(['all', 'unmappable'])]
            secact_agg = secact_sub.groupby('target')['spearman_rho'].median().dropna().reset_index()
        else:
            secact_agg = secact_sub[['target', 'spearman_rho']].dropna().drop_duplicates(subset=['target'], keep='first')

        if len(secact_agg) > 0:
            # Split into matched (shared with CytoSig) and extra
            matched = secact_agg[secact_agg['target'].isin(cyto_targets)].sort_values('spearman_rho', ascending=False)
            extra = secact_agg[~secact_agg['target'].isin(cyto_targets)].sort_values('spearman_rho', ascending=False)
            # Take matched + enough extra to reach n_cyto total
            n_extra = max(0, n_cyto - len(matched))
            selected = pd.concat([matched, extra.head(n_extra)], ignore_index=True)
            selected = selected.sort_values('spearman_rho', ascending=False)
            # Mark which are matched vs additional
            is_matched = [t in cyto_targets for t in selected['target']]
            result[label]['secact'] = {
                'targets': selected['target'].tolist(),
                'rhos': [round(r, 4) for r in selected['spearman_rho'].tolist()],
                'matched': is_matched,
                'median': round(float(selected['spearman_rho'].median()), 3),
                'n': int(len(selected)),
                'total_secact': int(len(secact_agg)),
            }
    return result


def prepare_scatter_data():
    """Prepare donor scatter data for key targets (CytoSig + SecAct)."""
    key_targets = ['IFNG', 'IL1B', 'TNFA', 'TGFB1', 'IL6', 'IL10', 'VEGFA', 'CD40L', 'TRAIL', 'HGF']
    # SecAct uses gene names that differ for some targets
    SECACT_NAME_MAP = {'CD40L': 'CD40LG', 'TRAIL': 'TNFSF10', 'TNFA': 'TNF'}
    atlas_files = {
        'GTEx': 'gtex_cytosig.json',
        'TCGA': 'tcga_cytosig.json',
        'CIMA': 'cima_cytosig.json',
        'Inflammation Main': 'inflammation_main_cytosig.json',
        'scAtlas (Normal)': 'scatlas_normal_cytosig.json',
        'scAtlas (Cancer)': 'scatlas_cancer_cytosig.json',
    }
    secact_files = {
        'GTEx': 'gtex_secact.json',
        'TCGA': 'tcga_secact.json',
        'CIMA': 'cima_secact.json',
        'Inflammation Main': 'inflammation_main_secact.json',
        'scAtlas (Normal)': 'scatlas_normal_secact.json',
        'scAtlas (Cancer)': 'scatlas_cancer_secact.json',
    }

    result = {}
    for atlas_label, filename in atlas_files.items():
        data = load_scatter(SCATTER_DONOR, filename)
        if data is None:
            continue
        result[atlas_label] = {}
        for target in key_targets:
            if target in data:
                d = data[target]
                pts = d.get('points', [])
                result[atlas_label][target] = {
                    'x': [round(p[0], 3) for p in pts],
                    'y': [round(p[1], 3) for p in pts],
                    'rho': round(d.get('rho', 0), 4),
                    'pval': d.get('pval', 1),
                    'n': d.get('n', len(pts)),
                }

        # Load SecAct scatter data
        secact_filename = secact_files.get(atlas_label)
        if secact_filename:
            secact_data = load_scatter(SCATTER_DONOR, secact_filename)
            if secact_data:
                for target in key_targets:
                    secact_name = SECACT_NAME_MAP.get(target, target)
                    if secact_name in secact_data:
                        d = secact_data[secact_name]
                        pts = d.get('points', [])
                        result[atlas_label][f'{target}_secact'] = {
                            'x': [round(p[0], 3) for p in pts],
                            'y': [round(p[1], 3) for p in pts],
                            'rho': round(d.get('rho', 0), 4),
                            'pval': d.get('pval', 1),
                            'n': d.get('n', len(pts)),
                        }
    return result


def prepare_celltype_comparison(df):
    """Prepare three-way celltype-level comparison: LinCytoSig vs CytoSig vs SecAct.

    Returns:
      - scatter: per-atlas LinCytoSig vs CytoSig scatter points (for interactive chart)
      - celltype_table: per-celltype win/loss summary (Δρ vs CytoSig and vs SecAct)
      - cytokine_table: per-cytokine win/loss summary
      - overall: aggregate win/loss counts and mean Δρ
    """
    mc = load_method_comparison()
    if not mc:
        return {}
    categories = mc.get('categories', [])
    matched = mc.get('matched_targets', {})
    cyto_rhos = mc.get('cytosig', {}).get('rhos', {})
    lin_rhos = mc.get('lincytosig', {}).get('rhos', {})
    sec_rhos = mc.get('secact', {}).get('rhos', {})

    cat_keys = [c['key'] for c in categories]

    # 1) Per-atlas scatter data (LinCytoSig vs CytoSig)
    scatter = {}
    for cat in categories:
        key = cat['key']
        label = cat['label'].replace('Inflammation Atlas', 'Inflammation Main')
        points = []
        for lin_target, mapping in matched.items():
            cyto_target = mapping.get('cytosig')
            if not cyto_target:
                continue
            cyto_val = cyto_rhos.get(cyto_target, {}).get(key)
            lin_val = lin_rhos.get(lin_target, {}).get(key)
            if cyto_val is not None and lin_val is not None:
                points.append({
                    'target': lin_target,
                    'cytosig': round(float(cyto_val), 4),
                    'lincytosig': round(float(lin_val), 4),
                })
        n_lin_win = sum(1 for p in points if p['lincytosig'] > p['cytosig'])
        n_cyto_win = sum(1 for p in points if p['cytosig'] > p['lincytosig'])
        scatter[label] = {
            'points': points,
            'n_lin_win': n_lin_win,
            'n_cyto_win': n_cyto_win,
            'n_tie': len(points) - n_lin_win - n_cyto_win,
        }

    # 2) Per-celltype table (aggregated across all atlases)
    from collections import defaultdict
    ct_vs_cyto = defaultdict(list)
    ct_vs_sec = defaultdict(list)
    for lin_target, mapping in matched.items():
        parts = lin_target.split('__')
        if len(parts) != 2:
            continue
        ct, cyto = parts
        cyto_target = mapping.get('cytosig')
        sec_target = mapping.get('secact')
        for key in cat_keys:
            lin_val = lin_rhos.get(lin_target, {}).get(key)
            if lin_val is not None and cyto_target and cyto_target in cyto_rhos:
                cyto_val = cyto_rhos[cyto_target].get(key)
                if cyto_val is not None:
                    ct_vs_cyto[ct].append(lin_val - cyto_val)
            if lin_val is not None and sec_target and sec_target in sec_rhos:
                sec_val = sec_rhos[sec_target].get(key)
                if sec_val is not None:
                    ct_vs_sec[ct].append(lin_val - sec_val)

    celltype_table = []
    for ct in sorted(ct_vs_cyto, key=lambda c: sum(1 for d in ct_vs_cyto[c] if d > 0) / max(len(ct_vs_cyto[c]), 1), reverse=True):
        diffs_cyto = ct_vs_cyto[ct]
        diffs_sec = ct_vs_sec.get(ct, [])
        sem_cyto = float(np.std(diffs_cyto) / np.sqrt(len(diffs_cyto))) if len(diffs_cyto) > 1 else 0.0
        sem_sec = float(np.std(diffs_sec) / np.sqrt(len(diffs_sec))) if len(diffs_sec) > 1 else 0.0
        celltype_table.append({
            'celltype': ct,
            'vs_cyto_mean': round(float(sum(diffs_cyto) / len(diffs_cyto)), 4) if diffs_cyto else None,
            'vs_cyto_sem': round(sem_cyto, 4),
            'vs_cyto_win': sum(1 for d in diffs_cyto if d > 0),
            'vs_cyto_loss': sum(1 for d in diffs_cyto if d < 0),
            'vs_cyto_n': len(diffs_cyto),
            'vs_sec_mean': round(float(sum(diffs_sec) / len(diffs_sec)), 4) if diffs_sec else None,
            'vs_sec_sem': round(sem_sec, 4),
            'vs_sec_win': sum(1 for d in diffs_sec if d > 0),
            'vs_sec_loss': sum(1 for d in diffs_sec if d < 0),
            'vs_sec_n': len(diffs_sec),
        })

    # 3) Per-cytokine table
    cyto_vs_cyto = defaultdict(list)
    cyto_vs_sec = defaultdict(list)
    for lin_target, mapping in matched.items():
        parts = lin_target.split('__')
        if len(parts) != 2:
            continue
        ct, cyto = parts
        cyto_target = mapping.get('cytosig')
        sec_target = mapping.get('secact')
        for key in cat_keys:
            lin_val = lin_rhos.get(lin_target, {}).get(key)
            if lin_val is not None and cyto_target and cyto_target in cyto_rhos:
                cyto_val = cyto_rhos[cyto_target].get(key)
                if cyto_val is not None:
                    cyto_vs_cyto[cyto].append(lin_val - cyto_val)
            if lin_val is not None and sec_target and sec_target in sec_rhos:
                sec_val = sec_rhos[sec_target].get(key)
                if sec_val is not None:
                    cyto_vs_sec[cyto].append(lin_val - sec_val)

    cytokine_table = []
    for cyto in sorted(cyto_vs_cyto, key=lambda c: sum(cyto_vs_cyto[c]) / max(len(cyto_vs_cyto[c]), 1), reverse=True):
        diffs_cyto = cyto_vs_cyto[cyto]
        diffs_sec = cyto_vs_sec.get(cyto, [])
        cytokine_table.append({
            'cytokine': cyto,
            'vs_cyto_mean': round(float(sum(diffs_cyto) / len(diffs_cyto)), 4) if diffs_cyto else None,
            'vs_cyto_win': sum(1 for d in diffs_cyto if d > 0),
            'vs_cyto_loss': sum(1 for d in diffs_cyto if d < 0),
            'vs_cyto_n': len(diffs_cyto),
            'vs_sec_mean': round(float(sum(diffs_sec) / len(diffs_sec)), 4) if diffs_sec else None,
            'vs_sec_win': sum(1 for d in diffs_sec if d > 0),
            'vs_sec_loss': sum(1 for d in diffs_sec if d < 0),
            'vs_sec_n': len(diffs_sec),
        })

    # 4) Per-atlas summary
    atlas_summary = []
    for cat in categories:
        key = cat['key']
        label = cat['label'].replace('Inflammation Atlas', 'Inflammation Main')
        diffs_cyto, diffs_sec = [], []
        for lin_target, mapping in matched.items():
            cyto_target = mapping.get('cytosig')
            sec_target = mapping.get('secact')
            lin_val = lin_rhos.get(lin_target, {}).get(key)
            if lin_val is not None and cyto_target and cyto_target in cyto_rhos:
                cyto_val = cyto_rhos[cyto_target].get(key)
                if cyto_val is not None:
                    diffs_cyto.append(lin_val - cyto_val)
            if lin_val is not None and sec_target and sec_target in sec_rhos:
                sec_val = sec_rhos[sec_target].get(key)
                if sec_val is not None:
                    diffs_sec.append(lin_val - sec_val)
        atlas_summary.append({
            'atlas': label,
            'vs_cyto_mean': round(float(sum(diffs_cyto) / len(diffs_cyto)), 4) if diffs_cyto else None,
            'vs_cyto_win': sum(1 for d in diffs_cyto if d > 0),
            'vs_cyto_loss': sum(1 for d in diffs_cyto if d < 0),
            'vs_cyto_n': len(diffs_cyto),
            'vs_sec_mean': round(float(sum(diffs_sec) / len(diffs_sec)), 4) if diffs_sec else None,
            'vs_sec_win': sum(1 for d in diffs_sec if d > 0),
            'vs_sec_loss': sum(1 for d in diffs_sec if d < 0),
            'vs_sec_n': len(diffs_sec),
        })

    # 5) Overall
    all_cyto = [d for diffs in ct_vs_cyto.values() for d in diffs]
    all_sec = [d for diffs in ct_vs_sec.values() for d in diffs]
    overall = {
        'vs_cyto_mean': round(float(sum(all_cyto) / len(all_cyto)), 4) if all_cyto else None,
        'vs_cyto_win': sum(1 for d in all_cyto if d > 0),
        'vs_cyto_loss': sum(1 for d in all_cyto if d < 0),
        'vs_cyto_n': len(all_cyto),
        'vs_sec_mean': round(float(sum(all_sec) / len(all_sec)), 4) if all_sec else None,
        'vs_sec_win': sum(1 for d in all_sec if d > 0),
        'vs_sec_loss': sum(1 for d in all_sec if d < 0),
        'vs_sec_n': len(all_sec),
    }

    return {
        'scatter': scatter,
        'celltypeTable': celltype_table,
        'cytokineTable': cytokine_table,
        'atlasSummary': atlas_summary,
        'overall': overall,
    }


def prepare_level_comparison_data():
    """Prepare data for section 5.2: Effect of Aggregation Level.

    Loads pre-computed three-way matched comparison at each celltype
    aggregation level. Downsamples rho arrays for boxplot rendering.
    """
    if not LEVEL_COMPARISON_PATH.exists():
        print(f'  WARNING: {LEVEL_COMPARISON_PATH} not found')
        return {}

    with open(LEVEL_COMPARISON_PATH) as f:
        raw = json.load(f)

    rng = np.random.default_rng(42)
    max_pts = 500

    result = {}
    for atlas_raw, levels in raw.items():
        atlas = atlas_raw.replace('Inflammation Atlas', 'Inflammation Main')
        atlas_data = []
        for lv in levels:
            entry = {
                'level': lv['level'],
                'n': lv['n_matched_3way'],
                'n_2way': lv['n_matched_2way'],
                'vs_cyto_win': lv['vs_cyto_win'],
                'vs_cyto_loss': lv['vs_cyto_loss'],
                'vs_cyto_mean': round(lv['vs_cyto_mean'], 4),
                'cyto_median': round(lv['cyto_median'], 4),
                'lin_median': round(lv['lin_median'], 4),
            }
            if lv.get('vs_sec_win') is not None:
                entry['vs_sec_win'] = lv['vs_sec_win']
                entry['vs_sec_loss'] = lv['vs_sec_loss']
                entry['vs_sec_mean'] = round(lv['vs_sec_mean'], 4)
                entry['sec_median'] = round(lv['sec_median'], 4)

            # Subsample rho arrays for boxplots (keep size manageable)
            for key, method in [('cyto_rhos', 'cytosig'), ('lin_rhos', 'lincytosig'), ('sec_rhos', 'secact')]:
                arr = lv.get(key, [])
                if len(arr) > max_pts:
                    idx = rng.choice(len(arr), max_pts, replace=False)
                    arr = [arr[int(i)] for i in sorted(idx)]
                entry[method] = [round(v, 4) for v in arr]

            atlas_data.append(entry)
        result[atlas] = atlas_data

    return result


def prepare_good_bad_data(df):
    """Prepare top/bottom correlated targets per atlas for CytoSig and SecAct.

    Data source:
        Merged correlation CSVs at donor level. 6 atlases: GTEx, TCGA, CIMA,
        Inflammation Main, scAtlas Normal, scAtlas Cancer. Two signature
        types: CytoSig and SecAct.
    Method:
        For each atlas x signature: deduplicate by target (keep first),
        sort by spearman_rho descending. Take top 15 and bottom 15.
        Return target name and rho (rounded to 4 decimal places).
    Output:
        Nested dict: {sig: {atlas_label: {top: [...], bottom: [...]}}}
        for interactive top/bottom bar charts (Section 4.2).
    """
    atlas_configs = [
        ('gtex', 'donor_only', 'GTEx'),
        ('tcga', 'donor_only', 'TCGA'),
        ('cima', 'donor_only', 'CIMA'),
        ('inflammation_main', 'donor_only', 'Inflammation Main'),
        ('scatlas_normal', 'donor_organ', 'scAtlas (Normal)'),
        ('scatlas_cancer', 'tumor_only', 'scAtlas (Cancer)'),
    ]
    result = {}
    for sig_type in ['cytosig', 'secact']:
        result[sig_type] = {}
        for atlas, level, label in atlas_configs:
            sub = df[(df['atlas'] == atlas) & (df['level'] == level) & (df['signature'] == sig_type)]
            # Drop duplicates by target (keep first = highest level match)
            sub = sub.drop_duplicates(subset=['target'], keep='first')
            sub = sub.sort_values('spearman_rho', ascending=False).reset_index(drop=True)
            if len(sub) == 0:
                continue
            top = sub.head(15)
            bottom = sub.tail(15).sort_values('spearman_rho', ascending=True).reset_index(drop=True)
            result[sig_type][label] = {
                'top': [{'target': r['target'], 'rho': round(float(r['spearman_rho']), 4)}
                        for _, r in top.iterrows()],
                'bottom': [{'target': r['target'], 'rho': round(float(r['spearman_rho']), 4)}
                           for _, r in bottom.iterrows()],
            }
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# HTML GENERATION
# ═══════════════════════════════════════════════════════════════════════════════

def generate_html(summary_table, boxplot_data, consistency_data, heatmap_data,
                  levels_data, bulk_data, scatter_data, good_bad_data,
                  method_boxplot_data, celltype_comparison_data,
                  level_comparison_data, stratified_data):

    # Serialize all data as JSON for embedding
    data_json = json.dumps({
        'summary': summary_table,
        'boxplot': boxplot_data,
        'stratified': stratified_data,
        'consistency': consistency_data,
        'heatmap': heatmap_data,
        'levels': levels_data,
        'bulk': bulk_data,
        'scatter': scatter_data,
        'goodBad': good_bad_data,
        'atlasLabels': ATLAS_LABELS,
        'sigColors': SIG_COLORS,
        'familyColors': FAMILY_COLORS,
        'methodBoxplot': method_boxplot_data,
        'celltypeComp': celltype_comparison_data,
        'levelComp': level_comparison_data,
    }, separators=(',', ':'))

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>CytoAtlas Report</title>
<script src="https://cdn.plot.ly/plotly-2.35.0.min.js" charset="utf-8"></script>
<style>
  :root {{
    --blue: #2563EB; --amber: #D97706; --emerald: #059669; --red: #DC2626;
    --gray-50: #F9FAFB; --gray-100: #F3F4F6; --gray-200: #E5E7EB;
    --gray-300: #D1D5DB; --gray-500: #6B7280; --gray-700: #374151; --gray-900: #111827;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{
    font-family: 'Segoe UI', system-ui, -apple-system, sans-serif;
    color: var(--gray-900); background: var(--gray-50);
    line-height: 1.7; font-size: 15px;
  }}
  .container {{ max-width: 1200px; margin: 0 auto; background: white; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
  .header {{
    background: linear-gradient(135deg, #1E3A5F 0%, #2563EB 50%, #059669 100%);
    color: white; padding: 50px 60px 40px;
  }}
  .header h1 {{ font-size: 32px; font-weight: 700; margin-bottom: 8px; }}
  .header .subtitle {{ font-size: 18px; opacity: 0.9; margin-bottom: 24px; }}
  .header .meta {{ font-size: 14px; opacity: 0.8; line-height: 1.8; }}
  .header .meta strong {{ opacity: 1; }}
  .content {{ padding: 40px 60px 60px; }}
  h2 {{
    font-size: 24px; font-weight: 700; color: var(--gray-900);
    border-bottom: 3px solid var(--blue); padding-bottom: 8px;
    margin: 48px 0 20px;
  }}
  h2:first-child {{ margin-top: 0; }}
  h3 {{ font-size: 18px; font-weight: 600; color: var(--gray-700); margin: 32px 0 12px; }}
  p {{ margin: 0 0 14px; }}
  ul, ol {{ margin: 0 0 14px 24px; }}
  li {{ margin-bottom: 6px; }}
  table {{ width: 100%; border-collapse: collapse; margin: 16px 0 20px; font-size: 13.5px; }}
  th {{
    background: var(--gray-700); color: white; padding: 10px 12px;
    text-align: left; font-weight: 600; font-size: 12.5px;
    text-transform: uppercase; letter-spacing: 0.3px;
  }}
  td {{ padding: 8px 12px; border-bottom: 1px solid var(--gray-200); vertical-align: top; }}
  tr:nth-child(even) td {{ background: var(--gray-50); }}
  tr:hover td {{ background: #EFF6FF; }}
  .figure {{
    margin: 28px 0; background: var(--gray-50);
    border: 1px solid var(--gray-200); border-radius: 8px; overflow: hidden;
  }}
  .figure img {{
    width: 100%; display: block; border-bottom: 1px solid var(--gray-200);
    cursor: zoom-in; transition: filter 0.2s;
  }}
  .figure img:hover {{ filter: brightness(0.93); }}
  .figure .caption {{ padding: 12px 16px; font-size: 13px; color: var(--gray-500); }}
  .figure .caption strong {{ color: var(--gray-700); }}
  .callout {{
    background: #EFF6FF; border-left: 4px solid var(--blue);
    padding: 16px 20px; margin: 20px 0; border-radius: 0 6px 6px 0;
  }}
  .callout.green {{ background: #ECFDF5; border-color: var(--emerald); }}
  .callout.amber {{ background: #FFFBEB; border-color: var(--amber); }}
  .callout.red {{ background: #FEF2F2; border-color: var(--red); }}
  .callout p:last-child {{ margin-bottom: 0; }}
  .stats-row {{ display: flex; gap: 16px; margin: 24px 0; flex-wrap: wrap; }}
  .stat-card {{
    flex: 1; min-width: 140px; background: white;
    border: 1px solid var(--gray-200); border-radius: 8px;
    padding: 16px; text-align: center;
  }}
  .stat-card .number {{ font-size: 28px; font-weight: 700; color: var(--blue); display: block; }}
  .stat-card .label {{
    font-size: 12px; color: var(--gray-500); text-transform: uppercase;
    letter-spacing: 0.5px; margin-top: 4px;
  }}
  hr {{ border: none; border-top: 2px solid var(--gray-200); margin: 40px 0; }}
  code {{
    background: var(--gray-100); padding: 2px 6px; border-radius: 4px;
    font-family: 'Consolas', 'Monaco', monospace; font-size: 13px;
  }}
  .badge {{
    display: inline-block; padding: 2px 8px; border-radius: 12px;
    font-size: 11px; font-weight: 600; text-transform: uppercase;
  }}
  .badge.blue {{ background: #DBEAFE; color: var(--blue); }}
  .badge.amber {{ background: #FEF3C7; color: var(--amber); }}
  .badge.green {{ background: #D1FAE5; color: var(--emerald); }}
  .win {{ color: var(--emerald); font-weight: 600; }}
  .lose {{ color: var(--red); font-weight: 600; }}
  .toc {{
    background: var(--gray-50); border: 1px solid var(--gray-200);
    border-radius: 8px; padding: 24px 28px; margin: 24px 0;
  }}
  .toc h3 {{ margin: 0 0 12px; font-size: 16px; }}
  .toc ol {{ margin: 0 0 0 20px; }}
  .toc li {{ margin-bottom: 4px; }}
  .toc a {{ color: var(--blue); text-decoration: none; }}
  .toc a:hover {{ text-decoration: underline; }}
  /* Plotly containers */
  .plotly-container {{
    margin: 20px 0; border: 1px solid var(--gray-200); border-radius: 8px;
    overflow: hidden; background: white;
  }}
  .plotly-container .caption {{
    padding: 10px 16px; font-size: 13px; color: var(--gray-500);
    border-top: 1px solid var(--gray-200); background: var(--gray-50);
  }}
  /* Tabs for heatmap/levels */
  .tab-bar {{
    display: flex; gap: 0; border-bottom: 2px solid var(--gray-200); margin-bottom: 0;
  }}
  .tab-btn {{
    padding: 10px 20px; border: none; background: var(--gray-50);
    cursor: pointer; font-size: 13px; font-weight: 600;
    color: var(--gray-500); border-bottom: 3px solid transparent;
    transition: all 0.2s;
  }}
  .tab-btn:hover {{ background: var(--gray-100); }}
  .tab-btn.active {{ color: var(--blue); border-bottom-color: var(--blue); background: white; }}
  .tab-btn.cytosig {{ color: var(--blue); }}
  .tab-btn.lincyto_best {{ color: #92400E; }}
  .tab-btn.lincytosig {{ color: var(--amber); }}
  .tab-btn.secact {{ color: var(--emerald); }}
  /* Dropdown controls */
  .controls {{ padding: 12px 16px; display: flex; gap: 16px; align-items: center; flex-wrap: wrap; }}
  .controls label {{ font-size: 13px; font-weight: 600; color: var(--gray-700); }}
  .controls select {{
    padding: 6px 12px; border: 1px solid var(--gray-300); border-radius: 6px;
    font-size: 13px; background: white; cursor: pointer;
  }}
  /* Lightbox */
  .lightbox-overlay {{
    position: fixed; top: 0; left: 0; width: 100%; height: 100%;
    background: rgba(0,0,0,0.92); z-index: 1000;
    display: flex; align-items: center; justify-content: center;
    opacity: 0; pointer-events: none; transition: opacity 0.3s ease;
  }}
  .lightbox-overlay.active {{ opacity: 1; pointer-events: auto; }}
  .lightbox-overlay img {{
    max-width: 92vw; max-height: 88vh; object-fit: contain;
    border-radius: 4px; box-shadow: 0 8px 40px rgba(0,0,0,0.6);
  }}
  .lightbox-close {{
    position: absolute; top: 16px; right: 20px; color: #fff; font-size: 32px;
    cursor: pointer; background: rgba(255,255,255,0.1); border: none; border-radius: 50%;
    width: 44px; height: 44px; display: flex; align-items: center; justify-content: center;
  }}
  .lightbox-close:hover {{ background: rgba(255,255,255,0.25); }}
  @media print {{
    body {{ background: white; }}
    .container {{ box-shadow: none; }}
    .header {{ -webkit-print-color-adjust: exact; print-color-adjust: exact; }}
    .figure {{ break-inside: avoid; }}
    h2 {{ break-after: avoid; }}
    .lightbox-overlay {{ display: none !important; }}
    .plotly-container {{ break-inside: avoid; }}
  }}
</style>
</head>
<body>

<div class="container">

<!-- HEADER -->
<div class="header">
  <h1>CytoAtlas</h1>
  <div class="subtitle">Pan-Disease Single-Cell Cytokine Activity Atlas</div>
  <div class="meta">
    <strong>Date:</strong> February 12, 2026
  </div>
</div>

<div class="content">

<!-- EXECUTIVE SUMMARY -->
<h2>Executive Summary</h2>

<p>CytoAtlas is a comprehensive computational resource that maps cytokine and secreted protein signaling activity across <strong>29 million human cells</strong> from four independent datasets spanning healthy donors, inflammatory diseases, cancers, and cytokine perturbations. The system uses <strong>linear ridge regression</strong> against experimentally derived signature matrices to infer activity &mdash; producing fully interpretable, conditional z-scores rather than black-box predictions.</p>

<div class="stats-row">
  <div class="stat-card"><span class="number">29M</span><span class="label">Total Cells</span></div>
  <div class="stat-card"><span class="number">4</span><span class="label">Datasets</span></div>
  <div class="stat-card"><span class="number">1,213</span><span class="label">Signatures</span></div>
  <div class="stat-card"><span class="number">6</span><span class="label">Validation Atlases</span></div>
  <div class="stat-card"><span class="number">262</span><span class="label">API Endpoints</span></div>
  <div class="stat-card"><span class="number">12</span><span class="label">Web Pages</span></div>
</div>

<div class="callout green">
<p><strong>Key results:</strong></p>
<ul>
  <li>1,213 signatures (43 CytoSig + 1,170 SecAct), plus 178 cell-type-specific LinCytoSig variants, validated across 4 independent atlases</li>
  <li>Spearman correlations reach &rho;=0.6&ndash;0.9 for well-characterized cytokines (IL1B, TNFA, VEGFA, TGFB family)</li>
  <li>Cross-atlas consistency demonstrates signatures generalize across CIMA, Inflammation Main, scAtlas, GTEx, and TCGA</li>
  <li>SecAct achieves the highest correlations in bulk &amp; organ-level analyses (median &rho;=0.40 in GTEx/TCGA)</li>
</ul>
</div>

<!-- TOC (inline only, no floating TOC — item 0) -->
<div class="toc">
  <h3>Table of Contents</h3>
  <ol>
    <li><a href="#sec1">System Architecture and Design Rationale</a></li>
    <li><a href="#sec2">Dataset Catalog</a></li>
    <li><a href="#sec3">Scientific Value Proposition</a></li>
    <li><a href="#sec4">Validation Results</a></li>
    <li><a href="#sec5">CytoSig vs LinCytoSig vs SecAct Comparison</a></li>
    <li><a href="#sec6">Key Takeaways for Scientific Discovery</a></li>
    <li><a href="#sec7">Appendix: Technical Specifications</a></li>
  </ol>
</div>

<hr>

<!-- SECTION 1 -->
<h2 id="sec1">1. System Architecture and Design Rationale</h2>

<h3>1.1 Why This Architecture?</h3>

<p>CytoAtlas was designed around three principles that distinguish it from typical bioinformatics databases:</p>

<p><strong>Principle 1: Linear interpretability over complex models.</strong><br>
Ridge regression (L2-regularized linear regression) was chosen deliberately over methods like autoencoders, graph neural networks, or foundation models. The resulting activity z-scores are <strong>conditional on the specific genes in the signature matrix</strong>, meaning every prediction can be traced to a weighted combination of known gene responses.</p>

<p><strong>Principle 2: Multi-level validation at every aggregation.</strong><br>
CytoAtlas validates at five levels: donor-level pseudobulk, donor &times; cell-type pseudobulk, single-cell, bulk RNA-seq (GTEx/TCGA), and bootstrap resampled with confidence intervals.</p>

<p><strong>Principle 3: Reproducibility through separation of concerns.</strong></p>

<table>
  <tr><th>Component</th><th>Technology</th><th>Purpose</th></tr>
  <tr><td><strong>Pipeline</strong></td><td>Python + CuPy (GPU)</td><td>Activity inference, 10&ndash;34x speedup</td></tr>
  <tr><td><strong>Storage</strong></td><td>DuckDB (3 databases, 68 tables)</td><td>Columnar analytics, no server needed</td></tr>
  <tr><td><strong>API</strong></td><td>FastAPI (262 endpoints)</td><td>RESTful data access, caching, auth</td></tr>
  <tr><td><strong>Frontend</strong></td><td>React 19 + TypeScript</td><td>Interactive exploration (12 pages)</td></tr>
</table>

<h3>1.2 Processing Scale</h3>

<table>
  <tr><th>Dataset</th><th>Cells/Samples</th><th>Processing Time</th><th>Hardware</th></tr>
  <tr><td>GTEx</td><td>19,788 bulk samples</td><td>~10min</td><td>A100 80GB</td></tr>
  <tr><td>TCGA</td><td>11,069 bulk samples</td><td>~10min</td><td>A100 80GB</td></tr>
  <tr><td>CIMA</td><td>6.5M cells</td><td>~2h</td><td>A100 80GB</td></tr>
  <tr><td>Inflammation Atlas</td><td>6.3M cells</td><td>~2h</td><td>A100 80GB</td></tr>
  <tr><td>scAtlas Normal</td><td>2.3M cells</td><td>~1h</td><td>A100 80GB</td></tr>
  <tr><td>scAtlas Cancer</td><td>4.1M cells</td><td>~1h</td><td>A100 80GB</td></tr>
  <tr><td>parse_10M</td><td>9.7M cells</td><td>~3h</td><td>A100 80GB</td></tr>
</table>
<p style="font-size:0.85em;color:#555;margin-top:4px;"><strong>Total:</strong> ~29M single cells + ~31K bulk RNA-seq samples, processed through ridge regression against 3 signature matrices (CytoSig, LinCytoSig, SecAct). <strong>Processing Time</strong> = wall-clock time for full activity inference on a single NVIDIA A100 GPU. See Section 2.1 for per-dataset details and cleaning considerations.</p>

<div class="figure">
  <img src="../figures/fig1_dataset_overview.png" alt="Figure 1: Dataset Overview">
  <div class="caption"><strong>Figure 1.</strong> CytoAtlas overview. (A) Cell counts across 4 datasets totaling 29M cells. (B) Three signature matrices. (C) Multi-level validation strategy.</div>
</div>

<hr>

<!-- SECTION 2: Dataset Catalog -->
<h2 id="sec2">2. Dataset Catalog</h2>

<h3>2.1 Datasets and Scale <a href="../reports/weekly/dataset_analytics.html" style="font-size:14px;font-weight:400;color:var(--blue);text-decoration:none;">[detailed analytics]</a></h3>

<!-- Item 1: Updated references -->
<table>
  <tr><th>#</th><th>Dataset</th><th>Type</th><th>Cells/Samples</th><th>Donors</th><th>Cell Types</th><th>Reference</th></tr>
  <tr><td>1</td><td><strong>GTEx</strong></td><td>Bulk RNA-seq</td><td>19,788 samples</td><td>946 donors</td><td>&mdash;</td><td>GTEx Consortium, v11</td></tr>
  <tr><td>2</td><td><strong>TCGA</strong></td><td>Bulk RNA-seq</td><td>11,069 samples</td><td>10,274 donors</td><td>&mdash;</td><td>TCGA PanCancer</td></tr>
  <tr><td>3</td><td><strong>CIMA</strong></td><td>scRNA-seq</td><td>6,484,974</td><td>421 donors</td><td>27 L2 / 100+ L3</td><td>J. Yin et al., <em>Science</em>, 2026</td></tr>
  <tr><td>4</td><td><strong>Inflammation Main</strong></td><td>scRNA-seq</td><td>4,918,140</td><td>817 samples</td><td>66+</td><td>Jimenez-Gracia et al., <em>Nature Medicine</em>, 2026</td></tr>
  <tr><td>5</td><td><strong>Inflammation Val</strong></td><td>scRNA-seq</td><td>849,922</td><td>144 samples</td><td>66+</td><td>Validation cohort</td></tr>
  <tr><td>6</td><td><strong>Inflammation Ext</strong></td><td>scRNA-seq</td><td>572,872</td><td>86 samples</td><td>66+</td><td>External cohort</td></tr>
  <tr><td>7</td><td><strong>scAtlas Normal</strong></td><td>scRNA-seq</td><td>2,293,951</td><td>317 donors</td><td>102 subCluster</td><td>Q. Shi et al., <em>Nature</em>, 2025</td></tr>
  <tr><td>8</td><td><strong>scAtlas Cancer</strong></td><td>scRNA-seq</td><td>4,146,975</td><td>717 donors (601 tumor-only)</td><td>162 cellType1</td><td>Q. Shi et al., <em>Nature</em>, 2025</td></tr>
  <tr><td>9</td><td><strong>parse_10M</strong></td><td>scRNA-seq</td><td>9,697,974</td><td>12 donors &times; 91 cytokines</td><td>18 PBMC types</td><td>Oesinghaus et al., <em>bioRxiv</em>, 2026</td></tr>
</table>
<p><strong>Grand total:</strong> ~29 million single cells + ~31K bulk samples across 9 datasets, 100+ cell types.</p>

<h3>2.2 Disease and Condition Categories</h3>

<p><strong>CIMA (421 healthy donors):</strong> Healthy population atlas with paired blood biochemistry (19 markers: ALT, AST, glucose, lipid panel, etc.) and plasma metabolomics (1,549 features). Enables age, BMI, sex, and smoking correlations with cytokine activity.</p>
<p><strong>Inflammation Atlas (20 diseases):</strong> RA, SLE, Sjogren's, PSA, Crohn's, UC, COVID-19, Sepsis, HIV, HBV, BRCA, CRC, HNSCC, NPC, COPD, Cirrhosis, MS, Asthma, Atopic Dermatitis</p>
<p><strong>scAtlas Normal (317 donors):</strong> 35 organs, 12 tissues with &ge;20 donors for per-organ stratification (Breast 124, Lung 97, Colon 65, Heart 52, Liver 43, etc.)</p>
<p><strong>scAtlas Cancer (717 donors, 601 tumor-only):</strong> 29 cancer types, 11 with &ge;20 tumor-only donors for per-cancer stratification (HCC 88, PAAD 58, CRC 51, ESCA 48, HNSC 39, LUAD 36, NPC 36, KIRC 31, BRCA 30, ICC 29, STAD 27)</p>
<!-- Item 2: parse_10M is NOT ground truth -->
<p><strong>parse_10M:</strong> 90 cytokines &times; 12 donors &mdash; independent <em>in vitro</em> perturbation dataset for comparison. A considerable portion of cytokines (~58%) are produced in <em>E. coli</em>, with the remainder from insect (Sf21, 12%) and mammalian (CHO, NS0, HEK293, ~30%) expression systems. Because exogenous perturbagens may induce effects differing from endogenously produced cytokines, parse_10M serves as an independent comparison rather than strict ground truth. CytoSig/SecAct has a potential advantage in this regard, as it infers relationships directly from physiologically relevant samples.</p>

<h3>2.3 Signature Matrices</h3>

<!-- Item 3: Updated SecAct reference -->
<table>
  <tr><th>Matrix</th><th>Targets</th><th>Construction</th><th>Reference</th></tr>
  <tr><td><span class="badge blue">CytoSig</span></td><td>43 cytokines</td><td>Median log2FC across all experimental bulk RNA-seq</td><td>Jiang et al., <em>Nature Methods</em>, 2021</td></tr>
  <tr><td><span class="badge amber">LinCytoSig</span></td><td>178 (45 cell types &times; 1&ndash;13 cytokines)</td><td>Cell-type-stratified median from CytoSig database (<a href="methodology.html">methodology</a>)</td><td>This work</td></tr>
  <tr><td><span class="badge green">SecAct</span></td><td>1,170 secreted proteins</td><td>Median global Moran's I across 1,000 Visium datasets</td><td>Ru et al., <em>Nature Methods</em>, 2026 (in press)</td></tr>
</table>

<hr>

<!-- SECTION 3: Scientific Value -->
<h2 id="sec3">3. Scientific Value Proposition</h2>

<h3>3.1 What Makes CytoAtlas Different from Deep Learning Approaches?</h3>

<p>Most single-cell analysis tools use complex models (VAEs, GNNs, transformers) that produce <strong>aggregated, non-linear representations</strong> difficult to interpret biologically. CytoAtlas takes the opposite approach:</p>

<table>
  <tr><th>Property</th><th>CytoAtlas (Ridge Regression)</th><th>Typical DL Approach</th></tr>
  <tr><td><strong>Model</strong></td><td>Linear (z = X&beta; + &epsilon;)</td><td>Non-linear (multi-layer NN)</td></tr>
  <tr><td><strong>Interpretability</strong></td><td>Every gene's contribution is a coefficient</td><td>Feature importance approximated post-hoc</td></tr>
  <tr><td><strong>Conditionality</strong></td><td>Activity conditional on specific gene set</td><td>Latent space mixes all features</td></tr>
  <tr><td><strong>Confidence</strong></td><td>Permutation-based z-scores with CI</td><td>Often point estimates only</td></tr>
  <tr><td><strong>Generalization</strong></td><td>Tested across 6 independent cohorts</td><td>Often held-out splits of same cohort</td></tr>
  <tr><td><strong>Bias</strong></td><td>Transparent &mdash; limited by signature matrix genes</td><td>Hidden in architecture and training data</td></tr>
</table>

<!-- Item 4: DNN-based → DL-based -->
<div class="callout">
<p><strong>The key insight:</strong> CytoAtlas is not trying to replace DL-based tools. It provides an <strong>orthogonal, complementary signal</strong> that a human scientist can directly inspect. When CytoAtlas says "IFNG activity is elevated in CD8+ T cells from RA patients," you can verify this by checking the IFNG signature genes in those cells.</p>
</div>

<h3>3.2 What Scientific Questions Does CytoAtlas Answer?</h3>
<ol>
  <li><strong>Which cytokines are active in which cell types across diseases?</strong> &mdash; IL1B/TNFA in monocytes/macrophages, IFNG in CD8+ T and NK cells, IL17A in Th17, VEGFA in endothelial/tumor cells, TGFB family in stromal cells &mdash; quantified across 20 diseases, 35 organs, and 15 cancer types.</li>
  <li><strong>Are cytokine activities consistent across independent cohorts?</strong> &mdash; Yes. IL1B, TNFA, VEGFA, and TGFB family show consistent positive correlations across all 6 validation atlases (Figure 7).</li>
  <li><strong>Does cell-type-specific biology matter for cytokine inference?</strong> &mdash; For select immune types, yes: LinCytoSig improves prediction for Basophils (+0.21 &Delta;&rho;), NK cells (+0.19), and DCs (+0.18), but global CytoSig wins overall (Figures 10&ndash;11).</li>
  <li><strong>Which secreted proteins beyond cytokines show validated activity?</strong> &mdash; SecAct (1,170 targets) achieves the highest correlations across all atlases (median &rho;=0.33&ndash;0.49), with novel validated targets like Activin A (&rho;=0.98), CXCL12 (&rho;=0.92), and BMP family (Figure 12).</li>
  <li><strong>Can we predict treatment response from cytokine activity?</strong> &mdash; We are incorporating cytokine-blocking therapy outcomes from bulk RNA-seq to test whether predicted cytokine activity associates with therapy response. Additionally, Inflammation Atlas responder/non-responder labels enable treatment response prediction using cytokine activity profiles as features.</li>
</ol>

<h3>3.3 Validation Philosophy</h3>
<p>CytoAtlas validates against a simple but powerful principle: <strong>if CytoSig predicts high IFNG activity for a sample, that sample should have high IFNG gene expression.</strong> This expression-activity correlation is computed via Spearman rank correlation across donors/samples.</p>
<p>This is a <em>conservative</em> validation &mdash; it only captures signatures where the target gene itself is expressed. Signatures that act through downstream effectors would not be captured, meaning our validation <strong>underestimates</strong> true accuracy.</p>

<hr>

<!-- SECTION 4: Validation Results -->
<h2 id="sec4">4. Validation Results</h2>

<!-- Item 5: Embedded summary table instead of figure -->
<h3>4.1 Overall Performance Summary</h3>

<div id="summary-table-container"></div>

<div class="callout">
<p><strong>PRIMARY independent level:</strong> The summary table above reports results at each dataset&rsquo;s PRIMARY independent level &mdash; the aggregation level where samples are fully independent (each donor counted once). This ensures correlation statistics are not inflated by donor duplication. See the &ldquo;Primary Level&rdquo; column for each dataset&rsquo;s level.</p>
<p><strong>How &ldquo;N Targets&rdquo; is determined:</strong> A target is included in the validation for a given atlas only if (1) the target&rsquo;s signature genes overlap sufficiently with the atlas gene expression matrix, and (2) the target gene itself is expressed in enough samples to compute a meaningful Spearman correlation. Targets whose gene is absent or not detected in a dataset are excluded. CytoSig defines 43 cytokines and SecAct defines 1,170 secreted proteins. Inflammation Main retains only 33 of 43 CytoSig targets and 805 of 1,170 SecAct targets because 10 cytokine genes (BDNF, BMP4, CXCL12, GCSF, IFN1, IL13, IL17A, IL36, IL4, WNT3A) are not sufficiently expressed in these blood/PBMC samples.</p>
<p><strong>Stratified levels</strong> (GTEx by_tissue, TCGA primary_by_cancer): Correlations are computed within each tissue/cancer type (ensuring independence), then summarized across groups. N Targets counts unique targets at the &ldquo;all&rdquo; aggregate level. Finer per-tissue or per-cancer breakdowns are available in Section 4.2 below.</p>
</div>

<!-- Per-tissue/per-cancer stratified validation -->
<h3>4.2 Per-Tissue and Per-Cancer Stratified Validation
    <a href="stats_section_4.3.html" style="font-size:14px;font-weight:400;color:var(--blue);text-decoration:none;">[Statistical Methods]</a></h3>
<div class="plotly-container">
  <div class="tab-bar" id="stratified-tabs">
    <button class="tab-btn active" onclick="renderStratified('total')">Total</button>
    <button class="tab-btn" onclick="renderStratified('matched')">Matched</button>
  </div>
  <div class="controls">
    <label>Dataset:</label>
    <select id="stratified-dataset-select" onchange="renderStratified()">
      <option value="GTEx">GTEx (by tissue)</option>
      <option value="TCGA">TCGA (by cancer type)</option>
    </select>
  </div>
  <div id="stratified-chart" style="height:700px;"></div>
  <div id="stratified-caption" class="caption">
    <strong>Figure 2.</strong> Per-tissue/cancer CytoSig vs SecAct median Spearman &rho; comparison. BH-FDR corrected significance: *** q&lt;0.001, ** q&lt;0.01, * q&lt;0.05.
  </div>
</div>

<div class="callout">
<p><strong>Stratified validation:</strong> Instead of aggregating tissues/cancers into a single median-of-medians, this view shows the CytoSig vs SecAct comparison <em>within each individual tissue</em> (GTEx) or <em>cancer type</em> (TCGA). Mann-Whitney U test (Total tab: all targets) and Wilcoxon signed-rank test (Matched tab: 32 shared targets) with BH-FDR correction across all strata within each dataset.</p>
</div>

<!-- Item 6: Interactive boxplot with Total / Matched tabs -->
<h3>4.3 Cross-Dataset Comparison: CytoSig vs SecAct</h3>
<div class="plotly-container">
  <div class="tab-bar" id="boxplot-tabs">
    <button class="tab-btn active" onclick="renderBoxplot('total')">Total</button>
    <button class="tab-btn" onclick="renderBoxplot('matched')">Matched</button>
  </div>
  <div id="boxplot-chart" style="height:500px;"></div>
  <div id="boxplot-caption" class="caption"><strong>Figure 3.</strong> Spearman &rho; distributions across atlases for CytoSig (43 targets) vs SecAct (1,170 targets). Independence-corrected: GTEx/TCGA use median-of-medians. Mann-Whitney U test p-values shown above each atlas.</div>
</div>

<div class="callout">
<p><strong>Why does SecAct appear to underperform CytoSig in the Inflammation Main atlas?</strong></p>
<p>This is a <strong>composition effect</strong>, not a genuine performance gap, confirmed by two complementary statistical tests:</p>
<p><strong>Total comparison (Mann&ndash;Whitney U test):</strong> Compares the full &rho; distributions of CytoSig (43 cytokine signatures) vs SecAct (~1,170 secreted protein signatures) using <strong>independence-corrected</strong> values. For GTEx/TCGA, each target&rsquo;s representative &rho; is the median across per-tissue/cancer values (median-of-medians); for other atlases, donor_only/tumor_only &rho; is used directly. SecAct achieves a significantly higher median &rho; in 5 of 6 atlases (GTEx: <em>p</em> = 4.75 &times; 10<sup>&minus;4</sup>; TCGA: <em>p</em> = 2.80 &times; 10<sup>&minus;3</sup>; CIMA: <em>p</em> = 3.18 &times; 10<sup>&minus;2</sup>; scAtlas Normal: <em>p</em> = 1.04 &times; 10<sup>&minus;4</sup>; scAtlas Cancer: <em>p</em> = 1.06 &times; 10<sup>&minus;5</sup>). Inflammation Main is the sole exception (U = 14,101, <em>p</em> = 0.548, not significant) and the only atlas where CytoSig&rsquo;s median &rho; (0.323) exceeds SecAct&rsquo;s (0.173).</p>
<p><strong>Matched comparison (Wilcoxon signed-rank test):</strong> Restricts to the 32 targets shared between both methods (22 direct + 10 alias-resolved), each target serving as its own control. SecAct&rsquo;s median &rho; is consistently higher across all 6 atlases, reaching significance in 5 (GTEx: <em>p</em> = 3.54 &times; 10<sup>&minus;5</sup>; TCGA: <em>p</em> = 3.24 &times; 10<sup>&minus;6</sup>; CIMA: <em>p</em> = 2.28 &times; 10<sup>&minus;2</sup>; scAtlas Normal: <em>p</em> = 3.54 &times; 10<sup>&minus;5</sup>; scAtlas Cancer: <em>p</em> = 3.54 &times; 10<sup>&minus;5</sup>). Inflammation Main is not significant (<em>p</em> = 0.141).</p>
<p>Inflammation Main is largely blood-derived, so many SecAct targets that perform well in multi-organ contexts contribute near-zero or negative correlations here. In fact, 99 SecAct targets are negative <em>only</em> in Inflammation Main but positive in all other atlases, reflecting tissue-specific expression limitations rather than inference failure. The &ldquo;Matched&rdquo; tab above demonstrates the fair comparison on equal footing.</p>
</div>

<h3>4.4 Best and Worst Correlated Targets</h3>
<div class="plotly-container">
  <div class="controls">
    <label>Signature:</label>
    <select id="goodbad-sig-select" onchange="updateGoodBad()">
      <option value="cytosig">CytoSig</option>
      <option value="secact">SecAct</option>
    </select>
    <label>Atlas:</label>
    <select id="goodbad-atlas-select" onchange="updateGoodBad()">
    </select>
  </div>
  <div id="goodbad-chart" style="height:700px;"></div>
  <div class="caption"><strong>Figure 4.</strong> Top 15 (best) and bottom 15 (worst) correlated targets. Select signature type and atlas from dropdowns.</div>
</div>

<p><strong>Consistently well-correlated targets (&rho; &gt; 0.3 across multiple atlases):</strong></p>
<ul>
  <li><strong>IL1B</strong> (&rho; = 0.67 CIMA, 0.68 Inflammation Main) &mdash; canonical inflammatory cytokine</li>
  <li><strong>TNFA</strong> (&rho; = 0.63 CIMA, 0.60 Inflammation Main) &mdash; master inflammatory regulator</li>
  <li><strong>VEGFA</strong> (&rho; = 0.79 Inflammation Main, 0.92 scAtlas) &mdash; angiogenesis factor</li>
  <li><strong>TGFB1/2/3</strong> (&rho; = 0.35&ndash;0.55 across atlases)</li>
  <li><strong>BMP2/4</strong> (&rho; = 0.26&ndash;0.92 depending on atlas)</li>
</ul>

<p><strong>Consistently poorly correlated targets (&rho; &lt; 0 in multiple atlases):</strong></p>
<ul>
  <li><strong>CD40L</strong> (&rho; = &minus;0.48 CIMA, &minus;0.56 Inflammation Main) &mdash; membrane-bound, not secreted</li>
  <li><strong>TRAIL</strong> (&rho; = &minus;0.46 CIMA, &minus;0.55 Inflammation Main) &mdash; apoptosis inducer</li>
  <li><strong>LTA</strong> (&rho; = &minus;0.33 CIMA), <strong>HGF</strong> (&rho; = &minus;0.25 CIMA)</li>
</ul>

<!-- Gene mapping verification and mechanistic analysis -->
<div class="callout amber">
<p><strong>Gene mapping verified:</strong> All four targets are correctly mapped (CD40L&rarr;CD40LG, TRAIL&rarr;TNFSF10, LTA&rarr;LTA, HGF&rarr;HGF). No gene ID confusion exists. The poor correlations reflect specific molecular mechanisms:</p>
</div>

<table>
  <tr><th>Target</th><th>Gene</th><th>Dominant Mechanism</th><th>Contributing Factors</th></tr>
  <tr>
    <td><strong>CD40L</strong></td><td>CD40LG</td>
    <td>Platelet-derived sCD40L invisible to scRNA-seq (~95% of circulating CD40L); ADAM10-mediated membrane shedding</td>
    <td>Unstable mRNA (3&prime;-UTR destabilizing element); transient expression kinetics (peak 6&ndash;8h post-activation); paracrine disconnect (T cell &rarr; B cell/DC)</td>
  </tr>
  <tr>
    <td><strong>TRAIL</strong></td><td>TNFSF10</td>
    <td>Three decoy receptors (DcR1/TNFRSF10C, DcR2/TNFRSF10D, OPG/TNFRSF11B) competitively sequester ligand without signaling</td>
    <td>Non-functional splice variants (TRAIL-beta, TRAIL-gamma lack exon 3) inflate mRNA counts; cathepsin E-mediated shedding; apoptosis-induced survival bias in scRNA-seq data</td>
  </tr>
  <tr>
    <td><strong>LTA</strong></td><td>LTA</td>
    <td>Obligate heteromeric complex with LTB: the dominant form (LT&alpha;1&beta;2) requires <em>LTB</em> co-expression and signals through LTBR, not TNFR1/2</td>
    <td>Mathematical collinearity with TNFA in ridge regression (LTA3 homotrimer binds the same TNFR1/2 receptors as TNF-&alpha;); 7 known splice variants; low/transient expression</td>
  </tr>
  <tr>
    <td><strong>HGF</strong></td><td>HGF</td>
    <td>Obligate mesenchymal-to-epithelial paracrine topology: HGF produced by fibroblasts/stellate cells, MET receptor on epithelial cells</td>
    <td>Secreted as inactive pro-HGF requiring proteolytic cleavage by HGFAC/uPA (post-translational activation is rate-limiting); ECM/heparin sequestration creates stored protein pool invisible to transcriptomics</td>
  </tr>
</table>

<div class="callout">
<p><strong>Key insight:</strong> None of these targets have isoforms or subunits mapping to different gene IDs that would cause gene ID confusion. The poor correlations are supposedly driven by <strong>post-translational regulation</strong> (membrane shedding, proteolytic activation, decoy receptor sequestration), <strong>paracrine signaling topology</strong> (producer and responder cells are different cell types), and <strong>heteromeric complex dependence</strong> (LTA requires LTB). These represent fundamental limitations of correlating ligand mRNA abundance and predicted activity as validation strategy of cytokine activity prediction model.</p>
</div>

<div class="callout amber">
<p><strong>However, SecAct rescues all four targets.</strong> The poor correlations above are <strong>CytoSig-specific</strong>, not universal. SecAct achieves strong positive correlations for every one of these targets (mean &rho; across atlases):</p>
<table style="margin:0.5em 0;">
  <tr><th>Target</th><th>CytoSig Gene</th><th>CytoSig Mean &rho;</th><th>SecAct Gene</th><th>SecAct Mean &rho;</th></tr>
  <tr><td>CD40L</td><td>CD40LG</td><td>&minus;0.006</td><td>CD40LG</td><td><strong>+0.420</strong></td></tr>
  <tr><td>TRAIL</td><td>TNFSF10</td><td>&minus;0.016</td><td>TNFSF10</td><td><strong>+0.418</strong></td></tr>
  <tr><td>LTA</td><td>LTA</td><td>&minus;0.019</td><td>LTA</td><td><strong>+0.474</strong></td></tr>
  <tr><td>HGF</td><td>HGF</td><td>+0.034</td><td>HGF</td><td><strong>+0.540</strong></td></tr>
</table>
<p>The key difference is <strong>how the signature matrices are constructed</strong>. CytoSig derives signatures from <strong>log2 fold-change in cytokine stimulation experiments</strong> (in vitro), which fails when the relationship between ligand mRNA and downstream activity is confounded by post-translational regulation, decoy receptors, or paracrine topology. SecAct derives signatures from <strong>spatial co-expression correlations</strong> (Moran&rsquo;s I across 1,000+ Visium spatial transcriptomics datasets), which captures the actual tissue-level gene&ndash;protein relationships regardless of whether the signaling mechanism involves membrane shedding, proteolytic activation, or cross-cell-type paracrine signaling. Select &ldquo;SecAct&rdquo; in the dropdown above to verify these correlations interactively.</p>
</div>

<!-- Item 10: Interactive consistency plot -->
<h3>4.5 Cross-Atlas Consistency</h3>
<div class="plotly-container">
  <div id="consistency-chart" style="height:550px;"></div>
  <div class="caption"><strong>Figure 5.</strong> Key cytokine target correlations tracked across 6 independent atlases (donor-level). Solid lines = CytoSig; dashed lines = SecAct. Lines colored by cytokine family. Click legend entries to show/hide targets.</div>
</div>

<!-- Item 11: Aggregation levels with Total / Matched tabs -->
<h3>4.6 Effect of Aggregation Level
    <a href="stats_section_4.6.html" style="font-size:14px;font-weight:400;color:var(--blue);text-decoration:none;">[Statistical Methods]</a></h3>
<div class="plotly-container">
  <div class="tab-bar" id="levels-tabs">
    <button class="tab-btn active" onclick="renderLevels('total')">Total</button>
    <button class="tab-btn" onclick="renderLevels('matched')">Matched</button>
  </div>
  <div class="controls">
    <label>Atlas:</label>
    <select id="levels-atlas-select" onchange="renderLevels()">
    </select>
  </div>
  <div id="levels-chart" style="height:500px;"></div>
  <div id="levels-caption" class="caption"><strong>Figure 6.</strong> Effect of cell-type annotation granularity on validation correlations. <strong>Total:</strong> CytoSig (43 targets) vs SecAct (1,170 targets). <strong>Matched:</strong> 32 shared targets only. Select atlas from dropdown.</div>
</div>

<div class="callout">
<p><strong>Aggregation levels explained:</strong> Pseudobulk profiles are aggregated at increasingly fine cell-type resolution. At coarser levels, each pseudobulk profile averages more cells, yielding smoother expression estimates but masking cell-type-specific signals. At finer levels, each profile is more cell-type-specific but based on fewer cells.</p>
</div>

<table>
  <tr><th>Atlas</th><th>Level</th><th>Description</th><th>N Cell Types</th></tr>
  <tr><td rowspan="5"><strong>CIMA</strong></td>
    <td>Donor Only</td><td>Whole-sample pseudobulk per donor</td><td>1 (all)</td></tr>
  <tr><td>Donor &times; L1</td><td>Broad lineages (B, CD4_T, CD8_T, Myeloid, NK, etc.)</td><td>7</td></tr>
  <tr><td>Donor &times; L2</td><td>Intermediate (CD4_memory, CD8_naive, DC, Mono, etc.)</td><td>28</td></tr>
  <tr><td>Donor &times; L3</td><td>Fine-grained (CD4_Tcm, cMono, Switched_Bm, etc.)</td><td>39</td></tr>
  <tr><td>Donor &times; L4</td><td>Finest marker-annotated (CD4_Th17-like_RORC, cMono_IL1B, etc.)</td><td>73</td></tr>
  <tr><td rowspan="3"><strong>Inflammation Main</strong></td>
    <td>Donor Only</td><td>Whole-sample pseudobulk per donor</td><td>1 (all)</td></tr>
  <tr><td>Donor &times; L1</td><td>Broad categories (B, DC, Mono, T_CD4/CD8 subsets, etc.)</td><td>18</td></tr>
  <tr><td>Donor &times; L2</td><td>Fine-grained (Th1, Th2, Tregs, NK_adaptive, etc.)</td><td>65</td></tr>
  <tr><td rowspan="3"><strong>scAtlas Normal</strong></td>
    <td>Donor &times; Organ</td><td>Per-organ pseudobulk (Bladder, Blood, Breast, Lung, etc.)</td><td>25 organs</td></tr>
  <tr><td>Donor &times; Organ &times; CT1</td><td>Broad cell types within each organ</td><td>191</td></tr>
  <tr><td>Donor &times; Organ &times; CT2</td><td>Fine cell types within each organ</td><td>356</td></tr>
  <tr><td rowspan="3"><strong>scAtlas Cancer</strong></td>
    <td>Tumor Only</td><td>Whole-sample pseudobulk per tumor donor</td><td>1 (all)</td></tr>
  <tr><td>Tumor &times; Cancer</td><td>Per-cancer type pseudobulk (HCC, PAAD, CRC, etc.)</td><td>29 types</td></tr>
  <tr><td>Tumor &times; Cancer &times; CT1</td><td>Broad cell types within each cancer type</td><td>~120</td></tr>
</table>

<h3>Representative Scatter Plots</h3>
<div class="plotly-container">
  <div class="controls">
    <label>Target:</label>
    <select id="scatter-target-select" onchange="updateScatter()">
    </select>
    <label>Atlas:</label>
    <select id="scatter-atlas-select" onchange="updateScatter()">
    </select>
    <label>Signature:</label>
    <select id="scatter-sig-select" onchange="updateScatter()">
      <option value="cytosig">CytoSig</option>
      <option value="secact">SecAct</option>
    </select>
  </div>
  <div id="scatter-chart" style="height:500px;"></div>
  <div class="caption"><strong>Figure 7.</strong> Donor-level expression vs predicted activity. Select target, atlas, and signature method from dropdowns.</div>
</div>

<!-- Item 12: Interactive heatmap with tabs -->
<h3>Biologically Important Targets Heatmap</h3>
<div class="plotly-container">
  <div class="tab-bar" id="heatmap-tabs">
    <button class="tab-btn cytosig active" onclick="switchHeatmapTab('cytosig')">CytoSig</button>
    <button class="tab-btn secact" onclick="switchHeatmapTab('secact')">SecAct</button>
  </div>
  <div id="heatmap-chart" style="min-height:500px;"></div>
  <div class="caption"><strong>Figure 8.</strong> Spearman &rho; heatmap for biologically important targets across all atlases. Switch between signature types. Hover over cells for details.</div>
</div>

<div class="callout">
<p><strong>How each correlation value is computed:</strong> For each (target, atlas) cell, we compute Spearman rank correlation between <em>predicted cytokine activity</em> (ridge regression z-score) and <em>target gene expression</em> across all donor-level pseudobulk samples. Specifically:</p>
<ol>
  <li><strong>Pseudobulk aggregation:</strong> For each atlas, gene expression is aggregated to the donor level (one profile per donor or donor &times; cell type).</li>
  <li><strong>Activity inference:</strong> Ridge regression (<code>secactpy.ridge</code>, &lambda;=5&times;10<sup>5</sup>) is applied using the signature matrix (CytoSig: 4,881 genes &times; 43 cytokines; SecAct: 7,919 genes &times; 1,170 targets) to predict activity z-scores for each pseudobulk sample.</li>
  <li><strong>Correlation:</strong> Spearman &rho; is computed between the predicted activity z-score and the original expression of the target gene across all donor-level samples within that atlas. A positive &rho; means higher predicted activity tracks with higher target gene expression.</li>
</ol>
<p>GTEx uses per-tissue pseudobulk (median-of-medians across 29 tissues); TCGA uses per-cancer type (median-of-medians across 33 cancers); CIMA/Inflammation Main use donor-only; scAtlas Normal uses donor-only; scAtlas Cancer uses tumor-only.</p>
</div>

<!-- Item 13: Interactive comprehensive validation -->
<h3>Comprehensive Validation Across All Datasets</h3>
<div class="plotly-container">
  <div class="controls">
    <label>Dataset:</label>
    <select id="bulk-dataset-select" onchange="updateBulk()">
      <option value="GTEx">GTEx</option>
      <option value="TCGA">TCGA</option>
      <option value="CIMA">CIMA</option>
      <option value="Inflammation Main">Inflammation Main</option>
      <option value="scAtlas (Normal)">scAtlas (Normal)</option>
      <option value="scAtlas (Cancer)">scAtlas (Cancer)</option>
    </select>
    <label>Signature:</label>
    <select id="bulk-sig-select" onchange="updateBulk()">
      <option value="cytosig">CytoSig</option>
      <option value="secact">SecAct (matched + top additional)</option>
    </select>
  </div>
  <div id="bulk-chart" style="height:500px;"></div>
  <div class="caption"><strong>Figure 9.</strong> Validation: targets ranked by Spearman &rho; across all datasets and signature types. Select dataset and signature from dropdowns.</div>
</div>

<hr>

<!-- SECTION 5 -->
<h2 id="sec5">5. CytoSig vs LinCytoSig vs SecAct Comparison</h2>

<h3>5.1 Method Overview</h3>

<table>
  <tr>
    <th>Method</th><th>Targets</th><th>Genes</th><th>Specificity</th><th>Selection</th>
  </tr>
  <tr>
    <td><span class="badge blue">CytoSig</span></td>
    <td>43 cytokines</td><td>4,881 curated</td><td>Global (all cell types)</td><td>&mdash;</td>
  </tr>
  <tr>
    <td><span class="badge amber">LinCytoSig (orig)</span></td>
    <td>178 (45 CT &times; cytokines)</td><td>All ~20K</td><td>Cell-type specific</td><td>Matched cell type</td>
  </tr>
  <tr>
    <td><span style="background:#F59E0B;color:#fff;padding:1px 8px;border-radius:3px;font-size:0.85em">LinCytoSig (gene-filtered)</span></td>
    <td>178</td><td>4,881 (CytoSig overlap)</td><td>Cell-type specific</td><td>Matched cell type</td>
  </tr>
  <tr>
    <td><span style="background:#B45309;color:#fff;padding:1px 8px;border-radius:3px;font-size:0.85em">LinCytoSig Best (combined)</span></td>
    <td>43 (1 per cytokine)</td><td>All ~20K</td><td>Best CT per cytokine</td><td>Max combined GTEx+TCGA &rho;</td>
  </tr>
  <tr>
    <td><span style="background:#92400E;color:#fff;padding:1px 8px;border-radius:3px;font-size:0.85em">LinCytoSig Best (comb+filt)</span></td>
    <td>43 (1 per cytokine)</td><td>4,881 (CytoSig overlap)</td><td>Best CT per cytokine</td><td>Max combined &rho; (filtered)</td>
  </tr>
  <tr>
    <td><span style="background:#7C3AED;color:#fff;padding:1px 8px;border-radius:3px;font-size:0.85em">LinCytoSig Best (GTEx)</span></td>
    <td>43 (1 per cytokine)</td><td>All ~20K</td><td>Best CT per cytokine</td><td>Max GTEx &rho;</td>
  </tr>
  <tr>
    <td><span style="background:#EC4899;color:#fff;padding:1px 8px;border-radius:3px;font-size:0.85em">LinCytoSig Best (TCGA)</span></td>
    <td>43 (1 per cytokine)</td><td>All ~20K</td><td>Best CT per cytokine</td><td>Max TCGA &rho;</td>
  </tr>
  <tr>
    <td><span style="background:#A78BFA;color:#fff;padding:1px 8px;border-radius:3px;font-size:0.85em">LinCytoSig Best (GTEx+filt)</span></td>
    <td>43 (1 per cytokine)</td><td>4,881 (CytoSig overlap)</td><td>Best CT per cytokine</td><td>Max GTEx &rho; (filtered)</td>
  </tr>
  <tr>
    <td><span style="background:#F9A8D4;color:#333;padding:1px 8px;border-radius:3px;font-size:0.85em">LinCytoSig Best (TCGA+filt)</span></td>
    <td>43 (1 per cytokine)</td><td>4,881 (CytoSig overlap)</td><td>Best CT per cytokine</td><td>Max TCGA &rho; (filtered)</td>
  </tr>
  <tr>
    <td><span class="badge green">SecAct</span></td>
    <td>1,170 secreted proteins</td><td>Spatial Moran&rsquo;s I</td><td>Global (all cell types)</td><td>&mdash;</td>
  </tr>
</table>
<p style="margin-top:0.5em;font-size:0.9em;color:#6B7280;">
  <strong>Gene filter:</strong> LinCytoSig signatures restricted from ~20K to CytoSig&rsquo;s 4,881 curated genes.
  <strong>Best selection:</strong> For each cytokine, test all cell-type-specific LinCytoSig signatures and select the one with the highest bulk RNA-seq correlation. &ldquo;Combined&rdquo; uses pooled GTEx+TCGA; &ldquo;GTEx&rdquo; and &ldquo;TCGA&rdquo; select independently per bulk dataset.
  &ldquo;+filt&rdquo; variants apply the same cell-type selection but restrict to CytoSig gene space.
  See <a href="methodology.html">LinCytoSig Methodology</a> for details.
</p>

<div class="plotly-container">
  <div class="controls">
    <label>View:</label>
    <select id="method-boxplot-view" onchange="updateMethodBoxplot()">
      <option value="all">All Atlases</option>
    </select>
  </div>
  <div id="method-boxplot-chart" style="height:600px;"></div>
  <div class="caption"><strong>Figure 10.</strong> Ten-way signature method comparison at matched (cell type, cytokine) pair level across 4 combined atlases. All 10 methods are evaluated on the <em>same set</em> of matched pairs per atlas (identical n). Use dropdown to view individual atlas boxplots. For LinCytoSig construction, see <a href="methodology.html">LinCytoSig Methodology</a>.</div>
</div>

<div class="callout">
<p><strong>Ten methods compared on identical matched pairs across 4 combined atlases:</strong></p>
<ol>
  <li><strong><span style="color:#2563EB">CytoSig</span></strong> &mdash; 43 cytokines, 4,881 curated genes, global (all cell types)</li>
  <li><strong><span style="color:#D97706">LinCytoSig (orig)</span></strong> &mdash; cell-type-matched signatures, all ~20K genes</li>
  <li><strong><span style="color:#F59E0B">LinCytoSig (gene-filtered)</span></strong> &mdash; cell-type-matched signatures, restricted to CytoSig&rsquo;s 4,881 genes</li>
  <li><strong><span style="color:#B45309">LinCytoSig Best (combined)</span></strong> &mdash; best cell-type signature per cytokine (selected by combined GTEx+TCGA bulk &rho;), all ~20K genes</li>
  <li><strong><span style="color:#92400E">LinCytoSig Best (comb+filt)</span></strong> &mdash; best combined bulk signature, restricted to 4,881 genes</li>
  <li><strong><span style="color:#7C3AED">LinCytoSig Best (GTEx)</span></strong> &mdash; best per cytokine selected by GTEx-only bulk &rho;, all ~20K genes</li>
  <li><strong><span style="color:#EC4899">LinCytoSig Best (TCGA)</span></strong> &mdash; best per cytokine selected by TCGA-only bulk &rho;, all ~20K genes</li>
  <li><strong><span style="color:#A78BFA">LinCytoSig Best (GTEx+filt)</span></strong> &mdash; GTEx-selected best, restricted to 4,881 genes</li>
  <li><strong><span style="color:#F9A8D4">LinCytoSig Best (TCGA+filt)</span></strong> &mdash; TCGA-selected best, restricted to 4,881 genes</li>
  <li><strong><span style="color:#059669">SecAct</span></strong> &mdash; 1,170 secreted proteins (Moran&rsquo;s I), subset matching CytoSig targets</li>
</ol>
<p><strong>Key findings:</strong></p>
<ul>
  <li><strong>SecAct achieves the highest median &rho;</strong> across all 4 combined atlases, benefiting from spatial-transcriptomics-derived signatures.</li>
  <li><strong>CytoSig outperforms most LinCytoSig variants</strong> at donor level, with one notable exception: scAtlas Normal Best-orig (0.298) exceeds CytoSig (0.216).</li>
  <li><strong>Gene filtering improves LinCytoSig</strong> in most atlases (CIMA +102%, Inflammation Main), confirming noise reduction from restricting the gene space.</li>
  <li><strong>GTEx-selected best</strong> performs comparably to combined-selected in most atlases but slightly better in scAtlas Cancer (0.300 vs 0.275). <strong>TCGA-selected best</strong> generally underperforms other selection strategies, suggesting GTEx&rsquo;s broader tissue coverage provides more generalizable selections.</li>
  <li><strong>Gene filtering of GTEx/TCGA-selected:</strong> GTEx+filt and TCGA+filt show mixed results &mdash; filtering sometimes improves (e.g., TCGA+filt in Inflammation Main: 0.260 vs TCGA-orig 0.168) but can also reduce performance, indicating the optimal gene space depends on both the selection dataset and atlas context.</li>
  <li><strong>General ranking:</strong> SecAct &gt; CytoSig &gt; LinCytoSig Best variants &gt; LinCytoSig (filt) &gt; LinCytoSig (orig), though atlas-specific exceptions exist.</li>
</ul>
</div>

<h3>5.2 Effect of Aggregation Level</h3>

<div class="callout amber">
<p><strong>Methodology:</strong> At each cell-type aggregation level (CIMA: L1&ndash;L4 = 7&ndash;73 cell types; Inflammation: L1&ndash;L2; scAtlas: CT1&ndash;CT2 = coarse/fine), we match CytoSig, LinCytoSig, and SecAct on identical (cytokine, cell type) pairs &mdash; using the <em>exact same pseudobulk samples</em> and <em>identical n</em> for all three methods. For each pair, Spearman &rho; measures agreement between predicted activity and target gene expression. If lineage-specific aggregation helps, LinCytoSig should increasingly outperform CytoSig as cell-type resolution increases (L1 &rarr; L4).</p>
</div>

<h4>5.2.1 Distribution at Each Level</h4>
<div class="plotly-container">
  <div class="controls">
    <label>Atlas:</label>
    <select id="level-comp-atlas" onchange="updateLevelComp()">
    </select>
  </div>
  <div id="level-comp-chart" style="height:550px;"></div>
  <div class="caption"><strong>Figure 11.</strong> Distribution of Spearman &rho; at each cell-type aggregation level. All three methods evaluated on identical matched pairs per level. Finer levels (more cell types) should theoretically favor lineage-specific methods.</div>
</div>

<h4>5.2.2 Summary</h4>
<div id="level-summary-container"></div>
<p style="font-size:0.9em;color:#6B7280;">n = number of three-way matched pairs. &Delta;&rho; = LinCytoSig &minus; competitor (negative = LinCytoSig underperforms).</p>

<!-- Per-celltype summary table (from method_comparison.json, finest level) -->
<h4>5.2.3 Which Cell Types Benefit?</h4>
<div id="celltype-table-container" style="max-height:500px;overflow-y:auto;"></div>
<p style="font-size:0.9em;color:#6B7280;">Aggregated across all atlases at finest celltype level. Green = LinCytoSig wins more; red = LinCytoSig loses more.</p>

<!-- Per-cytokine summary table -->
<h4>5.2.4 Which Cytokines Benefit?</h4>
<div id="cytokine-table-container" style="max-height:500px;overflow-y:auto;"></div>
<p style="font-size:0.9em;color:#6B7280;">Sorted by mean &Delta;&rho; vs CytoSig (best to worst).</p>

<div class="callout">
<p><strong>Key finding: Lineage-specific aggregation provides no systematic advantage at any level.</strong></p>
<ul>
  <li><strong>At every level, LinCytoSig underperforms CytoSig</strong> (mean &Delta;&rho; ranges from &minus;0.08 at coarse L1 to &minus;0.02 at fine L4 in CIMA). Finer cell types reduce the gap slightly but never close it.</li>
  <li><strong>SecAct wins at every level</strong> in CIMA and scAtlas. In Inflammation Main L2, LinCytoSig is nearly tied with SecAct (&Delta;&rho; = +0.01) but still loses to CytoSig.</li>
  <li><strong>Per cell type:</strong> Only 5 of 43 cell types show consistent LinCytoSig advantage vs CytoSig (NK Cell, Basophil, DC, Trophoblast, Arterial Endothelial). No cell type beats SecAct.</li>
  <li><strong>Interpretation:</strong> CytoSig&rsquo;s global signature, derived from median log2FC across all cell types, already captures the dominant transcriptional response. Restricting to a single cell type&rsquo;s response introduces noise from small sample sizes without gaining meaningful lineage specificity. The hypothesis that finer resolution should favor LinCytoSig is not supported by the data.</li>
</ul>
</div>

<h3>5.3 SecAct: Breadth Over Depth</h3>

<ul>
  <li><strong>Highest median &rho;</strong> in organ-level analyses (scAtlas normal: 0.307, cancer: 0.363)</li>
  <li><strong>Highest median &rho;</strong> in bulk RNA-seq (GTEx: 0.395, TCGA: 0.415)</li>
  <li><strong>97.1% positive correlation</strong> in TCGA</li>
  <li><strong>Wins decisively at celltype level</strong> against both CytoSig and LinCytoSig in scAtlas (19/3 wins vs CytoSig in scAtlas Normal, 20/2 in Cancer)</li>
</ul>

<div class="plotly-container">
  <div id="celltype-delta-rho-chart" style="height:900px;"></div>
  <div class="caption"><strong>Figure 12.</strong> Per-celltype mean &Delta;&rho; (LinCytoSig &minus; CytoSig) aggregated across 4 atlases at donor &times; celltype level. Orange = LinCytoSig advantage; blue = CytoSig advantage. Error bars show SEM.</div>
</div>

<hr>

<!-- SECTION 6 -->
<h2 id="sec6">6. Key Takeaways for Scientific Discovery</h2>

<h3>6.1 What CytoAtlas Enables</h3>
<ol>
  <li><strong>Quantitative cytokine activity per cell type per disease</strong> &mdash; 43 CytoSig cytokines + 1,170 SecAct secreted proteins across 29M cells</li>
  <li><strong>Cross-disease comparison</strong> &mdash; same signatures validated across 20 diseases, 35 organs, 15 cancer types</li>
  <li><strong>Independent perturbation comparison</strong> &mdash; parse_10M provides 90 cytokine perturbations &times; 12 donors &times; 18 cell types for independent comparison with CytoSig predictions</li>
  <li><strong>Multi-level validation</strong> &mdash; donor, donor &times; celltype, bulk RNA-seq (GTEx/TCGA), and resampled bootstrap validation across 6 atlases</li>
</ol>

<h3>6.2 Limitations</h3>
<ol>
  <li><strong>Linear model:</strong> Cannot capture non-linear cytokine interactions</li>
  <li><strong>Transcriptomics-only:</strong> Post-translational regulation invisible</li>
  <li><strong>Signature matrix bias:</strong> Underrepresented cell types have weaker signatures</li>
  <li><strong>Validation metric:</strong> Expression-activity correlation underestimates true accuracy (signatures acting through downstream effectors are not captured)</li>
</ol>

<h3>6.3 Future Directions</h3>
<ol>
  <li>scGPT cohort integration (~35M cells)</li>
  <li>cellxgene Census integration</li>
  <li>Classification of cytokine blocking therapy</li>
</ol>

<hr>

<!-- SECTION 7 -->
<h2 id="sec7">7. Appendix: Technical Specifications</h2>

<h3>A. Computational Infrastructure</h3>
<ul>
  <li><strong>GPU:</strong> NVIDIA A100 80GB (SLURM gpu partition)</li>
  <li><strong>Memory:</strong> 256&ndash;512GB host RAM per node</li>
  <li><strong>Pipeline:</strong> 24 Python scripts, 18 pipeline subpackages (~18.7K lines)</li>
  <li><strong>API:</strong> 262 REST endpoints across 17 routers</li>
  <li><strong>Frontend:</strong> 12 pages, 122 source files, 11.4K LOC</li>
</ul>

<h3>B. Statistical Methods</h3>
<ul>
  <li><strong>Activity inference:</strong> Ridge regression (&lambda;=5&times;10<sup>5</sup>, z-score normalization, permutation-based significance)</li>
  <li><strong>Correlation:</strong> Spearman rank correlation</li>
  <li><strong>Multiple testing:</strong> Benjamini-Hochberg FDR (q &lt; 0.05)</li>
  <li><strong>Bootstrap:</strong> 100&ndash;1000 resampling iterations</li>
  <li><strong>Differential:</strong> Wilcoxon rank-sum test with effect size</li>
</ul>

</div><!-- content -->
</div><!-- container -->

<!-- Lightbox for static images -->
<div class="lightbox-overlay" id="lightbox">
  <button class="lightbox-close" onclick="closeLightbox()">&times;</button>
  <img src="" alt="">
</div>

<script>
// ═══════════════════════════════════════════════════════════════════════════
// DATA
// ═══════════════════════════════════════════════════════════════════════════
var DATA = {data_json};

var PLOTLY_CONFIG = {{responsive: true, displayModeBar: true, modeBarButtonsToRemove: ['lasso2d','select2d']}};
var SIG_NAMES = {{'cytosig':'CytoSig','lincyto_best':'LinCytoSig Best (comb+filt)','lincytosig':'LinCytoSig','secact':'SecAct'}};

// ═══════════════════════════════════════════════════════════════════════════
// 4.1 SUMMARY TABLE (embedded HTML)
// ═══════════════════════════════════════════════════════════════════════════
(function() {{
  var rows = DATA.summary;
  var sigBg = {{'CYTOSIG':'#DBEAFE','SECACT':'#D1FAE5'}};
  var html = '<table><tr>';
  var cols = ['Atlas','Signature','Primary Level','N Targets','Median \\u03c1','Mean \\u03c1','Std \\u03c1','Min \\u03c1','Max \\u03c1','% Significant','% Positive'];
  cols.forEach(function(c) {{ html += '<th>' + c + '</th>'; }});
  html += '</tr>';
  rows.forEach(function(r) {{
    var bg = sigBg[r.signature] || '';
    html += '<tr>';
    html += '<td style="background:'+bg+'">' + r.atlas + '</td>';
    html += '<td style="background:'+bg+';font-weight:600">' + r.signature + '</td>';
    html += '<td style="background:'+bg+';font-size:0.85em">' + r.primary_level + '</td>';
    html += '<td style="background:'+bg+'">' + r.n_targets + '</td>';
    html += '<td style="background:'+bg+';font-weight:600">' + r.median_rho + '</td>';
    html += '<td style="background:'+bg+'">' + r.mean_rho + '</td>';
    html += '<td style="background:'+bg+'">' + r.std_rho + '</td>';
    html += '<td style="background:'+bg+'">' + r.min_rho + '</td>';
    html += '<td style="background:'+bg+'">' + r.max_rho + '</td>';
    html += '<td style="background:'+bg+'">' + r.pct_sig + '%</td>';
    html += '<td style="background:'+bg+'">' + r.pct_pos + '%</td>';
    html += '</tr>';
  }});
  html += '</table>';
  html += '<p style="font-size:0.85em;color:#555;margin-top:8px;"><strong>Column definitions:</strong> ';
  html += '<strong>Primary Level</strong> = aggregation level where samples are fully independent (each donor counted once). ';
  html += '<strong>Median/Mean \\u03c1</strong> = Spearman rank correlation between predicted activity and target gene expression across donors. ';
  html += '<strong>% Significant</strong> = fraction of targets with p &lt; 0.05. ';
  html += '<strong>% Positive</strong> = fraction of targets with \\u03c1 &gt; 0 (predicted activity positively correlates with target expression).</p>';
  document.getElementById('summary-table-container').innerHTML = html;
}})();

// ═══════════════════════════════════════════════════════════════════════════
// 4.3 BOXPLOT with Total / Matched tabs + statistical tests
// ═══════════════════════════════════════════════════════════════════════════
function formatPval(p) {{
  if (p === null || p === undefined) return '';
  if (p < 0.001) return 'p < 0.001';
  if (p < 0.01) return 'p = ' + p.toFixed(3);
  if (p < 0.05) return 'p = ' + p.toFixed(3);
  return 'p = ' + p.toFixed(2) + ' (ns)';
}}
function formatQval(q) {{
  if (q === null || q === undefined) return '';
  if (q < 0.001) return 'q < 0.001';
  if (q < 0.01) return 'q = ' + q.toFixed(3);
  if (q < 0.05) return 'q = ' + q.toFixed(3);
  return 'q = ' + q.toFixed(2) + ' (ns)';
}}
function sigStars(p) {{
  if (p === null || p === undefined) return '';
  if (p < 0.001) return '***';
  if (p < 0.01) return '**';
  if (p < 0.05) return '*';
  return 'ns';
}}
window.renderBoxplot = function(mode) {{
  // Update tab styling
  document.querySelectorAll('#boxplot-tabs .tab-btn').forEach(function(b) {{ b.classList.remove('active'); }});
  var btns = document.querySelectorAll('#boxplot-tabs .tab-btn');
  if (mode === 'total') btns[0].classList.add('active');
  else btns[1].classList.add('active');

  var bd = DATA.boxplot;
  var atlases = DATA.atlasLabels;
  var traces = [];
  var annotations = [];

  if (mode === 'total') {{
    // Total: CytoSig (all) vs SecAct (all)
    var cfgs = [
      {{key:'cytosig', name:'CytoSig', color:'#2563EB'}},
      {{key:'secact', name:'SecAct', color:'#059669'}},
    ];
    cfgs.forEach(function(cfg) {{
      var y = [], x = [];
      atlases.forEach(function(a) {{
        if (bd[a] && bd[a][cfg.key]) {{
          bd[a][cfg.key].forEach(function(v) {{ y.push(v); x.push(a); }});
        }}
      }});
      var n = bd[atlases[0]] && bd[atlases[0]][cfg.key] ? bd[atlases[0]][cfg.key].length : 0;
      traces.push({{
        type:'box', y:y, x:x,
        name: cfg.name + ' (n=' + (cfg.key === 'cytosig' ? '43' : '1,170') + ')',
        marker:{{color:cfg.color}}, boxpoints:false, line:{{width:1.5}},
      }});
    }});
    // Add significance annotations per atlas
    atlases.forEach(function(a, i) {{
      var st = bd[a] ? bd[a].stats_total : null;
      if (st && st.q_bh !== undefined) {{
        annotations.push({{
          x: a, y: 1.05, xref:'x', yref:'paper',
          text: '<b>' + sigStars(st.q_bh) + '</b><br><span style="font-size:9px">' + formatQval(st.q_bh) + '</span>',
          showarrow: false, font:{{size:11}}, align:'center',
        }});
      }}
    }});
    document.getElementById('boxplot-caption').innerHTML =
      '<strong>Figure 3.</strong> Spearman \\u03c1 distributions: CytoSig (43 targets) vs SecAct (1,170 targets) across atlases. ' +
      'Independence-corrected: GTEx/TCGA use median-of-medians (one \\u03c1 per target). ' +
      'Mann-Whitney U test, BH-FDR corrected across 6 atlases. Significance: *** q&lt;0.001, ** q&lt;0.01, * q&lt;0.05, ns = not significant.';
  }} else {{
    // Matched: CytoSig (32) vs SecAct (32) on shared targets
    var cfgs = [
      {{key:'cytosig_matched', name:'CytoSig (matched)', color:'#2563EB'}},
      {{key:'secact_matched', name:'SecAct (matched)', color:'#059669'}},
    ];
    cfgs.forEach(function(cfg) {{
      var y = [], x = [];
      atlases.forEach(function(a) {{
        if (bd[a] && bd[a][cfg.key]) {{
          bd[a][cfg.key].forEach(function(v) {{ y.push(v); x.push(a); }});
        }}
      }});
      traces.push({{
        type:'box', y:y, x:x,
        name: cfg.name + ' (n=32)',
        marker:{{color:cfg.color}}, boxpoints:false, line:{{width:1.5}},
      }});
    }});
    // Add significance annotations per atlas
    atlases.forEach(function(a, i) {{
      var sm = bd[a] ? bd[a].stats_matched : null;
      if (sm && sm.q_bh !== undefined) {{
        annotations.push({{
          x: a, y: 1.05, xref:'x', yref:'paper',
          text: '<b>' + sigStars(sm.q_bh) + '</b><br><span style="font-size:9px">' + formatQval(sm.q_bh) + '</span>',
          showarrow: false, font:{{size:11}}, align:'center',
        }});
      }}
    }});
    document.getElementById('boxplot-caption').innerHTML =
      '<strong>Figure 3.</strong> Spearman \\u03c1 distributions: CytoSig vs SecAct on 32 matched targets shared between both methods (22 direct + 10 alias-resolved). ' +
      'Wilcoxon signed-rank test (paired by target), BH-FDR corrected across 6 atlases. Significance: *** q&lt;0.001, ** q&lt;0.01, * q&lt;0.05, ns = not significant.';
  }}

  var title = mode === 'total'
    ? 'Total: CytoSig (43 targets) vs SecAct (1,170 targets)'
    : 'Matched: CytoSig vs SecAct (32 shared targets)';

  Plotly.newPlot('boxplot-chart', traces, {{
    title: title,
    yaxis: {{title: 'Spearman \\u03c1', zeroline: true, zerolinecolor: '#ccc'}},
    boxmode: 'group',
    legend: {{orientation:'h', y:1.15, x:0.5, xanchor:'center'}},
    margin: {{t:100, b:80}},
    annotations: annotations,
  }}, PLOTLY_CONFIG);
}};
renderBoxplot('total');

// ═══════════════════════════════════════════════════════════════════════════
// 4.2 STRATIFIED VALIDATION (per-tissue/per-cancer)
// ═══════════════════════════════════════════════════════════════════════════
var currentStratifiedMode = 'total';
window.renderStratified = function(mode) {{
  if (mode) currentStratifiedMode = mode;
  else mode = currentStratifiedMode;
  document.querySelectorAll('#stratified-tabs .tab-btn').forEach(function(b) {{ b.classList.remove('active'); }});
  var btns = document.querySelectorAll('#stratified-tabs .tab-btn');
  if (mode === 'total') btns[0].classList.add('active');
  else btns[1].classList.add('active');

  var dataset = document.getElementById('stratified-dataset-select').value;
  var sd = DATA.stratified[dataset];
  if (!sd) return;

  var strata = sd.strata;
  var d = sd.data;
  // Reverse so first stratum appears at top of horizontal chart
  var labels = strata.slice().reverse();

  // Select correct data arrays based on mode
  var cytoKey = mode === 'total' ? 'cytosig' : 'cytosig_matched';
  var secKey = mode === 'total' ? 'secact' : 'secact_matched';
  var statsKey = mode === 'total' ? 'stats_total' : 'stats_matched';

  // Build boxplot traces: repeat category label for each rho value
  var cytoY = [], cytoX = [], secY = [], secX = [];
  labels.forEach(function(s) {{
    var entry = d[s];
    var cytoArr = entry[cytoKey] || [];
    var secArr = entry[secKey] || [];
    cytoArr.forEach(function(v) {{ cytoY.push(s); cytoX.push(v); }});
    secArr.forEach(function(v) {{ secY.push(s); secX.push(v); }});
  }});

  var nLabel = mode === 'total' ? '43 vs ~1,170' : '32 matched';
  var traces = [
    {{
      type: 'box', orientation: 'h',
      y: cytoY, x: cytoX,
      name: 'CytoSig', marker: {{color: '#2563EB'}},
      line: {{width: 1.5}}, boxpoints: false,
    }},
    {{
      type: 'box', orientation: 'h',
      y: secY, x: secX,
      name: 'SecAct', marker: {{color: '#059669'}},
      line: {{width: 1.5}}, boxpoints: false,
    }},
  ];

  // Significance annotations — anchored to right margin to avoid overlap
  var annotations = [];
  labels.forEach(function(s) {{
    var entry = d[s];
    var st = entry[statsKey];
    if (st && st.q_bh !== undefined && st.q_bh !== null) {{
      annotations.push({{
        x: 1.01, y: s, xref: 'paper', yref: 'y',
        text: sigStars(st.q_bh), showarrow: false,
        font: {{size: 9, color: '#374151'}}, xanchor: 'left',
      }});
    }}
  }});

  var height = Math.max(500, strata.length * 28 + 140);
  document.getElementById('stratified-chart').style.height = height + 'px';

  Plotly.newPlot('stratified-chart', traces, {{
    xaxis: {{title: 'Spearman \\u03c1', zeroline: true, zerolinecolor: '#ccc'}},
    yaxis: {{automargin: true, tickfont: {{size: 10}}}},
    boxmode: 'group',
    legend: {{orientation: 'h', y: 1.02, x: 0.5, xanchor: 'center'}},
    margin: {{t: 30, l: 180, b: 50, r: 50}},
    annotations: annotations,
  }}, PLOTLY_CONFIG);

  var dsType = dataset === 'GTEx' ? 'tissue' : 'cancer type';
  var capMode = mode === 'total'
    ? 'All targets compared (' + nLabel + '; Mann-Whitney U test, BH-FDR corrected).'
    : 'Matched targets (' + nLabel + '; Wilcoxon signed-rank test, BH-FDR corrected).';
  document.getElementById('stratified-caption').innerHTML =
    '<strong>Figure 2.</strong> Per-' + dsType +
    ' Spearman \\u03c1 distributions: CytoSig vs SecAct. ' + capMode +
    ' Significance: *** q&lt;0.001, ** q&lt;0.01, * q&lt;0.05.';
}};
renderStratified('total');

// ═══════════════════════════════════════════════════════════════════════════
// 4.4 GOOD/BAD CORRELATIONS (interactive with signature + atlas dropdown)
// ═══════════════════════════════════════════════════════════════════════════
(function() {{
  var gb = DATA.goodBad;
  var atlasSel = document.getElementById('goodbad-atlas-select');
  var sigSel = document.getElementById('goodbad-sig-select');
  // Populate atlas dropdown in consistent order matching boxplot (Section 4.2)
  var firstSig = Object.keys(gb)[0];
  var availableAtlases = Object.keys(gb[firstSig]);
  DATA.atlasLabels.forEach(function(a) {{
    if (availableAtlases.indexOf(a) !== -1) {{
      var opt = document.createElement('option');
      opt.value = a; opt.textContent = a;
      atlasSel.appendChild(opt);
    }}
  }});
  window.updateGoodBad = function() {{
    var sig = sigSel.value;
    var atlas = atlasSel.value;
    var sigData = gb[sig];
    if (!sigData || !sigData[atlas]) return;
    var d = sigData[atlas];
    // Top 15 (reverse for horizontal bar layout)
    var topTargets = d.top.map(function(r){{return r.target;}}).reverse();
    var topRhos = d.top.map(function(r){{return r.rho;}}).reverse();
    // Bottom 15
    var botTargets = d.bottom.map(function(r){{return r.target;}}).reverse();
    var botRhos = d.bottom.map(function(r){{return r.rho;}}).reverse();
    var sigColor = DATA.sigColors[sig] || '#6B7280';
    var traces = [
      {{type:'bar', y:topTargets, x:topRhos, orientation:'h',
        marker:{{color:'#DC2626'}},
        xaxis:'x', yaxis:'y', text:topRhos.map(function(r){{return r.toFixed(3);}}), textposition:'outside',
        hovertemplate:'<b>%{{y}}</b><br>\\u03c1 = %{{x:.4f}}<extra></extra>'}},
      {{type:'bar', y:botTargets, x:botRhos, orientation:'h',
        marker:{{color:'#2563EB'}},
        xaxis:'x2', yaxis:'y2', text:botRhos.map(function(r){{return r.toFixed(3);}}), textposition:'outside',
        hovertemplate:'<b>%{{y}}</b><br>\\u03c1 = %{{x:.4f}}<extra></extra>'}},
    ];
    Plotly.newPlot('goodbad-chart', traces, {{
      title: atlas + ' — ' + SIG_NAMES[sig] + ': Best & Worst Correlated Targets',
      grid: {{rows:1, columns:2, pattern:'independent'}},
      xaxis: {{title:'Spearman \\u03c1', domain:[0, 0.45]}},
      yaxis: {{automargin: true}},
      xaxis2: {{title:'Spearman \\u03c1', domain:[0.55, 1]}},
      yaxis2: {{automargin: true}},
      annotations: [
        {{text:'<b>Top 15 (Best)</b>', xref:'paper', yref:'paper', x:0.22, y:1.06, showarrow:false, font:{{color:'#DC2626',size:14}}}},
        {{text:'<b>Bottom 15 (Worst)</b>', xref:'paper', yref:'paper', x:0.78, y:1.06, showarrow:false, font:{{color:'#2563EB',size:14}}}},
      ],
      showlegend: false, margin: {{t:80, l:120, r:20}},
    }}, PLOTLY_CONFIG);
  }};
  updateGoodBad();
}})();

// ═══════════════════════════════════════════════════════════════════════════
// 4.5 CROSS-ATLAS CONSISTENCY (interactive Plotly line chart)
// ═══════════════════════════════════════════════════════════════════════════
(function() {{
  var cd = DATA.consistency;
  var atlases = DATA.atlasLabels;
  var traces = [];
  Object.keys(cd).forEach(function(key) {{
    var d = cd[key];
    var isSecact = d.sig_type === 'secact';
    var displayName = isSecact ? key.replace('_secact','') + ' (SecAct)' : key;
    traces.push({{
      type: 'scatter', mode: 'lines+markers',
      x: atlases, y: d.rhos,
      name: displayName,
      line: {{color: d.color, width: isSecact ? 1.2 : 2, dash: isSecact ? 'dash' : 'solid'}},
      marker: {{size: isSecact ? 5 : 8, symbol: isSecact ? 'square' : 'circle'}},
      opacity: isSecact ? 0.6 : 0.9,
      legendgroup: isSecact ? key.replace('_secact','') : key,
      hovertemplate: '<b>' + displayName + '</b><br>%{{x}}: \\u03c1=%{{y:.3f}}<extra></extra>',
    }});
  }});
  Plotly.newPlot('consistency-chart', traces, {{
    title: 'Cross-Atlas Consistency: CytoSig (solid) vs SecAct (dashed)',
    yaxis: {{title: 'Spearman \\u03c1', zeroline: true, zerolinecolor: '#ccc'}},
    legend: {{font: {{size: 9}}}},
    hovermode: 'closest',
    margin: {{t:60, b:80, r:20}},
  }}, PLOTLY_CONFIG);
}})();

// ═══════════════════════════════════════════════════════════════════════════
// 4.6 AGGREGATION LEVELS (Total / Matched tabs, per atlas)
// ═══════════════════════════════════════════════════════════════════════════
(function() {{
  var ld = DATA.levels;
  var atlasNames = Object.keys(ld);
  var atlasSel = document.getElementById('levels-atlas-select');
  atlasNames.forEach(function(a) {{
    var opt = document.createElement('option');
    opt.value = a; opt.textContent = a;
    atlasSel.appendChild(opt);
  }});
  var currentMode = 'total';
  window.renderLevels = function(mode) {{
    if (mode) currentMode = mode;
    // Update tab active state
    var tabs = document.querySelectorAll('#levels-tabs .tab-btn');
    tabs.forEach(function(btn) {{
      btn.classList.remove('active');
      if ((currentMode === 'total' && btn.textContent === 'Total') ||
          (currentMode === 'matched' && btn.textContent === 'Matched'))
        btn.classList.add('active');
    }});
    var atlasName = atlasSel.value;
    var atlasData = ld[atlasName];
    if (!atlasData) return;
    // Select data keys based on mode
    var cytoKey = currentMode === 'total' ? 'cytosig' : 'cytosig_matched';
    var secKey = currentMode === 'total' ? 'secact' : 'secact_matched';
    var cytoLabel = currentMode === 'total' ? 'CytoSig (all targets)' : 'CytoSig (32 matched)';
    var secLabel = currentMode === 'total' ? 'SecAct (all targets)' : 'SecAct (32 matched)';
    var sigConfigs = [
      {{key: cytoKey, name: cytoLabel, color: '#2563EB'}},
      {{key: secKey, name: secLabel, color: '#059669'}},
    ];
    var statsKey = currentMode === 'total' ? 'total' : 'matched';
    // Collect all levels
    var allLevels = [];
    sigConfigs.forEach(function(cfg) {{
      if (atlasData[cfg.key]) {{
        Object.keys(atlasData[cfg.key]).forEach(function(lv) {{
          if (allLevels.indexOf(lv) === -1) allLevels.push(lv);
        }});
      }}
    }});
    var traces = [];
    sigConfigs.forEach(function(cfg) {{
      var sigData = atlasData[cfg.key];
      if (!sigData) return;
      var y = [], x = [];
      allLevels.forEach(function(level) {{
        if (sigData[level]) {{
          sigData[level].rhos.forEach(function(v) {{
            y.push(v);
            x.push(level);
          }});
        }}
      }});
      traces.push({{
        type: 'box', y: y, x: x,
        name: cfg.name,
        marker: {{color: cfg.color}},
        boxpoints: false, line: {{width: 1.5}},
      }});
    }});
    // Significance annotations above each level
    var annotations = [];
    var levelStats = atlasData.stats || {{}};
    allLevels.forEach(function(level) {{
      var ls = levelStats[level];
      if (ls && ls[statsKey] && ls[statsKey].q_bh !== undefined) {{
        var q = ls[statsKey].q_bh;
        annotations.push({{
          x: level, y: 1.05, xref: 'x', yref: 'paper',
          text: '<b>' + sigStars(q) + '</b><br><span style="font-size:9px">' + formatQval(q) + '</span>',
          showarrow: false, font: {{size: 11}}, align: 'center',
        }});
      }}
    }});
    var testLabel = currentMode === 'total' ? 'Mann-Whitney U test' : 'Wilcoxon signed-rank test';
    Plotly.newPlot('levels-chart', traces, {{
      yaxis: {{title: 'Spearman \\u03c1', zeroline: true, zerolinecolor: '#ccc'}},
      boxmode: 'group',
      legend: {{orientation:'h', y:1.12, x:0.5, xanchor:'center'}},
      margin: {{t:70, b:100}},
      xaxis: {{tickangle: -20}},
      annotations: annotations,
    }}, PLOTLY_CONFIG);
    // Update caption
    var cap = document.getElementById('levels-caption');
    if (currentMode === 'total') {{
      cap.innerHTML = '<strong>Figure 6.</strong> ' + atlasName + ' \u2014 CytoSig (all ~43 targets) vs SecAct (all ~1,170 targets) across aggregation levels. Mann-Whitney U test, BH-FDR corrected. Significance: *** q&lt;0.001, ** q&lt;0.01, * q&lt;0.05, ns = not significant.';
    }} else {{
      cap.innerHTML = '<strong>Figure 6.</strong> ' + atlasName + ' \u2014 CytoSig vs SecAct restricted to 32 shared targets across aggregation levels. Wilcoxon signed-rank test, BH-FDR corrected. Significance: *** q&lt;0.001, ** q&lt;0.01, * q&lt;0.05, ns = not significant.';
    }}
  }};
  renderLevels('total');
}})();

// ═══════════════════════════════════════════════════════════════════════════
// SCATTER PLOTS (interactive with dropdown)
// ═══════════════════════════════════════════════════════════════════════════
(function() {{
  var sd = DATA.scatter;
  var targetSel = document.getElementById('scatter-target-select');
  var atlasSel = document.getElementById('scatter-atlas-select');
  // Populate dropdowns
  var allTargets = new Set();
  var allAtlases = Object.keys(sd);
  allAtlases.forEach(function(a) {{
    Object.keys(sd[a]).forEach(function(t) {{ allTargets.add(t); }});
    var opt = document.createElement('option');
    opt.value = a; opt.textContent = a;
    atlasSel.appendChild(opt);
  }});
  ['IFNG','IL1B','TNFA','TGFB1','IL6','IL10','VEGFA','CD40L','TRAIL','HGF'].forEach(function(t) {{
    if (allTargets.has(t)) {{
      var opt = document.createElement('option');
      opt.value = t; opt.textContent = t;
      targetSel.appendChild(opt);
    }}
  }});
  window.updateScatter = function() {{
    var target = targetSel.value;
    var atlas = atlasSel.value;
    var sigSel = document.getElementById('scatter-sig-select');
    var sig = sigSel ? sigSel.value : 'cytosig';
    var dataKey = (sig === 'secact') ? target + '_secact' : target;
    var sigLabel = sig === 'secact' ? 'SecAct' : 'CytoSig';
    if (!sd[atlas] || !sd[atlas][dataKey]) {{
      Plotly.newPlot('scatter-chart', [], {{title: 'No data for ' + target + ' (' + sigLabel + ') in ' + atlas, margin:{{t:40}}}}, PLOTLY_CONFIG);
      return;
    }}
    var d = sd[atlas][dataKey];
    var color = d.rho > 0.2 ? '#059669' : (d.rho < 0 ? '#DC2626' : '#6B7280');
    var traces = [{{
      type: 'scatter', mode: 'markers',
      x: d.x, y: d.y,
      marker: {{color: color, size: 6, opacity: 0.5}},
      hovertemplate: 'Expr: %{{x:.2f}}<br>Activity: %{{y:.2f}}<extra></extra>',
    }}];
    // Regression line
    if (d.x.length > 2) {{
      var sorted = d.x.map(function(v,i) {{ return [v, d.y[i]]; }}).sort(function(a,b) {{ return a[0]-b[0]; }});
      var n = sorted.length, sx=0, sy=0, sxx=0, sxy=0;
      sorted.forEach(function(p) {{ sx+=p[0]; sy+=p[1]; sxx+=p[0]*p[0]; sxy+=p[0]*p[1]; }});
      var slope = (n*sxy - sx*sy) / (n*sxx - sx*sx);
      var intercept = (sy - slope*sx) / n;
      var x0 = sorted[0][0], x1 = sorted[sorted.length-1][0];
      traces.push({{type:'scatter', mode:'lines', x:[x0,x1], y:[slope*x0+intercept, slope*x1+intercept],
        line:{{color:'black', dash:'dash', width:2}}, showlegend:false}});
    }}
    var pstr = d.pval < 0.001 ? d.pval.toExponential(1) : d.pval.toFixed(4);
    Plotly.newPlot('scatter-chart', traces, {{
      title: target + ' — ' + atlas + ' (' + sigLabel + ')  (\\u03c1=' + d.rho.toFixed(3) + ', p=' + pstr + ', n=' + d.n + ')',
      xaxis: {{title: 'Mean Expression'}}, yaxis: {{title: 'Predicted Activity'}},
      showlegend: false, margin: {{t:60}},
    }}, PLOTLY_CONFIG);
  }};
  updateScatter();
}})();

// ═══════════════════════════════════════════════════════════════════════════
// HEATMAP (tabbed Plotly)
// ═══════════════════════════════════════════════════════════════════════════
var currentHeatmapTab = 'cytosig';
window.switchHeatmapTab = function(sig) {{
  currentHeatmapTab = sig;
  document.querySelectorAll('#heatmap-tabs .tab-btn').forEach(function(b){{b.classList.remove('active');}});
  document.querySelector('#heatmap-tabs .tab-btn.' + sig).classList.add('active');
  renderHeatmap(sig);
}};

function renderHeatmap(sig) {{
  var hd = DATA.heatmap[sig];
  if (!hd || hd.targets.length === 0) {{
    Plotly.newPlot('heatmap-chart', [], {{title: 'No data for ' + SIG_NAMES[sig]}}, PLOTLY_CONFIG);
    return;
  }}
  var h = Math.max(500, hd.targets.length * 18 + 100);
  document.getElementById('heatmap-chart').style.height = h + 'px';
  // Replace nulls with NaN for plotly
  var z = hd.matrix.map(function(row) {{
    return row.map(function(v) {{ return v === null ? NaN : v; }});
  }});
  var hovertext = hd.targets.map(function(t, i) {{
    return hd.atlases.map(function(a, j) {{
      var v = hd.matrix[i][j];
      return t + '<br>' + a + '<br>\\u03c1 = ' + (v !== null ? v.toFixed(3) : 'N/A');
    }});
  }});
  // Red-white-blue: high=red, middle=white, low=blue
  var rwb = [[0,'#2563EB'],[0.25,'#93C5FD'],[0.5,'#FFFFFF'],[0.75,'#FCA5A5'],[1,'#DC2626']];
  Plotly.newPlot('heatmap-chart', [{{
    type: 'heatmap', z: z, x: hd.atlases, y: hd.targets,
    colorscale: rwb, zmin: -0.5, zmax: 0.8,
    hovertext: hovertext, hoverinfo: 'text',
    colorbar: {{title: 'Spearman \\u03c1', len: 0.5}},
  }}], {{
    title: SIG_NAMES[sig] + ' — Biologically Important Targets',
    xaxis: {{side: 'bottom', tickangle: -30}},
    yaxis: {{autorange: 'reversed', tickfont: {{size: 9}}}},
    margin: {{t:50, l:160, r:60, b:100}},
  }}, PLOTLY_CONFIG);
}}
renderHeatmap('cytosig');

// ═══════════════════════════════════════════════════════════════════════════
// COMPREHENSIVE VALIDATION (CytoSig + SecAct, all datasets)
// ═══════════════════════════════════════════════════════════════════════════
window.updateBulk = function() {{
  var dataset = document.getElementById('bulk-dataset-select').value;
  var sig = document.getElementById('bulk-sig-select').value;
  var bd = DATA.bulk;
  if (!bd[dataset] || !bd[dataset][sig]) {{
    Plotly.newPlot('bulk-chart', [], {{title: 'No data for ' + dataset + ' / ' + SIG_NAMES[sig], margin:{{t:40}}}}, PLOTLY_CONFIG);
    return;
  }}
  var d = bd[dataset][sig];
  var n = d.targets.length;
  var colors;
  if (sig === 'secact' && d.matched) {{
    // Color matched targets with SecAct green, additional targets with gray
    colors = d.matched.map(function(m, i) {{
      return m ? '#059669' : '#94A3B8';
    }});
  }} else {{
    var sigColor = DATA.sigColors[sig] || '#6B7280';
    colors = d.rhos.map(function(r) {{
      return r > 0.2 ? sigColor : (r < -0.1 ? '#DC2626' : '#6B7280');
    }});
  }}
  var totalNote = (sig === 'secact' && d.total_secact) ?
    '  [showing ' + n + ' of ' + d.total_secact + ' total SecAct targets: matched CytoSig targets + top additional]' : '';
  var traces = [{{
    type: 'bar', x: d.targets, y: d.rhos,
    marker: {{color: colors}},
    hovertemplate: '<b>%{{x}}</b><br>\\u03c1 = %{{y:.3f}}<extra></extra>',
  }}];
  var title = dataset + ' — ' + SIG_NAMES[sig] + '  (n=' + n + ', median \\u03c1=' + d.median + ')';
  var layout = {{
    title: title,
    xaxis: {{title: 'Target (ranked by \\u03c1)', tickangle: -90, tickfont: {{size: 7}}}},
    yaxis: {{title: 'Spearman \\u03c1'}},
    showlegend: false, margin: {{t:60, b:120}},
  }};
  if (sig === 'secact' && d.matched) {{
    layout.annotations = [{{
      text: 'Green = matched CytoSig target, Gray = additional SecAct target',
      xref:'paper', yref:'paper', x:0.5, y:1.06, showarrow:false,
      font:{{size:11, color:'#6B7280'}},
    }}];
  }}
  Plotly.newPlot('bulk-chart', traces, layout, PLOTLY_CONFIG);
}};
updateBulk();

// ═══════════════════════════════════════════════════════════════════════════
// 5.1 METHOD COMPARISON BOXPLOT (10-way with dropdown for per-atlas view)
// ═══════════════════════════════════════════════════════════════════════════
(function() {{
  var mb = DATA.methodBoxplot;
  if (!mb) return;
  var allAtlases = Object.keys(mb);
  var sigConfigs = [
    {{key:'cytosig',                name:'CytoSig',                      color:'#2563EB'}},
    {{key:'lincyto_orig',           name:'LinCytoSig (orig)',             color:'#D97706'}},
    {{key:'lincyto_filt',           name:'LinCytoSig (gene-filtered)',    color:'#F59E0B'}},
    {{key:'lincyto_best_orig',      name:'LinCytoSig Best (combined)',    color:'#B45309'}},
    {{key:'lincyto_best_filt',      name:'LinCytoSig Best (comb+filt)',   color:'#92400E'}},
    {{key:'lincyto_best_gtex',      name:'LinCytoSig Best (GTEx)',        color:'#7C3AED'}},
    {{key:'lincyto_best_tcga',      name:'LinCytoSig Best (TCGA)',        color:'#EC4899'}},
    {{key:'lincyto_best_gtex_filt', name:'LinCytoSig Best (GTEx+filt)',   color:'#A78BFA'}},
    {{key:'lincyto_best_tcga_filt', name:'LinCytoSig Best (TCGA+filt)',   color:'#F9A8D4'}},
    {{key:'secact',                 name:'SecAct',                        color:'#059669'}},
  ];

  // Populate dropdown
  var sel = document.getElementById('method-boxplot-view');
  allAtlases.forEach(function(a) {{
    var opt = document.createElement('option');
    opt.value = a; opt.textContent = a;
    sel.appendChild(opt);
  }});

  window.updateMethodBoxplot = function() {{
    var view = sel.value;
    var atlasesToShow = (view === 'all') ? allAtlases : [view];
    var isSingle = (view !== 'all');
    var traces = [];

    sigConfigs.forEach(function(cfg) {{
      var y = [], x = [];
      atlasesToShow.forEach(function(a) {{
        if (mb[a] && mb[a][cfg.key] && mb[a][cfg.key].length > 0) {{
          mb[a][cfg.key].forEach(function(v) {{
            y.push(v);
            x.push(isSingle ? cfg.name : a);
          }});
        }}
      }});
      if (y.length > 0) {{
        var n = mb[atlasesToShow[0]][cfg.key] ? mb[atlasesToShow[0]][cfg.key].length : 0;
        traces.push({{
          type: 'box', y: y, x: x,
          name: cfg.name + (isSingle ? '' : ' (n=' + n + ')'),
          marker: {{color: cfg.color}},
          boxpoints: isSingle ? 'all' : false,
          jitter: 0.3, pointpos: 0,
          line: {{width: 1.5}},
          showlegend: !isSingle,
        }});
      }}
    }});

    var title = isSingle
      ? view + ' \u2014 10-Way Method Comparison (n=' + (mb[view].cytosig ? mb[view].cytosig.length : 0) + ' matched pairs)'
      : '10-Way Signature Method Comparison (Expression\u2013Activity Correlation)';

    Plotly.newPlot('method-boxplot-chart', traces, {{
      title: title,
      yaxis: {{title: 'Spearman \\u03c1', zeroline: true, zerolinecolor: '#ccc'}},
      boxmode: 'group',
      legend: {{orientation:'h', y:1.25, x:0.5, xanchor:'center', font:{{size:9}}}},
      margin: {{t:isSingle ? 80 : 140, b:isSingle ? 120 : 80}},
    }}, PLOTLY_CONFIG);
  }};
  updateMethodBoxplot();
}})();

// ═══════════════════════════════════════════════════════════════════════════
// 5.2 EFFECT OF AGGREGATION LEVEL
// ═══════════════════════════════════════════════════════════════════════════
(function() {{
  var lc = DATA.levelComp;
  if (!lc) return;

  var atlasNames = Object.keys(lc);
  var sel = document.getElementById('level-comp-atlas');
  atlasNames.forEach(function(a) {{
    var opt = document.createElement('option');
    opt.value = a; opt.textContent = a;
    sel.appendChild(opt);
  }});

  var methodCfg = [
    {{key: 'cytosig', name: 'CytoSig', color: '#2563EB'}},
    {{key: 'lincytosig', name: 'LinCytoSig', color: '#D97706'}},
    {{key: 'secact', name: 'SecAct', color: '#059669'}},
  ];

  window.updateLevelComp = function() {{
    var atlas = sel.value;
    var levels = lc[atlas];
    if (!levels || levels.length === 0) {{
      Plotly.newPlot('level-comp-chart', [], {{title: 'No data for ' + atlas}}, PLOTLY_CONFIG);
      return;
    }}
    var traces = [];
    methodCfg.forEach(function(cfg) {{
      var y = [], x = [];
      levels.forEach(function(lv) {{
        var arr = lv[cfg.key] || [];
        arr.forEach(function(v) {{
          y.push(v);
          x.push(lv.level + ' (n=' + lv.n + ')');
        }});
      }});
      if (y.length > 0) {{
        traces.push({{
          type: 'box', y: y, x: x,
          name: cfg.name,
          marker: {{color: cfg.color}},
          boxpoints: false,
          line: {{width: 1.5}},
        }});
      }}
    }});
    Plotly.newPlot('level-comp-chart', traces, {{
      title: atlas + ' \u2014 CytoSig vs LinCytoSig vs SecAct by Aggregation Level',
      yaxis: {{title: 'Spearman \\u03c1', zeroline: true, zerolinecolor: '#ccc'}},
      boxmode: 'group',
      legend: {{orientation: 'h', y: 1.15, x: 0.5, xanchor: 'center'}},
      margin: {{t: 80, b: 80}},
    }}, PLOTLY_CONFIG);
  }};
  updateLevelComp();

  // 5.2.2 Summary table (all atlases × all levels)
  var html = '<table><tr><th>Atlas</th><th>Level</th><th>n (3-way)</th>';
  html += '<th>CytoSig med &rho;</th><th>LinCytoSig med &rho;</th><th>SecAct med &rho;</th>';
  html += '<th>vs CytoSig (W/L)</th><th>&Delta;&rho;</th>';
  html += '<th>vs SecAct (W/L)</th><th>&Delta;&rho;</th></tr>';
  atlasNames.forEach(function(atlas) {{
    var levels = lc[atlas];
    levels.forEach(function(lv, i) {{
      var cBg = lv.vs_cyto_mean > 0 ? '#D1FAE5' : '#FEE2E2';
      var sBg = lv.vs_sec_mean !== undefined ? (lv.vs_sec_mean > 0 ? '#D1FAE5' : '#FEE2E2') : '';
      html += '<tr>';
      if (i === 0) html += '<td rowspan="' + levels.length + '"><strong>' + atlas + '</strong></td>';
      html += '<td>' + lv.level + '</td>';
      html += '<td>' + lv.n + '</td>';
      html += '<td>' + lv.cyto_median + '</td>';
      html += '<td>' + lv.lin_median + '</td>';
      html += '<td>' + (lv.sec_median !== undefined ? lv.sec_median : '\u2014') + '</td>';
      html += '<td style="background:' + cBg + '">' + lv.vs_cyto_win + '/' + lv.vs_cyto_loss + '</td>';
      html += '<td style="background:' + cBg + '">' + (lv.vs_cyto_mean > 0 ? '+' : '') + lv.vs_cyto_mean + '</td>';
      if (lv.vs_sec_win !== undefined) {{
        html += '<td style="background:' + sBg + '">' + lv.vs_sec_win + '/' + lv.vs_sec_loss + '</td>';
        html += '<td style="background:' + sBg + '">' + (lv.vs_sec_mean > 0 ? '+' : '') + lv.vs_sec_mean + '</td>';
      }} else {{
        html += '<td colspan="2">\u2014</td>';
      }}
      html += '</tr>';
    }});
  }});
  html += '</table>';
  document.getElementById('level-summary-container').innerHTML = html;
}})();

// 5.2.3-5.2.4 CELLTYPE & CYTOKINE DETAIL TABLES (from celltypeComp)
(function() {{
  var cc = DATA.celltypeComp;
  if (!cc) return;

  // 5.2.3 Celltype summary table
  var ct = cc.celltypeTable || [];
  var ctHtml = '<table><tr><th>Cell Type</th><th colspan="3">vs CytoSig</th><th colspan="3">vs SecAct</th></tr>';
  ctHtml += '<tr><th></th><th>Win/Loss</th><th>Mean &Delta;&rho;</th><th>n</th><th>Win/Loss</th><th>Mean &Delta;&rho;</th><th>n</th></tr>';
  ct.forEach(function(r) {{
    var cBg = r.vs_cyto_mean > 0 ? '#D1FAE5' : (r.vs_cyto_mean < -0.15 ? '#FEE2E2' : '');
    var sBg = r.vs_sec_mean !== null ? (r.vs_sec_mean > 0 ? '#D1FAE5' : (r.vs_sec_mean < -0.15 ? '#FEE2E2' : '')) : '';
    ctHtml += '<tr><td><strong>' + r.celltype.replace(/_/g, ' ') + '</strong></td>';
    ctHtml += '<td style="background:' + cBg + '">' + r.vs_cyto_win + '/' + r.vs_cyto_loss + '</td>';
    ctHtml += '<td style="background:' + cBg + '">' + (r.vs_cyto_mean > 0 ? '+' : '') + r.vs_cyto_mean + '</td>';
    ctHtml += '<td>' + r.vs_cyto_n + '</td>';
    if (r.vs_sec_n > 0) {{
      ctHtml += '<td style="background:' + sBg + '">' + r.vs_sec_win + '/' + r.vs_sec_loss + '</td>';
      ctHtml += '<td style="background:' + sBg + '">' + (r.vs_sec_mean > 0 ? '+' : '') + r.vs_sec_mean + '</td>';
      ctHtml += '<td>' + r.vs_sec_n + '</td>';
    }} else {{
      ctHtml += '<td colspan="3" style="color:#9CA3AF">No matched SecAct target</td>';
    }}
    ctHtml += '</tr>';
  }});
  ctHtml += '</table>';
  document.getElementById('celltype-table-container').innerHTML = ctHtml;

  // 5.2.4 Cytokine summary table
  var cyto = cc.cytokineTable || [];
  var cyHtml = '<table><tr><th>Cytokine</th><th colspan="3">vs CytoSig</th><th colspan="3">vs SecAct</th></tr>';
  cyHtml += '<tr><th></th><th>Win/Loss</th><th>Mean &Delta;&rho;</th><th>n</th><th>Win/Loss</th><th>Mean &Delta;&rho;</th><th>n</th></tr>';
  cyto.forEach(function(r) {{
    var cBg = r.vs_cyto_mean > 0 ? '#D1FAE5' : (r.vs_cyto_mean < -0.15 ? '#FEE2E2' : '');
    var sBg = r.vs_sec_mean !== null ? (r.vs_sec_mean > 0 ? '#D1FAE5' : (r.vs_sec_mean < -0.15 ? '#FEE2E2' : '')) : '';
    cyHtml += '<tr><td><strong>' + r.cytokine + '</strong></td>';
    cyHtml += '<td style="background:' + cBg + '">' + r.vs_cyto_win + '/' + r.vs_cyto_loss + '</td>';
    cyHtml += '<td style="background:' + cBg + '">' + (r.vs_cyto_mean > 0 ? '+' : '') + r.vs_cyto_mean + '</td>';
    cyHtml += '<td>' + r.vs_cyto_n + '</td>';
    if (r.vs_sec_n > 0) {{
      cyHtml += '<td style="background:' + sBg + '">' + r.vs_sec_win + '/' + r.vs_sec_loss + '</td>';
      cyHtml += '<td style="background:' + sBg + '">' + (r.vs_sec_mean > 0 ? '+' : '') + r.vs_sec_mean + '</td>';
      cyHtml += '<td>' + r.vs_sec_n + '</td>';
    }} else {{
      cyHtml += '<td colspan="3" style="color:#9CA3AF">No matched SecAct target</td>';
    }}
    cyHtml += '</tr>';
  }});
  cyHtml += '</table>';
  document.getElementById('cytokine-table-container').innerHTML = cyHtml;
}})();

// ═══════════════════════════════════════════════════════════════════════════
// 5.3 CELLTYPE Δρ BAR CHART (Figure 12)
// ═══════════════════════════════════════════════════════════════════════════
(function() {{
  var cc = DATA.celltypeComp;
  if (!cc || !cc.celltypeTable) return;

  // Sort by mean Δρ descending
  var ct = cc.celltypeTable.slice().sort(function(a, b) {{ return a.vs_cyto_mean - b.vs_cyto_mean; }});

  var names = ct.map(function(r) {{ return r.celltype.replace(/_/g, ' '); }});
  var means = ct.map(function(r) {{ return r.vs_cyto_mean; }});
  var sems = ct.map(function(r) {{ return r.vs_cyto_sem || 0; }});
  var colors = means.map(function(v) {{ return v > 0 ? '#D97706' : '#2563EB'; }});
  var hoverTexts = ct.map(function(r) {{
    return '<b>' + r.celltype.replace(/_/g, ' ') + '</b><br>' +
      'Mean \u0394\u03c1: ' + (r.vs_cyto_mean > 0 ? '+' : '') + r.vs_cyto_mean +
      ' \u00b1 ' + (r.vs_cyto_sem || 0) + '<br>' +
      'Win/Loss: ' + r.vs_cyto_win + '/' + r.vs_cyto_loss +
      ' (n=' + r.vs_cyto_n + ')';
  }});

  var nPos = means.filter(function(v) {{ return v > 0; }}).length;

  Plotly.newPlot('celltype-delta-rho-chart', [{{
    type: 'bar',
    y: names,
    x: means,
    orientation: 'h',
    marker: {{color: colors, opacity: 0.85}},
    error_x: {{type: 'data', array: sems, visible: true, thickness: 1, width: 3}},
    text: hoverTexts,
    hovertemplate: '%{{text}}<extra></extra>',
  }}], {{
    title: 'Celltype-Specific \u0394\u03c1: LinCytoSig vs CytoSig<br><span style="font-size:12px;color:#6B7280">Aggregated across 4 atlases at donor \u00d7 celltype level</span>',
    xaxis: {{title: 'Mean \u0394\u03c1 (LinCytoSig \u2212 CytoSig) \u00b1 SEM', zeroline: true, zerolinewidth: 2, zerolinecolor: '#000'}},
    yaxis: {{automargin: true, tickfont: {{size: 10}}}},
    margin: {{l: 180, t: 80, b: 60, r: 30}},
    annotations: [{{
      text: nPos + '/' + ct.length + ' cell types favor LinCytoSig<br>Orange = LinCytoSig better, Blue = CytoSig better',
      xref: 'paper', yref: 'paper', x: 0.98, y: 0.02,
      showarrow: false, font: {{size: 11}}, align: 'right',
      bgcolor: '#FEF3C7', borderpad: 4,
    }}],
  }}, PLOTLY_CONFIG);
}})();

// ═══════════════════════════════════════════════════════════════════════════
// LIGHTBOX for static images
// ═══════════════════════════════════════════════════════════════════════════
document.querySelectorAll('.figure img').forEach(function(img) {{
  img.addEventListener('click', function() {{
    var lb = document.getElementById('lightbox');
    lb.querySelector('img').src = img.src;
    lb.classList.add('active');
    document.body.style.overflow = 'hidden';
  }});
}});
function closeLightbox() {{
  var lb = document.getElementById('lightbox');
  lb.classList.remove('active');
  document.body.style.overflow = '';
}}
document.getElementById('lightbox').addEventListener('click', function(e) {{
  if (e.target === this) closeLightbox();
}});
document.addEventListener('keydown', function(e) {{
  if (e.key === 'Escape') closeLightbox();
}});

</script>

</body>
</html>"""

    return html


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    print('CytoAtlas PI Report — Interactive HTML Generation')
    print('=' * 60)

    print('\nLoading correlation data...')
    df = load_all_correlations()
    print(f'  Loaded {len(df)} correlation records')
    # Use inflammation_main directly (no merge of val/ext)
    print(f'  Using inflammation_main directly')

    print('\nPreparing data for interactive figures...')
    summary_table = prepare_summary_table(df)
    print(f'  Summary table: {len(summary_table)} rows')

    boxplot_data = prepare_boxplot_data(df)
    print(f'  Boxplot data: {len(boxplot_data)} atlases')

    stratified_data = prepare_stratified_data(df)
    print(f'  Stratified data: {len(stratified_data)} datasets')

    consistency_data = prepare_consistency_data(df)
    print(f'  Consistency data: {len(consistency_data)} targets')

    heatmap_data = prepare_heatmap_data(df)
    for sig, d in heatmap_data.items():
        print(f'  Heatmap {sig}: {len(d["targets"])} targets')

    levels_data = prepare_levels_data(df)
    print(f'  Levels data: {len(levels_data)} atlases x 2 signatures')

    bulk_data = prepare_bulk_validation_data(df)
    print(f'  Bulk validation: {len(bulk_data)} datasets')

    scatter_data = prepare_scatter_data()
    print(f'  Scatter data: {len(scatter_data)} atlases')

    good_bad_data = prepare_good_bad_data(df)
    print(f'  Good/bad data: {len(good_bad_data)} signature types')

    method_boxplot_data = prepare_method_comparison_boxplot(df)
    print(f'  Method comparison boxplot: {len(method_boxplot_data)} atlases')

    celltype_comparison_data = prepare_celltype_comparison(df)
    print(f'  Celltype comparison: {len(celltype_comparison_data.get("celltypeTable", []))} cell types')

    level_comparison_data = prepare_level_comparison_data()
    print(f'  Level comparison: {len(level_comparison_data)} atlases')

    print('\nGenerating HTML...')
    html = generate_html(summary_table, boxplot_data, consistency_data,
                         heatmap_data, levels_data, bulk_data, scatter_data,
                         good_bad_data, method_boxplot_data, celltype_comparison_data,
                         level_comparison_data, stratified_data)

    output_path = REPORT_DIR / 'REPORT.html'
    output_path.write_text(html, encoding='utf-8')
    size_kb = len(html) / 1024
    print(f'  Written to: {output_path}')
    print(f'  Size: {size_kb:.0f} KB')

    print('\n' + '=' * 60)
    print('Done!')


if __name__ == '__main__':
    main()
