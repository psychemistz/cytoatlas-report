#!/usr/bin/env python3
"""
16_residual_correlation.py — Residual Correlation After Cell-Fraction Adjustment
================================================================================

For each single-cell dataset with matched cell-fraction data (CIMA, Inflammation
Main, scAtlas Cancer), computes:

  1. Direct Spearman correlation between target gene expression (X) and
     predicted activity (A) across donors.
  2. Residual correlation: fits X = b0 + b1·F + b2·(A×F) + ε (where F =
     cell-type fractions, A×F = interaction terms), then computes
     Spearman(ε, A). This tests whether predicted activity tracks expression
     AFTER removing the effect of cell-type composition.

Note: scAtlas Normal is excluded because donor IDs in the donor_only pseudobulk
(284C, 290B, ...) do not match the celltype1 pseudobulk (01, 02, ...) — the two
files use incompatible identifier systems.

Output: /data/parks34/projects/2cytoatlas/visualization/data/validation/residual_correlation.json
"""

import json
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression

warnings.filterwarnings('ignore')

# ─── Paths ────────────────────────────────────────────────────────────────────
BASE = Path('/data/parks34/projects/2cytoatlas')
RESULTS = BASE / 'results' / 'cross_sample_validation'
VIZ_DIR = BASE / 'visualization' / 'data' / 'validation'

# ─── Dataset configs ──────────────────────────────────────────────────────────
DATASETS = {
    'cima': {
        'expression': RESULTS / 'cima' / 'cima_donor_only_pseudobulk.h5ad',
        'activity_cytosig': RESULTS / 'cima' / 'cima_donor_only_cytosig.h5ad',
        'activity_secact': RESULTS / 'cima' / 'cima_donor_only_secact.h5ad',
        'celltype_pb': RESULTS / 'cima' / 'cima_donor_l1_pseudobulk.h5ad',
        'donor_col': 'donor',
        'ct_col': 'cell_type_l1',
        'label': 'CIMA',
        'gene_id': 'symbol',  # var index is already gene symbols
    },
    'inflammation_main': {
        'expression': RESULTS / 'inflammation_main' / 'inflammation_main_donor_only_pseudobulk.h5ad',
        'activity_cytosig': RESULTS / 'inflammation_main' / 'inflammation_main_donor_only_cytosig.h5ad',
        'activity_secact': RESULTS / 'inflammation_main' / 'inflammation_main_donor_only_secact.h5ad',
        'celltype_pb': RESULTS / 'inflammation_main' / 'inflammation_main_donor_l1_pseudobulk.h5ad',
        'donor_col': 'donor',
        'ct_col': 'Level1',
        'label': 'Inflammation Atlas Main',
        'gene_id': 'ensembl',  # var index is ENSG IDs; symbol in var['symbol']
    },
    # scAtlas Normal excluded: donor IDs in donor_only (284C, 290B, ...) don't match
    # celltype1 pseudobulk (01, 02, ...) — incompatible identifier systems.
    'scatlas_cancer': {
        'expression': RESULTS / 'scatlas_cancer' / 'scatlas_cancer_tumor_only_pseudobulk.h5ad',
        'activity_cytosig': RESULTS / 'scatlas_cancer' / 'scatlas_cancer_tumor_only_cytosig.h5ad',
        'activity_secact': RESULTS / 'scatlas_cancer' / 'scatlas_cancer_tumor_only_secact.h5ad',
        'celltype_pb': RESULTS / 'scatlas_cancer' / 'scatlas_cancer_tumor_by_cancer_celltype1_pseudobulk.h5ad',
        'donor_col': 'donor',
        'ct_col': 'cellType1',
        'label': 'scAtlas (Cancer)',
        'gene_id': 'symbol',
    },
}

# CytoSig alias map: CytoSig common name → gene symbol (same as report generator)
ALIAS_MAP = {
    'Activin A': 'INHBA', 'CD40L': 'CD40LG', 'GCSF': 'CSF3', 'GMCSF': 'CSF2',
    'IFNL': 'IFNL1', 'IL36': 'IL36A', 'MCSF': 'CSF1', 'TNFA': 'TNF',
    'TRAIL': 'TNFSF10', 'TWEAK': 'TNFSF12',
}

MAX_CELLTYPES = 15        # Cap on cell types (rare ones grouped into "Other")
MIN_DONOR_FRAC = 0.10     # Cell type must be present in ≥10% of donors


def compute_cell_fractions(celltype_pb_path, donor_col, ct_col, expression_donors):
    """Compute donor × cell-type fraction matrix from cell-type pseudobulk.

    Returns:
        fractions: DataFrame (donors × cell types), rows sum to 1
        celltypes_used: list of cell type names in the fraction matrix
    """
    adata = ad.read_h5ad(celltype_pb_path)
    df = adata.obs[[donor_col, ct_col, 'n_cells']].copy()
    df[donor_col] = df[donor_col].astype(str)
    df[ct_col] = df[ct_col].astype(str)

    # Pivot: donor × celltype → n_cells
    pivot = df.pivot_table(index=donor_col, columns=ct_col, values='n_cells',
                           aggfunc='sum', fill_value=0)

    # Keep only donors that appear in the expression data
    common_donors = sorted(set(pivot.index) & set(expression_donors))
    if len(common_donors) == 0:
        raise ValueError(f'No overlapping donors between expression ({len(expression_donors)}) '
                         f'and celltype pseudobulk ({len(pivot)})')
    pivot = pivot.loc[common_donors]

    # Filter rare cell types: present in <10% of donors
    presence = (pivot > 0).mean(axis=0)
    keep = presence[presence >= MIN_DONOR_FRAC].index.tolist()

    # If too many cell types, keep top-K by total cell count
    if len(keep) > MAX_CELLTYPES:
        totals = pivot[keep].sum(axis=0).sort_values(ascending=False)
        top_k = totals.index[:MAX_CELLTYPES - 1].tolist()
        other_cols = [c for c in keep if c not in top_k]
        pivot['Other'] = pivot[other_cols].sum(axis=1)
        keep = top_k + ['Other']

    pivot = pivot[keep]

    # Normalize to fractions (rows sum to 1)
    row_sums = pivot.sum(axis=1)
    # Drop donors with zero total cells
    valid = row_sums > 0
    pivot = pivot[valid]
    row_sums = row_sums[valid]
    fractions = pivot.div(row_sums, axis=0)

    # Drop the most abundant cell type (reference category) to avoid
    # multicollinearity with the intercept
    ref_ct = fractions.mean(axis=0).idxmax()
    celltypes_used = [c for c in fractions.columns if c != ref_ct]
    fractions = fractions[celltypes_used]

    return fractions, celltypes_used


def build_gene_index(expr, gene_id_type):
    """Build a symbol → positional index mapping.

    For Ensembl-indexed files (Inflammation Main), maps via the 'symbol' column.
    For symbol-indexed files, maps directly from var index.
    Returns dict: gene_symbol → column_index.
    """
    if gene_id_type == 'ensembl':
        # Use the 'symbol' column to map gene symbols to positional indices
        symbols = expr.var['symbol'].astype(str)
        gene_index = {}
        for i, sym in enumerate(symbols):
            if sym and sym != 'nan':
                gene_index[sym] = i
        return gene_index
    else:
        return {g: i for i, g in enumerate(expr.var.index)}


def process_signature(expression_path, activity_path, fractions, sig_name, gene_id_type):
    """Compute direct and residual correlations for one dataset × signature.

    Returns dict with targets, direct_rho, residual_rho, etc.
    """
    expr = ad.read_h5ad(expression_path)
    act = ad.read_h5ad(activity_path)

    # Align donors across expression, activity, and fractions
    expr_donors = [str(d) for d in expr.obs.index]
    act_donors = [str(d) for d in act.obs.index]
    frac_donors = list(fractions.index)
    common = sorted(set(expr_donors) & set(act_donors) & set(frac_donors))

    if len(common) < 20:
        print(f'    WARNING: Only {len(common)} common donors for {sig_name}, skipping')
        return None

    # Reindex all to common donors
    expr_idx = [expr_donors.index(d) for d in common]
    act_idx = [act_donors.index(d) for d in common]

    X_mat = expr.X[expr_idx]
    if hasattr(X_mat, 'toarray'):
        X_mat = X_mat.toarray()
    X_mat = np.asarray(X_mat, dtype=np.float64)

    A_mat = act.X[act_idx]
    if hasattr(A_mat, 'toarray'):
        A_mat = A_mat.toarray()
    A_mat = np.asarray(A_mat, dtype=np.float64)

    F_mat = fractions.loc[common].values.astype(np.float64)  # (n_donors, n_celltypes)

    act_targets = list(act.var.index)

    # Build gene symbol → column index mapping
    gene_index = build_gene_index(expr, gene_id_type)

    targets_out = []
    direct_rhos = []
    residual_rhos = []
    direct_pvals = []
    residual_pvals = []

    for t_idx, target in enumerate(act_targets):
        # Resolve gene name: target → gene symbol → expression column index
        gene_name = None
        if target in gene_index:
            gene_name = target
        elif target in ALIAS_MAP and ALIAS_MAP[target] in gene_index:
            gene_name = ALIAS_MAP[target]
        else:
            # Reverse alias check (SecAct gene → CytoSig alias)
            reverse = {v: k for k, v in ALIAS_MAP.items()}
            if target in reverse and reverse[target] in gene_index:
                gene_name = reverse[target]

        if gene_name is None:
            continue

        g_idx = gene_index[gene_name]
        x = X_mat[:, g_idx]  # expression vector
        a = A_mat[:, t_idx]  # activity vector

        # Skip if NaN present or no variance
        if np.any(np.isnan(x)) or np.any(np.isnan(a)):
            continue
        if np.std(x) < 1e-10 or np.std(a) < 1e-10:
            continue

        # Direct correlation
        rho_direct, p_direct = spearmanr(x, a)

        # Build design matrix: [F, a⊙F] (fractions + interaction terms)
        aF = a[:, None] * F_mat  # (n_donors, n_celltypes)
        design = np.hstack([F_mat, aF])

        # Check for NaN/Inf in design matrix
        if np.any(~np.isfinite(design)):
            continue

        # Fit OLS: x ~ design (with intercept)
        model = LinearRegression(fit_intercept=True)
        model.fit(design, x)
        residuals = x - model.predict(design)

        # Residual correlation
        if np.std(residuals) < 1e-10:
            continue
        rho_resid, p_resid = spearmanr(residuals, a)

        targets_out.append(target)
        direct_rhos.append(round(float(rho_direct), 4))
        residual_rhos.append(round(float(rho_resid), 4))
        direct_pvals.append(float(p_direct))
        residual_pvals.append(float(p_resid))

    if not targets_out:
        return None

    # Summary statistics
    dr = np.array(direct_rhos)
    rr = np.array(residual_rhos)
    deltas = rr - dr
    n_pos_direct = np.sum(dr > 0)

    return {
        'targets': targets_out,
        'direct_rho': direct_rhos,
        'residual_rho': residual_rhos,
        'direct_pval': direct_pvals,
        'residual_pval': residual_pvals,
        'n_donors': len(common),
        'n_celltypes': F_mat.shape[1],
        'celltypes_used': list(fractions.columns),
        'model_params': 1 + 2 * F_mat.shape[1],  # intercept + F + a*F
        'stats': {
            'median_direct': round(float(np.median(dr)), 4),
            'median_residual': round(float(np.median(rr)), 4),
            'median_delta': round(float(np.median(deltas)), 4),
            'pct_retained_positive': round(float(
                100 * np.sum((dr > 0) & (rr > 0)) / n_pos_direct
            ), 1) if n_pos_direct > 0 else 0.0,
            'pct_sign_preserved': round(float(
                100 * np.mean(np.sign(dr) == np.sign(rr))
            ), 1),
        },
    }


def main():
    print('Residual Correlation — Cell-Fraction Adjustment')
    print('=' * 60)

    results = {}

    for ds_key, cfg in DATASETS.items():
        label = cfg['label']
        print(f'\n{label}')
        print('-' * 40)

        # Load expression donors for alignment
        expr = ad.read_h5ad(cfg['expression'], backed='r')
        expr_donors = [str(d) for d in expr.obs.index]

        # Compute cell-type fractions
        print(f'  Computing cell-type fractions from {cfg["celltype_pb"].name}...')
        fractions, ct_used = compute_cell_fractions(
            cfg['celltype_pb'], cfg['donor_col'], cfg['ct_col'], expr_donors)
        print(f'  {fractions.shape[0]} donors × {fractions.shape[1]} cell types')
        print(f'  Cell types: {ct_used}')

        results[label] = {}

        for sig_name, sig_key in [('CytoSig', 'cytosig'), ('SecAct', 'secact')]:
            act_path = cfg[f'activity_{sig_key}']
            print(f'  {sig_name}: processing...')
            result = process_signature(
                cfg['expression'], act_path, fractions, sig_name, cfg['gene_id'])
            if result is None:
                print(f'    No results for {sig_name}')
                continue
            results[label][sig_key] = result
            s = result['stats']
            print(f'    {len(result["targets"])} targets, '
                  f'median direct ρ={s["median_direct"]:.3f}, '
                  f'residual ρ={s["median_residual"]:.3f}, '
                  f'Δ={s["median_delta"]:.3f}, '
                  f'{s["pct_sign_preserved"]:.0f}% sign preserved')

    # Write output
    out_path = VIZ_DIR / 'residual_correlation.json'
    with open(out_path, 'w') as f:
        json.dump(results, f, separators=(',', ':'))
    print(f'\nWritten to {out_path}')
    print(f'Size: {out_path.stat().st_size / 1024:.1f} KB')


if __name__ == '__main__':
    main()
