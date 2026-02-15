#!/usr/bin/env python3
"""
CytoAtlas PI Report — Static Figure Generation
================================================
Generates all figures for the comprehensive report to PI (Peng Jiang).

Output: report/figures/
"""

import json
import warnings
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings('ignore')

# ─── Paths ───────────────────────────────────────────────────────────────────
BASE = Path('/data/parks34/projects/2cytoatlas')
CORR_DIR = BASE / 'results' / 'cross_sample_validation' / 'correlations'
VIZ_DIR = BASE / 'visualization' / 'data' / 'validation'
LINCYTO_FILT_DIR = BASE / 'results' / 'lincytosig_gene_filter'
BEST_SELECTION = VIZ_DIR / 'best_lincytosig_selection.json'
SCATTER_DONOR = VIZ_DIR / 'donor_scatter'
SCATTER_CT = VIZ_DIR / 'celltype_scatter'
SCATTER_RESAMP = VIZ_DIR / 'resampled_scatter'
CELLTYPE_SIG = BASE / 'results' / 'celltype_signatures'
FIG_DIR = BASE / 'report' / 'figures'
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ─── Style ───────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.facecolor': 'white',
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# Color palette
COLORS = {
    'cytosig': '#2563EB',      # blue
    'lincyto_best': '#92400E', # dark brown
    'lincytosig': '#D97706',   # amber
    'secact': '#059669',       # emerald
    'good': '#059669',         # green
    'bad': '#DC2626',          # red
    'neutral': '#6B7280',      # gray
}

ATLAS_COLORS = {
    'cima': '#3B82F6',
    'inflammation': '#EF4444',
    'scatlas_normal': '#10B981',
    'scatlas_cancer': '#8B5CF6',
    'gtex': '#6366F1',
    'tcga': '#EC4899',
}

# Biologically important target families
BIO_FAMILIES = {
    'Interferon': ['IFNG', 'IFN1', 'IFNL'],
    'TGF-β': ['TGFB1', 'TGFB2', 'TGFB3', 'BMP2', 'BMP4', 'BMP6', 'GDF11'],
    'Interleukin': ['IL1A', 'IL1B', 'IL2', 'IL4', 'IL6', 'IL10', 'IL12', 'IL13',
                    'IL15', 'IL17A', 'IL21', 'IL22', 'IL27', 'IL33'],
    'TNF': ['TNFA', 'LTA', 'TRAIL', 'TWEAK', 'CD40L'],
    'Growth Factor': ['EGF', 'FGF2', 'HGF', 'VEGFA', 'PDGFB', 'IGF1'],
    'Chemokine': ['CXCL12'],
    'Colony-Stimulating': ['GMCSF', 'GCSF', 'MCSF'],
}

# Flatten for lookup
TARGET_TO_FAMILY = {}
for fam, targets in BIO_FAMILIES.items():
    for t in targets:
        TARGET_TO_FAMILY[t] = fam

# CytoSig → SecAct alias map (common name → HGNC gene symbol)
ALIAS_MAP = {
    'Activin A': 'INHBA', 'CD40L': 'CD40LG', 'GCSF': 'CSF3', 'GMCSF': 'CSF2',
    'IFNL': 'IFNL1', 'IL36': 'IL36A', 'MCSF': 'CSF1', 'TNFA': 'TNF',
    'TRAIL': 'TNFSF10', 'TWEAK': 'TNFSF12',
}
DIRECT_MATCHED = [
    'BDNF', 'BMP2', 'BMP4', 'BMP6', 'CXCL12', 'FGF2', 'GDF11', 'HGF',
    'IFNG', 'IL10', 'IL15', 'IL1A', 'IL1B', 'IL21', 'IL27', 'IL6',
    'LIF', 'LTA', 'OSM', 'TGFB1', 'TGFB3', 'VEGFA',
]
ALL_32_CS = DIRECT_MATCHED + list(ALIAS_MAP.keys())
ALL_32_SA = DIRECT_MATCHED + list(ALIAS_MAP.values())

# Independent levels per dataset
INDEPENDENT_LEVELS = {
    'gtex': ('by_tissue', True),        # (level, needs_median_of_medians)
    'tcga': ('primary_by_cancer', True),
    'cima': ('donor_only', False),
    'inflammation_main': ('donor_only', False),
    'scatlas_normal': ('donor_only', False),
    'scatlas_cancer': ('tumor_only', False),
}

FAMILY_COLORS = {
    'Interferon': '#DC2626',
    'TGF-β': '#2563EB',
    'Interleukin': '#059669',
    'TNF': '#D97706',
    'Growth Factor': '#8B5CF6',
    'Chemokine': '#EC4899',
    'Colony-Stimulating': '#6366F1',
    'Other': '#9CA3AF',
}


def load_all_correlations():
    """Load and merge all per-atlas + bulk correlations."""
    dfs = []
    for f in CORR_DIR.glob('*_correlations.csv'):
        if f.name == 'all_correlations.csv' or 'resampled' in f.name or 'summary' in f.name:
            continue
        df = pd.read_csv(f, low_memory=False)
        dfs.append(df)
    # Also load bulk
    bulk = pd.read_csv(CORR_DIR / 'all_correlations.csv', low_memory=False)
    dfs.append(bulk)
    merged = pd.concat(dfs, ignore_index=True)
    merged = merged.drop_duplicates(subset=['target', 'celltype', 'atlas', 'level', 'signature'])
    return merged


def merge_inflammation_atlases(df):
    """Replace inflammation_main/val/ext with single 'inflammation' atlas."""
    mask = df['atlas'].isin(['inflammation_main', 'inflammation_val', 'inflammation_ext'])
    inflam = df[mask].copy()
    inflam['atlas'] = 'inflammation'
    return pd.concat([df[~mask], inflam], ignore_index=True)


def load_method_comparison():
    with open(VIZ_DIR / 'method_comparison.json') as f:
        return json.load(f)


def load_scatter(directory, filename):
    fp = directory / filename
    if fp.exists():
        with open(fp) as f:
            return json.load(f)
    return None


SIG_DISPLAY = {
    'cytosig': 'CytoSig',
    'lincyto_best': 'LinCytoSig Best\n(comb+filt)',
    'lincytosig': 'LinCytoSig',
    'secact': 'SecAct',
}


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
# FIGURE 0: Dataset Summary Bar Chart (for dataset_analytics.html)
# ═══════════════════════════════════════════════════════════════════════════════
def fig0_dataset_summary():
    """Bar-chart summary of datasets, signatures, and validation layers.

    Used in dataset_analytics.html. Separates single-cell and bulk datasets
    into two sub-panels with appropriate scales.
    Output:
        fig0_dataset_summary.png, fig0_dataset_summary.pdf
    """
    fig = plt.figure(figsize=(16, 5.5))
    gs = gridspec.GridSpec(1, 4, width_ratios=[3, 2, 3, 3], wspace=0.35)

    # Panel A1: Single-cell datasets (millions)
    ax1 = fig.add_subplot(gs[0])
    sc_names = ['parse_10M', 'scAtlas', 'CIMA', 'Inflammation\nAtlas']
    sc_vals = [9.7, 6.4, 6.5, 6.3]
    sc_colors = ['#F59E0B', '#10B981', '#3B82F6', '#EF4444']
    bars = ax1.barh(sc_names, sc_vals, color=sc_colors, edgecolor='white', linewidth=0.5)
    for bar, v in zip(bars, sc_vals):
        ax1.text(bar.get_width() + 0.15, bar.get_y() + bar.get_height() / 2,
                 f'{v}M', va='center', fontsize=9, fontweight='bold')
    ax1.set_xlabel('Cells (millions)')
    ax1.set_title('A. Single-Cell Datasets', fontweight='bold')
    ax1.set_xlim(0, 13)

    # Panel A2: Bulk RNA-seq datasets (thousands)
    ax2 = fig.add_subplot(gs[1])
    bulk_names = ['GTEx', 'TCGA']
    bulk_vals = [19.8, 11.1]
    bulk_colors = ['#6366F1', '#EC4899']
    bars = ax2.barh(bulk_names, bulk_vals, color=bulk_colors, edgecolor='white', linewidth=0.5)
    for bar, v in zip(bars, bulk_vals):
        ax2.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height() / 2,
                 f'{v}K', va='center', fontsize=9, fontweight='bold')
    ax2.set_xlabel('Samples (thousands)')
    ax2.set_title('B. Bulk RNA-seq', fontweight='bold')
    ax2.set_xlim(0, 28)

    # Panel B: Signature matrices
    ax3 = fig.add_subplot(gs[2])
    sig_types = ['CytoSig\n(43)', 'LinCytoSig\n(178)', 'SecAct\n(1,170)']
    counts = [43, 178, 1170]
    sig_colors = [COLORS['cytosig'], COLORS['lincytosig'], COLORS['secact']]
    bars = ax3.bar(sig_types, counts, color=sig_colors, edgecolor='white', linewidth=0.5)
    for bar, v in zip(bars, counts):
        ax3.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 20,
                 f'{v:,}', ha='center', fontsize=9, fontweight='bold')
    ax3.set_ylabel('Signatures')
    ax3.set_title('C. Signature Matrices (1,391 total)', fontweight='bold')

    # Panel C: Validation layers
    ax4 = fig.add_subplot(gs[3])
    layers = ['Bulk RNA-seq\n(GTEx/TCGA)', 'Pseudobulk\n(donor)', 'Pseudobulk\n(donor×celltype)',
              'Single-cell\n(per cell)']
    layer_counts = [2, 8, 8, 6]
    layer_colors = ['#8B5CF6', '#3B82F6', '#10B981', '#F59E0B']
    bars = ax4.barh(layers, layer_counts, color=layer_colors, edgecolor='white')
    for bar, v in zip(bars, layer_counts):
        ax4.text(bar.get_width() + 0.15, bar.get_y() + bar.get_height() / 2,
                 str(v), va='center', fontsize=9, fontweight='bold')
    ax4.set_xlabel('Dataset × Level Combinations')
    ax4.set_title('D. Validation Layers', fontweight='bold')
    ax4.set_xlim(0, 12)

    fig.suptitle('Dataset Overview: ~29M Single Cells + ~31K Bulk RNA-seq Samples',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'fig0_dataset_summary.png')
    fig.savefig(FIG_DIR / 'fig0_dataset_summary.pdf')
    plt.close(fig)
    print('  ✓ Figure 0: Dataset summary (bar chart)')


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 1: Dataset Overview
# ═══════════════════════════════════════════════════════════════════════════════
def fig1_dataset_overview():
    """Schematic flow diagram: Data Sources → Activity Inference → Validation.

    Data source:
        Hardcoded values from H5AD .n_obs metadata:
        CIMA 6.5M, Inflammation Atlas Main 6.3M, scAtlas 6.4M, parse_10M 9.7M.
        GTEx 19.8K, TCGA 11.1K bulk RNA-seq samples.
        Signature counts from secactpy: CytoSig 43, LinCytoSig 178, SecAct 1,170.
    Method:
        Schematic — no statistical computation.
    Output:
        fig1_dataset_overview.png, fig1_dataset_overview.pdf
    """
    fig, ax = plt.subplots(1, 1, figsize=(16, 9))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 12)
    ax.axis('off')

    # ── Helper: draw a rounded box (no double-draw, clean single pass) ──
    def draw_box(x, y, w, h, color, alpha=0.12, lw=1.5):
        bg = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.1',
                            facecolor=color, edgecolor='none', alpha=alpha, linewidth=0)
        ax.add_patch(bg)
        border = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.1',
                                facecolor='none', edgecolor=color, linewidth=lw, alpha=0.6)
        ax.add_patch(border)

    # ── Helper: draw a group container (lighter, larger pad) ──
    def draw_group(x, y, w, h, color):
        bg = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.15',
                            facecolor=color, edgecolor='none', alpha=0.04, linewidth=0)
        ax.add_patch(bg)
        border = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.15',
                                facecolor='none', edgecolor=color, linewidth=1.2, alpha=0.3,
                                linestyle='--')
        ax.add_patch(border)

    # ── Helper: draw an arrow between columns ──
    def draw_arrow(x1, y1, x2, y2):
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='->', color='#64748B', lw=2.5,
                                    connectionstyle='arc3,rad=0'))

    # Card height and vertical gap
    CH = 0.95    # card height
    CG = 0.2     # gap between cards
    CS = CH + CG  # card step
    SH = 0.8     # shorter card (signatures, no subtitle)
    SS = SH + CG  # signature step
    PSH = 0.85   # pipeline step height
    PSS = PSH + 0.15
    GROUP_GAP = 0.8   # gap between groups within a column
    GROUP_PAD = 0.15  # padding inside group box around cards

    TOP = 10.5  # top of content (all columns start here)

    # ── Pre-compute natural bottoms for each column ──
    # Col 1: SC (4 cards) + gap + Bulk (2 cards)
    sc_content_h = len([1,2,3,4]) * CS - CG
    bulk_content_h = 2 * CS - CG
    col1_natural_bottom = TOP - sc_content_h - GROUP_GAP - bulk_content_h

    # Col 2: Sigs (3 cards) + gap + Pipeline (4 steps)
    sig_content_h = 3 * SS - CG
    pipe_content_h = 4 * PSS - 0.15
    col2_natural_bottom = TOP - sig_content_h - GROUP_GAP - pipe_content_h

    # Col 3: Validation (7 cards)
    val_content_h = 7 * CS - CG

    # Global bottom = lowest bottom across columns, with padding for group box
    global_bottom = min(col1_natural_bottom, col2_natural_bottom,
                        TOP - val_content_h) - GROUP_PAD

    # ═══ COLUMN 1: Data Sources (x=0.3–4.7) ═══
    col1_x, col1_w = 0.3, 4.4

    ax.text(col1_x + col1_w / 2, 11.5, 'Data Sources', fontsize=14, fontweight='bold',
            ha='center', va='center', color='#1E293B')
    ax.text(col1_x + col1_w / 2, 11.05, '~29M cells  +  ~31K bulk samples',
            fontsize=9, ha='center', va='center', color='#64748B')

    # -- Single-cell group --
    sc_top = TOP
    sc_datasets = [
        ('CIMA',               '6.5M cells',  '421 donors, healthy population',   '#3B82F6'),
        ('Inflammation Atlas', '6.3M cells',  'main/val/ext, 20 diseases',        '#EF4444'),
        ('scAtlas',            '6.4M cells',  'Normal + Cancer, multi-organ',     '#10B981'),
        ('parse_10M',          '9.7M cells',  '12 donors × 90 cytokines (+PBS)',  '#F59E0B'),
    ]
    sc_cards_bottom = sc_top - len(sc_datasets) * CS + CG
    draw_group(col1_x, sc_cards_bottom - GROUP_PAD, col1_w,
               sc_top - sc_cards_bottom + GROUP_PAD + 0.4, '#3B82F6')
    ax.text(col1_x + 0.2, sc_top + 0.15, 'Single-Cell', fontsize=10.5,
            fontweight='bold', va='bottom', color='#1E40AF')

    for i, (name, count, desc, color) in enumerate(sc_datasets):
        by = sc_top - i * CS - CH
        draw_box(col1_x + 0.2, by, col1_w - 0.4, CH, color)
        ax.text(col1_x + 0.45, by + CH * 0.62, name, fontsize=10.5, fontweight='bold',
                va='center', color='#1E293B')
        ax.text(col1_x + col1_w - 0.45, by + CH * 0.65, count, fontsize=9.5, fontweight='bold',
                va='center', ha='right', color=color)
        ax.text(col1_x + 0.45, by + CH * 0.25, desc, fontsize=8, va='center', color='#64748B')

    # -- Bulk RNA-seq group (extends to global bottom) --
    bulk_top = sc_cards_bottom - GROUP_GAP
    bulk_datasets = [
        ('GTEx', '19.8K samples', '~50 normal tissues, TPM', '#6366F1'),
        ('TCGA', '11.1K samples', '33 cancer types, RSEM',   '#EC4899'),
    ]
    draw_group(col1_x, global_bottom, col1_w,
               bulk_top - global_bottom + 0.4, '#8B5CF6')
    ax.text(col1_x + 0.2, bulk_top + 0.15, 'Bulk RNA-seq', fontsize=10.5,
            fontweight='bold', va='bottom', color='#5B21B6')

    for i, (name, count, desc, color) in enumerate(bulk_datasets):
        by = bulk_top - i * CS - CH
        draw_box(col1_x + 0.2, by, col1_w - 0.4, CH, color)
        ax.text(col1_x + 0.45, by + CH * 0.62, name, fontsize=10.5, fontweight='bold',
                va='center', color='#1E293B')
        ax.text(col1_x + col1_w - 0.45, by + CH * 0.65, count, fontsize=9.5, fontweight='bold',
                va='center', ha='right', color=color)
        ax.text(col1_x + 0.45, by + CH * 0.25, desc, fontsize=8, va='center', color='#64748B')

    # ═══ ARROWS: Column 1 → Column 2 ═══
    mid_y = (TOP + global_bottom) / 2
    draw_arrow(5.0, mid_y, 5.7, mid_y)

    # ═══ COLUMN 2: Activity Inference (x=5.8–10.2) ═══
    col2_x, col2_w = 5.8, 4.4

    ax.text(col2_x + col2_w / 2, 11.5, 'Activity Inference', fontsize=14, fontweight='bold',
            ha='center', va='center', color='#1E293B')
    ax.text(col2_x + col2_w / 2, 11.05, 'Ridge regression  ·  GPU-accelerated',
            fontsize=9, ha='center', va='center', color='#64748B')

    # -- Signature matrices group --
    sig_top = TOP
    sigs = [
        ('CytoSig',    '43 cytokines',           COLORS['cytosig']),
        ('LinCytoSig', '178 cell-type specific',  COLORS['lincytosig']),
        ('SecAct',     '1,170 secreted proteins', COLORS['secact']),
    ]
    sig_cards_bottom = sig_top - len(sigs) * SS + CG
    draw_group(col2_x, sig_cards_bottom - GROUP_PAD, col2_w,
               sig_top - sig_cards_bottom + GROUP_PAD + 0.4, '#64748B')
    ax.text(col2_x + 0.2, sig_top + 0.15, '3 Signature Matrices', fontsize=10.5,
            fontweight='bold', va='bottom', color='#374151')

    for i, (name, desc, color) in enumerate(sigs):
        by = sig_top - i * SS - SH
        draw_box(col2_x + 0.2, by, col2_w - 0.4, SH, color)
        ax.text(col2_x + 0.45, by + SH / 2, name, fontsize=10.5, fontweight='bold',
                va='center', color='#1E293B')
        ax.text(col2_x + col2_w - 0.45, by + SH / 2, desc, fontsize=9,
                va='center', ha='right', color=color)

    # -- Pipeline group (extends to global bottom) --
    pipe_top = sig_cards_bottom - GROUP_GAP
    steps = [
        ('1.', 'Pseudobulk aggregation',    '(donor, donor×celltype)'),
        ('2.', 'Gene expression matrix',     '(genes × samples)'),
        ('3.', 'Ridge regression',           '(signature × samples → activity)'),
        ('4.', 'Cross-sample correlation',   '(predicted vs observed)'),
    ]
    draw_group(col2_x, global_bottom, col2_w,
               pipe_top - global_bottom + 0.4, '#64748B')
    ax.text(col2_x + 0.2, pipe_top + 0.15, 'Pipeline', fontsize=10.5,
            fontweight='bold', va='bottom', color='#374151')

    for i, (num, step, detail) in enumerate(steps):
        sy_center = pipe_top - i * PSS - PSH / 2
        ax.text(col2_x + 0.45, sy_center + 0.1, num + '  ' + step,
                fontsize=9.5, va='center', color='#1E293B')
        ax.text(col2_x + 0.75, sy_center - 0.2, detail,
                fontsize=8, va='center', color='#94A3B8')

    # ═══ ARROWS: Column 2 → Column 3 ═══
    draw_arrow(10.5, mid_y, 11.2, mid_y)

    # ═══ COLUMN 3: Validation (x=11.3–15.7) ═══
    col3_x, col3_w = 11.3, 4.4

    ax.text(col3_x + col3_w / 2, 11.5, 'Validation', fontsize=14, fontweight='bold',
            ha='center', va='center', color='#1E293B')
    ax.text(col3_x + col3_w / 2, 11.05, '6 datasets  ×  3 signatures',
            fontsize=9, ha='center', va='center', color='#64748B')

    val_top = TOP
    val_items = [
        ('Overall Performance',      'Spearman ρ across all targets',       '#2563EB', '§4.1'),
        ('Per-Tissue Stratified',    'GTEx tissues / TCGA cancer types',    '#7C3AED', '§4.2'),
        ('Cross-Dataset Comparison', 'CytoSig vs SecAct across 6 datasets', '#059669', '§4.3'),
        ('Best/Worst Targets',       'Top and bottom correlated targets',   '#DC2626', '§4.4'),
        ('Cross-Dataset Consistency','Same targets across datasets',        '#D97706', '§4.5'),
        ('Aggregation Levels',       'Donor → celltype → single-cell',     '#0891B2', '§4.6'),
        ('Bulk RNA-seq Validation',  'GTEx + TCGA concordance',             '#6366F1', '§4.7'),
    ]
    # Validation outline box (extends to global bottom)
    draw_group(col3_x, global_bottom, col3_w,
               val_top - global_bottom + 0.4, '#64748B')
    ax.text(col3_x + 0.2, val_top + 0.15, '7 Validation Analyses', fontsize=10.5,
            fontweight='bold', va='bottom', color='#374151')

    for i, (title, desc, color, sec) in enumerate(val_items):
        by = val_top - i * CS - CH
        draw_box(col3_x + 0.2, by, col3_w - 0.4, CH, color, alpha=0.10)
        ax.text(col3_x + 0.45, by + CH * 0.62, title, fontsize=9.5, fontweight='bold',
                va='center', color='#1E293B')
        ax.text(col3_x + 0.45, by + CH * 0.25, desc, fontsize=8, va='center', color='#64748B')
        ax.text(col3_x + col3_w - 0.45, by + CH * 0.5, sec, fontsize=8.5, fontweight='bold',
                va='center', ha='right', color=color)

    # ═══ Title ═══
    fig.suptitle('CytoAtlas: Pan-Disease Single-Cell Cytokine Activity Atlas',
                 fontsize=15, fontweight='bold', y=0.99, color='#0F172A')

    fig.savefig(FIG_DIR / 'fig1_dataset_overview.png')
    fig.savefig(FIG_DIR / 'fig1_dataset_overview.pdf')
    plt.close(fig)
    print('  ✓ Figure 1: Dataset overview')


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 2: Correlation Summary Across All Atlases (Boxplot)
# ═══════════════════════════════════════════════════════════════════════════════
def fig2_correlation_summary(df):
    """Boxplot of Spearman rho across all atlas × level × signature combinations.

    Data source:
        correlations/{atlas}_correlations.csv (single-cell atlases),
        correlations/all_correlations.csv (GTEx, TCGA),
        lincytosig_gene_filter/*_lincyto_filt_correlations.csv (LinCytoSig Best).
    Method:
        Filter to donor-level (donor_only or donor_organ) per atlas.
        Group by atlas within each signature panel (2x2 grid).
        Box = Q1–Q3, whiskers = 1.5×IQR, outliers hidden.
        Each data point is one target's Spearman rho.
    Output:
        fig2_correlation_summary_boxplot.png, fig2_correlation_summary_boxplot.pdf
    """
    sig_types = ['cytosig', 'lincyto_best', 'lincytosig', 'secact']
    fig, axes = plt.subplots(2, 2, figsize=(18, 12), sharey=True)
    axes_flat = axes.flatten()

    for idx, sig_type in enumerate(sig_types):
        ax = axes_flat[idx]
        sub = df[df['signature'] == sig_type].copy()

        # Get donor-level only for consistent comparison
        donor_levels = sub[sub['level'].str.contains('donor', case=False)]

        # Group by atlas
        atlas_order = ['cima', 'inflammation', 'scatlas_normal', 'scatlas_cancer', 'gtex', 'tcga']
        atlas_labels = ['CIMA', 'Inflammation\nAtlas', 'scAtlas\n(Normal)', 'scAtlas\n(Cancer)', 'GTEx', 'TCGA']

        data_boxes = []
        positions = []
        box_colors = []
        used_labels = []
        for i, (atlas, label) in enumerate(zip(atlas_order, atlas_labels)):
            subset = donor_levels[donor_levels['atlas'] == atlas]
            if len(subset) > 0:
                data_boxes.append(subset['spearman_rho'].dropna().values)
                positions.append(i)
                box_colors.append(ATLAS_COLORS.get(atlas, '#6B7280'))
                used_labels.append(label)

        if data_boxes:
            bp = ax.boxplot(data_boxes, positions=positions, widths=0.6, patch_artist=True,
                           showfliers=False, medianprops=dict(color='black', linewidth=2))
            for patch, color in zip(bp['boxes'], box_colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)

            ax.set_xticks(positions)
            ax.set_xticklabels(used_labels, fontsize=8, rotation=0)

        ax.axhline(0, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
        ax.set_title(SIG_DISPLAY.get(sig_type, sig_type.upper()), fontsize=13,
                     fontweight='bold', color=COLORS[sig_type])
        if idx % 2 == 0:
            ax.set_ylabel('Spearman ρ (expression vs activity)')

    fig.suptitle('Cross-Sample Validation: Spearman Correlation Distributions\n(Donor-Level Pseudobulk)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'fig2_correlation_summary_boxplot.png')
    fig.savefig(FIG_DIR / 'fig2_correlation_summary_boxplot.pdf')
    plt.close(fig)
    print('  ✓ Figure 2: Correlation summary boxplot')


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 3: Top / Bottom Targets — Good vs Bad Correlations
# ═══════════════════════════════════════════════════════════════════════════════
def fig3_good_bad_correlations(df):
    """Bar charts of best and worst correlated targets for CytoSig across atlases.

    Data source:
        correlations/{atlas}_correlations.csv filtered to signature='cytosig'.
        Atlases: CIMA (donor_only), Inflammation Main (donor_only),
        scAtlas Normal (donor_organ).
    Method:
        Sort targets by spearman_rho descending. Take top 15 (good) and
        bottom 15 (bad). Color bars by cytokine family membership
        (Interferon, TGF-beta, Interleukin, TNF, Growth Factor, Chemokine,
        Colony-Stimulating). Annotate rho values on each bar.
    Output:
        fig3_good_bad_correlations_cytosig.png, fig3_good_bad_correlations_cytosig.pdf
    """
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))

    atlas_configs = [
        ('cima', 'donor_only', 'CIMA (Donor)'),
        ('inflammation_main', 'donor_only', 'Inflammation Main (Donor)'),
        ('scatlas_normal', 'donor_organ', 'scAtlas Normal (Organ)'),
    ]

    for col, (atlas, level, title) in enumerate(atlas_configs):
        sub = df[(df['atlas'] == atlas) & (df['level'] == level) & (df['signature'] == 'cytosig')]
        sub = sub.sort_values('spearman_rho', ascending=False)

        # Top 15 (good)
        ax = axes[0, col]
        top = sub.head(15)
        colors_top = [FAMILY_COLORS.get(TARGET_TO_FAMILY.get(t, 'Other'), '#9CA3AF')
                      for t in top['target']]
        bars = ax.barh(range(len(top)), top['spearman_rho'].values, color=colors_top)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top['target'].values, fontsize=9)
        ax.invert_yaxis()
        ax.set_xlabel('Spearman ρ')
        ax.set_title(f'{title}\nTop 15 (Best Correlation)', fontweight='bold')
        ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
        # Annotate rho values
        for i, (_, row) in enumerate(top.iterrows()):
            ax.text(row['spearman_rho'] + 0.01, i, f'{row["spearman_rho"]:.3f}',
                    va='center', fontsize=8, color='black')

        # Bottom 15 (bad)
        ax = axes[1, col]
        bottom = sub.tail(15).iloc[::-1]
        colors_bot = [FAMILY_COLORS.get(TARGET_TO_FAMILY.get(t, 'Other'), '#9CA3AF')
                      for t in bottom['target']]
        bars = ax.barh(range(len(bottom)), bottom['spearman_rho'].values, color=colors_bot)
        ax.set_yticks(range(len(bottom)))
        ax.set_yticklabels(bottom['target'].values, fontsize=9)
        ax.invert_yaxis()
        ax.set_xlabel('Spearman ρ')
        ax.set_title(f'{title}\nBottom 15 (Worst Correlation)', fontweight='bold')
        ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
        for i, (_, row) in enumerate(bottom.iterrows()):
            ax.text(max(row['spearman_rho'] - 0.15, -0.5), i, f'{row["spearman_rho"]:.3f}',
                    va='center', fontsize=8, color='black')

    # Legend for families
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c) for fam, c in FAMILY_COLORS.items()]
    fig.legend(handles, list(FAMILY_COLORS.keys()), loc='lower center', ncol=4, fontsize=9,
               frameon=True, title='Cytokine Family')

    fig.suptitle('CytoSig Validation: Best & Worst Correlated Targets by Atlas',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0.05, 1, 0.96])
    fig.savefig(FIG_DIR / 'fig3_good_bad_correlations_cytosig.png')
    fig.savefig(FIG_DIR / 'fig3_good_bad_correlations_cytosig.pdf')
    plt.close(fig)
    print('  ✓ Figure 3: Good/bad correlations (CytoSig)')


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 4: Biologically Important Targets Heatmap
# ═══════════════════════════════════════════════════════════════════════════════
def fig4_bio_targets_heatmap(df):
    """Heatmap of Spearman rho for biologically important targets across atlases.

    Data source:
        Merged correlation CSVs at donor level per atlas (6 atlases).
        Three panels: CytoSig (all 43 targets), LinCytoSig Best (43),
        SecAct (matched targets from BIO_FAMILIES).
    Method:
        Build matrix (targets × atlases) of donor-level spearman_rho.
        Display as imshow with RdBu_r colormap, range [-0.5, 0.8].
        Annotate cells with rho values; white text for |rho| > 0.4.
    Output:
        fig4_bio_targets_heatmap.png, fig4_bio_targets_heatmap.pdf
    """
    fig, axes = plt.subplots(1, 3, figsize=(22, 10))

    for idx, sig_type in enumerate(['cytosig', 'lincyto_best', 'secact']):
        ax = axes[idx]

        # Filter to donor-level (simplest aggregation)
        sub = df[df['signature'] == sig_type].copy()
        # Use lowest aggregation level per atlas
        level_map = {
            'cima': 'donor_only', 'inflammation': 'donor_only',
            'scatlas_normal': 'donor_organ', 'scatlas_cancer': 'donor_organ',
            'gtex': 'donor_only', 'tcga': 'donor_only',
        }
        rows = []
        for atlas, level in level_map.items():
            atlas_sub = sub[(sub['atlas'] == atlas) & (sub['level'] == level)]
            for _, row in atlas_sub.iterrows():
                rows.append(row)
        sub = pd.DataFrame(rows)

        if sig_type in ('cytosig', 'lincyto_best'):
            # All targets are cytokine names (43 targets)
            targets = sorted(sub['target'].unique())
        else:
            # SecAct: use same cytokine names as CytoSig for comparison
            bio_targets = list(TARGET_TO_FAMILY.keys())
            targets = [t for t in sub['target'].unique() if t in bio_targets]
            targets = sorted(targets)

        if not targets:
            ax.text(0.5, 0.5, 'No matching targets', ha='center', va='center', transform=ax.transAxes)
            continue

        atlas_order = ['cima', 'inflammation', 'scatlas_normal', 'scatlas_cancer', 'gtex', 'tcga']
        atlas_labels = ['CIMA', 'Inflam Atlas', 'scAtlas (N)', 'scAtlas (C)', 'GTEx', 'TCGA']

        # Build heatmap matrix
        matrix = np.full((len(targets), len(atlas_order)), np.nan)
        for i, target in enumerate(targets):
            for j, atlas in enumerate(atlas_order):
                level = level_map[atlas]
                match = sub[(sub['target'] == target) & (sub['atlas'] == atlas) & (sub['level'] == level)]
                if len(match) > 0:
                    matrix[i, j] = match['spearman_rho'].values[0]

        im = ax.imshow(matrix, aspect='auto', cmap='RdBu_r', vmin=-0.5, vmax=0.8,
                       interpolation='nearest')
        ax.set_xticks(range(len(atlas_labels)))
        ax.set_xticklabels(atlas_labels, rotation=45, ha='right', fontsize=8)
        ax.set_yticks(range(len(targets)))
        ax.set_yticklabels(targets, fontsize=7)
        ax.set_title(SIG_DISPLAY.get(sig_type, sig_type.upper()), fontsize=13,
                     fontweight='bold', color=COLORS[sig_type])

        # Add value annotations for significant ones
        for i in range(len(targets)):
            for j in range(len(atlas_order)):
                if not np.isnan(matrix[i, j]):
                    color = 'white' if abs(matrix[i, j]) > 0.4 else 'black'
                    ax.text(j, i, f'{matrix[i, j]:.2f}', ha='center', va='center',
                            fontsize=5, color=color)

    plt.colorbar(im, ax=axes, shrink=0.6, label='Spearman ρ')
    fig.suptitle('Cross-Sample Validation: Biologically Important Targets\n(Donor-Level Spearman Correlation)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 0.95, 0.95])
    fig.savefig(FIG_DIR / 'fig4_bio_targets_heatmap.png')
    fig.savefig(FIG_DIR / 'fig4_bio_targets_heatmap.pdf')
    plt.close(fig)
    print('  ✓ Figure 4: Biologically important targets heatmap')


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 5: Scatter Plots — Representative Good & Bad Correlations
# ═══════════════════════════════════════════════════════════════════════════════
def fig5_representative_scatter(df):
    """Representative scatter plots showing good and bad correlations.

    Data source:
        visualization/data/validation/donor_scatter/cima_cytosig.json
        for scatter points. Target selection from cima_correlations.csv
        sorted by spearman_rho.
    Method:
        Top 4 targets by rho (> 0.3) and bottom 4 targets. Each scatter
        shows donor-level mean expression (x) vs predicted activity (y)
        with scipy.stats.linregress trend line. Colored green (rho > 0.2),
        red (rho < 0), or gray.
    Output:
        fig5_representative_scatter_cima.png, fig5_representative_scatter_cima.pdf
    """
    # Find best/worst CytoSig targets in CIMA
    cima_donor = df[(df['atlas'] == 'cima') & (df['level'] == 'donor_only') & (df['signature'] == 'cytosig')]
    cima_donor_sorted = cima_donor.sort_values('spearman_rho', ascending=False)

    # Key biologically important targets: pick good and bad ones
    good_targets = cima_donor_sorted[cima_donor_sorted['spearman_rho'] > 0.3].head(4)['target'].tolist()
    # Pick the MOST negative targets (sort ascending, take first 4)
    bad_targets = cima_donor_sorted.nsmallest(4, 'spearman_rho')['target'].tolist()

    targets_to_plot = good_targets + bad_targets
    n_targets = len(targets_to_plot)

    if n_targets == 0:
        print('  ⚠ No scatter targets found for figure 5')
        return

    # Load CIMA donor scatter data
    scatter_data = load_scatter(SCATTER_DONOR, 'cima_cytosig.json')
    if scatter_data is None:
        print('  ⚠ cima_cytosig.json not found')
        return

    fig, axes = plt.subplots(2, 4, figsize=(20, 10))

    for idx, target in enumerate(targets_to_plot[:8]):
        row = idx // 4
        col = idx % 4
        ax = axes[row, col]

        if target in scatter_data:
            data = scatter_data[target]
            points = data.get('points', [])
            if points:
                x = [p[0] for p in points]
                y = [p[1] for p in points]
                rho = data.get('rho', np.nan)
                pval = data.get('pval', np.nan)
                n = data.get('n', len(points))

                color = COLORS['good'] if rho > 0.2 else (COLORS['bad'] if rho < 0 else COLORS['neutral'])
                ax.scatter(x, y, alpha=0.4, s=15, c=color, edgecolors='none')

                # Add regression line
                if len(x) > 2:
                    slope, intercept, _, _, _ = stats.linregress(x, y)
                    x_line = np.linspace(min(x), max(x), 100)
                    ax.plot(x_line, slope * x_line + intercept, 'k--', linewidth=1.5, alpha=0.7)

                pval_str = f'{pval:.1e}' if pval < 0.001 else f'{pval:.4f}'
                ax.set_title(f'{target}\nρ={rho:.3f}, p={pval_str}, n={n}', fontsize=10,
                            fontweight='bold', color=color)
            else:
                ax.text(0.5, 0.5, 'No points', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(target, fontsize=10)
        else:
            ax.text(0.5, 0.5, 'Not found', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(target, fontsize=10)

        ax.set_xlabel('Mean Expression')
        ax.set_ylabel('Predicted Activity')

    axes[0, 0].annotate('GOOD CORRELATIONS', xy=(0, 1.15), xycoords='axes fraction',
                        fontsize=12, fontweight='bold', color=COLORS['good'])
    axes[1, 0].annotate('POOR CORRELATIONS', xy=(0, 1.15), xycoords='axes fraction',
                        fontsize=12, fontweight='bold', color=COLORS['bad'])

    fig.suptitle('CIMA: Donor-Level Expression vs CytoSig Predicted Activity',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(FIG_DIR / 'fig5_representative_scatter_cima.png')
    fig.savefig(FIG_DIR / 'fig5_representative_scatter_cima.pdf')
    plt.close(fig)
    print('  ✓ Figure 5: Representative scatter plots')


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 6: Cross-Atlas Consistency for Key Targets
# ═══════════════════════════════════════════════════════════════════════════════
def fig6_cross_atlas_consistency(df):
    """Show how key target correlations are consistent across atlases.

    Data source:
        Merged correlation CSVs (CytoSig) and gene-filtered LinCytoSig
        Best CSVs at donor level per atlas. 14 key targets: IFNG, IL1B,
        TNFA, TGFB1, IL6, IL10, IL17A, IL4, BMP2, EGF, HGF, VEGFA,
        CXCL12, GMCSF.
    Method:
        For each target: extract donor-level rho at each of 6 atlases.
        Solid lines = CytoSig, dashed lines = LinCytoSig Best (comb+filt).
        Lines colored by cytokine family (FAMILY_COLORS).
    Output:
        fig6_cross_atlas_consistency.png, fig6_cross_atlas_consistency.pdf
    """
    key_targets = ['IFNG', 'IL1B', 'TNFA', 'TGFB1', 'IL6', 'IL10', 'IL17A',
                   'IL4', 'BMP2', 'EGF', 'HGF', 'VEGFA', 'CXCL12', 'GMCSF']

    cytosig = df[(df['signature'] == 'cytosig')].copy()
    lincyto_best = df[(df['signature'] == 'lincyto_best')].copy()

    # Use donor-level only
    level_map = {
        'cima': 'donor_only', 'inflammation': 'donor_only',
        'scatlas_normal': 'donor_organ', 'scatlas_cancer': 'donor_organ',
        'gtex': 'donor_only', 'tcga': 'donor_only',
    }

    fig, ax = plt.subplots(figsize=(14, 8))

    atlas_order = list(level_map.keys())
    atlas_labels = ['CIMA', 'Inflammation\nAtlas', 'scAtlas\n(Normal)', 'scAtlas\n(Cancer)', 'GTEx', 'TCGA']

    # CytoSig: solid lines
    for i, target in enumerate(key_targets):
        rhos = []
        for atlas in atlas_order:
            level = level_map[atlas]
            match = cytosig[(cytosig['target'] == target) & (cytosig['atlas'] == atlas) & (cytosig['level'] == level)]
            if len(match) > 0:
                rhos.append(match['spearman_rho'].values[0])
            else:
                rhos.append(np.nan)

        color = FAMILY_COLORS.get(TARGET_TO_FAMILY.get(target, 'Other'), '#9CA3AF')
        ax.plot(range(len(atlas_order)), rhos, 'o-', label=target, color=color,
                markersize=6, linewidth=1.5, alpha=0.8)

    # LinCytoSig Best: dashed lines (same targets, same colors, smaller markers)
    for i, target in enumerate(key_targets):
        rhos = []
        for atlas in atlas_order:
            level = level_map[atlas]
            match = lincyto_best[(lincyto_best['target'] == target) & (lincyto_best['atlas'] == atlas) & (lincyto_best['level'] == level)]
            if len(match) > 0:
                rhos.append(match['spearman_rho'].values[0])
            else:
                rhos.append(np.nan)

        color = FAMILY_COLORS.get(TARGET_TO_FAMILY.get(target, 'Other'), '#9CA3AF')
        label_str = f'{target} (Best)' if i == 0 else None
        ax.plot(range(len(atlas_order)), rhos, 's--', color=color,
                markersize=4, linewidth=1.0, alpha=0.5,
                label='_nolegend_')

    # Add a single legend entry for the dashed style
    from matplotlib.lines import Line2D
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Line2D([0], [0], color='black', linestyle='-', linewidth=1.5, marker='o', markersize=6, label='CytoSig'))
    handles.append(Line2D([0], [0], color='black', linestyle='--', linewidth=1.0, marker='s', markersize=4, alpha=0.5, label='LinCytoSig Best'))

    ax.set_xticks(range(len(atlas_labels)))
    ax.set_xticklabels(atlas_labels, rotation=30, ha='right')
    ax.set_ylabel('Spearman ρ')
    ax.set_title('Cross-Atlas Consistency of Key Cytokine Targets\n(Solid = CytoSig, Dashed = LinCytoSig Best; Donor-Level)',
                 fontsize=13, fontweight='bold')
    ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
    ax.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, frameon=True)

    plt.tight_layout()
    fig.savefig(FIG_DIR / 'fig6_cross_atlas_consistency.png')
    fig.savefig(FIG_DIR / 'fig6_cross_atlas_consistency.pdf')
    plt.close(fig)
    print('  ✓ Figure 6: Cross-atlas consistency')


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 7: Validation Levels Comparison (Donor vs Celltype vs Bulk)
# ═══════════════════════════════════════════════════════════════════════════════
def fig7_validation_levels(df):
    """Compare correlations across different aggregation levels for all atlases.

    Data source:
        Merged correlation CSVs filtered by atlas and aggregation level.
        LinCytoSig Best from gene-filtered CSVs via best_lincytosig_selection.json.
    Method:
        For each atlas x level x signature: extract rho array. Plot paired
        boxplots (CytoSig blue, LinCytoSig Best brown). Median labels above.
        CIMA: 5 levels (donor, L1-L4). Inflammation: 3 (donor, L1, L2).
        scAtlas: 3 (donor x organ, celltype1, celltype2).
        GTEx/TCGA: donor-only level (not shown — bulk RNA-seq).
    Output:
        fig7_validation_levels.png, fig7_validation_levels.pdf
    """
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes_flat = axes.flatten()

    atlas_configs = [
        ('cima', ['donor_only', 'donor_l1', 'donor_l2', 'donor_l3', 'donor_l4'], 'CIMA'),
        ('inflammation', ['donor_only', 'donor_l1', 'donor_l2'], 'Inflammation Atlas'),
        ('scatlas_normal', ['donor_organ', 'donor_organ_celltype1', 'donor_organ_celltype2'], 'scAtlas Normal'),
        ('scatlas_cancer', ['donor_organ', 'donor_organ_celltype1', 'donor_organ_celltype2'], 'scAtlas Cancer'),
    ]

    level_labels = {
        'donor_only': 'Donor\nOnly',
        'donor_l1': 'Donor\n× L1',
        'donor_l2': 'Donor\n× L2',
        'donor_l3': 'Donor\n× L3',
        'donor_l4': 'Donor\n× L4',
        'donor_organ': 'Donor\n× Organ',
        'donor_organ_celltype1': 'Donor × Organ\n× Celltype1',
        'donor_organ_celltype2': 'Donor × Organ\n× Celltype2',
    }

    sig_configs = [
        ('cytosig', COLORS['cytosig'], 'CytoSig'),
        ('lincyto_best', COLORS['lincyto_best'], 'LinCytoSig Best'),
    ]

    from matplotlib.patches import Patch

    for idx, (atlas, levels, title) in enumerate(atlas_configs):
        ax = axes_flat[idx]

        # Compute valid levels (where at least one sig has data)
        valid_levels = []
        for level in levels:
            for sig_type, _, _ in sig_configs:
                sub = df[(df['atlas'] == atlas) & (df['signature'] == sig_type)]
                level_data = sub[sub['level'] == level]['spearman_rho'].dropna()
                if len(level_data) > 0:
                    if level not in valid_levels:
                        valid_levels.append(level)
                    break

        n_levels = len(valid_levels)
        n_sigs = len(sig_configs)
        width = 0.35
        labels = [level_labels.get(lv, lv) for lv in valid_levels]

        for sig_idx, (sig_type, sig_color, sig_name) in enumerate(sig_configs):
            sub = df[(df['atlas'] == atlas) & (df['signature'] == sig_type)]
            data_boxes = []
            positions = []
            for lv_idx, level in enumerate(valid_levels):
                level_data = sub[sub['level'] == level]['spearman_rho'].dropna()
                if len(level_data) > 0:
                    data_boxes.append(level_data.values)
                    positions.append(lv_idx + 1 + (sig_idx - 0.5) * width)
                else:
                    data_boxes.append([])
                    positions.append(lv_idx + 1 + (sig_idx - 0.5) * width)

            non_empty = [(d, p) for d, p in zip(data_boxes, positions) if len(d) > 0]
            if non_empty:
                bp = ax.boxplot([d for d, _ in non_empty],
                               positions=[p for _, p in non_empty],
                               widths=width * 0.85, patch_artist=True, showfliers=False,
                               medianprops=dict(color='black', linewidth=2))
                for patch in bp['boxes']:
                    patch.set_facecolor(sig_color)
                    patch.set_alpha(0.7)

                for d, p in non_empty:
                    med = np.median(d)
                    ax.text(p, med + 0.02, f'{med:.3f}', ha='center', fontsize=6, fontweight='bold',
                            color=sig_color)

        ax.set_xticks(range(1, n_levels + 1))
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylabel('Spearman ρ' if idx % 2 == 0 else '')
        ax.set_title(f'{title}', fontweight='bold', fontsize=11)
        ax.axhline(0, color='gray', linestyle='--', alpha=0.5)

        if idx == 0:
            ax.legend(handles=[Patch(facecolor=c, alpha=0.7, label=n) for _, c, n in sig_configs],
                      fontsize=8, loc='upper right')

    fig.suptitle('Effect of Aggregation Level on Validation Correlations\n'
                 'CytoSig (blue) vs LinCytoSig Best (brown)\n'
                 '(Finer cell-type stratification → more data points but lower per-target correlations)',
                 fontsize=13, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(FIG_DIR / 'fig7_validation_levels.png')
    fig.savefig(FIG_DIR / 'fig7_validation_levels.pdf')
    plt.close(fig)
    print('  ✓ Figure 7: Validation levels comparison (4 atlases)')


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 8: 6-Way Method Comparison (CytoSig, 4× LinCytoSig, SecAct)
# ═══════════════════════════════════════════════════════════════════════════════
def fig8_method_comparison():
    """10-way comparison: CytoSig, 8x LinCytoSig variants, SecAct.

    Data source:
        visualization/data/validation/method_comparison_8way_all.json
        (pre-computed by compute_lincytosig_gene_filter.py).
    Method:
        Load rho arrays for ~20 matched cytokines per method per atlas.
        Display as boxplots with outliers. Median labels annotated.
        Uses donor-level (CIMA/Inflammation) or donor x organ (scAtlas)
        pseudobulk correlations for matched cytokines across 4 combined
        atlases: CIMA, Inflammation Atlas, scAtlas Normal, scAtlas Cancer.
    Output:
        fig8_method_comparison.png, fig8_method_comparison.pdf
    """
    eightway_path = VIZ_DIR / 'method_comparison_8way_all.json'
    if not eightway_path.exists():
        eightway_path = VIZ_DIR / 'method_comparison_6way_all.json'
    with open(eightway_path) as f:
        data = json.load(f)

    methods = [
        ('cytosig', 'CytoSig', '#2563EB'),
        ('lincyto_orig', 'LinCytoSig\n(no filter)', '#D97706'),
        ('lincyto_filt', 'LinCytoSig\n(gene filter)', '#F59E0B'),
        ('lincyto_best_orig', 'Best\n(combined)', '#B45309'),
        ('lincyto_best_filt', 'Best\n(comb+filt)', '#92400E'),
        ('lincyto_best_gtex', 'Best\n(GTEx)', '#7C3AED'),
        ('lincyto_best_tcga', 'Best\n(TCGA)', '#EC4899'),
        ('lincyto_best_gtex_filt', 'Best\n(GTEx+filt)', '#A78BFA'),
        ('lincyto_best_tcga_filt', 'Best\n(TCGA+filt)', '#F9A8D4'),
        ('secact', 'SecAct', '#059669'),
    ]
    atlas_order = ['CIMA', 'Inflammation Atlas', 'scAtlas Normal', 'scAtlas Cancer']

    fig, axes = plt.subplots(2, 2, figsize=(22, 14))
    axes_flat = axes.flatten()

    for ax_idx, atlas_name in enumerate(atlas_order):
        ax = axes_flat[ax_idx]
        if atlas_name not in data:
            ax.set_visible(False)
            continue

        atlas_data = data[atlas_name]
        n_cytokines = len(atlas_data['cytokines'])

        box_data = []
        box_labels = []
        box_colors = []
        for method_key, method_label, method_color in methods:
            vals = atlas_data.get(method_key, [])
            vals = [v for v in vals if v is not None]
            if vals:
                box_data.append(vals)
                box_labels.append(method_label)
                box_colors.append(method_color)

        if not box_data:
            ax.set_visible(False)
            continue

        bp = ax.boxplot(box_data, patch_artist=True, showfliers=True,
                       flierprops=dict(marker='o', markersize=3, alpha=0.3),
                       medianprops=dict(color='black', linewidth=2),
                       widths=0.55)
        for patch, color in zip(bp['boxes'], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)

        ax.set_xticklabels(box_labels, fontsize=5.5, rotation=40, ha='right')
        ax.set_ylabel('Spearman ρ')
        ax.set_title(f'{atlas_name} (n={n_cytokines} cytokines)', fontweight='bold', fontsize=11)
        ax.axhline(0, color='gray', linestyle='--', alpha=0.5)

        # Median labels
        for i, vals in enumerate(box_data):
            if vals:
                med = np.median(vals)
                ax.text(i + 1, med + 0.02, f'{med:.3f}', ha='center', fontsize=5,
                        fontweight='bold', color=box_colors[i])

    fig.suptitle('10-Way Method Comparison: CytoSig vs LinCytoSig Variants vs SecAct\n'
                 '(Donor-level pseudobulk validation per matched cytokine)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(FIG_DIR / 'fig8_method_comparison.png')
    fig.savefig(FIG_DIR / 'fig8_method_comparison.pdf')
    plt.close(fig)
    print('  ✓ Figure 8: 10-way method comparison (4 combined atlases)')


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 9: LinCytoSig vs CytoSig — Matched Target Comparison
# ═══════════════════════════════════════════════════════════════════════════════
def fig9_lincytosig_vs_cytosig():
    """Scatter: LinCytoSig vs CytoSig matched correlations, identify winners.

    Data source:
        visualization/data/validation/method_comparison.json.
        Uses matched_targets mapping and per-category rhos.
    Method:
        For each matched LinCytoSig-CytoSig pair per atlas category:
        plot CytoSig rho (x) vs LinCytoSig rho (y). Diagonal = equal.
        Win threshold: |diff| > 0.05. Color: amber (LinCytoSig wins),
        blue (CytoSig wins), gray (tie). Top 3 winners/losers labeled.
    Output:
        fig9_lincytosig_vs_cytosig_scatter.png, fig9_lincytosig_vs_cytosig_scatter.pdf
    """
    mc = load_method_comparison()
    matched = mc['matched_targets']  # dict: lincytosig_target -> {cytosig: str, secact: str}

    categories = mc['categories']
    cat_keys = [c['key'] for c in categories]
    cat_labels = [c['label'] for c in categories]

    n_cats = len(cat_keys)
    n_cols = 2
    n_rows = (n_cats + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6.5 * n_cols, 5.5 * n_rows))
    if n_cats > 1:
        axes = axes.flatten()
    else:
        axes = [axes]

    for cat_idx, (cat_key, cat_label) in enumerate(zip(cat_keys, cat_labels)):
        ax = axes[cat_idx]

        cytosig_rhos = []
        lincytosig_rhos = []
        labels = []
        celltype_labels = []

        for lin_target, match_info in matched.items():
            cyto_target = match_info.get('cytosig') if isinstance(match_info, dict) else match_info
            if cyto_target is None:
                continue
            lin_rho = mc['lincytosig']['rhos'].get(lin_target, {}).get(cat_key)
            cyto_rho = mc['cytosig']['rhos'].get(cyto_target, {}).get(cat_key)

            if lin_rho is not None and cyto_rho is not None:
                if not (isinstance(lin_rho, float) and np.isnan(lin_rho)):
                    if not (isinstance(cyto_rho, float) and np.isnan(cyto_rho)):
                        cytosig_rhos.append(cyto_rho)
                        lincytosig_rhos.append(lin_rho)
                        labels.append(lin_target)
                        # Extract celltype
                        parts = lin_target.rsplit('__', 1)
                        celltype_labels.append(parts[0] if len(parts) == 2 else 'Unknown')

        if not cytosig_rhos:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(cat_label)
            continue

        cytosig_arr = np.array(cytosig_rhos)
        lincytosig_arr = np.array(lincytosig_rhos)

        # Color by who wins
        colors = []
        for c, l in zip(cytosig_arr, lincytosig_arr):
            if l > c + 0.05:
                colors.append(COLORS['lincytosig'])  # LinCytoSig wins
            elif c > l + 0.05:
                colors.append(COLORS['cytosig'])  # CytoSig wins
            else:
                colors.append(COLORS['neutral'])  # Tie

        ax.scatter(cytosig_arr, lincytosig_arr, c=colors, alpha=0.6, s=40, edgecolors='white', linewidths=0.5)

        # Diagonal
        lim = [-0.5, 1.0]
        ax.plot(lim, lim, 'k--', alpha=0.4, linewidth=1)
        ax.set_xlim(lim)
        ax.set_ylim(lim)

        # Count winners
        n_lin_wins = np.sum(lincytosig_arr > cytosig_arr + 0.05)
        n_cyto_wins = np.sum(cytosig_arr > lincytosig_arr + 0.05)
        n_tie = len(cytosig_arr) - n_lin_wins - n_cyto_wins

        ax.set_xlabel('CytoSig ρ')
        ax.set_ylabel('LinCytoSig ρ')
        # Put win counts inside axes instead of title to avoid overlap with suptitle
        ax.set_title(f'{cat_label}', fontsize=10, fontweight='bold', pad=8)
        ax.text(0.02, 0.98,
                f'LinCytoSig\u2191: {n_lin_wins}   CytoSig\u2191: {n_cyto_wins}   Tie: {n_tie}',
                transform=ax.transAxes, fontsize=7.5, fontweight='bold',
                va='top', ha='left',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8, edgecolor='#ccc'))

        # Label notable points with widely spread offsets to avoid overlap
        diffs = lincytosig_arr - cytosig_arr
        top_indices = np.argsort(diffs)[-3:]  # Top 3 LinCytoSig winners
        bottom_indices = np.argsort(diffs)[:3]  # Top 3 CytoSig winners
        # Offsets: top winners go upper-left, bottom winners go lower-right
        top_offsets = [(-40, 20), (-35, 35), (-50, 8)]
        bottom_offsets = [(20, -30), (30, -20), (15, -40)]
        all_indices = list(top_indices) + list(bottom_indices)
        all_offsets = top_offsets + bottom_offsets
        for j, i in enumerate(all_indices):
            if abs(diffs[i]) > 0.1:
                ofs = all_offsets[j]
                ax.annotate(labels[i].replace('__', '\n'), (cytosig_arr[i], lincytosig_arr[i]),
                           fontsize=5, alpha=0.8, textcoords='offset points', xytext=ofs,
                           arrowprops=dict(arrowstyle='->', color='gray', lw=0.5, alpha=0.5))

    # Hide unused axes
    for i in range(cat_idx + 1, len(axes)):
        axes[i].set_visible(False)

    fig.suptitle('LinCytoSig vs CytoSig: Matched Target Correlation Comparison\n(Points above diagonal = LinCytoSig better)',
                 fontsize=13, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0.02, 1, 0.95], h_pad=3.0, w_pad=2.0)
    fig.savefig(FIG_DIR / 'fig9_lincytosig_vs_cytosig_scatter.png')
    fig.savefig(FIG_DIR / 'fig9_lincytosig_vs_cytosig_scatter.pdf')
    plt.close(fig)
    print('  ✓ Figure 9: LinCytoSig vs CytoSig matched comparison')


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 10: LinCytoSig Advantage Analysis — Which cell types benefit?
# ═══════════════════════════════════════════════════════════════════════════════
def fig10_lincytosig_advantage():
    """Analyze when and why LinCytoSig outperforms CytoSig by cell type.

    Data source:
        visualization/data/validation/method_comparison.json.
        Cell type extracted from LinCytoSig target names (CellType__Cytokine).
    Method:
        For each matched pair across all atlas categories: compute
        lin_rho - cyto_rho. Aggregate by cell type: mean difference,
        win count (diff > 0.05), loss count (diff < -0.05).
        Panel A: horizontal bar of mean delta-rho per cell type.
        Panel B: win/loss rate (%) per cell type.
    Output:
        fig10_lincytosig_advantage_by_celltype.png, fig10_lincytosig_advantage_by_celltype.pdf
    """
    mc = load_method_comparison()
    matched = mc['matched_targets']

    categories = mc['categories']
    cat_keys = [c['key'] for c in categories]

    # Collect per-celltype differences across all categories
    celltype_advantages = {}
    for lin_target, match_info in matched.items():
        cyto_target = match_info.get('cytosig') if isinstance(match_info, dict) else match_info
        if cyto_target is None:
            continue
        parts = lin_target.rsplit('__', 1)
        if len(parts) != 2:
            continue
        celltype, cytokine = parts

        if celltype not in celltype_advantages:
            celltype_advantages[celltype] = {'diffs': [], 'n_better': 0, 'n_worse': 0, 'n_total': 0}

        for cat_key in cat_keys:
            lin_rho = mc['lincytosig']['rhos'].get(lin_target, {}).get(cat_key)
            cyto_rho = mc['cytosig']['rhos'].get(cyto_target, {}).get(cat_key)
            if lin_rho is not None and cyto_rho is not None:
                if not (isinstance(lin_rho, float) and np.isnan(lin_rho)):
                    if not (isinstance(cyto_rho, float) and np.isnan(cyto_rho)):
                        diff = lin_rho - cyto_rho
                        celltype_advantages[celltype]['diffs'].append(diff)
                        celltype_advantages[celltype]['n_total'] += 1
                        if diff > 0.05:
                            celltype_advantages[celltype]['n_better'] += 1
                        elif diff < -0.05:
                            celltype_advantages[celltype]['n_worse'] += 1

    # Compute mean advantage per celltype
    for ct in celltype_advantages:
        diffs = celltype_advantages[ct]['diffs']
        celltype_advantages[ct]['mean_diff'] = np.mean(diffs) if diffs else 0
        celltype_advantages[ct]['median_diff'] = np.median(diffs) if diffs else 0

    # Sort by mean diff
    sorted_cts = sorted(celltype_advantages.items(), key=lambda x: x[1]['mean_diff'], reverse=True)

    fig, axes = plt.subplots(1, 2, figsize=(16, 10))

    # Panel A: Mean difference by celltype
    ax = axes[0]
    ct_names = [ct for ct, _ in sorted_cts]
    mean_diffs = [v['mean_diff'] for _, v in sorted_cts]
    colors = [COLORS['lincytosig'] if d > 0 else COLORS['cytosig'] for d in mean_diffs]

    bars = ax.barh(range(len(ct_names)), mean_diffs, color=colors, alpha=0.8)
    ax.set_yticks(range(len(ct_names)))
    ax.set_yticklabels(ct_names, fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, color='black', linewidth=1)
    ax.set_xlabel('Mean Δρ (LinCytoSig − CytoSig)')
    ax.set_title('A. Mean Advantage by Cell Type\n(positive = LinCytoSig better)', fontweight='bold')

    # Annotate values
    for i, (ct, v) in enumerate(sorted_cts):
        ax.text(v['mean_diff'] + 0.005 * np.sign(v['mean_diff']),
                i, f'{v["mean_diff"]:+.3f}', va='center', fontsize=7)

    # Panel B: Win rate by celltype
    ax = axes[1]
    win_rates = []
    for ct, v in sorted_cts:
        if v['n_total'] > 0:
            win_rates.append(100 * v['n_better'] / v['n_total'])
        else:
            win_rates.append(0)

    lose_rates = []
    for ct, v in sorted_cts:
        if v['n_total'] > 0:
            lose_rates.append(100 * v['n_worse'] / v['n_total'])
        else:
            lose_rates.append(0)

    ax.barh(range(len(ct_names)), win_rates, color=COLORS['lincytosig'], alpha=0.7, label='LinCytoSig wins')
    ax.barh(range(len(ct_names)), [-r for r in lose_rates], color=COLORS['cytosig'], alpha=0.7, label='CytoSig wins')
    ax.set_yticks(range(len(ct_names)))
    ax.set_yticklabels(ct_names, fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, color='black', linewidth=1)
    ax.set_xlabel('% Comparisons Won')
    ax.set_title('B. Win Rate by Cell Type', fontweight='bold')
    ax.legend(fontsize=9, loc='lower right')

    fig.suptitle('LinCytoSig Advantage Analysis:\nWhen Does Cell-Type Specificity Help?',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(FIG_DIR / 'fig10_lincytosig_advantage_by_celltype.png')
    fig.savefig(FIG_DIR / 'fig10_lincytosig_advantage_by_celltype.pdf')
    plt.close(fig)
    print('  ✓ Figure 10: LinCytoSig advantage by cell type')


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 11: Celltype-Specific Δρ (LinCytoSig − CytoSig)
# ═══════════════════════════════════════════════════════════════════════════════
def fig11_celltype_delta_rho():
    """Per-celltype mean delta-rho (LinCytoSig - CytoSig) across all atlases.

    Data source:
        visualization/data/validation/method_comparison.json
        (same as Fig 10 but different aggregation and visualization).
    Method:
        Same cell-type extraction as Fig 10. Compute per-cell-type:
        mean delta-rho, SEM (std / sqrt(n)), win count, loss count.
        Sort by mean delta-rho descending. Horizontal bar chart with
        SEM error bars. Summary annotation: overall mean, count of
        cell types favoring each method.
    Output:
        fig11_celltype_delta_rho.png, fig11_celltype_delta_rho.pdf
    """
    mc = load_method_comparison()
    matched = mc['matched_targets']
    categories = mc['categories']
    cat_keys = [c['key'] for c in categories]

    # Collect per-celltype Δρ across all atlases
    from collections import defaultdict
    ct_diffs = defaultdict(list)
    for lin_target, match_info in matched.items():
        cyto_target = match_info.get('cytosig') if isinstance(match_info, dict) else match_info
        if cyto_target is None:
            continue
        parts = lin_target.rsplit('__', 1)
        if len(parts) != 2:
            continue
        celltype = parts[0]

        for cat_key in cat_keys:
            lin_rho = mc['lincytosig']['rhos'].get(lin_target, {}).get(cat_key)
            cyto_rho = mc['cytosig']['rhos'].get(cyto_target, {}).get(cat_key)
            if lin_rho is not None and cyto_rho is not None:
                if not (isinstance(lin_rho, float) and np.isnan(lin_rho)):
                    if not (isinstance(cyto_rho, float) and np.isnan(cyto_rho)):
                        ct_diffs[celltype].append(lin_rho - cyto_rho)

    # Compute stats per celltype
    ct_stats = []
    for ct, diffs in ct_diffs.items():
        arr = np.array(diffs)
        ct_stats.append({
            'celltype': ct.replace('_', ' '),
            'mean_diff': np.mean(arr),
            'median_diff': np.median(arr),
            'std_diff': np.std(arr),
            'n': len(arr),
            'n_win': int((arr > 0).sum()),
            'n_loss': int((arr < 0).sum()),
        })
    ct_stats.sort(key=lambda x: x['mean_diff'], reverse=True)

    fig, ax = plt.subplots(figsize=(10, 12))

    names = [s['celltype'] for s in ct_stats]
    means = [s['mean_diff'] for s in ct_stats]
    stds = [s['std_diff'] / np.sqrt(s['n']) for s in ct_stats]  # SEM
    colors = [COLORS['lincytosig'] if d > 0 else COLORS['cytosig'] for d in means]

    y_pos = range(len(names))
    ax.barh(y_pos, means, xerr=stds, color=colors, alpha=0.8,
            edgecolor='white', linewidth=0.5, capsize=2, error_kw={'linewidth': 0.8})
    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, color='black', linewidth=1)
    ax.set_xlabel('Mean Δρ (LinCytoSig − CytoSig) ± SEM', fontsize=11)
    ax.set_title('Celltype-Specific Δρ: LinCytoSig vs CytoSig\n(Aggregated Across 4 Atlases at Donor × Celltype Level)',
                 fontsize=12, fontweight='bold')

    # Annotate mean Δρ and win/loss
    for i, s in enumerate(ct_stats):
        sign = '+' if s['mean_diff'] > 0 else ''
        label = f'{sign}{s["mean_diff"]:.3f}  ({s["n_win"]}W/{s["n_loss"]}L)'
        x_offset = 0.005 if s['mean_diff'] >= 0 else -0.005
        ha = 'left' if s['mean_diff'] >= 0 else 'right'
        ax.text(s['mean_diff'] + x_offset + stds[i] * (1 if s['mean_diff'] >= 0 else -1),
                i, label, va='center', ha=ha, fontsize=6.5, color='#374151')

    # Summary annotation
    all_diffs = [d for diffs in ct_diffs.values() for d in diffs]
    n_ct_pos = sum(1 for s in ct_stats if s['mean_diff'] > 0)
    ax.text(0.98, 0.98,
            f'Overall: mean Δρ = {np.mean(all_diffs):+.3f}\n'
            f'{n_ct_pos}/{len(ct_stats)} cell types favor LinCytoSig\n'
            f'Orange = LinCytoSig better, Blue = CytoSig better',
            transform=ax.transAxes, va='top', ha='right', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='#FEF3C7', alpha=0.9))

    plt.tight_layout()
    fig.savefig(FIG_DIR / 'fig11_celltype_delta_rho.png', dpi=150)
    fig.savefig(FIG_DIR / 'fig11_celltype_delta_rho.pdf')
    plt.close(fig)
    print('  ✓ Figure 11: Celltype-specific Δρ vs CytoSig')


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 12: Bulk RNA-seq Validation (GTEx & TCGA)
# ═══════════════════════════════════════════════════════════════════════════════
def fig12_bulk_validation(df):
    """GTEx and TCGA bulk RNA-seq validation results.

    Data source:
        correlations/all_correlations.csv filtered to atlas in (gtex, tcga),
        level='donor_only', per signature type (cytosig, lincytosig, secact).
    Method:
        2x3 grid: rows = GTEx, TCGA; columns = CytoSig, LinCytoSig, SecAct.
        Small target sets (<=60): ranked bar chart sorted by rho descending.
        Color: green (rho > 0.2), red (rho < -0.1), gray otherwise.
        Large sets (>60, SecAct): histogram with 40 bins + median line.
    Output:
        fig12_bulk_validation.png, fig12_bulk_validation.pdf
    """
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    for row_idx, (atlas, title) in enumerate([('gtex', 'GTEx'), ('tcga', 'TCGA')]):
        for col_idx, sig_type in enumerate(['cytosig', 'lincytosig', 'secact']):
            ax = axes[row_idx, col_idx]
            # Use donor_only level for consistent comparison
            sub = df[(df['atlas'] == atlas) & (df['signature'] == sig_type) & (df['level'] == 'donor_only')]
            rhos = sub['spearman_rho'].dropna().values

            if len(rhos) == 0:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f'{title} — {sig_type.upper()}', fontweight='bold')
                continue

            # Ranked bar chart for all
            sorted_rhos = np.sort(rhos)[::-1]
            n = len(sorted_rhos)
            bar_colors = [COLORS['good'] if r > 0.2 else (COLORS['bad'] if r < -0.1 else COLORS['neutral'])
                          for r in sorted_rhos]

            if n <= 60:
                ax.bar(range(n), sorted_rhos, color=bar_colors, alpha=0.8, width=0.8)
                ax.set_xlabel('Target (ranked by ρ)')
            else:
                # For large sets use histogram
                ax.hist(rhos, bins=40, color=COLORS[sig_type], alpha=0.7, edgecolor='white')
                ax.set_xlabel('Spearman ρ')
                ax.set_ylabel('Count')
                # Add vertical lines for median
                med = np.median(rhos)
                ax.axvline(med, color='black', linewidth=2, linestyle='-', label=f'median={med:.3f}')
                ax.legend(fontsize=9)

            ax.axhline(0, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
            n_sig = np.sum(np.abs(rhos) > 0)  # just for display
            med = np.median(rhos)
            ax.set_title(f'{title} — {sig_type.upper()}\n(n={n}, median ρ={med:.3f})',
                         fontweight='bold', color=COLORS[sig_type])
            if col_idx == 0 and n <= 60:
                ax.set_ylabel('Spearman ρ')

    fig.suptitle('Bulk RNA-seq Validation: Donor-Level Spearman Correlation\n(CytoSig/LinCytoSig/SecAct predicted activity vs target gene expression)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(FIG_DIR / 'fig12_bulk_validation.png')
    fig.savefig(FIG_DIR / 'fig12_bulk_validation.pdf')
    plt.close(fig)
    print('  ✓ Figure 12: Bulk RNA-seq validation')


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 13: LinCytoSig Cell-Type Specificity Deep Dive
# ═══════════════════════════════════════════════════════════════════════════════
def fig13_lincytosig_specificity():
    """Show that LinCytoSig captures cell-type-specific biology CytoSig misses.

    Data source:
        visualization/data/validation/method_comparison.json.
        Iterates all matched targets across all atlas categories.
    Method:
        Compute diff = lin_rho - cyto_rho for every matched pair x category.
        Rank by diff. Top 20 = LinCytoSig wins, bottom 20 = CytoSig wins.
        Display as overlapping horizontal bars (amber = LinCytoSig,
        blue = CytoSig) with delta labels annotated.
    Output:
        fig13_lincytosig_specificity.png, fig13_lincytosig_specificity.pdf
    """
    mc = load_method_comparison()
    matched = mc['matched_targets']

    # Find cases where LinCytoSig dramatically outperforms CytoSig
    examples = []
    for lin_target, match_info in matched.items():
        cyto_target = match_info.get('cytosig') if isinstance(match_info, dict) else match_info
        if cyto_target is None:
            continue
        parts = lin_target.rsplit('__', 1)
        if len(parts) != 2:
            continue
        celltype, cytokine = parts

        for cat_key in [c['key'] for c in mc['categories']]:
            lin_rho = mc['lincytosig']['rhos'].get(lin_target, {}).get(cat_key)
            cyto_rho = mc['cytosig']['rhos'].get(cyto_target, {}).get(cat_key)
            if lin_rho is not None and cyto_rho is not None:
                if not (isinstance(lin_rho, float) and np.isnan(lin_rho)):
                    if not (isinstance(cyto_rho, float) and np.isnan(cyto_rho)):
                        examples.append({
                            'lin_target': lin_target,
                            'cyto_target': cyto_target,
                            'celltype': celltype,
                            'cytokine': cytokine,
                            'category': cat_key,
                            'lin_rho': lin_rho,
                            'cyto_rho': cyto_rho,
                            'diff': lin_rho - cyto_rho,
                        })

    if not examples:
        print('  ⚠ No examples for figure 13')
        return

    ex_df = pd.DataFrame(examples)

    # Top LinCytoSig winners
    top_lin = ex_df.nlargest(20, 'diff')
    # Top CytoSig winners (LinCytoSig worst)
    top_cyto = ex_df.nsmallest(20, 'diff')

    fig, axes = plt.subplots(1, 2, figsize=(18, 10))

    # Panel A: LinCytoSig wins
    ax = axes[0]
    labels = [f'{r["celltype"]}__{r["cytokine"]}\n({r["category"].split("_")[0]})' for _, r in top_lin.iterrows()]
    x_pos = range(len(top_lin))
    ax.barh(x_pos, top_lin['lin_rho'].values, color=COLORS['lincytosig'], alpha=0.8, label='LinCytoSig')
    ax.barh(x_pos, top_lin['cyto_rho'].values, color=COLORS['cytosig'], alpha=0.5, label='CytoSig')
    ax.set_yticks(x_pos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel('Spearman ρ')
    ax.set_title('A. LinCytoSig Best Advantages\n(Δρ up to +1.6)', fontweight='bold', color=COLORS['lincytosig'])
    ax.legend(fontsize=9)
    for i, (_, row) in enumerate(top_lin.iterrows()):
        ax.text(max(row['lin_rho'], row['cyto_rho']) + 0.02, i,
                f'Δ={row["diff"]:+.3f}', va='center', fontsize=7, fontweight='bold')

    # Panel B: CytoSig wins
    ax = axes[1]
    labels = [f'{r["celltype"]}__{r["cytokine"]}\n({r["category"].split("_")[0]})' for _, r in top_cyto.iterrows()]
    x_pos = range(len(top_cyto))
    ax.barh(x_pos, top_cyto['cyto_rho'].values, color=COLORS['cytosig'], alpha=0.8, label='CytoSig')
    ax.barh(x_pos, top_cyto['lin_rho'].values, color=COLORS['lincytosig'], alpha=0.5, label='LinCytoSig')
    ax.set_yticks(x_pos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel('Spearman ρ')
    ax.set_title('B. CytoSig Best Advantages\n(Δρ up to -1.5)', fontweight='bold', color=COLORS['cytosig'])
    ax.legend(fontsize=9)
    for i, (_, row) in enumerate(top_cyto.iterrows()):
        ax.text(max(row['lin_rho'], row['cyto_rho']) + 0.02, i,
                f'Δ={row["diff"]:+.3f}', va='center', fontsize=7, fontweight='bold')

    fig.suptitle('Cell-Type Specificity: When LinCytoSig Wins vs CytoSig Wins',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIG_DIR / 'fig13_lincytosig_specificity.png')
    fig.savefig(FIG_DIR / 'fig13_lincytosig_specificity.pdf')
    plt.close(fig)
    print('  ✓ Figure 13: LinCytoSig specificity deep dive')


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 14: Celltype-Level Scatter Examples (Good vs Bad)
# ═══════════════════════════════════════════════════════════════════════════════
def fig14_celltype_scatter_examples(df):
    """Show celltype-level scatter plots for selected targets.

    Data source:
        visualization/data/validation/celltype_scatter/{atlas}_{level}_cytosig.json.
        4 key targets: IFNG, TGFB1, IL1B, IL6.
        3 atlases: CIMA (donor_l1), Inflammation Main (donor_l1),
        scAtlas Normal (donor_organ_celltype1).
    Method:
        Load scatter JSON for each atlas x level. Each point = cell-type x
        donor pseudobulk (mean expression, predicted activity). Points
        optionally colored by cell type if celltype index is available.
        4x3 grid layout.
    Output:
        fig14_celltype_scatter_examples.png, fig14_celltype_scatter_examples.pdf
    """
    # Pick a few key targets
    key_targets = ['IFNG', 'TGFB1', 'IL1B', 'IL6']

    fig, axes = plt.subplots(len(key_targets), 3, figsize=(18, 5 * len(key_targets)))

    for row_idx, target in enumerate(key_targets):
        for col_idx, (atlas, level, label) in enumerate([
            ('cima', 'donor_l1', 'CIMA (L1)'),
            ('inflammation_main', 'donor_l1', 'Inflammation Main (L1)'),
            ('scatlas_normal', 'donor_organ_celltype1', 'scAtlas Normal (Celltype1)'),
        ]):
            ax = axes[row_idx, col_idx] if len(key_targets) > 1 else axes[col_idx]

            # Load scatter data
            filename = f'{atlas}_{level}_cytosig.json'
            scatter_data = load_scatter(SCATTER_CT, filename)

            if scatter_data is None:
                # Try alternate filename patterns
                alt_names = [f'{atlas}_l1_cytosig.json', f'{atlas}_celltype_cytosig.json']
                for alt in alt_names:
                    scatter_data = load_scatter(SCATTER_CT, alt)
                    if scatter_data is not None:
                        break

            if scatter_data and target in scatter_data:
                data = scatter_data[target]
                points = data.get('points', [])
                rho = data.get('rho', np.nan)

                if points:
                    x = [p[0] for p in points]
                    y = [p[1] for p in points]

                    # Color by celltype if available
                    celltypes = data.get('celltypes', [])
                    if celltypes and len(points) > 0 and len(points[0]) > 2:
                        ct_indices = [p[2] if len(p) > 2 else 0 for p in points]
                        unique_cts = sorted(set(ct_indices))
                        cmap = plt.cm.Set3(np.linspace(0, 1, max(len(unique_cts), 1)))
                        for ct_idx in unique_cts:
                            mask = [i for i, c in enumerate(ct_indices) if c == ct_idx]
                            ct_name = celltypes[int(ct_idx)] if int(ct_idx) < len(celltypes) else f'CT{ct_idx}'
                            ax.scatter([x[i] for i in mask], [y[i] for i in mask],
                                      alpha=0.5, s=15, label=ct_name[:15], color=cmap[int(ct_idx) % len(cmap)])
                    else:
                        color = COLORS['good'] if rho > 0.2 else COLORS['bad'] if rho < 0 else COLORS['neutral']
                        ax.scatter(x, y, alpha=0.4, s=15, c=color)

                ax.set_title(f'{target} — {label}\nρ = {rho:.3f}', fontweight='bold', fontsize=10)
            else:
                ax.text(0.5, 0.5, f'{target}\nNo data', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f'{target} — {label}', fontsize=10)

            ax.set_xlabel('Mean Expression')
            ax.set_ylabel('Predicted Activity')

    fig.suptitle('Cell-Type-Level Validation: Key Cytokine Targets\n(Each point = celltype × donor pseudobulk)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG_DIR / 'fig14_celltype_scatter_examples.png')
    fig.savefig(FIG_DIR / 'fig14_celltype_scatter_examples.pdf')
    plt.close(fig)
    print('  ✓ Figure 14: Celltype-level scatter examples')


# ═══════════════════════════════════════════════════════════════════════════════
# HELPER: Independence-corrected ρ per target
# ═══════════════════════════════════════════════════════════════════════════════
def get_independent_rhos(df, atlas, sig_type):
    """Get one ρ per target at the independence-corrected level.

    For GTEx/TCGA: median across tissues/cancers (median-of-medians).
    For others: direct value at donor_only/tumor_only.
    """
    level, needs_mom = INDEPENDENT_LEVELS.get(atlas, ('donor_only', False))
    sub = df[(df['atlas'] == atlas) & (df['level'] == level) & (df['signature'] == sig_type)]
    if needs_mom:
        sub = sub[~sub['celltype'].isin(['all', 'unmappable'])]
        return sub.groupby('target')['spearman_rho'].median()
    else:
        return sub.set_index('target')['spearman_rho']


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 4a: GTEx Per-Tissue Stratified Validation
# ═══════════════════════════════════════════════════════════════════════════════
def fig4a_gtex_per_tissue(df):
    """Horizontal dot plot: CytoSig vs SecAct median ρ per GTEx tissue.

    Data source:
        gtex_correlations.csv at by_tissue level.
    Method:
        For each tissue: compute median ρ across all CytoSig targets and
        all SecAct targets. Sort by SecAct median ρ descending. Show as
        paired dots (CytoSig blue, SecAct green) per tissue row.
    Output:
        fig4a_gtex_per_tissue.png, fig4a_gtex_per_tissue.pdf
    """
    gtex = df[(df['atlas'] == 'gtex') & (df['level'] == 'by_tissue')]
    gtex = gtex[~gtex['celltype'].isin(['all', 'unmappable'])]
    tissues = sorted(gtex['celltype'].unique())

    rows = []
    for tissue in tissues:
        t_data = gtex[gtex['celltype'] == tissue]
        for sig_type in ['cytosig', 'secact']:
            sub = t_data[t_data['signature'] == sig_type]
            if len(sub) > 0:
                rhos = sub['spearman_rho'].dropna()
                n_samples = sub['n_samples'].iloc[0] if 'n_samples' in sub.columns else 0
                rows.append({
                    'tissue': tissue, 'signature': sig_type,
                    'median_rho': rhos.median(), 'n_targets': len(rhos),
                    'n_samples': n_samples,
                })
    plot_df = pd.DataFrame(rows)

    # Sort by SecAct median
    secact_order = plot_df[plot_df['signature'] == 'secact'].sort_values('median_rho', ascending=True)
    tissue_order = secact_order['tissue'].tolist()

    fig, ax = plt.subplots(figsize=(10, max(8, len(tissue_order) * 0.35)))
    y_pos = {t: i for i, t in enumerate(tissue_order)}

    for sig_type, color, marker, label in [
        ('cytosig', COLORS['cytosig'], 'o', 'CytoSig'),
        ('secact', COLORS['secact'], 's', 'SecAct'),
    ]:
        sub = plot_df[plot_df['signature'] == sig_type]
        for _, row in sub.iterrows():
            if row['tissue'] in y_pos:
                ax.scatter(row['median_rho'], y_pos[row['tissue']],
                          color=color, marker=marker, s=50, alpha=0.8,
                          zorder=3, label=label if row['tissue'] == tissue_order[0] else '')

    ax.set_yticks(range(len(tissue_order)))
    ax.set_yticklabels(tissue_order, fontsize=8)
    ax.set_xlabel('Median Spearman ρ (across targets)')
    ax.set_title('GTEx Per-Tissue Validation: CytoSig vs SecAct\n(29 tissues, by_tissue level)',
                 fontsize=12, fontweight='bold')
    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax.legend(fontsize=9, loc='lower right')

    plt.tight_layout()
    fig.savefig(FIG_DIR / 'fig4a_gtex_per_tissue.png')
    fig.savefig(FIG_DIR / 'fig4a_gtex_per_tissue.pdf')
    plt.close(fig)
    print('  ✓ Figure 4a: GTEx per-tissue stratified validation')


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 4b: TCGA Per-Cancer Stratified Validation
# ═══════════════════════════════════════════════════════════════════════════════
def fig4b_tcga_per_cancer(df):
    """Horizontal dot plot: CytoSig vs SecAct median ρ per TCGA cancer type.

    Data source:
        tcga_correlations.csv at primary_by_cancer level.
    Method:
        Same as fig4a but for TCGA. 33 cancer types sorted by SecAct median.
    Output:
        fig4b_tcga_per_cancer.png, fig4b_tcga_per_cancer.pdf
    """
    tcga = df[(df['atlas'] == 'tcga') & (df['level'] == 'primary_by_cancer')]
    tcga = tcga[~tcga['celltype'].isin(['all', 'unmappable'])]
    cancers = sorted(tcga['celltype'].unique())

    rows = []
    for cancer in cancers:
        c_data = tcga[tcga['celltype'] == cancer]
        for sig_type in ['cytosig', 'secact']:
            sub = c_data[c_data['signature'] == sig_type]
            if len(sub) > 0:
                rhos = sub['spearman_rho'].dropna()
                n_samples = sub['n_samples'].iloc[0] if 'n_samples' in sub.columns else 0
                rows.append({
                    'cancer': cancer, 'signature': sig_type,
                    'median_rho': rhos.median(), 'n_targets': len(rhos),
                    'n_samples': n_samples,
                })
    plot_df = pd.DataFrame(rows)

    secact_order = plot_df[plot_df['signature'] == 'secact'].sort_values('median_rho', ascending=True)
    cancer_order = secact_order['cancer'].tolist()

    fig, ax = plt.subplots(figsize=(10, max(8, len(cancer_order) * 0.35)))
    y_pos = {c: i for i, c in enumerate(cancer_order)}

    for sig_type, color, marker, label in [
        ('cytosig', COLORS['cytosig'], 'o', 'CytoSig'),
        ('secact', COLORS['secact'], 's', 'SecAct'),
    ]:
        sub = plot_df[plot_df['signature'] == sig_type]
        for _, row in sub.iterrows():
            if row['cancer'] in y_pos:
                ax.scatter(row['median_rho'], y_pos[row['cancer']],
                          color=color, marker=marker, s=50, alpha=0.8,
                          zorder=3, label=label if row['cancer'] == cancer_order[0] else '')

    ax.set_yticks(range(len(cancer_order)))
    ax.set_yticklabels(cancer_order, fontsize=7)
    ax.set_xlabel('Median Spearman ρ (across targets)')
    ax.set_title('TCGA Per-Cancer Validation: CytoSig vs SecAct\n(33 cancer types, primary_by_cancer level)',
                 fontsize=12, fontweight='bold')
    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax.legend(fontsize=9, loc='lower right')

    plt.tight_layout()
    fig.savefig(FIG_DIR / 'fig4b_tcga_per_cancer.png')
    fig.savefig(FIG_DIR / 'fig4b_tcga_per_cancer.pdf')
    plt.close(fig)
    print('  ✓ Figure 4b: TCGA per-cancer stratified validation')


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 4c: Cross-Platform Concordance
# ═══════════════════════════════════════════════════════════════════════════════
def fig4c_cross_platform_concordance(df):
    """Scatter: GTEx vs scAtlas Normal per-tissue ρ concordance.

    Data source:
        gtex_correlations.csv (by_tissue), scatlas_normal_correlations.csv
        (donor_organ). Tissues matched by name mapping.
    Method:
        For matching tissues in GTEx and scAtlas Normal, compute per-tissue
        median CytoSig ρ. Scatter GTEx ρ (x) vs scAtlas ρ (y). Compute
        Spearman correlation of the ρ values across tissues.
    Output:
        fig4c_cross_platform_concordance.png, fig4c_cross_platform_concordance.pdf
    """
    # Tissue/organ mappings (GTEx tissue → scAtlas organ)
    tissue_map_normal = {
        'Blood': 'Blood', 'Breast': 'Breast', 'Colon': 'Colon',
        'Esophagus': 'Esophagus', 'Heart': 'Heart', 'Kidney': 'Kidney',
        'Liver': 'Liver', 'Lung': 'Lung', 'Ovary': 'Ovary',
        'Skin': 'Skin', 'Spleen': 'Spleen',
    }

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: GTEx vs scAtlas Normal
    ax = axes[0]
    gtex_bt = df[(df['atlas'] == 'gtex') & (df['level'] == 'by_tissue') & (df['signature'] == 'cytosig')]
    sn_do = df[(df['atlas'] == 'scatlas_normal') & (df['level'] == 'donor_organ') & (df['signature'] == 'cytosig')]

    gtex_rhos, sn_rhos, labels = [], [], []
    for gtex_tissue, sn_organ in tissue_map_normal.items():
        g_sub = gtex_bt[gtex_bt['celltype'] == gtex_tissue]
        s_sub = sn_do[sn_do['celltype'] == sn_organ]
        if len(g_sub) > 0 and len(s_sub) > 0:
            gtex_rhos.append(g_sub['spearman_rho'].median())
            sn_rhos.append(s_sub['spearman_rho'].median())
            labels.append(gtex_tissue)

    if len(gtex_rhos) >= 3:
        rho_concordance, p_conc = stats.spearmanr(gtex_rhos, sn_rhos)
        ax.scatter(gtex_rhos, sn_rhos, c=COLORS['cytosig'], s=60, alpha=0.8, zorder=3)
        for i, label in enumerate(labels):
            ax.annotate(label, (gtex_rhos[i], sn_rhos[i]),
                       fontsize=7, xytext=(5, 5), textcoords='offset points')
        lim = [min(min(gtex_rhos), min(sn_rhos)) - 0.05,
               max(max(gtex_rhos), max(sn_rhos)) + 0.05]
        ax.plot(lim, lim, 'k--', alpha=0.4)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_xlabel('GTEx per-tissue median ρ (CytoSig)')
        ax.set_ylabel('scAtlas Normal per-organ median ρ (CytoSig)')
        ax.set_title(f'GTEx vs scAtlas Normal\n(n={len(labels)} tissues, '
                     f'ρ_concordance={rho_concordance:.3f}, p={p_conc:.3f})',
                     fontsize=10, fontweight='bold')
    else:
        ax.text(0.5, 0.5, 'Insufficient matching tissues', ha='center',
                va='center', transform=ax.transAxes)

    # Panel B: TCGA vs scAtlas Cancer
    ax = axes[1]
    cancer_map = {
        'Breast Invasive Carcinoma': 'BRCA',
        'Colon Adenocarcinoma': 'CRC',
        'Head & Neck Squamous Cell Carcinoma': 'HNSC',
        'Kidney Clear Cell Carcinoma': 'KIRC',
        'Liver Hepatocellular Carcinoma': 'HCC',
        'Lung Adenocarcinoma': 'LUAD',
        'Ovarian Serous Cystadenocarcinoma': 'OV',
        'Pancreatic Adenocarcinoma': 'PAAD',
    }

    tcga_pbc = df[(df['atlas'] == 'tcga') & (df['level'] == 'primary_by_cancer') & (df['signature'] == 'cytosig')]
    sc_tbc = df[(df['atlas'] == 'scatlas_cancer') & (df['level'] == 'tumor_by_cancer') & (df['signature'] == 'cytosig')]

    tcga_rhos, sc_rhos, cancer_labels = [], [], []
    for tcga_cancer, sc_cancer in cancer_map.items():
        t_sub = tcga_pbc[tcga_pbc['celltype'] == tcga_cancer]
        s_sub = sc_tbc[sc_tbc['celltype'] == sc_cancer]
        if len(t_sub) > 0 and len(s_sub) > 0:
            tcga_rhos.append(t_sub['spearman_rho'].median())
            sc_rhos.append(s_sub['spearman_rho'].median())
            cancer_labels.append(sc_cancer)

    if len(tcga_rhos) >= 3:
        rho_concordance, p_conc = stats.spearmanr(tcga_rhos, sc_rhos)
        ax.scatter(tcga_rhos, sc_rhos, c=COLORS['secact'], s=60, alpha=0.8, zorder=3)
        for i, label in enumerate(cancer_labels):
            ax.annotate(label, (tcga_rhos[i], sc_rhos[i]),
                       fontsize=7, xytext=(5, 5), textcoords='offset points')
        lim = [min(min(tcga_rhos), min(sc_rhos)) - 0.05,
               max(max(tcga_rhos), max(sc_rhos)) + 0.05]
        ax.plot(lim, lim, 'k--', alpha=0.4)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_xlabel('TCGA per-cancer median ρ (CytoSig)')
        ax.set_ylabel('scAtlas Cancer per-cancer median ρ (CytoSig)')
        ax.set_title(f'TCGA vs scAtlas Cancer\n(n={len(cancer_labels)} cancers, '
                     f'ρ_concordance={rho_concordance:.3f}, p={p_conc:.3f})',
                     fontsize=10, fontweight='bold')
    else:
        ax.text(0.5, 0.5, 'Insufficient matching cancers', ha='center',
                va='center', transform=ax.transAxes)

    fig.suptitle('Cross-Platform Concordance: Bulk RNA-seq vs Single-Cell Pseudobulk',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'fig4c_cross_platform_concordance.png')
    fig.savefig(FIG_DIR / 'fig4c_cross_platform_concordance.pdf')
    plt.close(fig)
    print('  ✓ Figure 4c: Cross-platform concordance')


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 (was 15): Summary Statistics Table
# ═══════════════════════════════════════════════════════════════════════════════
def fig15_summary_table(df):
    """Publication-ready summary statistics table.

    Data source:
        Merged correlation CSVs at donor level per atlas (6 atlases x
        4 signature types: CytoSig, LinCytoSig Best, LinCytoSig, SecAct).
    Method:
        For each atlas x signature: compute median, mean, std, min, max of
        spearman_rho. Count targets with spearman_pval < 0.05 for significance
        rate. Count rho > 0 for positive rate. Render as matplotlib table
        with colored rows by signature type. Also export as CSV.
    Output:
        fig15_summary_table.png, fig15_summary_table.pdf, summary_statistics.csv
    """
    level_map = {
        'cima': 'donor_only', 'inflammation': 'donor_only',
        'scatlas_normal': 'donor_organ', 'scatlas_cancer': 'donor_organ',
        'gtex': 'donor_only', 'tcga': 'donor_only',
    }

    rows = []
    for atlas in level_map:
        level = level_map[atlas]
        for sig_type in ['cytosig', 'lincyto_best', 'lincytosig', 'secact']:
            sub = df[(df['atlas'] == atlas) & (df['level'] == level) & (df['signature'] == sig_type)]
            if len(sub) == 0:
                continue
            rhos = sub['spearman_rho'].dropna()
            n_sig = (sub['spearman_pval'] < 0.05).sum()
            n_pos = (rhos > 0).sum()
            rows.append({
                'Atlas': atlas.replace('_', ' ').title(),
                'Signature': SIG_DISPLAY.get(sig_type, sig_type.upper()).replace('\n', ' '),
                'N Targets': len(rhos),
                'Median ρ': f'{rhos.median():.3f}',
                'Mean ρ': f'{rhos.mean():.3f}',
                'Std ρ': f'{rhos.std():.3f}',
                'Min ρ': f'{rhos.min():.3f}',
                'Max ρ': f'{rhos.max():.3f}',
                '% Significant': f'{100 * n_sig / len(sub):.1f}%',
                '% Positive': f'{100 * n_pos / len(rhos):.1f}%',
            })

    table_df = pd.DataFrame(rows)

    # Tight figure height: header + rows + padding
    row_height = 0.3
    fig_height = (len(rows) + 1) * row_height + 1.5
    fig, ax = plt.subplots(figsize=(18, fig_height))
    ax.axis('off')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    table = ax.table(cellText=table_df.values, colLabels=table_df.columns,
                     cellLoc='center', loc='upper center')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.4)  # Make rows taller for readability
    table.auto_set_column_width(col=list(range(len(table_df.columns))))

    # Color header
    for j in range(len(table_df.columns)):
        table[0, j].set_facecolor('#374151')
        table[0, j].set_text_props(color='white', fontweight='bold')

    # Color rows by signature type
    for i in range(len(rows)):
        sig = rows[i]['Signature']
        if sig == 'CytoSig':
            bg = '#DBEAFE'
        elif sig == 'LinCytoSig Best (comb+filt)':
            bg = '#FDE8D0'
        elif sig == 'LinCytoSig':
            bg = '#FEF3C7'
        elif sig == 'SecAct':
            bg = '#D1FAE5'
        else:
            bg = '#F3F4F6'
        for j in range(len(table_df.columns)):
            table[i + 1, j].set_facecolor(bg)

    ax.set_title('Cross-Sample Validation Summary Statistics\n(Donor-Level Pseudobulk Spearman Correlation)',
                 fontsize=14, fontweight='bold', pad=20)
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'fig2_summary_table.png')
    fig.savefig(FIG_DIR / 'fig2_summary_table.pdf')
    plt.close(fig)

    # Also save as CSV
    table_df.to_csv(FIG_DIR / 'summary_statistics.csv', index=False)
    print('  ✓ Figure 2: Summary statistics table')


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════
def main():
    print('CytoAtlas PI Report — Figure Generation')
    print('=' * 60)

    print('\nLoading correlation data...')
    df = load_all_correlations()
    print(f'  Loaded {len(df)} correlation records')
    df = merge_inflammation_atlases(df)
    print(f'  After inflammation merge: {len(df)} records')

    lincyto_best = load_lincyto_best_data()
    if len(lincyto_best) > 0:
        lincyto_best = merge_inflammation_atlases(lincyto_best)
        df = pd.concat([df, lincyto_best], ignore_index=True)
        df = df.drop_duplicates(subset=['target', 'celltype', 'atlas', 'level', 'signature'])
        print(f'  After lincyto_best merge: {len(df)} records')

    print(f'  Atlases: {df["atlas"].nunique()} | Signatures: {df["signature"].nunique()} | Targets: {df["target"].nunique()}')

    print('\nGenerating figures...')
    # Section 1
    fig1_dataset_overview()
    # Section 4: Validation Results
    fig15_summary_table(df)           # Now Fig 2: Summary statistics table
    fig2_correlation_summary(df)      # Now Fig 3: Cross-dataset boxplot
    fig4a_gtex_per_tissue(df)         # Fig 4a: GTEx per-tissue (NEW)
    fig4b_tcga_per_cancer(df)         # Fig 4b: TCGA per-cancer (NEW)
    fig4c_cross_platform_concordance(df)  # Fig 4c: Cross-platform (NEW)
    fig3_good_bad_correlations(df)    # Now Fig 5: Best/worst targets
    fig6_cross_atlas_consistency(df)  # Fig 6: Cross-atlas consistency
    fig7_validation_levels(df)        # Fig 7: Aggregation levels
    # Section 5: Method Comparison
    fig4_bio_targets_heatmap(df)      # Now Fig 8: Bio targets heatmap
    fig5_representative_scatter(df)   # Now Fig 9: Representative scatter
    fig8_method_comparison()          # Now Fig 10: 10-way comparison
    fig9_lincytosig_vs_cytosig()      # Now Fig 11: Lin vs Cyto scatter
    fig10_lincytosig_advantage()      # Now Fig 12: LinCytoSig advantage
    fig11_celltype_delta_rho()        # Now Fig 13: Celltype delta rho
    fig12_bulk_validation(df)         # Now Fig 14: Bulk validation
    fig13_lincytosig_specificity()    # Now Fig 15: LinCytoSig specificity
    fig14_celltype_scatter_examples(df)  # Now Fig 16: Celltype scatter

    print('\n' + '=' * 60)
    print(f'All figures saved to: {FIG_DIR}/')
    print('Done!')


if __name__ == '__main__':
    main()
