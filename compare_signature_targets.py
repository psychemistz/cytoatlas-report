#!/usr/bin/env python3
"""
Compare CytoSig and SecAct signature target lists,
then check their presence in GTEx and TCGA expression datasets.
"""

import h5py
import sys

# -- Paths ----------------------------------------------------------------
GTEX_DIR  = "/data/parks34/projects/2cytoatlas/results/cross_sample_validation/gtex"
TCGA_DIR  = "/data/parks34/projects/2cytoatlas/results/cross_sample_validation/tcga"

# -- Helper ----------------------------------------------------------------
def read_var_index(path):
    """Read var/_index from an h5ad file and return as a sorted set of strings."""
    with h5py.File(path, "r") as f:
        raw = f["var"]["_index"][:]
        return sorted(set(x.decode() if isinstance(x, bytes) else x for x in raw))

# -- 1. Load signature target lists ----------------------------------------
cytosig_targets = read_var_index(f"{GTEX_DIR}/gtex_by_tissue_cytosig.h5ad")
secact_targets  = read_var_index(f"{GTEX_DIR}/gtex_by_tissue_secact.h5ad")

# CytoSig target -> gene-symbol mapping (for expression lookup)
CYTOSIG_GENE_MAP = {
    "Activin A": "INHBA", "BDNF": "BDNF", "BMP2": "BMP2", "BMP4": "BMP4",
    "BMP6": "BMP6", "CD40L": "CD40LG", "CXCL12": "CXCL12", "EGF": "EGF",
    "FGF2": "FGF2", "GCSF": "CSF3", "GDF11": "GDF11", "GMCSF": "CSF2",
    "HGF": "HGF", "IFN1": "IFNA1", "IFNG": "IFNG", "IFNL": "IFNL1",
    "IL10": "IL10", "IL12": "IL12A", "IL13": "IL13", "IL15": "IL15",
    "IL17A": "IL17A", "IL1A": "IL1A", "IL1B": "IL1B", "IL2": "IL2",
    "IL21": "IL21", "IL22": "IL22", "IL27": "IL27", "IL3": "IL3",
    "IL36": "IL36A", "IL4": "IL4", "IL6": "IL6", "LIF": "LIF",
    "LTA": "LTA", "MCSF": "CSF1", "NO": "NOS1", "OSM": "OSM",
    "TGFB1": "TGFB1", "TGFB3": "TGFB3", "TNFA": "TNF", "TRAIL": "TNFSF10",
    "TWEAK": "TNFSF12", "VEGFA": "VEGFA", "WNT3A": "WNT3A",
}

cytosig_genes = set(CYTOSIG_GENE_MAP.values())   # gene symbols for the 43 targets
secact_genes  = set(secact_targets)               # already gene symbols

# -- 2. Load expression gene lists -----------------------------------------
gtex_genes = set(read_var_index(f"{GTEX_DIR}/gtex_by_tissue_expression.h5ad"))
tcga_genes = set(read_var_index(f"{TCGA_DIR}/tcga_by_cancer_expression.h5ad"))

# -- 3. Compute overlaps ---------------------------------------------------
shared_genes = cytosig_genes & secact_genes   # gene-symbol level overlap

cytosig_in_gtex = cytosig_genes & gtex_genes
cytosig_in_tcga = cytosig_genes & tcga_genes
secact_in_gtex  = secact_genes  & gtex_genes
secact_in_tcga  = secact_genes  & tcga_genes
shared_in_gtex  = shared_genes  & gtex_genes
shared_in_tcga  = shared_genes  & tcga_genes

# -- 4. Report --------------------------------------------------------------
sep = "=" * 72

print(sep)
print("SIGNATURE TARGET COMPARISON: CytoSig vs SecAct")
print(sep)

print(f"\n{'Metric':<50} {'Count':>6}")
print("-" * 58)
print(f"{'CytoSig targets (cytokine names)':<50} {len(cytosig_targets):>6}")
print(f"{'CytoSig unique gene symbols':<50} {len(cytosig_genes):>6}")
print(f"{'SecAct targets (gene symbols)':<50} {len(secact_genes):>6}")
print(f"{'GTEx expression genes':<50} {len(gtex_genes):>6}")
print(f"{'TCGA expression genes':<50} {len(tcga_genes):>6}")

print(f"\n{sep}")
print("OVERLAP: CytoSig (gene symbols) intersect SecAct")
print(sep)
print(f"{'Shared targets (CytoSig genes  intersect  SecAct genes)':<50} {len(shared_genes):>6}")
print(f"\nShared targets (sorted):")
for g in sorted(shared_genes):
    # find the CytoSig name(s) for this gene
    cs_names = [k for k, v in CYTOSIG_GENE_MAP.items() if v == g]
    label = f" (CytoSig: {', '.join(cs_names)})" if cs_names and cs_names[0] != g else ""
    print(f"  {g}{label}")

print(f"\n{sep}")
print("CytoSig-only targets (NOT in SecAct)")
print(sep)
cytosig_only = cytosig_genes - secact_genes
print(f"Count: {len(cytosig_only)}")
for g in sorted(cytosig_only):
    cs_names = [k for k, v in CYTOSIG_GENE_MAP.items() if v == g]
    label = f" (CytoSig: {', '.join(cs_names)})" if cs_names and cs_names[0] != g else ""
    print(f"  {g}{label}")

print(f"\n{sep}")
print("PRESENCE IN BULK EXPRESSION DATASETS")
print(sep)
print(f"\n{'Metric':<50} {'Count':>6} {'  / Total':>10} {'   %':>8}")
print("-" * 76)
print(f"{'CytoSig gene symbols in GTEx':<50} {len(cytosig_in_gtex):>6} {'  / ' + str(len(cytosig_genes)):>10} {100*len(cytosig_in_gtex)/len(cytosig_genes):>7.1f}%")
print(f"{'CytoSig gene symbols in TCGA':<50} {len(cytosig_in_tcga):>6} {'  / ' + str(len(cytosig_genes)):>10} {100*len(cytosig_in_tcga)/len(cytosig_genes):>7.1f}%")
print(f"{'SecAct targets in GTEx':<50} {len(secact_in_gtex):>6} {'  / ' + str(len(secact_genes)):>10} {100*len(secact_in_gtex)/len(secact_genes):>7.1f}%")
print(f"{'SecAct targets in TCGA':<50} {len(secact_in_tcga):>6} {'  / ' + str(len(secact_genes)):>10} {100*len(secact_in_tcga)/len(secact_genes):>7.1f}%")
print(f"{'Shared (CytoSig  intersect  SecAct) in GTEx':<50} {len(shared_in_gtex):>6} {'  / ' + str(len(shared_genes)):>10} {100*len(shared_in_gtex)/len(shared_genes):>7.1f}%")
print(f"{'Shared (CytoSig  intersect  SecAct) in TCGA':<50} {len(shared_in_tcga):>6} {'  / ' + str(len(shared_genes)):>10} {100*len(shared_in_tcga)/len(shared_genes):>7.1f}%")

# -- 5. Missing genes detail -----------------------------------------------
print(f"\n{sep}")
print("MISSING CytoSig gene symbols (not in expression data)")
print(sep)
cs_missing_gtex = cytosig_genes - gtex_genes
cs_missing_tcga = cytosig_genes - tcga_genes
if cs_missing_gtex:
    print(f"\n  Missing from GTEx ({len(cs_missing_gtex)}):")
    for g in sorted(cs_missing_gtex):
        cs_names = [k for k, v in CYTOSIG_GENE_MAP.items() if v == g]
        label = f" (CytoSig: {', '.join(cs_names)})" if cs_names and cs_names[0] != g else ""
        print(f"    {g}{label}")
else:
    print("\n  All CytoSig genes found in GTEx.")
if cs_missing_tcga:
    print(f"\n  Missing from TCGA ({len(cs_missing_tcga)}):")
    for g in sorted(cs_missing_tcga):
        cs_names = [k for k, v in CYTOSIG_GENE_MAP.items() if v == g]
        label = f" (CytoSig: {', '.join(cs_names)})" if cs_names and cs_names[0] != g else ""
        print(f"    {g}{label}")
else:
    print("\n  All CytoSig genes found in TCGA.")

print(f"\n{sep}")
print("MISSING SecAct targets (not in expression data)")
print(sep)
sa_missing_gtex = secact_genes - gtex_genes
sa_missing_tcga = secact_genes - tcga_genes
if sa_missing_gtex:
    print(f"\n  Missing from GTEx ({len(sa_missing_gtex)}):")
    for g in sorted(sa_missing_gtex):
        print(f"    {g}")
else:
    print("\n  All SecAct genes found in GTEx.")
if sa_missing_tcga:
    print(f"\n  Missing from TCGA ({len(sa_missing_tcga)}):")
    for g in sorted(sa_missing_tcga):
        print(f"    {g}")
else:
    print("\n  All SecAct genes found in TCGA.")

# -- 6. Shared missing detail ----------------------------------------------
shared_missing_gtex = shared_genes - gtex_genes
shared_missing_tcga = shared_genes - tcga_genes
if shared_missing_gtex or shared_missing_tcga:
    print(f"\n{sep}")
    print("MISSING shared (CytoSig intersect SecAct) targets in expression data")
    print(sep)
    if shared_missing_gtex:
        print(f"\n  Missing from GTEx ({len(shared_missing_gtex)}):")
        for g in sorted(shared_missing_gtex):
            print(f"    {g}")
    if shared_missing_tcga:
        print(f"\n  Missing from TCGA ({len(shared_missing_tcga)}):")
        for g in sorted(shared_missing_tcga):
            print(f"    {g}")

print(f"\n{sep}")
print("Done.")
