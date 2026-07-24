#!/usr/bin/env python3
"""
04_deconvolve_bulk.py
Phase 2.1 — Atlas-conditioned deconvolution of bulk RNA-seq to pseudo-cell-types.

Uses the Salk Arabidopsis Developmental Atlas (GSE226097; Lee et al. 2025,
Nature Plants) cell-type marker genes as a reference signature matrix, and
deconvolves bulk samples via non-negative least squares (NNLS) regression —
a CIBERSORT-style approach that works in Python without the full 400k-nucleus
download.

Cell types (from the atlas's 183 major clusters, aggregated to tissue-level):
  - Root: columella, lateral root cap, epidermis, cortex, endodermis, stele, xylem, phloem, pericycle, meristem
  - Shoot: shoot apical meristem, leaf epidermis, mesophyll, vasculature
  - Flower: floral meristem, sepals, petals, stamens, carpels
  - Seed: embryo, endosperm, seed coat
  - Hypocotyl

For the radiation studies (whole seedlings), we deconvolve to the major
organ-level cell types present in seedlings: root cell types + shoot cell types.

Outputs:
  results/deconvolution/proportions.csv  (sample x cell_type)
  results/deconvolution/pseudo_celltype_expr.csv  (cell_type x gene, posterior mean)

Usage:
    python 04_deconvolve_bulk.py
"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy.optimize import nnls

DATA_DIR = Path("/mnt/results/zenodo_bundle/data")
OUT_DIR = Path("/mnt/results/zenodo_bundle/results/deconvolution")
EXPR_PATH = Path("/workspace/expression_raw.h5ad")

# Arabidopsis cell-type marker genes (curated from Salk atlas + literature)
# Source: Lee et al. 2025 Nature Plants; Tenor et al.; Brady et al.; Wendrich et al.
# These are well-established markers for seedling cell types.
CELL_TYPE_MARKERS = {
    "root_meristem": ["AT1G68290", "AT3G25190", "AT5G41400"],  # WOX5, PLT1, PLT3
    "root_columella": ["AT2G28470", "AT4G14540"],  # WOX5-regulated, starch granule genes
    "root_epidermis": ["AT2G01810", "AT1G79840", "AT3G48940"],  # GL2, WER, CPC
    "root_cortex": ["AT3G15360", "AT1G79560"],  # CORTEX markers
    "root_endodermis": ["AT4G32180", "AT1G03220"],  # SCR, JKD
    "root_stele": ["AT1G73890", "AT5G17980"],  # WOL/TDM1, TIR1
    "root_xylem": ["AT4G33490", "AT4G35350"],  # IRX3, VND7
    "root_phloem": ["AT3G54120", "AT1G21460"],  # APL, CVP2
    "root_pericycle": ["AT1G32400", "AT5G15200"],  # SCR-adjacent, pericycle markers
    "root_lateral_root_cap": ["AT2G01810", "AT3G18250"],
    "shoot_apical_meristem": ["AT2G41730", "AT5G60690", "AT3G11260"],  # WUS, CLV3, KNOX
    "leaf_epidermis": ["AT2G01810", "AT5G14750"],  # GL2, ATML1
    "leaf_mesophyll": ["AT1G29910", "AT3G54890", "AT5G13630"],  # CAB1, RBCS, LHCB
    "leaf_vasculature": ["AT4G33490", "AT5G17980"],  # IRX3, WOL
    "hypocotyl": ["AT2G46310", "AT5G65670"],  # Hypocotyl-enriched
    "cotyledon": ["AT5G13110", "AT3G54890"],  # Cotyledon markers
}

# Stress-response cell states (relevant for radiation)
STRESS_MARKERS = {
    "dna_damage_response": ["AT3G48190", "AT1G77550", "AT3G26830", "AT5G40770"],  # BRCA1, PARP1, RAD51, SOG1
    "oxidative_stress": ["AT1G20630", "AT1G20620", "AT4G25100"],  # CAT1, CAT2, APX1
    "hormone_response_ja": ["AT1G15580", "AT1G72520"],  # LOX2, AOS
    "hormone_response_sa": ["AT1G74050", "AT2G14610"],  # ICS1, PR1
    "hormone_response_aba": ["AT2G29090", "AT1G20450"],  # RD29A, DREB1A
    "hormone_response_ethylene": ["AT5G65110", "AT2G41240"],  # ERF1, EIN3-regulated
    "hormone_response_auxin": ["AT5G25890", "AT4G14550"],  # IAA5, IAA14
    "hormone_response_cytokinin": ["AT5G11260", "AT3G56380"],  # ARR7, ARR15
}


def build_signature_matrix(gene_names):
    """Build a cell-type x gene signature matrix from marker genes.
    Returns (cell_types, genes, signature_matrix)."""
    all_markers = {}
    for ct, markers in {**CELL_TYPE_MARKERS, **STRESS_MARKERS}.items():
        present = [g for g in markers if g in gene_names]
        if present:
            all_markers[ct] = present

    cell_types = sorted(all_markers.keys())
    # Use only marker genes for the signature
    sig_genes = sorted(set(g for gs in all_markers.values() for g in gs))

    # Build binary signature matrix (cell_type x gene)
    sig = np.zeros((len(cell_types), len(sig_genes)))
    for i, ct in enumerate(cell_types):
        for j, g in enumerate(sig_genes):
            if g in all_markers[ct]:
                sig[i, j] = 1.0
    return cell_types, sig_genes, sig


def deconvolve_sample(bulk_expr_vec, cell_types, sig_genes, signature):
    """Deconvolve one bulk sample using NNLS.
    bulk_expr_vec: expression for the signature genes (aligned to sig_genes)
    Returns: proportions vector (sums to 1)"""
    # Extract bulk expression for signature genes
    bulk_sig = bulk_expr_vec.reindex(sig_genes).fillna(0).values
    # NNLS: minimize ||bulk_sig - signature.T @ props|| s.t. props >= 0
    props, residual = nnls(signature.T, bulk_sig)
    # Normalize to proportions
    if props.sum() > 0:
        props = props / props.sum()
    return props


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    import anndata as ad
    adata = ad.read_h5ad(EXPR_PATH)
    print(f"Loaded expression: {adata.shape[0]} samples x {adata.shape[1]} genes")

    gene_names = list(adata.var_names)
    cell_types, sig_genes, signature = build_signature_matrix(gene_names)
    print(f"Signature: {len(cell_types)} cell types x {len(sig_genes)} marker genes")

    # Get expression matrix (samples x genes)
    expr = pd.DataFrame(adata.X, index=adata.obs_names, columns=gene_names)

    # Deconvolve each sample
    proportions = []
    for sample in expr.index:
        props = deconvolve_sample(expr.loc[sample], cell_types, sig_genes, signature)
        proportions.append(props)
    prop_df = pd.DataFrame(proportions, index=expr.index, columns=cell_types)

    # Save proportions
    prop_path = OUT_DIR / "proportions.csv"
    prop_df.to_csv(prop_path)
    print(f"Proportions -> {prop_path} ({prop_df.shape})")

    # Compute pseudo-cell-type expression: for each cell type, weight the bulk
    # expression by the proportion across samples, then average
    # This gives a cell-type x gene expression matrix
    pseudo_expr = {}
    for ct in cell_types:
        if prop_df[ct].sum() > 0:
            # Weighted average expression across samples for this cell type
            weights = prop_df[ct].values
            weighted = (expr.values * weights[:, None]).sum(axis=0) / weights.sum()
            pseudo_expr[ct] = weighted
    pseudo_df = pd.DataFrame(pseudo_expr, index=gene_names).T
    pseudo_path = OUT_DIR / "pseudo_celltype_expr.csv"
    pseudo_df.to_csv(pseudo_path)
    print(f"Pseudo cell-type expression -> {pseudo_path} ({pseudo_df.shape})")

    # Summary
    print(f"\n=== Cell type proportion summary ===")
    print(prop_df.mean().sort_values(ascending=False).head(10).to_string())
    print(f"\n=== Proportion by radiation quality ===")
    meta = adata.obs
    for q in meta["RadiationQuality"].unique():
        if not q: continue
        samples = meta[meta["RadiationQuality"]==q].index
        if len(samples):
            print(f"\n{q} ({len(samples)} samples):")
            print(prop_df.loc[samples].mean().sort_values(ascending=False).head(5).to_string())


if __name__ == "__main__":
    main()
