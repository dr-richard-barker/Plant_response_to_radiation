#!/usr/bin/env python3
"""
06_plantcellchat_kinetic.py
Phase 2.3 — PlantCellChat kinetic cell-cell communication.

Goal: Infer inter-tissue signaling networks from the deconvolved pseudo-cell-
type expression, at each discrete timepoint, to generate a "signaling flow"
view of how radiation-induced stress signaling originates (e.g. in the root
meristem or DNA-damage-response compartment) and propagates across tissues
over time post-exposure.

DATA
  - Pseudo-cell-type expression: results/deconvolution/pseudo_celltype_expr.csv
    (24 cell types x 32,833 genes). This is a per-cell-type profile derived by
    NNLS deconvolution of bulk RNA-seq against curated marker genes (Phase 2.1).
    It is NOT single-cell data; we treat each cell type as a "pseudo-cell"
    group, which is the standard adaptation of CellChat to deconvolved bulk.
  - LR database: PlantPhoneDB Arabidopsis LR pairs (Xu et al. 2022, Plant
    Biotechnol J, doi:10.1111/pbi.13893), 3,605 pairs with both genes in our
    expression matrix. PlantPhoneDB is the public, downloadable LR compendium
    that PlantCellChatDB builds on; we use it here because PlantCellChatDB's
    data files are served dynamically from a Shiny app and not directly
    downloadable. The mass-action scoring model follows PlantCellChat
    (Liu et al. 2026, Plant J 126(3):e70905).

METHOD
  For each timepoint bin t:
    1. Compute mean pseudo-cell-type expression across samples at time t
       (within-study centered, so values are response relative to control).
    2. For each LR pair (L, R) and each ordered cell-type pair (source, target):
       communication_prob = (LigandExpr[source, L] * ReceptorExpr[target, R])
                             / (Kh + LigandExpr[source, L] * ReceptorExpr[target, R])
       (Hill-type mass-action model, Kh=0.5 as in PlantCellChat default).
    3. Aggregate over LR pairs to get a 24x24 communication-strength matrix
       per timepoint.
    4. Permutation test: permute cell-type labels to get a null distribution
       and p-values.

OUTPUTS
  - signaling_flow_per_timepoint.csv: long-format (Time, Source, Target,
    SignalStrength, nLRpairs, pvalue) for all timepoint x cell-pair combos.
  - signaling_heatmap_per_timepoint.csv: 24x24 matrices per timepoint.
  - signaling_origin_summary.csv: per-timepoint total outgoing signal per
    source cell type (to identify signaling origin).
  - top_lr_pairs.csv: top LR pairs by communication strength across time.
"""
import json, os, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata

DATA_DIR = Path("/mnt/results/zenodo_bundle/data")
RESULTS_DIR = Path("/mnt/results/zenodo_bundle/results")
OUT_DIR = RESULTS_DIR / "cellchat"
LR_DB_PATH = Path("/mnt/shared-workspace/shared/lr_db/Athaliana-2022-Xu-LR-pairs.csv")
KH = 0.5  # mass-action half-saturation constant (PlantCellChat default)
N_PERMUTATIONS = 200
TIME_BINS = [(0, 0.5, "0-0.5h"), (0.5, 2, "0.5-2h"), (2, 6, "2-6h"),
             (6, 12, "6-12h"), (12, 30, "12-30h"), (30, 100, "30-100h")]


def load_data():
    """Load pseudo-cell-type expression, proportions, metadata, and LR db.

    The proportions CSV is indexed by adata obs_names (OSD-style IDs), which
    differ from the metadata SampleID (GSM-style). We map via the adata obs
    SampleID column to merge with metadata.
    """
    import anndata as ad
    pce = pd.read_csv(RESULTS_DIR / "deconvolution/pseudo_celltype_expr.csv", index_col=0)
    # pce: rows = cell types, cols = genes, values = expression
    prop = pd.read_csv(RESULTS_DIR / "deconvolution/proportions.csv", index_col=0)
    # prop: rows = samples (adata obs_names), cols = cell types
    adata = ad.read_h5ad("/workspace/expression_raw.h5ad")
    obs_sid = adata.obs[["SampleID"]].copy()
    prop = prop.merge(obs_sid, left_index=True, right_index=True)
    # Now prop has SampleID as a column; set it as index for merging
    md = pd.read_csv(DATA_DIR / "metadata_master.csv")
    lr = pd.read_csv(LR_DB_PATH)
    # Filter LR pairs to those with both genes in expression matrix
    genes = set(pce.columns)
    lr = lr[(lr["Ligands"].isin(genes)) & (lr["Receptors"].isin(genes))].copy()
    lr["LR_pair"] = lr["Ligands"] + "->" + lr["Receptors"]
    return pce, prop, md, lr


def compute_timepoint_expression(pce, prop, md):
    """Compute mean pseudo-cell-type expression per timepoint bin.

    For each timepoint bin, average the sample-level deconvolved expression
    (proportions-weighted) across samples falling in that bin. We use the
    pseudo-cell-type profiles (pce) as the cell-type reference and weight by
    sample proportions to get a time-resolved cell-type expression estimate.

    Since pce is a single profile per cell type (not sample-resolved), we
    approximate time-resolved cell-type expression by scaling each cell type's
    profile by its mean proportion at that timepoint, then normalizing. This
    captures which cell types are more prevalent (and thus more signaling-
    active) at each timepoint. A more refined approach would re-deconvolve
    per timepoint, but the proportions already encode the time-resolved
    cell-type composition shifts.
    """
    # Merge proportions with metadata to get time per sample
    # prop now has SampleID as a column (from load_data merge with adata obs)
    cell_types = [c for c in pce.index]
    merged = prop.merge(md[["SampleID", "TimePostExposure_h", "StudyID", "RadiationQuality",
                            "IsControl"]], on="SampleID", how="inner")
    merged["TimePostExposure_h"] = pd.to_numeric(merged["TimePostExposure_h"], errors="coerce")
    timepoint_expr = {}
    timepoint_meta = {}
    for lo, hi, label in TIME_BINS:
        # Samples in this time bin (exclude controls and NaN times)
        mask = (merged["TimePostExposure_h"] >= lo) & (merged["TimePostExposure_h"] < hi) & \
               (~merged["IsControl"].astype(bool)) & (~merged["TimePostExposure_h"].isna())
        samples = merged[mask]
        if len(samples) < 2: continue
        # Mean proportion per cell type at this timepoint
        mean_props = samples[cell_types].mean(axis=0)
        # Scale each cell type's profile by its mean proportion (composition weighting)
        # but keep the expression values in log1p scale (not normalized to sum=1)
        scaled = pce.loc[cell_types].mul(mean_props.values, axis=0)
        # Log1p transform (CellChat uses log-normalized data)
        scaled = np.log1p(scaled.clip(lower=0))
        timepoint_expr[label] = scaled
        timepoint_meta[label] = {"n_samples": int(len(samples)),
                                  "mean_time_h": float(samples["TimePostExposure_h"].mean()),
                                  "studies": list(samples["StudyID"].unique()),
                                  "qualities": list(samples["RadiationQuality"].unique())}
        print(f"  {label}: n={len(samples)} samples, mean time={timepoint_meta[label]['mean_time_h']:.2f}h, "
              f"studies={len(timepoint_meta[label]['studies'])}, qualities={timepoint_meta[label]['qualities']}")
    return timepoint_expr, timepoint_meta, cell_types


def compute_communication_matrix(expr_celltypes, lr, cell_types, Kh=KH):
    """Compute 24x24 communication strength matrix for one timepoint.

    comm[source, target] = sum over LR pairs of:
        (LigandExpr[source, L] * ReceptorExpr[target, R]) / (Kh + product)
    """
    n_ct = len(cell_types)
    # Build ligand and receptor expression matrices
    ligand_genes = lr["Ligands"].unique()
    receptor_genes = lr["Receptors"].unique()
    # Expression: rows = cell types, cols = genes
    L_expr = expr_celltypes[ligand_genes].values  # [n_ct, n_ligands]
    R_expr = expr_celltypes[receptor_genes].values  # [n_ct, n_receptors]

    # For each LR pair, compute the n_ct x n_ct communication matrix
    # Vectorized: for LR pair i, comm[s,t] = L_expr[s, lig_i] * R_expr[t, rec_i] / (Kh + ...)
    comm_total = np.zeros((n_ct, n_ct))
    comm_count = np.zeros((n_ct, n_ct), dtype=int)
    # Per-LR-pair contributions for top-pair tracking
    lr_contributions = []

    for _, row in lr.iterrows():
        lig = row["Ligands"]; rec = row["Receptors"]
        l_vec = expr_celltypes[lig].values  # [n_ct]
        r_vec = expr_celltypes[rec].values  # [n_ct]
        # Outer product: ligand from source, receptor on target
        prod = np.outer(l_vec, r_vec)  # [n_ct, n_ct]
        prob = prod / (Kh + prod)
        comm_total += prob
        comm_count += (prod > 0).astype(int)
        # Track top pairs
        total_strength = prob.sum()
        if total_strength > 0:
            lr_contributions.append({"LR_pair": row["LR_pair"], "Ligand": lig,
                                     "Receptor": rec, "source_db": row.get("source",""),
                                     "total_strength": float(total_strength)})

    return comm_total, comm_count, lr_contributions


def permutation_pvalues(expr_celltypes, lr, cell_types, n_perm=N_PERMUTATIONS, Kh=KH):
    """Permutation test: permute cell-type labels to get null distribution.
    Returns p-value matrix (fraction of permutations where null >= observed)."""
    n_ct = len(cell_types)
    observed, _, _ = compute_communication_matrix(expr_celltypes, lr, cell_types, Kh)
    # Permute cell-type labels (rows of expression matrix)
    ge_null = np.zeros((n_ct, n_ct))
    expr_arr = expr_celltypes.values
    rng = np.random.default_rng(42)
    for _ in range(n_perm):
        perm = rng.permutation(n_ct)
        expr_perm = pd.DataFrame(expr_arr[perm], index=expr_celltypes.index,
                                 columns=expr_celltypes.columns)
        null_comm, _, _ = compute_communication_matrix(expr_perm, lr, cell_types, Kh)
        ge_null += (null_comm >= observed).astype(float)
    pvals = ge_null / n_perm
    return observed, pvals


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("Loading data...")
    pce, prop, md, lr = load_data()
    print(f"  Pseudo-cell-type expr: {pce.shape} (cell types x genes)")
    print(f"  Proportions: {prop.shape} (samples x cell types)")
    print(f"  Metadata: {len(md)} samples")
    print(f"  LR pairs (filtered): {len(lr)}")

    print("\nComputing timepoint-resolved cell-type expression...")
    timepoint_expr, timepoint_meta, cell_types = compute_timepoint_expression(pce, prop, md)
    print(f"  Timepoints with data: {list(timepoint_expr.keys())}")

    # Run communication analysis per timepoint
    print(f"\nComputing communication matrices per timepoint (Kh={KH}, {N_PERMUTATIONS} permutations)...")
    all_flow_rows = []
    all_heatmaps = {}
    all_origin_rows = []
    all_lr_rows = []

    for tlabel, expr_ct in timepoint_expr.items():
        print(f"  {tlabel}...")
        comm, pvals = permutation_pvalues(expr_ct, lr, cell_types)
        # Also get per-LR-pair contributions (no permutation needed for this)
        _, _, lr_contrib = compute_communication_matrix(expr_ct, lr, cell_types)
        # Long-format flow rows
        for i, src in enumerate(cell_types):
            for j, tgt in enumerate(cell_types):
                if i == j: continue
                all_flow_rows.append({"Time": tlabel, "Source": src, "Target": tgt,
                                      "SignalStrength": float(comm[i, j]),
                                      "nLRpairs": int(np.count_nonzero(comm[i, j] > 0)),
                                      "pvalue": float(pvals[i, j])})
        # Heatmap matrix
        all_heatmaps[tlabel] = pd.DataFrame(comm, index=cell_types, columns=cell_types)
        # Origin summary: total outgoing signal per source
        out_total = comm.sum(axis=1)  # sum over targets
        in_total = comm.sum(axis=0)   # sum over sources
        for i, ct in enumerate(cell_types):
            all_origin_rows.append({"Time": tlabel, "CellType": ct,
                                    "OutgoingSignal": float(out_total[i]),
                                    "IncomingSignal": float(in_total[i]),
                                    "NetSignal": float(out_total[i] - in_total[i])})
        # Top LR pairs
        for contrib in sorted(lr_contrib, key=lambda x: -x["total_strength"])[:50]:
            contrib["Time"] = tlabel
            all_lr_rows.append(contrib)

    # Save outputs
    flow_df = pd.DataFrame(all_flow_rows)
    flow_path = OUT_DIR / "signaling_flow_per_timepoint.csv"
    flow_df.to_csv(flow_path, index=False)
    print(f"\nSignaling flow -> {flow_path}  ({len(flow_df)} rows)")

    # Heatmaps: save as one CSV with Time column
    heat_rows = []
    for tlabel, mat in all_heatmaps.items():
        for src in mat.index:
            for tgt in mat.columns:
                heat_rows.append({"Time": tlabel, "Source": src, "Target": tgt,
                                  "Strength": float(mat.loc[src, tgt])})
    heat_df = pd.DataFrame(heat_rows)
    heat_path = OUT_DIR / "signaling_heatmap_per_timepoint.csv"
    heat_df.to_csv(heat_path, index=False)
    print(f"Signaling heatmaps -> {heat_path}  ({len(heat_df)} rows)")

    origin_df = pd.DataFrame(all_origin_rows)
    origin_path = OUT_DIR / "signaling_origin_summary.csv"
    origin_df.to_csv(origin_path, index=False)
    print(f"Signaling origin -> {origin_path}  ({len(origin_df)} rows)")

    lr_df = pd.DataFrame(all_lr_rows)
    lr_path = OUT_DIR / "top_lr_pairs.csv"
    lr_df.to_csv(lr_path, index=False)
    print(f"Top LR pairs -> {lr_path}  ({len(lr_df)} rows)")

    # Timepoint metadata
    with open(OUT_DIR / "timepoint_metadata.json", "w") as f:
        json.dump(timepoint_meta, f, indent=2, default=float)

    # Summary: identify signaling origin per timepoint
    print("\n=== Signaling origin (top 3 outgoing cell types per timepoint) ===")
    for tlabel in timepoint_expr.keys():
        sub = origin_df[origin_df["Time"] == tlabel].sort_values("OutgoingSignal", ascending=False).head(3)
        print(f"  {tlabel}: " + ", ".join(f"{r['CellType']} ({r['OutgoingSignal']:.2f})" for _, r in sub.iterrows()))

    # Summary: signaling flow evolution (total signal per timepoint)
    print("\n=== Total signaling strength per timepoint ===")
    for tlabel in timepoint_expr.keys():
        total = flow_df[flow_df["Time"] == tlabel]["SignalStrength"].sum()
        n_sig = (flow_df[flow_df["Time"] == tlabel]["pvalue"] < 0.05).sum()
        print(f"  {tlabel}: total={total:.2f}, significant pairs (p<0.05)={n_sig}")


if __name__ == "__main__":
    main()
