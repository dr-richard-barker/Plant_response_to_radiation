#!/usr/bin/env python3
"""
08_radiation_resilience_index.py
Phase 3.2 — Radiation Resilience Index (RRI): a dynamic, time-resolved
composite metric quantifying how tissue stability degrades or recovers after
radiation exposure.

DESIGN (per the user's "you decide" directive)
  RRI is a composite of three components, each normalized to [0,1] where
  higher = more resilient (closer to unperturbed state):

  1. Latent-distance component (PRIMARY, weight 0.50):
     Distance of each sample's GP-AE latent embedding from the control
     centroid in latent space. Samples close to the control manifold are
     resilient; samples far away are disrupted.
     RRI_latent = exp(-||z_sample - z_control_mean||^2 / sigma^2)
     where sigma = median pairwise distance among controls.

  2. Pathway-balance component (SECONDARY, weight 0.25):
     Balance between stress-response pathway activations. A balanced
     response (DNA repair, oxidative stress, hormone signaling all
     moderately activated) is more resilient than a single pathway
     dominating. We compute the Shannon evenness of pathway scores.
     RRI_pathway = evenness(pathway_scores) = H / log(n_pathways)
     where H = -sum(p_i * log(p_i)) over pathway activation proportions.

  3. Module-preservation component (SECONDARY, weight 0.25):
     How well-preserved are the sustained (housekeeping) WGCNA modules
     under radiation. High preservation = resilient; low = disrupted.
     Measured as the correlation between the sample's expression of
     sustained-module genes and the control mean.
     RRI_module = max(0, cor(expr_sustained[sample], ctrl_mean_sustained))

  RRI = 0.50 * RRI_latent + 0.25 * RRI_pathway + 0.25 * RRI_module

  The RRI is computed per sample and then averaged per timepoint bin to
  show the dynamic degradation/recovery trajectory.

OUTPUTS
  - rri_per_sample.csv: SampleID, RRI, RRI_latent, RRI_pathway, RRI_module,
    TimePostExposure_h, StudyID, RadiationQuality, AbsorbedDose_Gy
  - rri_per_timepoint.csv: time bin, mean RRI, std, n, per-quality breakdown
  - rri_summary.json: summary statistics and component weights
"""
import json, os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DATA_DIR = Path("/mnt/results/zenodo_bundle/data")
RESULTS_DIR = Path("/mnt/results/zenodo_bundle/results")
OUT_DIR = RESULTS_DIR / "rri"
WGCNA_DIR = RESULTS_DIR / "wgcna"
TRAJ_DIR = RESULTS_DIR / "trajectories"
DEVICE = None  # CPU-only for this script

# Component weights
W_LATENT = 0.50
W_PATHWAY = 0.25
W_MODULE = 0.25

# Pathway gene sets (curated Arabidopsis stress-response pathways)
# These are well-established marker genes for each pathway from the literature.
PATHWAY_GENES = {
    "dna_repair": [
        "AT1G07290",  # BRCA1
        "AT3G05120",  # PARP2
        "AT5G40820",  # RAD51
        "AT5G20850",  # BRCA2
        "AT3G48190",  # MRE11
        "AT1G65470",  # ATM
        "AT3G02680",  # ATR
        "AT4G19730",  # PARP1
        "AT5G24280",  # KU70
        "AT2G32750",  # KU80
    ],
    "oxidative_stress": [
        "AT1G02920",  # GST1
        "AT4G25130",  # SOD1 (CSD1)
        "AT5G18100",  # GPX6
        "AT3G10920",  # CAT2 (catalase)
        "AT1G20630",  # CAT1
        "AT2G28190",  # CSD2 (SOD)
        "AT3G26060",  # PER5 (peroxidase)
        "AT1G32940",  # GST6
        "AT5G03490",  # FSD1 (Fe SOD)
        "AT4G35090",  # APX1 (ascorbate peroxidase)
    ],
    "hormone_signaling": [
        "AT2G14610",  # PR1 (SA)
        "AT1G19670",  # ABA1
        "AT1G02820",  # BGL2 (PDF1.2-related)
        "AT5G44420",  # PDF1.2 (JA/ETH)
        "AT2G39300",  # JAZ1 (JA)
        "AT3G04720",  # JAZ3
        "AT1G13220",  # JAZ10
        "AT5G57050",  # ABI2 (ABA)
        "AT3G24220",  # ABA2
        "AT1G75800",  # ASA1 (auxin)
    ],
}


def load_latent_embeddings():
    """Load GP-AE latent embeddings with metadata."""
    z_df = pd.read_csv(TRAJ_DIR / "latent_embeddings_annotated.csv", index_col=0)
    return z_df


def compute_latent_component(z_df):
    """RRI_latent = exp(-||z - z_control_mean||^2 / sigma^2)
    where sigma = median pairwise distance among controls."""
    latent_cols = [c for c in z_df.columns if c.startswith("z")]
    Z = z_df[latent_cols].values

    # The annotated latent CSV already has SampleID, StudyID, RadiationQuality,
    # AbsorbedDose_Gy, TimePostExposure_h, LET_keV_um as columns.
    # We need IsControl from metadata.
    md = pd.read_csv(DATA_DIR / "metadata_master.csv")
    # z_df has SampleID as index. Merge needs SampleID as column, so reset index,
    # merge, then restore SampleID as index.
    z_df = z_df.reset_index()
    if "IsControl" not in z_df.columns:
        z_df = z_df.merge(md[["SampleID", "IsControl"]], on="SampleID", how="left")
    z_df = z_df.set_index("SampleID")
    # Fix IsControl (categorical with '', 'False', 'True')
    ic_str = z_df["IsControl"].astype(str)
    z_df["IsControl"] = ic_str.str.lower().isin(["true", "1", "yes"])
    z_df["TimePostExposure_h"] = pd.to_numeric(z_df["TimePostExposure_h"], errors="coerce")

    ctrl_mask = z_df["IsControl"].values
    if ctrl_mask.sum() == 0:
        # Fallback: use time=0 samples as controls
        ctrl_mask = (z_df["TimePostExposure_h"] == 0).values
    if ctrl_mask.sum() == 0:
        # Fallback: use samples closest to latent centroid
        centroid = Z.mean(axis=0)
        dists = np.linalg.norm(Z - centroid, axis=1)
        ctrl_mask = dists < np.median(dists)

    Z_ctrl = Z[ctrl_mask]
    z_ctrl_mean = Z_ctrl.mean(axis=0)

    # sigma = median pairwise distance among controls
    from scipy.spatial.distance import pdist
    if len(Z_ctrl) > 1:
        ctrl_dists = pdist(Z_ctrl)
        sigma = np.median(ctrl_dists)
    else:
        sigma = np.std(Z, axis=0).mean()
    sigma = max(sigma, 1e-6)

    # RRI_latent for all samples
    dists = np.linalg.norm(Z - z_ctrl_mean, axis=1)
    rri_latent = np.exp(-dists**2 / sigma**2)

    z_df["RRI_latent"] = rri_latent
    z_df["latent_distance"] = dists
    return z_df, sigma


def compute_pathway_component(z_df):
    """RRI_pathway = Shannon evenness of pathway activation scores.
    Evenness = H / log(n_pathways), where H = -sum(p_i * log(p_i))."""
    import anndata as ad
    from collections import defaultdict
    adata = ad.read_h5ad("/workspace/expression_raw.h5ad")
    expr = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
    obs_map = adata.obs["SampleID"].to_dict()
    sid_to_obs = defaultdict(list)
    for obs_name, sid in obs_map.items():
        sid_to_obs[sid].append(obs_name)

    # Compute pathway scores aligned to z_df index
    pathway_names = list(PATHWAY_GENES.keys())
    ps_matrix = np.zeros((len(z_df), len(pathway_names)))
    for i, sid in enumerate(z_df.index):
        obs_names = sid_to_obs.get(sid, [])
        if not obs_names or obs_names[0] not in expr.index:
            continue
        sample_expr = expr.loc[obs_names[0]]
        for j, pname in enumerate(pathway_names):
            genes = PATHWAY_GENES[pname]
            available = [g for g in genes if g in expr.columns]
            if len(available) > 0:
                ps_matrix[i, j] = float(np.log1p(sample_expr[available].values).mean())

    # Normalize each row to sum=1 (proportions)
    row_sums = ps_matrix.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    ps_props = ps_matrix / row_sums

    # Shannon evenness: H / log(n)
    n_pathways = ps_props.shape[1]
    with np.errstate(divide='ignore', invalid='ignore'):
        H = -np.sum(ps_props * np.log(ps_props), axis=1)
    H = np.nan_to_num(H, nan=0.0)
    evenness = H / np.log(n_pathways)

    z_df["RRI_pathway"] = evenness
    return z_df


def compute_module_component(z_df):
    """RRI_module = correlation of sample's sustained-module gene expression
    with the control mean. High correlation = preserved = resilient."""
    import anndata as ad
    from collections import defaultdict
    mods = pd.read_csv(WGCNA_DIR / "modules.csv")
    cls = pd.read_csv(WGCNA_DIR / "module_classification.csv")
    sustained_modules = cls[cls["Classification"] == "sustained"]["Module"].tolist()
    sustained_genes = mods[mods["Module"].isin(sustained_modules)]["Gene"].tolist()
    adata = ad.read_h5ad("/workspace/expression_raw.h5ad")
    expr = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
    sustained_genes = [g for g in sustained_genes if g in expr.columns]
    print(f"  Sustained-module genes for preservation: {len(sustained_genes)}")

    if len(sustained_genes) == 0:
        z_df["RRI_module"] = 0.5
        return z_df

    expr_sub = np.log1p(expr[sustained_genes].values.astype(float))
    obs = adata.obs
    ic_str = obs["IsControl"].astype(str)
    is_ctrl = ic_str.str.lower().isin(["true", "1", "yes"]).values
    if is_ctrl.sum() == 0:
        times = pd.to_numeric(obs["TimePostExposure_h"], errors="coerce").values
        is_ctrl = (times == 0)
    ctrl_mean = expr_sub[is_ctrl].mean(axis=0) if is_ctrl.sum() > 0 else expr_sub.mean(axis=0)

    # Per-obs_name correlation with control mean
    rri_by_obs = np.zeros(expr_sub.shape[0])
    for i in range(expr_sub.shape[0]):
        if np.std(expr_sub[i]) < 1e-8 or np.std(ctrl_mean) < 1e-8:
            rri_by_obs[i] = 0.5
        else:
            r = np.corrcoef(expr_sub[i], ctrl_mean)[0, 1]
            rri_by_obs[i] = max(0, r)
    rri_series = pd.Series(rri_by_obs, index=expr.index)

    # Map by SampleID (z_df index) -> obs_name -> RRI
    obs_map = adata.obs["SampleID"].to_dict()
    sid_to_obs = defaultdict(list)
    for obs_name, sid in obs_map.items():
        sid_to_obs[sid].append(obs_name)

    rri_module_mapped = []
    for sid in z_df.index:
        obs_names = sid_to_obs.get(sid, [])
        if obs_names and obs_names[0] in rri_series.index:
            rri_module_mapped.append(float(rri_series[obs_names[0]]))
        else:
            rri_module_mapped.append(0.5)
    z_df["RRI_module"] = rri_module_mapped
    return z_df


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("Loading GP-AE latent embeddings...")
    z_df = load_latent_embeddings()
    print(f"  Latent embeddings: {z_df.shape}")

    print("\nComputing RRI latent-distance component (weight={})...".format(W_LATENT))
    z_df, sigma = compute_latent_component(z_df)
    print(f"  Control sigma (latent distance scale): {sigma:.3f}")
    print(f"  RRI_latent range: [{z_df['RRI_latent'].min():.3f}, {z_df['RRI_latent'].max():.3f}]")

    print("\nComputing RRI pathway-balance component (weight={})...".format(W_PATHWAY))
    z_df = compute_pathway_component(z_df)
    print(f"  RRI_pathway range: [{z_df['RRI_pathway'].min():.3f}, {z_df['RRI_pathway'].max():.3f}]")

    print("\nComputing RRI module-preservation component (weight={})...".format(W_MODULE))
    z_df = compute_module_component(z_df)
    print(f"  RRI_module range: [{z_df['RRI_module'].min():.3f}, {z_df['RRI_module'].max():.3f}]")

    # Composite RRI
    z_df["RRI"] = (W_LATENT * z_df["RRI_latent"] +
                   W_PATHWAY * z_df["RRI_pathway"] +
                   W_MODULE * z_df["RRI_module"])

    # Save per-sample RRI
    out_cols = ["SampleID", "RRI", "RRI_latent", "RRI_pathway", "RRI_module",
                "latent_distance", "TimePostExposure_h", "StudyID",
                "RadiationQuality", "AbsorbedDose_Gy", "IsControl"]
    out_cols = [c for c in out_cols if c in z_df.columns]
    rri_per_sample = z_df[out_cols].copy()
    rri_path = OUT_DIR / "rri_per_sample.csv"
    rri_per_sample.to_csv(rri_path, index=False)
    print(f"\nRRI per sample -> {rri_path}  ({len(rri_per_sample)} samples)")

    # Per-timepoint aggregation (non-control samples only)
    rri_per_sample["TimePostExposure_h"] = pd.to_numeric(rri_per_sample["TimePostExposure_h"], errors="coerce")
    rri_per_sample["IsControl"] = rri_per_sample["IsControl"].astype(str).str.lower().isin(["true","1","yes"])
    non_ctrl = rri_per_sample[~rri_per_sample["IsControl"] & rri_per_sample["TimePostExposure_h"].notna()].copy()

    # Time bins
    bins = [0, 0.5, 2, 6, 12, 30, 100]
    labels = ["0-0.5h", "0.5-2h", "2-6h", "6-12h", "12-30h", "30-100h"]
    non_ctrl["time_bin"] = pd.cut(non_ctrl["TimePostExposure_h"], bins=bins, labels=labels, right=False)

    timepoint_rows = []
    for tbin in labels:
        sub = non_ctrl[non_ctrl["time_bin"] == tbin]
        if len(sub) == 0: continue
        row = {"time_bin": tbin, "n": len(sub),
               "RRI_mean": sub["RRI"].mean(), "RRI_std": sub["RRI"].std(),
               "RRI_latent_mean": sub["RRI_latent"].mean(),
               "RRI_pathway_mean": sub["RRI_pathway"].mean(),
               "RRI_module_mean": sub["RRI_module"].mean()}
        # Per-quality breakdown
        for q in sub["RadiationQuality"].unique():
            qsub = sub[sub["RadiationQuality"] == q]
            row[f"RRI_{q}_mean"] = qsub["RRI"].mean()
            row[f"n_{q}"] = len(qsub)
        timepoint_rows.append(row)
    tp_df = pd.DataFrame(timepoint_rows)
    tp_path = OUT_DIR / "rri_per_timepoint.csv"
    tp_df.to_csv(tp_path, index=False)
    print(f"RRI per timepoint -> {tp_path}")

    # Summary
    print("\n=== RRI trajectory (mean RRI per timepoint) ===")
    print(tp_df[["time_bin", "n", "RRI_mean", "RRI_std",
                 "RRI_latent_mean", "RRI_pathway_mean", "RRI_module_mean"]].to_string(index=False))

    # Count sustained genes
    mods_df = pd.read_csv(WGCNA_DIR / "modules.csv")
    cls_df = pd.read_csv(WGCNA_DIR / "module_classification.csv")
    sustained_mods = cls_df[cls_df["Classification"] == "sustained"]["Module"].tolist()
    n_sustained_genes = int(mods_df[mods_df["Module"].isin(sustained_mods)].shape[0])

    ctrl_rri = float(rri_per_sample[rri_per_sample["IsControl"]]["RRI"].mean()) if rri_per_sample["IsControl"].any() else None
    summary = {
        "weights": {"latent": W_LATENT, "pathway": W_PATHWAY, "module": W_MODULE},
        "n_samples": int(len(rri_per_sample)),
        "n_non_control_with_time": int(len(non_ctrl)),
        "rri_overall_mean": float(non_ctrl["RRI"].mean()),
        "rri_overall_std": float(non_ctrl["RRI"].std()),
        "rri_control_mean": ctrl_rri,
        "latent_sigma": float(sigma),
        "n_pathways": len(PATHWAY_GENES),
        "n_sustained_genes": n_sustained_genes
    }
    with open(OUT_DIR / "rri_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print(f"\nSummary -> {OUT_DIR / 'rri_summary.json'}")


if __name__ == "__main__":
    main()
