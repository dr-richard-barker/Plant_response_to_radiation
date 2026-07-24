#!/usr/bin/env python3
"""
12_lr_pair_specificity.py
Phase 5 — Ligand-receptor pair specificity across radiation qualities.

Goal: Extend the CellChat signaling analysis to compute per-radiation-quality,
per-timepoint signaling networks for ALL 5 qualities (gamma, GCR, spaceflight-LEO,
HZE-Fe, UV-B), then test which LR pairs show quality-specific signaling patterns
using Kruskal-Wallis differential signaling tests.

KEY DISCOVERY: HZE-Fe (OSD-46/OSD-320) and UV-B (OSD-296) microarray data is
available on disk, already log-normalized, with sample IDs matching metadata.
This enables EMPIRICAL CellChat on all 5 qualities — no GP-AE decoder needed.

DATA
  RNA-seq qualities (in expression_arabidopsis_space.h5ad):
    - gamma: 108 samples, 5 timepoint bins
    - GCR: 14 samples, 1 timepoint (2-6h)
    - spaceflight-LEO: 53 samples, static (no time)

  Microarray qualities (in /mnt/shared-workspace/shared/raw/counts/):
    - HZE-Fe: OSD-46 (18 HZE-Fe samples) + OSD-320 (18 HZE-Fe samples), 4 timepoint bins
    - UV-B: OSD-296 (36 samples), 3 timepoint bins

METHOD
  1. Deconvolve microarray samples (NNLS, same markers as script 04)
  2. For each (quality, time_bin) with >=2 non-control samples:
     a. Compute mean pseudo-cell-type expression (proportions-weighted)
     b. Run mass-action CellChat scoring (Kh=0.5, 200 permutations)
     c. Track per-LR-pair signaling strength
  3. At shared timepoints, test LR pair specificity:
     a. Kruskal-Wallis across qualities per LR pair
     b. BH-FDR correction
     c. Classify: quality-specific (FDR<0.05, max/min>2) vs shared (FDR>0.1)
     d. Post-hoc pairwise Mann-Whitney U for enriched quality

OUTPUTS
  results/deconvolution/proportions_microarray.csv
  results/deconvolution/pseudo_celltype_expr_microarray.csv
  results/cellchat/signaling_flow_per_quality_timepoint.csv
  results/cellchat/lr_pair_strength_per_quality.csv
  results/cellchat/lr_pair_specificity.csv
  results/cellchat/quality_timepoint_metadata.json
"""
import json, os, warnings, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import kruskal, mannwhitneyu
from itertools import combinations

warnings.filterwarnings("ignore", category=FutureWarning)

DATA_DIR = Path("/mnt/results/zenodo_bundle/data")
RESULTS_DIR = Path("/mnt/results/zenodo_bundle/results")
RAW_DIR = Path("/mnt/shared-workspace/shared/raw/counts")
OUT_DIR = RESULTS_DIR / "cellchat"
DEV_DIR = RESULTS_DIR / "deconvolution"
LR_DB_PATH = Path("/mnt/shared-workspace/shared/lr_db/Athaliana-2022-Xu-LR-pairs.csv")

KH = 0.5
N_PERMUTATIONS = 200
TIME_BINS = [(0, 0.5, "0-0.5h"), (0.5, 2, "0.5-2h"), (2, 6, "2-6h"),
             (6, 12, "6-12h"), (12, 30, "12-30h"), (30, 100, "30h+")]

# ---- Marker genes (from script 04, with microarray-compatible fallbacks) ----
CELL_TYPE_MARKERS = {
    "root_meristem": ["AT1G68290", "AT3G25190", "AT5G41400"],
    "root_columella": ["AT2G28470", "AT4G14540"],
    "root_epidermis": ["AT2G01810", "AT1G79840", "AT3G48940"],
    "root_cortex": ["AT3G15360", "AT1G79560"],
    "root_endodermis": ["AT4G32180", "AT1G03220", "AT4G28400"],  # SCR, JKD, +MYB36 fallback
    "root_stele": ["AT1G73890", "AT5G17980", "AT3G62940"],  # WOL, TIR1, +TIR1 alt
    "root_xylem": ["AT4G33490", "AT4G35350"],
    "root_phloem": ["AT3G54120", "AT1G21460"],
    "root_pericycle": ["AT1G32400", "AT5G15200"],
    "root_lateral_root_cap": ["AT2G01810", "AT3G18250"],
    "shoot_apical_meristem": ["AT2G41730", "AT5G60690", "AT3G11260"],
    "leaf_epidermis": ["AT2G01810", "AT5G14750"],
    "leaf_mesophyll": ["AT1G29910", "AT3G54890", "AT5G13630", "AT1G61550"],  # CAB1, RBCS, LHCB, +RBCS alt
    "leaf_vasculature": ["AT4G33490", "AT5G17980"],
    "hypocotyl": ["AT2G46310", "AT5G65670"],
    "cotyledon": ["AT5G13110", "AT3G54890"],
}

STRESS_MARKERS = {
    "dna_damage_response": ["AT3G48190", "AT1G77550", "AT3G26830", "AT5G40770"],
    "oxidative_stress": ["AT1G20630", "AT1G20620", "AT4G25100"],
    "hormone_response_ja": ["AT1G15580", "AT1G72520"],
    "hormone_response_sa": ["AT1G74050", "AT2G14610"],
    "hormone_response_aba": ["AT2G29090", "AT1G20450"],
    "hormone_response_ethylene": ["AT5G65110", "AT2G41240"],
    "hormone_response_auxin": ["AT5G25890", "AT4G14550"],
    "hormone_response_cytokinin": ["AT5G11260", "AT3G56380"],
}


def build_signature_matrix(gene_names):
    """Build cell-type x gene signature matrix from marker genes."""
    all_markers = {}
    for ct, markers in {**CELL_TYPE_MARKERS, **STRESS_MARKERS}.items():
        present = [g for g in markers if g in gene_names]
        if len(present) >= 1:
            all_markers[ct] = present
    cell_types = sorted(all_markers.keys())
    sig_genes = sorted(set(g for gs in all_markers.values() for g in gs))
    sig = np.zeros((len(cell_types), len(sig_genes)))
    for i, ct in enumerate(cell_types):
        for j, g in enumerate(sig_genes):
            if g in all_markers[ct]:
                sig[i, j] = 1.0
    return cell_types, sig_genes, sig


def deconvolve_sample(bulk_sig_vec, signature):
    """NNLS deconvolution of one sample. Returns proportions (sums to 1)."""
    props, residual = nnls(signature.T, bulk_sig_vec)
    if props.sum() > 0:
        props = props / props.sum()
    return props


# ============================================================
# Part 1: Microarray deconvolution (HZE-Fe + UV-B)
# ============================================================

def load_microarray(study_id):
    """Load a microarray probeset expression file and return (genes, samples, expr_df)."""
    study_dir = RAW_DIR / study_id
    # Find the probeset expression file
    files = list(study_dir.glob("*_array_normalized_expression_probeset*.csv"))
    if not files:
        raise FileNotFoundError(f"No probeset expression file in {study_dir}")
    df = pd.read_csv(files[0])
    # First column is TAIR gene ID
    gene_col = df.columns[0]
    # Sample columns: everything after the metadata columns
    meta_cols = {"TAIR", "SYMBOL", "GENENAME", "REFSEQ", "ENTREZID",
                 "STRING_id", "GOSLIM_IDS", "ProbesetID", "count_ENSEMBL_mappings"}
    sample_cols = [c for c in df.columns if c not in meta_cols]
    expr = df[["TAIR"] + sample_cols].copy()
    expr = expr.dropna(subset=["TAIR"])
    # Remove duplicate genes (keep max expression)
    expr = expr.groupby("TAIR")[sample_cols].max()
    return list(expr.index), sample_cols, expr


def deconvolve_microarray():
    """Deconvolve HZE-Fe and UV-B microarray samples."""
    print("\n=== Part 1: Microarray deconvolution ===")
    md = pd.read_csv(DATA_DIR / "metadata_master.csv")

    all_props = []
    all_expr_weighted = []

    for study_id in ["OSD-46", "OSD-320", "OSD-296"]:
        print(f"\n  Loading {study_id}...")
        genes, samples, expr = load_microarray(study_id)
        print(f"    {len(genes)} genes x {len(samples)} samples")

        # Build signature matrix
        cell_types, sig_genes, signature = build_signature_matrix(genes)
        print(f"    Signature: {len(cell_types)} cell types x {len(sig_genes)} markers")

        # Get metadata for this study
        study_md = md[md["StudyID"] == study_id].copy()
        study_md["IsControl_bool"] = study_md["IsControl"].astype(str).str.lower().isin(["true", "1", "yes"])

        # Deconvolve each sample
        for sample_col in samples:
            # Map sample column to metadata SampleID
            sample_md = study_md[study_md["SampleID"] == sample_col]
            if len(sample_md) == 0:
                continue

            bulk_sig = expr.loc[sig_genes, sample_col].fillna(0).values
            props = deconvolve_sample(bulk_sig, signature)

            row = {"SampleID": sample_col, "StudyID": study_id}
            for i, ct in enumerate(cell_types):
                row[ct] = props[i]
            row["RadiationQuality"] = sample_md["RadiationQuality"].iloc[0]
            row["TimePostExposure_h"] = sample_md["TimePostExposure_h"].iloc[0]
            row["IsControl"] = sample_md["IsControl_bool"].iloc[0]
            all_props.append(row)

            # Weighted expression for pseudo-cell-type profile
            # Weight each gene's expression by the cell-type proportion
            sample_expr = expr[sample_col].fillna(0).values
            for i, ct in enumerate(cell_types):
                if props[i] > 0:
                    all_expr_weighted.append({
                        "cell_type": ct, "sample": sample_col,
                        "expr": sample_expr * props[i]
                    })

    prop_df = pd.DataFrame(all_props)
    # Ensure all cell types are columns
    all_ct = sorted(set(CELL_TYPE_MARKERS.keys()) | set(STRESS_MARKERS.keys()))
    for ct in all_ct:
        if ct not in prop_df.columns:
            prop_df[ct] = 0.0

    prop_path = DEV_DIR / "proportions_microarray.csv"
    prop_df.to_csv(prop_path, index=False)
    print(f"\n  Proportions -> {prop_path} ({prop_df.shape})")

    # Compute pseudo-cell-type expression (weighted average)
    print("  Computing pseudo-cell-type expression...")
    gene_index = None
    ct_expr = {}
    for ct in all_ct:
        ct_samples = [r for r in all_expr_weighted if r["cell_type"] == ct]
        if not ct_samples:
            continue
        if gene_index is None:
            gene_index = genes
        total_weight = sum(1 for s in ct_samples)  # unweighted mean of weighted expr
        summed = np.zeros(len(gene_index))
        for s in ct_samples:
            summed += s["expr"]
        ct_expr[ct] = summed / total_weight if total_weight > 0 else summed

    pseudo_df = pd.DataFrame(ct_expr, index=gene_index).T
    pseudo_path = DEV_DIR / "pseudo_celltype_expr_microarray.csv"
    pseudo_df.to_csv(pseudo_path)
    print(f"  Pseudo cell-type expr -> {pseudo_path} ({pseudo_df.shape})")

    # Print proportion summary by quality
    print("\n  === Proportion summary by quality ===")
    for q in prop_df["RadiationQuality"].unique():
        if not q:
            continue
        q_df = prop_df[prop_df["RadiationQuality"] == q]
        print(f"\n  {q} ({len(q_df)} samples):")
        means = q_df[all_ct].mean().sort_values(ascending=False).head(5)
        for ct, val in means.items():
            if val > 0:
                print(f"    {ct}: {val:.3f}")

    return prop_df, pseudo_df, all_ct


# ============================================================
# Part 2: Per-quality, per-timepoint CellChat
# ============================================================

def load_rnaseq_data():
    """Load RNA-seq expression, proportions, and metadata."""
    import anndata as ad
    adata = ad.read_h5ad(DATA_DIR / "expression_arabidopsis_space.h5ad")
    pce = pd.read_csv(RESULTS_DIR / "deconvolution/pseudo_celltype_expr.csv", index_col=0)
    prop = pd.read_csv(RESULTS_DIR / "deconvolution/proportions.csv", index_col=0)
    # Merge proportions with metadata via SampleID
    obs_sid = adata.obs[["SampleID"]].copy()
    prop = prop.merge(obs_sid, left_index=True, right_index=True)
    md = pd.read_csv(DATA_DIR / "metadata_master.csv")
    return adata, pce, prop, md


def compute_quality_timepoint_expression(pce, prop, md, cell_types, platform="rnaseq"):
    """Compute mean pseudo-cell-type expression per (quality, timepoint) bin.

    For RNA-seq: uses proportions from deconvolution to weight pseudo-cell-type profiles.
    For microarray: uses microarray proportions and pseudo-cell-type profiles.
    """
    # Merge proportions with metadata
    merged = prop.merge(md[["SampleID", "TimePostExposure_h", "StudyID",
                            "RadiationQuality", "IsControl"]], on="SampleID", how="inner")
    merged["TimePostExposure_h"] = pd.to_numeric(merged["TimePostExposure_h"], errors="coerce")
    merged["IsControl_bool"] = merged["IsControl"].astype(str).str.lower().isin(["true", "1", "yes"])

    qt_expr = {}
    qt_meta = {}

    # Per-quality, per-timepoint
    for quality in merged["RadiationQuality"].unique():
        if not quality or quality == "":
            continue
        q_df = merged[merged["RadiationQuality"] == quality]

        for lo, hi, label in TIME_BINS:
            mask = (q_df["TimePostExposure_h"] >= lo) & \
                   (q_df["TimePostExposure_h"] < hi) & \
                   (~q_df["IsControl_bool"]) & \
                   (~q_df["TimePostExposure_h"].isna())
            samples = q_df[mask]
            if len(samples) < 2:
                continue

            # Mean proportion per cell type
            available_cts = [ct for ct in cell_types if ct in samples.columns]
            mean_props = samples[available_cts].mean(axis=0)

            # Scale pseudo-cell-type profiles by mean proportions
            if quality in pce.index.get_level_values(0) if isinstance(pce.index, pd.MultiIndex) else False:
                # This shouldn't happen with our data structure
                pass

            # pce: rows = cell types, cols = genes
            # We need to scale each cell type's profile by its mean proportion
            ct_in_pce = [ct for ct in available_cts if ct in pce.index]
            if not ct_in_pce:
                continue

            scaled = pce.loc[ct_in_pce].mul(mean_props[ct_in_pce].values, axis=0)
            scaled = np.log1p(scaled.clip(lower=0))

            key = f"{quality}__{label}"
            qt_expr[key] = scaled
            qt_meta[key] = {
                "quality": quality,
                "time_bin": label,
                "n_samples": int(len(samples)),
                "mean_time_h": float(samples["TimePostExposure_h"].mean()),
                "studies": list(samples["StudyID"].unique()),
                "platform": platform,
            }
            print(f"    {key}: n={len(samples)}, mean_time={qt_meta[key]['mean_time_h']:.2f}h")

    # Also add static profiles for spaceflight-LEO (no time data)
    sfl = merged[(merged["RadiationQuality"] == "spaceflight-LEO") & (~merged["IsControl_bool"])]
    if len(sfl) >= 2:
        available_cts = [ct for ct in cell_types if ct in sfl.columns]
        mean_props = sfl[available_cts].mean(axis=0)
        ct_in_pce = [ct for ct in available_cts if ct in pce.index]
        if ct_in_pce:
            scaled = pce.loc[ct_in_pce].mul(mean_props[ct_in_pce].values, axis=0)
            scaled = np.log1p(scaled.clip(lower=0))
            key = "spaceflight-LEO__static"
            qt_expr[key] = scaled
            qt_meta[key] = {
                "quality": "spaceflight-LEO",
                "time_bin": "static",
                "n_samples": int(len(sfl)),
                "mean_time_h": float("nan"),
                "studies": list(sfl["StudyID"].unique()),
                "platform": platform,
            }
            print(f"    {key}: n={len(sfl)} (static, chronic)")

    return qt_expr, qt_meta


def compute_communication_with_lr(expr_celltypes, lr, cell_types, Kh=KH):
    """Compute NxN communication matrix AND per-LR-pair contributions.

    Vectorized: builds ligand and receptor expression matrices, then computes
    all LR pair contributions via batched outer products.

    Returns: (comm_total, lr_strength_dict)
    lr_strength_dict: {LR_pair: total_signaling_strength}
    """
    # Use actual cell types present in the expression matrix
    actual_cts = list(expr_celltypes.index)
    n_ct = len(actual_cts)
    # Filter LR pairs to available genes
    available_genes = set(expr_celltypes.columns)
    lr_avail = lr[(lr["Ligands"].isin(available_genes)) & (lr["Receptors"].isin(available_genes))].copy()
    if len(lr_avail) == 0:
        return np.zeros((n_ct, n_ct)), {}

    # Build expression matrices: [n_ct, n_lr_pairs]
    L_mat = expr_celltypes.loc[:, lr_avail["Ligands"].values].values  # [n_ct, n_lr]
    R_mat = expr_celltypes.loc[:, lr_avail["Receptors"].values].values  # [n_ct, n_lr]

    comm_total = np.zeros((n_ct, n_ct))
    lr_strength = {}

    batch_size = 500
    for start in range(0, len(lr_avail), batch_size):
        end = min(start + batch_size, len(lr_avail))
        L_batch = L_mat[:, start:end]  # [n_ct, batch]
        R_batch = R_mat[:, start:end]  # [n_ct, batch]
        prod = L_batch[:, None, :] * R_batch[None, :, :]  # [n_ct, n_ct, batch]
        prob = prod / (Kh + prod)
        comm_total += prob.sum(axis=2)
        pair_strengths = prob.sum(axis=(0, 1))  # [batch]
        for idx_in_batch, k in enumerate(range(start, end)):
            if pair_strengths[idx_in_batch] > 0:
                lr_strength[lr_avail["LR_pair"].iloc[k]] = float(pair_strengths[idx_in_batch])

    return comm_total, lr_strength


def permutation_pvalues_fast(expr_celltypes, lr, cell_types, n_perm=N_PERMUTATIONS, Kh=KH):
    """Permutation test for communication matrix significance.
    Vectorized version using precomputed expression matrices."""
    actual_cts = list(expr_celltypes.index)
    n_ct = len(actual_cts)
    observed, _ = compute_communication_with_lr(expr_celltypes, lr, cell_types, Kh)

    # Precompute ligand and receptor expression matrices
    available_genes = set(expr_celltypes.columns)
    lr_avail = lr[(lr["Ligands"].isin(available_genes)) & (lr["Receptors"].isin(available_genes))].copy()
    if len(lr_avail) == 0:
        return observed, np.ones((n_ct, n_ct))

    L_mat = expr_celltypes.loc[:, lr_avail["Ligands"].values].values  # [n_ct, n_lr]
    R_mat = expr_celltypes.loc[:, lr_avail["Receptors"].values].values  # [n_ct, n_lr]
    n_lr = len(lr_avail)

    ge_null = np.zeros((n_ct, n_ct))
    rng = np.random.default_rng(42)
    batch_size = 500

    for _ in range(n_perm):
        perm = rng.permutation(n_ct)
        L_perm = L_mat[perm, :]  # permute cell-type rows
        R_perm = R_mat[perm, :]
        null_comm = np.zeros((n_ct, n_ct))
        for start in range(0, n_lr, batch_size):
            end = min(start + batch_size, n_lr)
            prod = L_perm[:, None, start:end] * R_perm[None, :, start:end]
            prob = prod / (Kh + prod)
            null_comm += prob.sum(axis=2)
        ge_null += (null_comm >= observed).astype(float)

    pvals = ge_null / n_perm
    return observed, pvals


def run_per_quality_cellchat(qt_expr, qt_meta, lr_rnaseq, lr_microarray, cell_types):
    """Run CellChat for each quality-timepoint combination."""
    print("\n=== Part 2: Per-quality CellChat ===")
    all_flow = []
    all_lr_strength = []

    for key, expr_ct in qt_expr.items():
        meta = qt_meta[key]
        quality = meta["quality"]
        platform = meta["platform"]

        # Use appropriate LR database
        lr = lr_rnaseq if platform == "rnaseq" else lr_microarray

        # Filter LR pairs to genes in expression
        available_genes = set(expr_ct.columns)
        lr_avail = lr[(lr["Ligands"].isin(available_genes)) & (lr["Receptors"].isin(available_genes))].copy()
        if len(lr_avail) < 10:
            print(f"  Skipping {key}: only {len(lr_avail)} LR pairs available")
            continue

        print(f"  {key}: {len(lr_avail)} LR pairs, {expr_ct.shape[0]} cell types...")

        # Communication matrix + per-LR-pair strength
        comm, lr_strength = compute_communication_with_lr(expr_ct, lr_avail, cell_types)

        # Permutation p-values
        _, pvals = permutation_pvalues_fast(expr_ct, lr_avail, cell_types)

        # Flow rows (use actual cell types from expression matrix)
        actual_cts = list(expr_ct.index)
        for i, src in enumerate(actual_cts):
            for j, tgt in enumerate(actual_cts):
                if i == j:
                    continue
                all_flow.append({
                    "Quality": quality,
                    "Time": meta["time_bin"],
                    "Source": src,
                    "Target": tgt,
                    "SignalStrength": float(comm[i, j]),
                    "pvalue": float(pvals[i, j]),
                })

        # Per-LR-pair strength
        for lr_pair, strength in lr_strength.items():
            all_lr_strength.append({
                "Quality": quality,
                "Time": meta["time_bin"],
                "LR_pair": lr_pair,
                "SignalStrength": strength,
                "n_samples": meta["n_samples"],
                "platform": platform,
            })

    flow_df = pd.DataFrame(all_flow)
    lr_str_df = pd.DataFrame(all_lr_strength)

    flow_path = OUT_DIR / "signaling_flow_per_quality_timepoint.csv"
    flow_df.to_csv(flow_path, index=False)
    print(f"\n  Signaling flow -> {flow_path} ({len(flow_df)} rows)")

    lr_str_path = OUT_DIR / "lr_pair_strength_per_quality.csv"
    lr_str_df.to_csv(lr_str_path, index=False)
    print(f"  LR pair strength -> {lr_str_path} ({len(lr_str_df)} rows)")

    return flow_df, lr_str_df


def compute_per_sample_lr_strength(qt_expr, qt_meta, lr_rnaseq, lr_microarray,
                                    pce_rnaseq, prop_rnaseq, md_rnaseq,
                                    pseudo_ma, prop_ma, md_ma, cell_types):
    """Compute per-sample LR pair signaling strengths for specificity testing.

    For each sample in each quality-timepoint group, scale pseudo-cell-type
    profiles by that sample's proportions, then compute per-LR-pair signaling.
    This gives n_samples values per quality per LR pair for Kruskal-Wallis.
    """
    print("\n  Computing per-sample LR pair strengths...")

    # Build per-sample expression for RNA-seq qualities
    # prop_rnaseq has SampleID column, cell types as other columns
    merged_rna = prop_rnaseq.merge(
        md_rnaseq[["SampleID", "TimePostExposure_h", "RadiationQuality", "IsControl"]],
        on="SampleID", how="inner"
    )
    merged_rna["IsControl_bool"] = merged_rna["IsControl"].astype(str).str.lower().isin(["true", "1", "yes"])
    merged_rna["TimePostExposure_h"] = pd.to_numeric(merged_rna["TimePostExposure_h"], errors="coerce")

    # Build per-sample expression for microarray qualities
    # prop_ma already has SampleID, RadiationQuality, TimePostExposure_h, IsControl
    merged_ma = prop_ma.copy()
    merged_ma["IsControl_bool"] = merged_ma["IsControl"].astype(str).str.lower().isin(["true", "1", "yes"])
    merged_ma["TimePostExposure_h"] = pd.to_numeric(merged_ma["TimePostExposure_h"], errors="coerce")

    all_sample_lr = []

    for lo, hi, label in TIME_BINS:
        # RNA-seq samples at this timepoint
        for platform, merged_df, pce_df, lr_db in [
            ("rnaseq", merged_rna, pce_rnaseq, lr_rnaseq),
            ("microarray", merged_ma, pseudo_ma, lr_microarray),
        ]:
            mask = (merged_df["TimePostExposure_h"] >= lo) & \
                   (merged_df["TimePostExposure_h"] < hi) & \
                   (~merged_df["IsControl_bool"]) & \
                   (~merged_df["TimePostExposure_h"].isna())
            samples = merged_df[mask]
            if len(samples) < 2:
                continue

            for _, sample_row in samples.iterrows():
                sample_id = sample_row["SampleID"]
                quality = sample_row["RadiationQuality"]
                if not quality:
                    continue

                # Get this sample's proportions
                ct_cols = [ct for ct in cell_types if ct in sample_row.index]
                sample_props = sample_row[ct_cols].values.astype(float)

                # Scale pseudo-cell-type profiles by this sample's proportions
                ct_in_pce = [ct for ct in ct_cols if ct in pce_df.index]
                if not ct_in_pce:
                    continue
                prop_vals = np.array([sample_row[ct] for ct in ct_in_pce])
                scaled = pce_df.loc[ct_in_pce].mul(prop_vals, axis=0)
                scaled = np.log1p(scaled.clip(lower=0))

                # Compute per-LR-pair strength for this sample
                available_genes = set(scaled.columns)
                lr_avail = lr_db[(lr_db["Ligands"].isin(available_genes)) &
                                 (lr_db["Receptors"].isin(available_genes))].copy()
                if len(lr_avail) < 10:
                    continue

                L_mat = scaled.loc[:, lr_avail["Ligands"].values].values
                R_mat = scaled.loc[:, lr_avail["Receptors"].values].values
                n_lr = len(lr_avail)

                batch_size = 500
                for start in range(0, n_lr, batch_size):
                    end = min(start + batch_size, n_lr)
                    prod = L_mat[:, None, start:end] * R_mat[None, :, start:end]
                    prob = prod / (KH + prod)
                    pair_strengths = prob.sum(axis=(0, 1))
                    for idx, k in enumerate(range(start, end)):
                        if pair_strengths[idx] > 0:
                            all_sample_lr.append({
                                "Quality": quality,
                                "Time": label,
                                "SampleID": sample_id,
                                "LR_pair": lr_avail["LR_pair"].iloc[k],
                                "SignalStrength": float(pair_strengths[idx]),
                                "platform": platform,
                            })

    # Also handle spaceflight-LEO static
    sfl_rna = merged_rna[(merged_rna["RadiationQuality"] == "spaceflight-LEO") & (~merged_rna["IsControl_bool"])]
    for _, sample_row in sfl_rna.iterrows():
        sample_id = sample_row["SampleID"]
        ct_cols = [ct for ct in cell_types if ct in sample_row.index]
        ct_in_pce = [ct for ct in ct_cols if ct in pce_rnaseq.index]
        if not ct_in_pce:
            continue
        prop_vals = np.array([sample_row[ct] for ct in ct_in_pce])
        scaled = pce_rnaseq.loc[ct_in_pce].mul(prop_vals, axis=0)
        scaled = np.log1p(scaled.clip(lower=0))
        available_genes = set(scaled.columns)
        lr_avail = lr_rnaseq[(lr_rnaseq["Ligands"].isin(available_genes)) &
                             (lr_rnaseq["Receptors"].isin(available_genes))].copy()
        if len(lr_avail) < 10:
            continue
        L_mat = scaled.loc[:, lr_avail["Ligands"].values].values
        R_mat = scaled.loc[:, lr_avail["Receptors"].values].values
        n_lr = len(lr_avail)
        batch_size = 500
        for start in range(0, n_lr, batch_size):
            end = min(start + batch_size, n_lr)
            prod = L_mat[:, None, start:end] * R_mat[None, :, start:end]
            prob = prod / (KH + prod)
            pair_strengths = prob.sum(axis=(0, 1))
            for idx, k in enumerate(range(start, end)):
                if pair_strengths[idx] > 0:
                    all_sample_lr.append({
                        "Quality": "spaceflight-LEO",
                        "Time": "static",
                        "SampleID": sample_id,
                        "LR_pair": lr_avail["LR_pair"].iloc[k],
                        "SignalStrength": float(pair_strengths[idx]),
                        "platform": "rnaseq",
                    })

    sample_lr_df = pd.DataFrame(all_sample_lr)
    sample_lr_path = OUT_DIR / "lr_pair_strength_per_sample.csv"
    sample_lr_df.to_csv(sample_lr_path, index=False)
    print(f"  Per-sample LR strength -> {sample_lr_path} ({len(sample_lr_df)} rows)")
    return sample_lr_df


# ============================================================
# Part 3: LR pair specificity analysis
# ============================================================

def compute_specificity(sample_lr_df):
    """Test LR pair specificity across radiation qualities using per-sample strengths.

    Uses Kruskal-Wallis test (non-parametric) across qualities at shared timepoints,
    with BH-FDR correction. This is valid because we have multiple samples per quality.
    """
    print("\n=== Part 3: LR pair specificity analysis (per-sample) ===")

    # Find timepoints with >=2 qualities
    timepoint_qualities = {}
    for t in sample_lr_df["Time"].unique():
        if t == "static":
            continue
        quals = sample_lr_df[sample_lr_df["Time"] == t]["Quality"].unique()
        if len(quals) >= 2:
            timepoint_qualities[t] = sorted(list(quals))

    print(f"  Shared timepoints: {timepoint_qualities}")

    results = []
    for timepoint, qualities in timepoint_qualities.items():
        print(f"\n  Testing at {timepoint}: {qualities}")

        # Get per-sample strengths for this timepoint
        tp_df = sample_lr_df[sample_lr_df["Time"] == timepoint]

        # Group by LR pair and quality
        # For each LR pair, collect sample strengths per quality
        lr_pivot = tp_df.pivot_table(index="LR_pair", columns="Quality",
                                      values="SignalStrength", aggfunc=list)

        # Only keep LR pairs present in all qualities
        lr_pivot = lr_pivot.dropna()
        print(f"    Common LR pairs (in all qualities): {len(lr_pivot)}")

        if len(lr_pivot) < 10:
            continue

        # For each LR pair, run Kruskal-Wallis across qualities
        pvalues = []
        for lr_pair in lr_pivot.index:
            groups = [lr_pivot.loc[lr_pair, q] for q in qualities]
            # Ensure each group has >=2 samples
            if any(len(g) < 2 for g in groups):
                pvalues.append(1.0)
                continue
            try:
                stat, pval = kruskal(*groups)
                pvalues.append(pval)
            except Exception:
                pvalues.append(1.0)

        pvalues = np.array(pvalues)

        # BH-FDR correction
        n = len(pvalues)
        sorted_idx = np.argsort(pvalues)
        sorted_pvals = pvalues[sorted_idx]
        sorted_adj = sorted_pvals * n / np.arange(1, n + 1)
        sorted_adj = np.minimum.accumulate(sorted_adj[::-1])[::-1]
        padj = np.zeros(n)
        padj[sorted_idx] = np.minimum(sorted_adj, 1.0)

        # Build results
        for i, lr_pair in enumerate(lr_pivot.index):
            groups = [lr_pivot.loc[lr_pair, q] for q in qualities]
            group_means = [np.mean(g) for g in groups]
            group_stds = [np.std(g) for g in groups]

            max_idx = np.argmax(group_means)
            min_idx = np.argmin(group_means)
            max_qual = qualities[max_idx]
            min_qual = qualities[min_idx]
            max_mean = group_means[max_idx]
            min_mean = group_means[min_idx]
            max_min_ratio = max_mean / min_mean if min_mean > 0 else float("inf")

            # Specificity classification
            if padj[i] < 0.05 and max_min_ratio > 2.0:
                specificity = "quality-specific"
            elif padj[i] > 0.1 or max_min_ratio < 1.5:
                specificity = "quality-shared"
            else:
                specificity = "intermediate"

            row = {
                "LR_pair": lr_pair,
                "timepoint": timepoint,
                "n_qualities": len(qualities),
                "max_quality": max_qual,
                "min_quality": min_qual,
                "max_mean_strength": max_mean,
                "min_mean_strength": min_mean,
                "max_min_ratio": max_min_ratio,
                "pvalue": pvalues[i],
                "padj": padj[i],
                "specificity": specificity,
            }
            for j, q in enumerate(qualities):
                row[f"mean_{q}"] = group_means[j]
                row[f"std_{q}"] = group_stds[j]
                row[f"n_{q}"] = len(groups[j])
            results.append(row)

    spec_df = pd.DataFrame(results)
    spec_path = OUT_DIR / "lr_pair_specificity.csv"
    spec_df.to_csv(spec_path, index=False)
    print(f"\n  LR pair specificity -> {spec_path} ({len(spec_df)} rows)")

    # Summary
    if len(spec_df) > 0:
        print("\n  === Specificity summary ===")
        for tp in spec_df["timepoint"].unique():
            tp_df = spec_df[spec_df["timepoint"] == tp]
            counts = tp_df["specificity"].value_counts()
            print(f"  {tp}: {dict(counts)}")
            specific = tp_df[tp_df["specificity"] == "quality-specific"]
            if len(specific) > 0:
                print(f"    Top 5 specific pairs:")
                for _, r in specific.nlargest(5, "max_min_ratio").iterrows():
                    print(f"      {r['LR_pair']}: ratio={r['max_min_ratio']:.2f}, "
                          f"enriched in {r['max_quality']}, padj={r['padj']:.4f}")

    return spec_df


# ============================================================
# Main
# ============================================================

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DEV_DIR.mkdir(parents=True, exist_ok=True)

    # Load LR databases
    lr_full = pd.read_csv(LR_DB_PATH)
    lr_full["LR_pair"] = lr_full["Ligands"] + "->" + lr_full["Receptors"]

    # ---- Part 1: Microarray deconvolution ----
    prop_ma, pseudo_ma, cell_types_ma = deconvolve_microarray()

    # ---- Part 2: Per-quality CellChat ----
    # Load RNA-seq data
    print("\n=== Loading RNA-seq data ===")
    adata, pce_rnaseq, prop_rnaseq, md = load_rnaseq_data()
    print(f"  RNA-seq expression: {adata.shape}")
    print(f"  RNA-seq pseudo-cell-type expr: {pce_rnaseq.shape}")
    print(f"  RNA-seq proportions: {prop_rnaseq.shape}")

    # Filter LR pairs for RNA-seq
    rnaseq_genes = set(pce_rnaseq.columns)
    lr_rnaseq = lr_full[(lr_full["Ligands"].isin(rnaseq_genes)) &
                        (lr_full["Receptors"].isin(rnaseq_genes))].copy()
    print(f"  LR pairs (RNA-seq): {len(lr_rnaseq)}")

    # Filter LR pairs for microarray
    ma_genes = set(pseudo_ma.columns)
    lr_microarray = lr_full[(lr_full["Ligands"].isin(ma_genes)) &
                            (lr_full["Receptors"].isin(ma_genes))].copy()
    print(f"  LR pairs (microarray): {len(lr_microarray)}")

    # All cell types (union)
    all_cell_types = sorted(set(CELL_TYPE_MARKERS.keys()) | set(STRESS_MARKERS.keys()))
    print(f"  All cell types: {len(all_cell_types)}")

    # Compute quality-timepoint expression for RNA-seq qualities
    print("\n  Computing RNA-seq quality-timepoint expression...")
    qt_expr_rna, qt_meta_rna = compute_quality_timepoint_expression(
        pce_rnaseq, prop_rnaseq, md, all_cell_types, platform="rnaseq"
    )

    # Compute quality-timepoint expression for microarray qualities
    print("\n  Computing microarray quality-timepoint expression...")
    # Prepare microarray proportions in the same format as RNA-seq
    prop_ma_for_merge = prop_ma.copy()
    prop_ma_for_merge["SampleID"] = prop_ma["SampleID"]
    # Need to merge with metadata for time/quality
    md_ma = md[["SampleID", "TimePostExposure_h", "StudyID", "RadiationQuality", "IsControl"]].copy()
    prop_ma_merged = prop_ma_for_merge.merge(md_ma, on="SampleID", how="inner")

    # Build a proportions DataFrame with SampleID as a column (for merge)
    ct_cols = [c for c in prop_ma_merged.columns if c in all_cell_types]
    prop_ma_for_qt = prop_ma_merged[["SampleID"] + ct_cols].copy()

    # For microarray, we need to use the microarray pseudo-cell-type expression
    qt_expr_ma, qt_meta_ma = compute_quality_timepoint_expression(
        pseudo_ma, prop_ma_for_qt, md, all_cell_types, platform="microarray"
    )

    # Merge all quality-timepoint expressions
    qt_expr_all = {**qt_expr_rna, **qt_expr_ma}
    qt_meta_all = {**qt_meta_rna, **qt_meta_ma}

    print(f"\n  Total quality-timepoint combinations: {len(qt_expr_all)}")
    for key in sorted(qt_expr_all.keys()):
        m = qt_meta_all[key]
        print(f"    {key}: n={m['n_samples']}, platform={m['platform']}")

    # Run CellChat
    flow_df, lr_str_df = run_per_quality_cellchat(
        qt_expr_all, qt_meta_all, lr_rnaseq, lr_microarray, all_cell_types
    )

    # Save quality-timepoint metadata
    with open(OUT_DIR / "quality_timepoint_metadata.json", "w") as f:
        json.dump(qt_meta_all, f, indent=2, default=lambda x: None if np.isnan(x) else float(x))

    # ---- Part 3: Per-sample LR pair strengths for specificity testing ----
    # Prepare microarray proportions for per-sample computation
    md_ma = md[["SampleID", "TimePostExposure_h", "StudyID", "RadiationQuality", "IsControl"]].copy()
    prop_ma_for_sample = prop_ma.copy()

    sample_lr_df = compute_per_sample_lr_strength(
        qt_expr_all, qt_meta_all, lr_rnaseq, lr_microarray,
        pce_rnaseq, prop_rnaseq, md,
        pseudo_ma, prop_ma_for_sample, md_ma, all_cell_types
    )

    # ---- Part 4: LR pair specificity (Kruskal-Wallis) ----
    spec_df = compute_specificity(sample_lr_df)

    print("\n=== Done ===")


if __name__ == "__main__":
    main()
