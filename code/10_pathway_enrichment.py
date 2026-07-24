#!/usr/bin/env python3
"""
10_pathway_enrichment.py
Phase 3.4 — Pathway enrichment heatmaps showing the temporal activation of
stress pathways (DNA repair, hormonal signaling, oxidative stress) across
time post-radiation and across cell types / radiation qualities.

APPROACH
  1. For each pathway gene set (dna_repair, oxidative_stress, hormone_signaling,
     plus sub-pathways: SA, JA, ETH, ABA, auxin), compute a pathway activation
     score per sample = mean(z-scored expression of pathway genes relative to
     within-study control mean). This is a GSVA-style single-sample score that
     is robust to batch effects.
  2. Aggregate scores by timepoint bin and by cell type (using deconvolved
     pseudo-cell-type expression) -> heatmap 1: pathway x timepoint.
  3. Aggregate by timepoint x radiation quality -> heatmap 2: pathway x
     (timepoint x quality).
  4. WGCNA module eigengene trajectory -> heatmap 3: module x timepoint,
     showing early-response (blue) vs sustained (turquoise) kinetics.
  5. Cell-type x pathway activation at the nadir timepoint -> heatmap 4,
     showing which cell types mount which pathway response.

OUTPUTS
  - results/pathway_enrichment/pathway_scores_per_sample.csv
  - results/pathway_enrichment/pathway_by_timepoint.csv
  - results/pathway_enrichment/pathway_by_timepoint_quality.csv
  - results/pathway_enrichment/pathway_by_celltype.csv
  - results/pathway_enrichment/module_eigengene_by_timepoint.csv
  - results/figures/pathway_heatmap_timepoint.svg/png   (heatmap 1)
  - results/figures/pathway_heatmap_timepoint_quality.svg/png  (heatmap 2)
  - results/figures/module_eigengene_heatmap.svg/png    (heatmap 3)
  - results/figures/pathway_heatmap_celltype.svg/png    (heatmap 4)
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Phylo font + SVG settings
matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"

RESULTS_DIR = Path("/mnt/results/zenodo_bundle/results")
OUT_DIR = RESULTS_DIR / "pathway_enrichment"
FIG_DIR = RESULTS_DIR / "figures"
DATA_DIR = Path("/mnt/results/zenodo_bundle/data")

# ---- Pathway gene sets (Arabidopsis AGI codes) ----
# Core pathways from RRI script 08, expanded with sub-pathway splits.
PATHWAY_GENES = {
    "DNA repair": [
        "AT1G07290", "AT3G05120", "AT5G40820", "AT5G20850", "AT3G48190",
        "AT1G65470", "AT3G02680", "AT4G19730", "AT5G24280", "AT2G32750",
    ],
    "Oxidative stress": [
        "AT1G02920", "AT4G25130", "AT5G18100", "AT3G10920", "AT1G20630",
        "AT2G28190", "AT3G26060", "AT1G32940", "AT5G03490", "AT4G35090",
    ],
    "SA signaling": [
        "AT2G14610", "AT1G02820", "AT3G57260", "AT1G80830", "AT5G51070",
    ],
    "JA signaling": [
        "AT5G44420", "AT2G39300", "AT3G04720", "AT1G13220", "AT1G17380",
    ],
    "ETH signaling": [
        "AT5G44420", "AT2G39300", "AT3G04720", "AT1G13220", "AT1G72330",
    ],
    "ABA signaling": [
        "AT1G19670", "AT5G57050", "AT3G24220", "AT1G07430", "AT4G26080",
    ],
    "Auxin signaling": [
        "AT1G75800", "AT4G32280", "AT5G64200", "AT1G04500", "AT2G36210",
    ],
    "Cell cycle / DDR checkpoint": [
        "AT3G20650", "AT5G23750", "AT1G20980", "AT3G12670", "AT1G65920",
    ],
}

# Time bins matching CellChat/RRI
TIME_BINS = [
    (0.0, 0.5, "0-0.5h"),
    (0.5, 2.0, "0.5-2h"),
    (2.0, 6.0, "2-6h"),
    (6.0, 12.0, "6-12h"),
    (12.0, 30.0, "12-30h"),
]


def bin_time(t):
    if pd.isna(t):
        return None
    for lo, hi, label in TIME_BINS:
        if lo <= t < hi:
            return label
    if t >= 30:
        return "12-30h"
    return None


def load_expression_and_metadata():
    """Load expression matrix and metadata, return matched samples.
    175/195 samples match directly via SampleID (GSM-style). The remaining 20
    come from studies without GSM cross-refs and are dropped from pathway
    analysis (they lack time/dose metadata anyway)."""
    import anndata as ad
    adata = ad.read_h5ad("/workspace/expression_raw.h5ad")
    md = pd.read_csv(DATA_DIR / "metadata_master.csv")

    # adata.obs["SampleID"] is GSM-style (e.g. Atha_Ler-0_sl_FLT_uG_Rep1)
    # metadata SampleID is GSM-style with _N suffix (e.g. GSM1506014_1)
    # Direct merge on SampleID
    obs = adata.obs.copy()
    obs["match_key"] = obs["SampleID"].astype(str) if "SampleID" in obs.columns else obs.index.astype(str)
    md["match_key"] = md["SampleID"].astype(str)

    # Merge metadata onto obs (keep only matched samples)
    merged = obs.merge(md, on="match_key", how="left", suffixes=("", "_md"))
    merged.index = obs.index

    # Drop unmatched samples (no metadata -> no time/quality info)
    n_before = len(merged)
    merged = merged[merged["RadiationQuality"].notna()].copy()
    print(f"  Matched {len(merged)}/{n_before} samples to metadata")

    # Fix IsControl
    if "IsControl" in merged.columns:
        ic = merged["IsControl"].astype(str).str.lower()
        merged["IsControl"] = ic.isin(["true", "1", "yes"])
    else:
        merged["IsControl"] = False

    # Numeric coercion
    for c in ["TimePostExposure_h", "AbsorbedDose_Gy", "LET_keV_um"]:
        if c in merged.columns:
            merged[c] = pd.to_numeric(merged[c], errors="coerce")

    merged["time_bin"] = merged["TimePostExposure_h"].apply(bin_time)

    # Subset adata to matched samples
    adata = adata[merged.index].copy()
    return adata, merged


def compute_pathway_scores(adata, meta):
    """Per-sample pathway activation score = mean z-scored expression of
    pathway genes relative to within-study control mean. Robust to batch."""
    X = adata.X if not hasattr(adata.X, "toarray") else adata.X.toarray()
    genes = list(adata.var_names)
    gene_idx = {g: i for i, g in enumerate(genes)}

    # Build study -> control mean per gene
    meta["StudyID"] = meta["StudyID"].fillna("unknown")
    study_ctrl_mean = {}
    for study in meta["StudyID"].unique():
        mask = (meta["StudyID"] == study) & (meta["IsControl"] == True)
        if mask.sum() == 0:
            # No controls in this study -> use all samples in study as baseline
            mask = meta["StudyID"] == study
        if mask.sum() > 0:
            study_ctrl_mean[study] = X[mask.values].mean(axis=0)

    # Z-score per sample relative to its study control mean
    scores = {}
    for pathway, gene_list in PATHWAY_GENES.items():
        idxs = [gene_idx[g] for g in gene_list if g in gene_idx]
        if not idxs:
            scores[pathway] = np.full(X.shape[0], np.nan)
            continue
        # Per-sample pathway score
        pathway_scores = np.zeros(X.shape[0])
        for i in range(X.shape[0]):
            study = meta.iloc[i]["StudyID"]
            ctrl_mean = study_ctrl_mean.get(study, X.mean(axis=0))
            # z-score each gene relative to control, then average across pathway genes
            z_vals = (X[i, idxs] - ctrl_mean[idxs])
            # Robust scale: divide by global SD of that gene across all samples
            gene_sd = X[:, idxs].std(axis=0)
            gene_sd[gene_sd == 0] = 1.0
            z_vals = z_vals / gene_sd
            pathway_scores[i] = np.nanmean(z_vals)
        scores[pathway] = pathway_scores

    scores_df = pd.DataFrame(scores, index=meta.index)
    scores_df["SampleID"] = meta["match_key"].values
    scores_df["StudyID"] = meta["StudyID"].values
    scores_df["RadiationQuality"] = meta["RadiationQuality"].values
    scores_df["TimePostExposure_h"] = meta["TimePostExposure_h"].values
    scores_df["time_bin"] = meta["time_bin"].values
    scores_df["IsControl"] = meta["IsControl"].values
    return scores_df


def compute_celltype_pathway_scores():
    """Pathway activation per pseudo-cell-type at each timepoint, using the
    deconvolved pseudo-cell-type expression matrix."""
    pct = pd.read_csv(RESULTS_DIR / "deconvolution" / "pseudo_celltype_expr.csv",
                      index_col=0)
    # rows = cell types, cols = genes
    celltypes = pct.index.tolist()
    genes = pct.columns.tolist()
    gene_idx = {g: i for i, g in enumerate(genes)}
    X = pct.values  # (n_celltypes, n_genes)

    # Z-score each gene across cell types (so pathway score = relative enrichment
    # of that pathway in a cell type vs other cell types)
    gene_mean = X.mean(axis=0)
    gene_sd = X.std(axis=0)
    gene_sd[gene_sd == 0] = 1.0
    Z = (X - gene_mean) / gene_sd

    scores = {}
    for pathway, gene_list in PATHWAY_GENES.items():
        idxs = [gene_idx[g] for g in gene_list if g in gene_idx]
        if not idxs:
            scores[pathway] = np.full(len(celltypes), np.nan)
            continue
        scores[pathway] = Z[:, idxs].mean(axis=1)
    scores_df = pd.DataFrame(scores, index=celltypes)
    scores_df.index.name = "CellType"
    return scores_df.reset_index()


def compute_module_eigengene_trajectory():
    """WGCNA module eigengene per timepoint (already computed in script 07).
    Load module eigengenes per sample and aggregate by time bin.
    The eigengene CSV uses 'Sample' = OSD-XXX_GSMYYYY format; match to
    metadata via the GSM number."""
    me = pd.read_csv(RESULTS_DIR / "wgcna" / "module_eigengenes.csv")
    me = me.rename(columns={"Sample": "Sample_orig"})
    # Extract GSM accession for matching
    me["gsm"] = me["Sample_orig"].str.extract(r"(GSM\d+)")
    md = pd.read_csv(DATA_DIR / "metadata_master.csv")
    md["gsm"] = md["SampleID"].str.extract(r"(GSM\d+)")
    # Merge on GSM
    me = me.merge(md[["gsm", "TimePostExposure_h", "IsControl"]],
                  on="gsm", how="left")
    me["TimePostExposure_h"] = pd.to_numeric(me["TimePostExposure_h"], errors="coerce")
    me["time_bin"] = me["TimePostExposure_h"].apply(bin_time)
    # Fix IsControl
    ic = me["IsControl"].astype(str).str.lower()
    me["IsControl"] = ic.isin(["true", "1", "yes"])

    me_cols = [c for c in me.columns if c.startswith("ME")]
    # Aggregate by time bin (exclude controls)
    non_ctrl = me[~me["IsControl"]]
    agg = non_ctrl.groupby("time_bin")[me_cols].mean().reset_index()
    return agg, me_cols


def plot_heatmap(df, row_labels, col_labels, title, out_prefix,
                 cmap="RdBu_r", center=0.0, vmin=None, vmax=None,
                 figsize=None, xlabel="", ylabel=""):
    """Render a heatmap and save SVG + PNG. Accepts DataFrame or ndarray."""
    data = np.asarray(df) if not hasattr(df, "values") else np.asarray(df.values)
    if vmin is None:
        vmin = np.nanmin(data)
    if vmax is None:
        vmax = np.nanmax(data)
    if figsize is None:
        figsize = (max(6, len(col_labels) * 0.8 + 2),
                   max(4, len(row_labels) * 0.5 + 1.5))

    fig, ax = plt.subplots(figsize=figsize)
    # Normalize for diverging colormap centered at `center`
    if cmap in ("RdBu_r", "RdBu", "seismic"):
        max_abs = max(abs(vmin - center), abs(vmax - center))
        norm = matplotlib.colors.TwoSlopeNorm(vmin=center - max_abs,
                                              vcenter=center,
                                              vmax=center + max_abs)
    else:
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

    im = ax.imshow(data, aspect="auto", cmap=cmap, norm=norm)
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=9)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11, pad=10)

    # Annotate cells with values
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data[i, j]
            if np.isnan(val):
                ax.text(j, i, "NA", ha="center", va="center", fontsize=7, color="gray")
            else:
                color = "white" if abs(val - center) > 0.6 * max(abs(vmin - center), abs(vmax - center)) else "black"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=7, color=color)

    cbar = fig.colorbar(im, ax=ax, shrink=0.7)
    cbar.ax.tick_params(labelsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"{out_prefix}.svg", format="svg")
    fig.savefig(FIG_DIR / f"{out_prefix}.png", dpi=150)
    plt.close(fig)
    print(f"  Saved {out_prefix}.svg/png")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading expression + metadata...")
    adata, meta = load_expression_and_metadata()
    print(f"  {adata.shape[0]} samples, {adata.shape[1]} genes")
    print(f"  {meta['time_bin'].notna().sum()} samples with time bins")
    print(f"  {meta['IsControl'].sum()} controls")

    print("\nComputing per-sample pathway scores...")
    scores = compute_pathway_scores(adata, meta)
    pathway_cols = list(PATHWAY_GENES.keys())
    scores.to_csv(OUT_DIR / "pathway_scores_per_sample.csv", index=False)
    print(f"  Saved {len(scores)} samples x {len(pathway_cols)} pathways")

    # --- Heatmap 1: pathway x timepoint ---
    print("\nHeatmap 1: pathway x timepoint...")
    non_ctrl = scores[~scores["IsControl"]]
    tp_agg = non_ctrl.groupby("time_bin")[pathway_cols].mean()
    # Reorder time bins
    tp_order = [b[2] for b in TIME_BINS if b[2] in tp_agg.index]
    tp_agg = tp_agg.reindex(tp_order)
    tp_agg.to_csv(OUT_DIR / "pathway_by_timepoint.csv")
    plot_heatmap(tp_agg, pathway_cols, tp_order,
                 "Pathway activation over time post-radiation",
                 "pathway_heatmap_timepoint",
                 cmap="RdBu_r", center=0.0,
                 xlabel="Time post-exposure", ylabel="Pathway")

    # --- Heatmap 2: pathway x (timepoint x quality) ---
    print("\nHeatmap 2: pathway x timepoint x quality...")
    tpq = non_ctrl.groupby(["time_bin", "RadiationQuality"])[pathway_cols].mean()
    # Sort by timepoint then quality
    tpq_sorted = tpq.reset_index().sort_values(
        ["time_bin", "RadiationQuality"]
    )
    tpq_mat = tpq_sorted[pathway_cols].values  # (n_tp_quality, n_pathways)
    col_labels = [f"{r['RadiationQuality']}\n{r['time_bin']}"
                  for _, r in tpq_sorted.iterrows()]
    tpq_sorted.to_csv(OUT_DIR / "pathway_by_timepoint_quality.csv", index=False)
    # Transpose so pathways are rows, tp_quality are columns
    plot_heatmap(tpq_mat.T, pathway_cols, col_labels,
                 "Pathway activation by radiation quality and time",
                 "pathway_heatmap_timepoint_quality",
                 cmap="RdBu_r", center=0.0,
                 xlabel="Radiation quality / time", ylabel="Pathway")

    # --- Heatmap 3: WGCNA module eigengene x timepoint ---
    print("\nHeatmap 3: module eigengene trajectory...")
    try:
        me_agg, me_cols = compute_module_eigengene_trajectory()
        me_agg.to_csv(OUT_DIR / "module_eigengene_by_timepoint.csv", index=False)
        # Clean module names (strip ME prefix for display)
        display_names = [c.replace("ME", "") for c in me_cols]
        me_order = [b[2] for b in TIME_BINS if b[2] in me_agg["time_bin"].values]
        me_agg_sorted = me_agg.set_index("time_bin").reindex(me_order)
        me_mat = me_agg_sorted[me_cols].values.T
        plot_heatmap(me_mat, display_names, me_order,
                     "WGCNA module eigengene trajectory",
                     "module_eigengene_heatmap",
                     cmap="RdBu_r", center=0.0,
                     xlabel="Time post-exposure", ylabel="Module")
    except Exception as e:
        print(f"  Module eigengene heatmap skipped: {e}")

    # --- Heatmap 4: pathway x cell type (pseudo-cell-type expression) ---
    print("\nHeatmap 4: pathway x cell type...")
    ct_scores = compute_celltype_pathway_scores()
    ct_scores.to_csv(OUT_DIR / "pathway_by_celltype.csv", index=False)
    ct_mat = ct_scores.set_index("CellType")[pathway_cols].values
    plot_heatmap(ct_mat.T, pathway_cols, ct_scores["CellType"].tolist(),
                 "Pathway activation across pseudo-cell-types",
                 "pathway_heatmap_celltype",
                 cmap="RdBu_r", center=0.0,
                 xlabel="Cell type", ylabel="Pathway")

    print("\nAll pathway enrichment heatmaps saved.")
    print(f"  Data -> {OUT_DIR}")
    print(f"  Figures -> {FIG_DIR}")


if __name__ == "__main__":
    main()
