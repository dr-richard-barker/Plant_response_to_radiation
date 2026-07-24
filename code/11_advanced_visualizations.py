#!/usr/bin/env python3
"""
11_advanced_visualizations.py
Advanced visualizations for the plant ionizing-radiation kinetic landscape:

  1A. 3D RRI surface over dose-time grid (GP-AE decoder sampling, 3 qualities)
  1B. Per-quality empirical 2D RRI plots (ground-truth for the 3D surface)
  1C. Animated signaling flow GIF (network graph, 5 timepoint frames)
  1D. 3D latent space scatter (colored by quality + time, + rotation GIF)
  1E. RRI component waterfall chart (stacked bar of latent/pathway/module)
  1F. Radiation quality prioritization table + bar chart

All figures: SVG + PNG, Liberation Sans font, Phylo color palette.
"""
import json, os, io
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.mplot3d import Axes3D
import torch, torch.nn as nn
import networkx as nx
import imageio.v2 as imageio
from scipy.spatial.distance import pdist

matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"

RESULTS_DIR = Path("/mnt/results/zenodo_bundle/results")
FIG_DIR = RESULTS_DIR / "figures"
DATA_DIR = Path("/mnt/results/zenodo_bundle/data")
LATENT_DIM = 16

# Phylo palette
PHYLO = {
    "blue": "#0279EE", "orange": "#FF9400", "green": "#75A025",
    "pink": "#FD9BED", "yellow": "#E9ED4C", "cream": "#FAF9F3",
    "light": "#ECE9E2", "black": "#000000", "gray": "#4D4D4D",
}
QUALITY_COLORS = {
    "gamma": "#0279EE", "GCR": "#FF9400", "HZE-Fe": "#FD9BED",
    "spaceflight-LEO": "#75A025", "UV-B": "#E9ED4C",
    "control": "#4D4D4D", "unknown": "#AAAAAA",
}

# ---- GP-AE model architecture (must match checkpoint) ----
class Encoder(nn.Module):
    def __init__(self, n_genes, latent_dim, dropout=0.2):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_genes, 256), nn.LayerNorm(256), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(256, 128), nn.LayerNorm(128), nn.ReLU(),
        )
        self.fc_mu = nn.Linear(128, latent_dim)
        self.fc_logvar = nn.Linear(128, latent_dim)
    def forward(self, x):
        return self.fc_mu(self.net(x)), self.fc_logvar(self.net(x))

class Decoder(nn.Module):
    def __init__(self, latent_dim, n_covariates, n_genes, dropout=0.15):
        super().__init__()
        self.cov_proj = nn.Sequential(
            nn.Linear(n_covariates, 32), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(32, latent_dim), nn.ReLU(),
        )
        self.net = nn.Sequential(
            nn.Linear(latent_dim * 2, 256), nn.LayerNorm(256), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(256, 256), nn.LayerNorm(256), nn.ReLU(),
            nn.Linear(256, n_genes),
        )
    def forward(self, z, covariates):
        return self.net(torch.cat([z, self.cov_proj(covariates)], dim=-1))

class CovariatePredictor(nn.Module):
    def __init__(self, latent_dim, n_covariates):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(latent_dim, 32), nn.ReLU(), nn.Linear(32, n_covariates))

class GPAutoencoder(nn.Module):
    def __init__(self, n_genes, latent_dim, n_covariates=3):
        super().__init__()
        self.encoder = Encoder(n_genes, latent_dim)
        self.decoder = Decoder(latent_dim, n_covariates, n_genes)
        self.cov_predictor = CovariatePredictor(latent_dim, n_covariates)
        self.latent_dim = latent_dim

# ---- Pathway gene sets (from script 08) ----
PATHWAY_GENES = {
    "dna_repair": ["AT1G07290","AT3G05120","AT5G40820","AT5G20850","AT3G48190",
                   "AT1G65470","AT3G02680","AT4G19730","AT5G24280","AT2G32750"],
    "oxidative_stress": ["AT1G02920","AT4G25130","AT5G18100","AT3G10920","AT1G20630",
                         "AT2G28190","AT3G26060","AT1G32940","AT5G03490","AT4G35090"],
    "hormone_signaling": ["AT2G14610","AT1G19670","AT1G02820","AT5G44420","AT2G39300",
                          "AT3G04720","AT1G13220","AT5G57050","AT3G24220","AT1G75800"],
}


def load_model_and_data():
    """Load GP-AE model, latent embeddings, RRI data, and expression."""
    ckpt = torch.load(DATA_DIR / "gp_ae_checkpoint.pt", map_location="cpu", weights_only=False)
    n_genes = len(ckpt["gene_names"])
    model = GPAutoencoder(n_genes, LATENT_DIM, n_covariates=3)
    model.load_state_dict(ckpt["model_state"])
    model.eval()

    gene_names = list(ckpt["gene_names"])
    cov_mean = np.array(ckpt["cov_mean"])
    cov_std = np.array(ckpt["cov_std"])
    gene_mean = np.array(ckpt["gene_mean"])
    gene_std = np.array(ckpt["gene_std"])

    # Latent embeddings
    le = pd.read_csv(RESULTS_DIR / "trajectories" / "latent_embeddings_annotated.csv", index_col=0)
    md = pd.read_csv(DATA_DIR / "metadata_master.csv")
    if "IsControl" not in le.columns:
        le = le.reset_index().merge(md[["SampleID", "IsControl"]], on="SampleID", how="left").set_index("SampleID")
    ic = le["IsControl"].astype(str).str.lower()
    le["IsControl"] = ic.isin(["true", "1", "yes"])
    for c in ["TimePostExposure_h", "AbsorbedDose_Gy", "LET_keV_um"]:
        if c in le.columns:
            le[c] = pd.to_numeric(le[c], errors="coerce")
    latent_cols = [c for c in le.columns if c.startswith("z")]

    # Control centroid
    z_ctrl = le.loc[le["IsControl"], latent_cols].values.mean(axis=0)
    # sigma for RRI_latent
    ctrl_z = le.loc[le["IsControl"], latent_cols].values
    if len(ctrl_z) > 1:
        sigma = np.median(pdist(ctrl_z))
    else:
        sigma = np.std(ctrl_z, axis=0).mean() + 1e-8

    # RRI per sample
    rri = pd.read_csv(RESULTS_DIR / "rri" / "rri_per_sample.csv")
    rri_tp = pd.read_csv(RESULTS_DIR / "rri" / "rri_per_timepoint.csv")

    # WGCNA sustained module genes
    modules = pd.read_csv(RESULTS_DIR / "wgcna" / "modules.csv")
    sustained_genes = modules[modules["Module"].isin(["turquoise", "grey"])]["Gene"].tolist()
    sustained_genes = [g for g in sustained_genes if g in gene_names]

    # Expression for control mean (module preservation)
    import anndata as ad
    adata = ad.read_h5ad("/workspace/expression_raw.h5ad")
    if "SampleID" in adata.obs.columns:
        obs_sids = adata.obs["SampleID"].astype(str).tolist()
    else:
        obs_sids = adata.obs_names.tolist()
    # Match controls
    ctrl_sids = set(md[md["IsControl"].astype(str).str.lower().isin(["true","1","yes"])]["SampleID"].astype(str))
    ctrl_indices = [i for i, s in enumerate(obs_sids) if s in ctrl_sids]
    if not ctrl_indices:
        ctrl_indices = list(range(min(10, len(obs_sids))))
    X = adata.X if not hasattr(adata.X, "toarray") else adata.X.toarray()
    ctrl_expr_mean = X[ctrl_indices].mean(axis=0)
    gene_to_idx = {g: i for i, g in enumerate(adata.var_names)}

    return {
        "model": model, "gene_names": gene_names, "cov_mean": cov_mean, "cov_std": cov_std,
        "gene_mean": gene_mean, "gene_std": gene_std,
        "le": le, "latent_cols": latent_cols, "z_ctrl": z_ctrl, "sigma": sigma,
        "rri": rri, "rri_tp": rri_tp,
        "sustained_genes": sustained_genes,
        "ctrl_expr_mean": ctrl_expr_mean, "gene_to_idx": gene_to_idx,
        "adata": adata, "X": X,
    }


# ---- 1A. 3D RRI Surface ----
def compute_rri_from_decoded(x_hat, data):
    """Compute 3-component RRI from a decoded expression vector."""
    model = data["model"]
    gene_names = data["gene_names"]
    z_ctrl = data["z_ctrl"]
    sigma = data["sigma"]

    # RRI_latent: re-encode decoded expression, compute distance from control centroid
    x_tensor = torch.tensor(x_hat, dtype=torch.float32).unsqueeze(0)
    with torch.no_grad():
        mu, _ = model.encoder(x_tensor)
    z_hat = mu.numpy()[0]
    rri_latent = np.exp(-np.sum((z_hat - z_ctrl) ** 2) / (sigma ** 2))

    # RRI_pathway: Shannon evenness of pathway gene expression
    gene_idx = {g: i for i, g in enumerate(gene_names)}
    pathway_scores = []
    for pname, genes in PATHWAY_GENES.items():
        idxs = [gene_idx[g] for g in genes if g in gene_idx]
        if idxs:
            pathway_scores.append(np.mean(x_hat[idxs]))
        else:
            pathway_scores.append(0.0)
    pathway_scores = np.array(pathway_scores)
    # Shift to positive, normalize to proportions
    ps_pos = pathway_scores - pathway_scores.min() + 1e-8
    props = ps_pos / ps_pos.sum()
    H = -np.sum(props * np.log(props))
    rri_pathway = H / np.log(len(props)) if len(props) > 1 else 1.0

    # RRI_module: correlation of sustained module genes with control mean
    sustained_idx = [gene_idx[g] for g in data["sustained_genes"] if g in gene_idx]
    if len(sustained_idx) > 2:
        decoded_sustained = x_hat[sustained_idx]
        # Map sustained genes to adata gene indices for control mean
        ctrl_sustained = []
        for g in data["sustained_genes"]:
            if g in data["gene_to_idx"] and g in gene_idx:
                ctrl_sustained.append(data["ctrl_expr_mean"][data["gene_to_idx"][g]])
        if len(ctrl_sustained) == len(sustained_idx):
            rri_module = max(0, np.corrcoef(decoded_sustained, np.array(ctrl_sustained))[0, 1])
        else:
            rri_module = 0.8  # fallback
    else:
        rri_module = 0.8

    rri = 0.50 * rri_latent + 0.25 * rri_pathway + 0.25 * rri_module
    return float(np.clip(rri, 0, 1)), float(rri_latent), float(rri_pathway), float(rri_module)


def plot_3d_rri_surface(data):
    """Generate 3D RRI surface for 3 radiation qualities via GP-AE decoder."""
    model = data["model"]
    z_ctrl = data["z_ctrl"]
    cov_mean = data["cov_mean"]
    cov_std = data["cov_std"]

    # Dose-time grid (log-spaced)
    dose_grid = np.logspace(-1, 2, 30)  # 0.1 to 100 Gy
    time_grid = np.logspace(-1, np.log10(24), 30)  # 0.1 to 24h

    # LET values per quality
    qualities = [
        ("gamma", 0.2, PHYLO["blue"]),
        ("HZE-Fe", 200.0, PHYLO["pink"]),
        ("GCR", 50.0, PHYLO["orange"]),  # mixed LET, use mid-value
    ]

    z_ctrl_tensor = torch.tensor(z_ctrl, dtype=torch.float32).unsqueeze(0)

    for qname, let_val, color in qualities:
        print(f"  Sampling RRI surface for {qname} (LET={let_val})...")
        rri_grid = np.zeros((len(dose_grid), len(time_grid)))
        rri_latent_grid = np.zeros_like(rri_grid)

        for i, dose in enumerate(dose_grid):
            for j, time_h in enumerate(time_grid):
                # Construct covariates: [log(dose), log(time), log(LET)]
                cov_raw = np.array([np.log(dose), np.log(time_h), np.log(let_val)])
                cov_stdized = (cov_raw - cov_mean) / cov_std
                cov_tensor = torch.tensor(cov_stdized, dtype=torch.float32).unsqueeze(0)

                with torch.no_grad():
                    x_hat = model.decoder(z_ctrl_tensor, cov_tensor).numpy()[0]

                rri_val, rri_lat, _, _ = compute_rri_from_decoded(x_hat, data)
                rri_grid[i, j] = rri_val
                rri_latent_grid[i, j] = rri_lat

        # Plot 3D surface
        D, T = np.meshgrid(dose_grid, time_grid, indexing="ij")
        fig = plt.figure(figsize=(9, 7))
        ax = fig.add_subplot(111, projection="3d")

        # Color map: Phylo blue (low RRI) -> cream (mid) -> green (high RRI)
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            "phylo_rri", [PHYLO["orange"], PHYLO["cream"], PHYLO["green"]]
        )
        norm = matplotlib.colors.Normalize(vmin=0.5, vmax=1.0)
        surf = ax.plot_surface(np.log10(D), np.log10(T), rri_grid, cmap=cmap, norm=norm,
                               alpha=0.9, edgecolor="none")

        # Overlay measured points
        rri_df = data["rri"]
        q_mask = rri_df["RadiationQuality"] == qname
        q_data = rri_df[q_mask & ~rri_df["IsControl"]] if "IsControl" in rri_df.columns else rri_df[q_mask]
        if len(q_data) > 0:
            for _, row in q_data.iterrows():
                d = row.get("AbsorbedDose_Gy", np.nan)
                t = row.get("TimePostExposure_h", np.nan)
                r = row.get("RRI", np.nan)
                if not np.isnan(d) and not np.isnan(t) and not np.isnan(r) and d > 0 and t > 0:
                    ax.scatter(np.log10(d), np.log10(t), r, c="black", s=30,
                              edgecolors="white", linewidths=0.5, zorder=5)

        ax.set_xlabel("log10(Dose Gy)", fontsize=9, labelpad=8)
        ax.set_ylabel("log10(Time h)", fontsize=9, labelpad=8)
        ax.set_zlabel("RRI", fontsize=9, labelpad=4)
        ax.set_zlim(0.4, 1.0)
        ax.set_title(f"RRI Surface — {qname} (LET={let_val} keV/μm)\n"
                     f"GP-AE model prediction (black dots = measured data)",
                     fontsize=10, pad=15)
        fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, label="RRI")
        ax.view_init(elev=25, azim=225)
        fig.tight_layout()
        fig.savefig(FIG_DIR / f"rri_surface_3d_{qname.lower().replace('-','_')}.svg", format="svg")
        fig.savefig(FIG_DIR / f"rri_surface_3d_{qname.lower().replace('-','_')}.png", dpi=150)
        plt.close(fig)
        print(f"    Saved rri_surface_3d_{qname.lower().replace('-','_')}.svg/png")

        # Save data
        surface_df = pd.DataFrame({
            "dose_Gy": D.flatten(), "time_h": T.flatten(),
            "RRI": rri_grid.flatten(), "RRI_latent": rri_latent_grid.flatten(),
            "quality": qname, "LET": let_val,
        })
        surface_df.to_csv(RESULTS_DIR / "rri" / f"rri_surface_{qname.lower().replace('-','_')}.csv", index=False)


# ---- 1B. Empirical 2D RRI Plots ----
def plot_empirical_2d(data):
    """Side-by-side RRI vs time (gamma) and RRI vs dose (GCR)."""
    rri = data["rri"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Left: RRI vs time for gamma
    gamma = rri[(rri["RadiationQuality"] == "gamma") & (~rri["IsControl"])]
    gamma = gamma[gamma["TimePostExposure_h"] > 0]
    if len(gamma) > 0:
        ax1.scatter(gamma["TimePostExposure_h"], gamma["RRI"], c=PHYLO["blue"],
                   alpha=0.6, s=40, edgecolors="white", linewidths=0.5, zorder=3)
        # Binned mean
        for tbin in [(0, 0.5), (0.5, 2), (2, 6), (6, 12), (12, 30)]:
            mask = (gamma["TimePostExposure_h"] >= tbin[0]) & (gamma["TimePostExposure_h"] < tbin[1])
            if mask.sum() > 0:
                ax1.errorbar(np.mean([tbin[0], tbin[1]]), gamma.loc[mask, "RRI"].mean(),
                           yerr=gamma.loc[mask, "RRI"].std(), fmt="o-", c=PHYLO["black"],
                           markersize=8, linewidth=2, capsize=4, zorder=4)
        ax1.set_xscale("log")
        ax1.set_xlabel("Time post-exposure (h)", fontsize=10)
        ax1.set_ylabel("RRI", fontsize=10)
        ax1.set_title("Gamma (100 Gy) — RRI vs time", fontsize=11)
        ax1.axhline(y=0.77, color=PHYLO["green"], linestyle="--", alpha=0.5, label="control mean")
        ax1.legend(fontsize=8)
        ax1.set_ylim(0.4, 1.0)

    # Right: RRI vs dose for GCR
    gcr = rri[(rri["RadiationQuality"] == "GCR") & (~rri["IsControl"])]
    if len(gcr) > 0:
        ax2.scatter(gcr["AbsorbedDose_Gy"], gcr["RRI"], c=PHYLO["orange"],
                   alpha=0.7, s=60, edgecolors="white", linewidths=0.5, zorder=3)
        for dose in sorted(gcr["AbsorbedDose_Gy"].unique()):
            mask = gcr["AbsorbedDose_Gy"] == dose
            if mask.sum() > 0:
                ax2.errorbar(dose, gcr.loc[mask, "RRI"].mean(),
                           yerr=gcr.loc[mask, "RRI"].std(), fmt="s-", c=PHYLO["black"],
                           markersize=10, linewidth=2, capsize=4, zorder=4)
        ax2.set_xlabel("Absorbed dose (Gy)", fontsize=10)
        ax2.set_ylabel("RRI", fontsize=10)
        ax2.set_title("GCR (3h post-exposure) — RRI vs dose", fontsize=11)
        ax2.axhline(y=0.77, color=PHYLO["green"], linestyle="--", alpha=0.5, label="control mean")
        ax2.legend(fontsize=8)
        ax2.set_ylim(0.4, 1.0)

    fig.suptitle("Empirical RRI: measured data (not model interpolation)", fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "rri_empirical_2d.svg", format="svg")
    fig.savefig(FIG_DIR / "rri_empirical_2d.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  Saved rri_empirical_2d.svg/png")


# ---- 1C. Animated Signaling Flow GIF ----
def create_signaling_gif():
    """Network graph animation of CellChat signaling across 5 timepoints."""
    sf = pd.read_csv(RESULTS_DIR / "cellchat" / "signaling_flow_per_timepoint.csv")
    origin = pd.read_csv(RESULTS_DIR / "cellchat" / "signaling_origin_summary.csv")

    cell_types = sorted(origin["CellType"].unique())
    timepoints = sorted(origin["Time"].unique())

    # Fixed node positions: root at bottom, leaf/shoot at top, stress in center
    pos = {}
    root_types = [c for c in cell_types if c.startswith("root_")]
    leaf_types = [c for c in cell_types if c.startswith("leaf_") or c == "shoot_apical_meristem"]
    stress_types = [c for c in cell_types if c.startswith("hormone_response_") or
                    c in ("dna_damage_response", "oxidative_stress")]
    other_types = [c for c in cell_types if c not in root_types + leaf_types + stress_types]

    # Arrange in clusters
    for i, ct in enumerate(root_types):
        angle = np.pi * (i / max(1, len(root_types) - 1))  # 0 to pi
        pos[ct] = (np.cos(angle) * 2, -3 + np.sin(angle) * 1.5)
    for i, ct in enumerate(leaf_types):
        angle = np.pi * (i / max(1, len(leaf_types) - 1))
        pos[ct] = (np.cos(angle) * 2, 3 + np.sin(angle) * 1.5)
    for i, ct in enumerate(stress_types):
        angle = 2 * np.pi * (i / max(1, len(stress_types)))
        pos[ct] = (np.cos(angle) * 2.5, np.sin(angle) * 2.5)
    for i, ct in enumerate(other_types):
        angle = 2 * np.pi * (i / max(1, len(other_types)))
        pos[ct] = (np.cos(angle) * 1.5, np.sin(angle) * 1.0)

    frames = []
    for tp in timepoints:
        tp_sf = sf[sf["Time"] == tp].copy()
        tp_origin = origin[origin["Time"] == tp].copy()

        # Normalize outgoing signal for node sizing
        max_out = tp_origin["OutgoingSignal"].max()
        tp_origin["node_size"] = 100 + 500 * (tp_origin["OutgoingSignal"] / max_out)

        # Top 50 edges by signal strength
        tp_sf = tp_sf.sort_values("SignalStrength", ascending=False).head(50)
        max_str = tp_sf["SignalStrength"].max()

        fig, ax = plt.subplots(figsize=(10, 10))
        G = nx.DiGraph()
        G.add_nodes_from(cell_types)

        # Draw edges
        for _, row in tp_sf.iterrows():
            src, tgt = row["Source"], row["Target"]
            if src in pos and tgt in pos:
                G.add_edge(src, tgt, weight=row["SignalStrength"])
                width = 0.5 + 4 * (row["SignalStrength"] / max_str)
                alpha = 0.2 + 0.6 * max(0, 1 - row["pvalue"])
                ax.annotate("", xy=pos[tgt], xytext=pos[src],
                           arrowprops=dict(arrowstyle="-|>", color=PHYLO["gray"],
                                          alpha=alpha, lw=width, connectionstyle="arc3,rad=0.1"))

        # Draw nodes
        for _, row in tp_origin.iterrows():
            ct = row["CellType"]
            if ct in pos:
                size = row["node_size"]
                # Color by outgoing signal: low=blue, high=orange
                norm_val = row["OutgoingSignal"] / max_out
                color = PHYLO["blue"] if norm_val < 0.5 else PHYLO["orange"]
                ax.scatter(pos[ct][0], pos[ct][1], s=size, c=color,
                          alpha=0.8, edgecolors="white", linewidths=1.0, zorder=5)
                # Label
                label = ct.replace("hormone_response_", "hr_").replace("root_", "r_").replace("leaf_", "l_")
                ax.annotate(label, pos[ct], fontsize=6, ha="center", va="bottom",
                           xytext=(0, 8), textcoords="offset points", zorder=6)

        ax.set_xlim(-4, 4)
        ax.set_ylim(-5, 5)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(f"Signaling Flow — {tp}\n(node size = outgoing signal, edge width = flow strength)",
                    fontsize=12, pad=15)

        # Save frame
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
        buf.seek(0)
        frames.append(imageio.imread(buf))
        fig.savefig(FIG_DIR / f"signaling_frame_{tp.replace('-','_').replace('.','_')}.png",
                   dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Frame {tp} rendered")

    # Save GIF
    gif_path = FIG_DIR / "signaling_flow_animation.gif"
    imageio.mimsave(str(gif_path), frames, duration=1.5, loop=0)
    print(f"  Saved signaling_flow_animation.gif ({len(frames)} frames)")


# ---- 1D. 3D Latent Space Scatter ----
def plot_latent_3d(data):
    """3D scatter of samples in first 3 latent dimensions."""
    le = data["le"]
    latent_cols = data["latent_cols"]

    # Use first 3 latent dims
    z0, z1, z2 = le[latent_cols[0]], le[latent_cols[1]], le[latent_cols[2]]

    # View A: colored by radiation quality
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    for qname in le["RadiationQuality"].dropna().unique():
        mask = le["RadiationQuality"] == qname
        color = QUALITY_COLORS.get(qname, PHYLO["gray"])
        ax.scatter(z0[mask], z1[mask], z2[mask], c=color, s=40, alpha=0.7,
                  edgecolors="white", linewidths=0.3, label=qname)
    # Controls
    ctrl_mask = le["IsControl"]
    ax.scatter(z0[ctrl_mask], z1[ctrl_mask], z2[ctrl_mask], c=PHYLO["gray"],
              s=30, alpha=0.4, marker="x", label="control")

    ax.set_xlabel(f"{latent_cols[0]}", fontsize=9)
    ax.set_ylabel(f"{latent_cols[1]}", fontsize=9)
    ax.set_zlabel(f"{latent_cols[2]}", fontsize=9)
    ax.set_title("GP-AE latent space (first 3 dimensions)\ncolored by radiation quality", fontsize=11)
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(0, 1))
    ax.view_init(elev=20, azim=30)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "latent_3d_scatter.svg", format="svg")
    fig.savefig(FIG_DIR / "latent_3d_scatter.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  Saved latent_3d_scatter.svg/png")

    # Rotation GIF — fixed figure size (no bbox_inches='tight' to keep frames uniform)
    frames = []
    fig_rot = plt.figure(figsize=(8, 7))
    ax_rot = fig_rot.add_subplot(111, projection="3d")
    for qname in le["RadiationQuality"].dropna().unique():
        mask = le["RadiationQuality"] == qname
        color = QUALITY_COLORS.get(qname, PHYLO["gray"])
        ax_rot.scatter(z0[mask], z1[mask], z2[mask], c=color, s=30, alpha=0.7,
                       edgecolors="white", linewidths=0.3)
    ax_rot.scatter(z0[ctrl_mask], z1[ctrl_mask], z2[ctrl_mask], c=PHYLO["gray"],
                   s=20, alpha=0.3, marker="x")
    ax_rot.set_xlabel(latent_cols[0], fontsize=8)
    ax_rot.set_ylabel(latent_cols[1], fontsize=8)
    ax_rot.set_zlabel(latent_cols[2], fontsize=8)

    for angle in range(0, 360, 10):
        ax_rot.view_init(elev=20, azim=angle)
        ax_rot.set_title(f"Latent space rotation ({angle}°)", fontsize=10)
        buf = io.BytesIO()
        fig_rot.savefig(buf, format="png", dpi=80)
        buf.seek(0)
        frames.append(imageio.imread(buf))

    plt.close(fig_rot)
    imageio.mimsave(str(FIG_DIR / "latent_3d_rotation.gif"), frames, duration=0.3, loop=0)
    print("  Saved latent_3d_rotation.gif")


# ---- 1E. RRI Component Waterfall ----
def plot_rri_waterfall(data):
    """Stacked bar showing RRI component contributions per timepoint."""
    rri_tp = data["rri_tp"]

    time_bins = ["0-0.5h", "0.5-2h", "2-6h", "6-12h", "12-30h"]
    rri_tp = rri_tp.set_index("time_bin").reindex(time_bins).reset_index()

    fig, ax = plt.subplots(figsize=(9, 5))

    latent_contrib = 0.50 * rri_tp["RRI_latent_mean"]
    pathway_contrib = 0.25 * rri_tp["RRI_pathway_mean"]
    module_contrib = 0.25 * rri_tp["RRI_module_mean"]

    x = range(len(time_bins))
    ax.bar(x, latent_contrib, color=PHYLO["blue"], label="Latent (0.50×)", width=0.6)
    ax.bar(x, pathway_contrib, bottom=latent_contrib, color=PHYLO["green"],
          label="Pathway (0.25×)", width=0.6)
    ax.bar(x, module_contrib, bottom=latent_contrib + pathway_contrib,
          color=PHYLO["orange"], label="Module (0.25×)", width=0.6)

    # Total RRI line
    totals = latent_contrib + pathway_contrib + module_contrib
    ax.plot(x, totals, "o-", color=PHYLO["black"], linewidth=2, markersize=8,
           label="Composite RRI", zorder=5)

    # Control mean reference
    ax.axhline(y=0.77, color=PHYLO["green"], linestyle="--", alpha=0.5, label="Control mean")

    # Annotate nadir
    nadir_idx = totals.idxmin()
    ax.annotate("nadir", xy=(nadir_idx, totals[nadir_idx]),
               xytext=(nadir_idx + 0.5, totals[nadir_idx] - 0.08),
               arrowprops=dict(arrowstyle="->", color=PHYLO["orange"]),
               fontsize=9, color=PHYLO["orange"])

    ax.set_xticks(list(x))
    ax.set_xticklabels(time_bins, fontsize=9)
    ax.set_ylabel("RRI contribution", fontsize=10)
    ax.set_title("RRI component waterfall — what drives the nadir at 2-6h", fontsize=11)
    ax.legend(fontsize=8, loc="lower right")
    ax.set_ylim(0, 1.05)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "rri_component_waterfall.svg", format="svg")
    fig.savefig(FIG_DIR / "rri_component_waterfall.png", dpi=150)
    plt.close(fig)
    print("  Saved rri_component_waterfall.svg/png")


# ---- 1F. Radiation Quality Prioritization ----
def compute_prioritization():
    """Score and rank 11 missing radiation qualities by data availability."""
    # Data from PLAN.md OSDR taxonomy
    qualities = [
        {"quality": "proton", "osdr_studies": 12, "let_class": "low",
         "let_gap_score": 0.25, "mission_score": 1.00,
         "note": "Dominant GCR component; NSRL + cyclotron. Fills low-LET particle gap."},
        {"quality": "cosmic-mixed", "osdr_studies": 10, "let_class": "mixed",
         "let_gap_score": 0.75, "mission_score": 0.80,
         "note": "Mixed-LET; spaceflight-relevant. Usually maps to GCR+spaceflight."},
        {"quality": "helium", "osdr_studies": 6, "let_class": "low",
         "let_gap_score": 0.25, "mission_score": 0.70,
         "note": "Alpha/He ions; GCR component. Low-LET particle."},
        {"quality": "beta", "osdr_studies": 5, "let_class": "low",
         "let_gap_score": 0.25, "mission_score": 0.20,
         "note": "Beta emitters; often dosimetry context. Minimal spaceflight relevance."},
        {"quality": "neutron", "osdr_studies": 2, "let_class": "mixed",
         "let_gap_score": 0.75, "mission_score": 0.60,
         "note": "Mixed-LET; secondary radiation hazard. Fills LET gap."},
        {"quality": "HZE-Si", "osdr_studies": 1, "let_class": "high",
         "let_gap_score": 0.50, "mission_score": 0.90,
         "note": "Silicon ions; major GCR component. High-LET. Reprocess OSD-658."},
        {"quality": "HZE-O", "osdr_studies": 1, "let_class": "high",
         "let_gap_score": 0.50, "mission_score": 0.85,
         "note": "Oxygen ions; GCR component. High-LET. Reprocess OSD-658."},
        {"quality": "X-ray", "osdr_studies": 2, "let_class": "low",
         "let_gap_score": 0.25, "mission_score": 0.10,
         "note": "Ground-based; not spaceflight-relevant. Low-LET photon."},
        {"quality": "UV-A", "osdr_studies": 2, "let_class": "n/a",
         "let_gap_score": 0.00, "mission_score": 0.30,
         "note": "ISS UV exposure; non-ionizing. No LET."},
        {"quality": "solar-particle-event", "osdr_studies": 0, "let_class": "mixed",
         "let_gap_score": 0.75, "mission_score": 1.00,
         "note": "Acute mission hazard; 0 OSDR studies. Needs dedicated NSRL campaign."},
        {"quality": "UV-C", "osdr_studies": 0, "let_class": "n/a",
         "let_gap_score": 0.00, "mission_score": 0.10,
         "note": "Germicidal; not spaceflight-relevant. 0 OSDR studies. GEO only."},
    ]

    max_osdr = 12  # proton has the most
    for q in qualities:
        q["availability_score"] = min(1.0, q["osdr_studies"] / max_osdr)
        q["composite_score"] = (0.50 * q["availability_score"] +
                                0.25 * q["let_gap_score"] +
                                0.25 * q["mission_score"])

    # Sort by composite score
    qualities.sort(key=lambda x: x["composite_score"], reverse=True)
    for i, q in enumerate(qualities):
        q["rank"] = i + 1
        if q["composite_score"] >= 0.6:
            q["tier"] = "Tier 1 (acquire now)"
        elif q["composite_score"] >= 0.35:
            q["tier"] = "Tier 2 (acquire if feasible)"
        elif q["composite_score"] >= 0.15:
            q["tier"] = "Tier 3 (opportunistic)"
        else:
            q["tier"] = "Tier 4 (requires new experiments)"

    df = pd.DataFrame(qualities)
    df.to_csv(RESULTS_DIR / "pathway_enrichment" / "radiation_quality_prioritization.csv", index=False)

    # Bar chart
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = []
    for s in df["composite_score"]:
        if s >= 0.6:
            colors.append(PHYLO["green"])
        elif s >= 0.35:
            colors.append(PHYLO["yellow"])
        elif s >= 0.15:
            colors.append(PHYLO["orange"])
        else:
            colors.append(PHYLO["gray"])

    y_pos = range(len(df))
    ax.barh(y_pos, df["composite_score"], color=colors, edgecolor="white", height=0.6)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(df["quality"], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Composite priority score", fontsize=10)
    ax.set_title("Radiation quality prioritization for future data acquisition\n"
                 "(50% data availability + 25% LET gap-filling + 25% mission relevance)",
                 fontsize=11)
    ax.set_xlim(0, 1.0)

    # Annotate scores
    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(row["composite_score"] + 0.02, i,
               f"{row['composite_score']:.2f} ({row['tier']})", fontsize=7, va="center")

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=PHYLO["green"], label="Tier 1 (acquire now)"),
        Patch(facecolor=PHYLO["yellow"], label="Tier 2 (if feasible)"),
        Patch(facecolor=PHYLO["orange"], label="Tier 3 (opportunistic)"),
        Patch(facecolor=PHYLO["gray"], label="Tier 4 (new experiments)"),
    ]
    ax.legend(handles=legend_elements, fontsize=7, loc="lower right")

    fig.tight_layout()
    fig.savefig(FIG_DIR / "radiation_quality_prioritization.svg", format="svg")
    fig.savefig(FIG_DIR / "radiation_quality_prioritization.png", dpi=150)
    plt.close(fig)
    print("  Saved radiation_quality_prioritization.svg/png")
    print(f"  Prioritization table: {len(df)} qualities ranked")
    for _, row in df.iterrows():
        print(f"    #{row['rank']} {row['quality']}: {row['composite_score']:.2f} — {row['tier']}")


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    print("Loading model and data...")
    data = load_model_and_data()
    print(f"  Model loaded, {len(data['gene_names'])} genes, {data['le'].shape[0]} samples")
    print(f"  Control centroid computed from {data['le']['IsControl'].sum()} controls")

    print("\n1A. 3D RRI surface (GP-AE decoder sampling)...")
    plot_3d_rri_surface(data)

    print("\n1B. Empirical 2D RRI plots...")
    plot_empirical_2d(data)

    print("\n1C. Animated signaling flow GIF...")
    create_signaling_gif()

    print("\n1D. 3D latent space scatter...")
    plot_latent_3d(data)

    print("\n1E. RRI component waterfall...")
    plot_rri_waterfall(data)

    print("\n1F. Radiation quality prioritization...")
    compute_prioritization()

    print("\nAll advanced visualizations saved.")
    print(f"  Figures -> {FIG_DIR}")


if __name__ == "__main__":
    main()
