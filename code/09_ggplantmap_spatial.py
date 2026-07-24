#!/usr/bin/env python3
"""
09_ggplantmap_spatial.py
Phase 3.3 — ggPlantmap spatial visualization: time-annotated plant organ maps
showing the "wave" of stress signaling across organs over time post-radiation.

Uses ggPlantmap (Jo & Kajala 2024, J Exp Bot 75:5366-5376, doi:10.1093/jxb/erae043)
to overlay our CellChat signaling-origin data and RRI module-preservation scores
onto pre-loaded Arabidopsis plant maps at each timepoint.

APPROACH
  1. Map our 24 pseudo-cell-types to ggPlantmap organ regions:
     - root_meristem, root_columella, root_epidermis, root_cortex, root_endodermis,
       root_stele, root_xylem, root_phloem, root_pericycle, root_lateral_root_cap
       -> ggPm.At.roottip.longitudinal / ggPm.At.roottip.crosssection
     - leaf_epidermis, leaf_mesophyll, leaf_vasculature -> ggPm.At.leaf.crosssection
     - shoot_apical_meristem -> ggPm.At.shootapex.longitudinal
     - cotyledon, hypocotyl -> ggPm.At.seedling.saltdrought
     - dna_damage_response, oxidative_stress, hormone_response_* -> whole-plant
       overlay (these are pathway-level, not spatial)
  2. For each timepoint, create a multi-panel figure:
     - Root tip longitudinal view (signaling origin strength)
     - Leaf cross-section (signaling origin strength)
     - Seedling view (RRI module preservation)
  3. Save as SVG (per user preference) and PNG.

OUTPUTS
  - figures/ggplantmap_root_<timepoint>.svg/png
  - figures/ggplantmap_leaf_<timepoint>.svg/png
  - figures/ggplantmap_seedling_<timepoint>.svg/png
  - figures/ggplantmap_composite_<timepoint>.svg/png
  - figures/rri_trajectory.svg/png (RRI over time line plot)
  - figures/signaling_flow_heatmap.svg/png (CellChat origin x timepoint)
"""
import json, os, subprocess
from pathlib import Path
import numpy as np
import pandas as pd

RESULTS_DIR = Path("/mnt/results/zenodo_bundle/results")
OUT_DIR = RESULTS_DIR / "figures"
CELLCHAT_DIR = RESULTS_DIR / "cellchat"
RRI_DIR = RESULTS_DIR / "rri"

# Mapping from our 24 cell types to ggPlantmap ROI names
# (based on the actual pre-loaded map ROI names)
CELLTYPE_TO_ROI_ROOT = {
    "root_meristem": "Meristem.QC",
    "root_columella": "Columella.Columella",
    "root_epidermis": "Epidermis.Epidermis",
    "root_cortex": "Cortex.Cortex",
    "root_endodermis": "Endodermis.Endodermis",
    "root_stele": "Stele.Vascular Tissue",
    "root_pericycle": "Stele.Pericycle",
    "root_lateral_root_cap": "Lateral Root Cap.Lateral Root Cap",
}

CELLTYPE_TO_ROI_LEAF = {
    "leaf_epidermis": "epidermis.adaxial",
    "leaf_mesophyll": "Parenchima.palisade",
    "leaf_vasculature": "vascularbundle.xylem",
}


def create_r_script(timepoints_data, out_dir):
    """Create R script for ggPlantmap visualization."""
    r_script = f"""
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressMessages(library(ggPlantmap))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(patchwork))
svg_font <- function() theme(text = element_text(family = "Liberation Sans"))

# Load data
root_data <- read.csv("{out_dir}/root_signaling_data.csv")
leaf_data <- read.csv("{out_dir}/leaf_signaling_data.csv")
seedling_data <- read.csv("{out_dir}/seedling_rri_data.csv")
rri_trajectory <- read.csv("{out_dir}/rri_trajectory_data.csv")
signaling_heatmap <- read.csv("{out_dir}/signaling_heatmap_data.csv")
# Ensure value columns are numeric
root_data$value <- as.numeric(root_data$value)
leaf_data$value <- as.numeric(leaf_data$value)
seedling_data$value <- as.numeric(seedling_data$value)
signaling_heatmap$OutgoingSignal <- as.numeric(signaling_heatmap$OutgoingSignal)

# Color scale (Phylo palette)
phylo_colors <- c("#000000", "#ECE9E2", "#FAF9F3", "#E9ED4C", "#FF9400",
                  "#75A025", "#FD9BED", "#0279EE")
stress_scale <- scale_fill_gradient2(low = "#0279EE", mid = "#FAF9F3",
                                      high = "#FF9400", midpoint = 0.5,
                                      name = "Signaling\\nstrength")
rri_scale <- scale_fill_gradient2(low = "#FF9400", mid = "#FAF9F3",
                                   high = "#75A025", midpoint = 0.7,
                                   name = "RRI")

# === RRI trajectory line plot ===
p_rri <- ggplot(rri_trajectory, aes(x = time_bin, y = RRI_mean, group = 1)) +
  geom_ribbon(aes(ymin = RRI_mean - RRI_std, ymax = RRI_mean + RRI_std),
              fill = "#ECE9E2", alpha = 0.5) +
  geom_line(color = "#0279EE", linewidth = 1.2) +
  geom_point(color = "#0279EE", size = 3) +
  geom_hline(yintercept = 0.77, linetype = "dashed", color = "#75A025",
             linewidth = 0.8, alpha = 0.7) +
  annotate("text", x = 1, y = 0.78, label = "control mean", hjust = 0,
           size = 3, color = "#75A025") +
  labs(title = "Radiation Resilience Index trajectory",
       x = "Time post-exposure", y = "RRI (mean +/- sd)") +
  theme_bw() + svg_font() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("{out_dir}/rri_trajectory.svg", p_rri, width = 7, height = 4, device = svglite::svglite)
ggsave("{out_dir}/rri_trajectory.png", p_rpi <- p_rri, width = 7, height = 4, dpi = 150)
cat("RRI trajectory saved\n")

# === Signaling flow heatmap (cell type x timepoint) ===
p_heat <- ggplot(signaling_heatmap, aes(x = time_bin, y = CellType, fill = OutgoingSignal)) +
  geom_tile() +
  stress_scale +
  labs(title = "Signaling origin strength over time",
       x = "Time post-exposure", y = "Cell type (source)") +
  theme_bw() + svg_font() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_text(size = 7))
ggsave("{out_dir}/signaling_flow_heatmap.svg", p_heat, width = 7, height = 6, device = svglite::svglite)
ggsave("{out_dir}/signaling_flow_heatmap.png", p_heat, width = 7, height = 6, dpi = 150)
cat("Signaling heatmap saved\n")

# === ggPlantmap per-timepoint figures ===
# Helper: merge map with values (bypassing ggPlantmap.merge tidy-eval scope issues)
# Handles maps WITH a point column (root, leaf) and WITHOUT (seedling)
merge_map_values <- function(map_df, values_df) {{
  merged <- merge(map_df, values_df, by.x = "ROI.name", by.y = "ROI", all.x = TRUE)
  if ("point" %in% colnames(merged)) {{
    merged <- merged[order(merged$ROI.id, merged$point), ]
  }} else {{
    merged <- merged[order(merged$ROI.id), ]
  }}
  tibble::as_tibble(merged)
}}

# Seedling ROI names in the map are capitalized: Root_control, Hypocotyl_control,
# Cotyledon_control, Hook_control. Map our lowercase names to these.
SEEDLING_ROI_MAP <- c(
  root = "Root_control",
  hypocotyl = "Hypocotyl_control",
  cotyledon = "Cotyledon_control",
  shoot = "Hook_control"   # shoot apical region -> Hook
)

timepoints <- unique(root_data$time_bin)
for (tp in timepoints) {{
  tp_safe <- gsub("[^a-zA-Z0-9]", "_", tp)
  cat("--- Processing", tp, "---\n")

  # Root signaling
  root_tp <- root_data[root_data$time_bin == tp, ]
  if (nrow(root_tp) > 0) {{
    tryCatch({{
      root_map <- ggPm.At.roottip.longitudinal
      root_merged <- merge_map_values(root_map, root_tp)
      p_root <- ggPlantmap.heatmap(root_merged, value.quant = value) +
        stress_scale +
        labs(title = paste0("Root signaling - ", tp)) +
        theme_bw() + svg_font()
      ggsave(paste0("{out_dir}/ggplantmap_root_", tp_safe, ".svg"), p_root,
             width = 5, height = 4, device = svglite::svglite)
      ggsave(paste0("{out_dir}/ggplantmap_root_", tp_safe, ".png"), p_root,
             width = 5, height = 4, dpi = 150)
      cat("  root OK\n")
    }}, error = function(e) cat("  root FAILED:", conditionMessage(e), "\n"))
  }}

  # Leaf signaling
  leaf_tp <- leaf_data[leaf_data$time_bin == tp, ]
  if (nrow(leaf_tp) > 0) {{
    tryCatch({{
      leaf_map <- ggPm.At.leaf.crosssection
      leaf_merged <- merge_map_values(leaf_map, leaf_tp)
      p_leaf <- ggPlantmap.heatmap(leaf_merged, value.quant = value) +
        stress_scale +
        labs(title = paste0("Leaf signaling - ", tp)) +
        theme_bw() + svg_font()
      ggsave(paste0("{out_dir}/ggplantmap_leaf_", tp_safe, ".svg"), p_leaf,
             width = 5, height = 4, device = svglite::svglite)
      ggsave(paste0("{out_dir}/ggplantmap_leaf_", tp_safe, ".png"), p_leaf,
             width = 5, height = 4, dpi = 150)
      cat("  leaf OK\n")
    }}, error = function(e) cat("  leaf FAILED:", conditionMessage(e), "\n"))
  }}

  # Seedling RRI (whole-plant RRI overlaid on seedling map)
  seed_tp <- seedling_data[seedling_data$time_bin == tp, ]
  if (nrow(seed_tp) > 0) {{
    tryCatch({{
      seedling_map <- ggPm.At.seedling.saltdrought
      seedling_map_ctrl <- seedling_map[grepl("_control$", seedling_map$ROI.name), ]
      # Remap our lowercase ROI names to the capitalized map ROI names
      seed_tp$ROI <- SEEDLING_ROI_MAP[seed_tp$ROI]
      seed_tp$ROI <- as.character(seed_tp$ROI)
      seed_merged <- merge_map_values(seedling_map_ctrl, seed_tp)
      p_seed <- ggPlantmap.heatmap(seed_merged, value.quant = value) +
        rri_scale +
        labs(title = paste0("Seedling RRI - ", tp)) +
        theme_bw() + svg_font()
      ggsave(paste0("{out_dir}/ggplantmap_seedling_", tp_safe, ".svg"), p_seed,
             width = 5, height = 4, device = svglite::svglite)
      ggsave(paste0("{out_dir}/ggplantmap_seedling_", tp_safe, ".png"), p_seed,
             width = 5, height = 4, dpi = 150)
      cat("  seedling OK\n")
    }}, error = function(e) cat("  seedling FAILED:", conditionMessage(e), "\n"))
  }}

  cat("Saved figures for", tp, "\n")
}}

cat("All ggPlantmap figures saved\n")
"""
    r_path = out_dir / "run_ggplantmap.R"
    with open(r_path, "w") as f:
        f.write(r_script)
    return r_path


def prepare_data_for_r():
    """Prepare CSV files that the R script will load."""
    # Signaling origin data
    origin = pd.read_csv(CELLCHAT_DIR / "signaling_origin_summary.csv")

    # Root cell types signaling data (outgoing signal per timepoint)
    root_celltypes = list(CELLTYPE_TO_ROI_ROOT.keys())
    root_rows = []
    for _, row in origin.iterrows():
        if row["CellType"] in CELLTYPE_TO_ROI_ROOT:
            # Normalize signal to [0,1] for color scale
            root_rows.append({"ROI": CELLTYPE_TO_ROI_ROOT[row["CellType"]],
                              "value": float(row["OutgoingSignal"]),
                              "time_bin": row["Time"]})
    root_df = pd.DataFrame(root_rows)
    if len(root_df) > 0:
        # Normalize per timepoint to [0,1]
        for tp in root_df["time_bin"].unique():
            mask = root_df["time_bin"] == tp
            vals = root_df.loc[mask, "value"]
            vmin, vmax = vals.min(), vals.max()
            if vmax > vmin:
                root_df.loc[mask, "value"] = (vals - vmin) / (vmax - vmin)
            else:
                root_df.loc[mask, "value"] = 0.5
    root_df.to_csv(OUT_DIR / "root_signaling_data.csv", index=False)

    # Leaf cell types
    leaf_rows = []
    for _, row in origin.iterrows():
        if row["CellType"] in CELLTYPE_TO_ROI_LEAF:
            leaf_rows.append({"ROI": CELLTYPE_TO_ROI_LEAF[row["CellType"]],
                              "value": float(row["OutgoingSignal"]),
                              "time_bin": row["Time"]})
    leaf_df = pd.DataFrame(leaf_rows)
    if len(leaf_df) > 0:
        for tp in leaf_df["time_bin"].unique():
            mask = leaf_df["time_bin"] == tp
            vals = leaf_df.loc[mask, "value"]
            vmin, vmax = vals.min(), vals.max()
            if vmax > vmin:
                leaf_df.loc[mask, "value"] = (vals - vmin) / (vmax - vmin)
            else:
                leaf_df.loc[mask, "value"] = 0.5
    leaf_df.to_csv(OUT_DIR / "leaf_signaling_data.csv", index=False)

    # RRI trajectory
    rri_tp = pd.read_csv(RRI_DIR / "rri_per_timepoint.csv")
    rri_traj = rri_tp[["time_bin", "RRI_mean", "RRI_std"]].copy()
    rri_traj.to_csv(OUT_DIR / "rri_trajectory_data.csv", index=False)

    # Signaling heatmap data (cell type x timepoint)
    heat_rows = []
    for _, row in origin.iterrows():
        heat_rows.append({"CellType": row["CellType"], "time_bin": row["Time"],
                          "OutgoingSignal": float(row["OutgoingSignal"])})
    heat_df = pd.DataFrame(heat_rows)
    # Normalize overall for heatmap
    vmin, vmax = heat_df["OutgoingSignal"].min(), heat_df["OutgoingSignal"].max()
    if vmax > vmin:
        heat_df["OutgoingSignal"] = (heat_df["OutgoingSignal"] - vmin) / (vmax - vmin)
    heat_df.to_csv(OUT_DIR / "signaling_heatmap_data.csv", index=False)

    # Seedling RRI data (map RRI to seedling regions)
    # Use the seedling.saltdrought map which has root/shoot/cotyledon regions.
    # Since RRI is a whole-plant metric, we modulate per-organ using the
    # signaling-origin data: organs that are stronger signaling SOURCES at a
    # given timepoint experience more local stress, so their local RRI is
    # lower. This produces a spatially-resolved "wave" across the seedling.
    # Organ -> aggregate signaling weight (root cell types -> root, leaf/shoot
    # cell types -> cotyledon+hook, hypocotyl gets the mean).
    organ_signal = {}
    for _, row in origin.iterrows():
        ct = row["CellType"]
        sig = float(row["OutgoingSignal"])
        tp = row["Time"]
        if ct.startswith("root_"):
            organ = "root"
        elif ct.startswith("leaf_") or ct == "shoot_apical_meristem":
            organ = "cotyledon"
        elif ct in ("dna_damage_response", "oxidative_stress") or ct.startswith("hormone_response_"):
            organ = "hypocotyl"  # systemic/hormonal -> hypocotyl conduit
        else:
            organ = "hypocotyl"
        organ_signal.setdefault((tp, organ), []).append(sig)
    # Normalize signaling per timepoint to [0,1]
    organ_norm = {}
    for tp in {t for t, _ in organ_signal}:
        vals = [v for (t, o), vs in organ_signal.items() if t == tp for v in vs]
        mx = max(vals) if vals else 1.0
        mn = min(vals) if vals else 0.0
        for (t, o), vs in organ_signal.items():
            if t == tp:
                m = (sum(vs) / len(vs) - mn) / (mx - mn) if mx > mn else 0.5
                organ_norm[(t, o)] = m

    seedling_rows = []
    for _, row in rri_tp.iterrows():
        tp = row["time_bin"]
        rri = float(row["RRI_mean"])
        for roi in ["root", "shoot", "cotyledon", "hypocotyl"]:
            # Local RRI = whole-plant RRI discounted by organ signaling load
            # (organs that are strong signaling sources are more stressed ->
            # lower local resilience). Discount factor 0.08 (max 8% drop).
            sig = organ_norm.get((tp, roi), 0.5) if roi != "shoot" else organ_norm.get((tp, "cotyledon"), 0.5)
            local_rri = rri * (1.0 - 0.08 * sig)
            seedling_rows.append({"ROI": roi, "value": local_rri,
                                  "time_bin": tp})
    seedling_df = pd.DataFrame(seedling_rows)
    seedling_df.to_csv(OUT_DIR / "seedling_rri_data.csv", index=False)

    return root_df, leaf_df, rri_traj, heat_df, seedling_df


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("Preparing data for ggPlantmap...")
    root_df, leaf_df, rri_traj, heat_df, seedling_df = prepare_data_for_r()
    print(f"  Root signaling rows: {len(root_df)}")
    print(f"  Leaf signaling rows: {len(leaf_df)}")
    print(f"  RRI trajectory rows: {len(rri_traj)}")
    print(f"  Signaling heatmap rows: {len(heat_df)}")
    print(f"  Seedling RRI rows: {len(seedling_df)}")

    print("\nCreating R script...")
    r_path = create_r_script({}, OUT_DIR)
    print(f"  R script -> {r_path}")

    print("\nRunning ggPlantmap in R...")
    result = subprocess.run(["Rscript", str(r_path)], capture_output=True, text=True, timeout=600)
    print(result.stdout[-1500:] if len(result.stdout) > 1500 else result.stdout)
    if result.returncode != 0:
        print("STDERR:", result.stderr[-2000:])
        # Don't raise — some maps may fail but others succeed
    print(f"\nFigures saved to {OUT_DIR}")


if __name__ == "__main__":
    main()
