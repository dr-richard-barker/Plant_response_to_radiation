
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressMessages(library(ggPlantmap))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(patchwork))
svg_font <- function() theme(text = element_text(family = "Liberation Sans"))

# Load data
root_data <- read.csv("/mnt/results/zenodo_bundle/results/figures/root_signaling_data.csv")
leaf_data <- read.csv("/mnt/results/zenodo_bundle/results/figures/leaf_signaling_data.csv")
seedling_data <- read.csv("/mnt/results/zenodo_bundle/results/figures/seedling_rri_data.csv")
rri_trajectory <- read.csv("/mnt/results/zenodo_bundle/results/figures/rri_trajectory_data.csv")
signaling_heatmap <- read.csv("/mnt/results/zenodo_bundle/results/figures/signaling_heatmap_data.csv")
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
                                      name = "Signaling\nstrength")
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
ggsave("/mnt/results/zenodo_bundle/results/figures/rri_trajectory.svg", p_rri, width = 7, height = 4, device = svglite::svglite)
ggsave("/mnt/results/zenodo_bundle/results/figures/rri_trajectory.png", p_rpi <- p_rri, width = 7, height = 4, dpi = 150)
cat("RRI trajectory saved
")

# === Signaling flow heatmap (cell type x timepoint) ===
p_heat <- ggplot(signaling_heatmap, aes(x = time_bin, y = CellType, fill = OutgoingSignal)) +
  geom_tile() +
  stress_scale +
  labs(title = "Signaling origin strength over time",
       x = "Time post-exposure", y = "Cell type (source)") +
  theme_bw() + svg_font() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_text(size = 7))
ggsave("/mnt/results/zenodo_bundle/results/figures/signaling_flow_heatmap.svg", p_heat, width = 7, height = 6, device = svglite::svglite)
ggsave("/mnt/results/zenodo_bundle/results/figures/signaling_flow_heatmap.png", p_heat, width = 7, height = 6, dpi = 150)
cat("Signaling heatmap saved
")

# === ggPlantmap per-timepoint figures ===
# Helper: merge map with values (bypassing ggPlantmap.merge tidy-eval scope issues)
# Handles maps WITH a point column (root, leaf) and WITHOUT (seedling)
merge_map_values <- function(map_df, values_df) {
  merged <- merge(map_df, values_df, by.x = "ROI.name", by.y = "ROI", all.x = TRUE)
  if ("point" %in% colnames(merged)) {
    merged <- merged[order(merged$ROI.id, merged$point), ]
  } else {
    merged <- merged[order(merged$ROI.id), ]
  }
  tibble::as_tibble(merged)
}

# Seedling ROI names in the map are capitalized: Root_control, Hypocotyl_control,
# Cotyledon_control, Hook_control. Map our lowercase names to these.
SEEDLING_ROI_MAP <- c(
  root = "Root_control",
  hypocotyl = "Hypocotyl_control",
  cotyledon = "Cotyledon_control",
  shoot = "Hook_control"   # shoot apical region -> Hook
)

timepoints <- unique(root_data$time_bin)
for (tp in timepoints) {
  tp_safe <- gsub("[^a-zA-Z0-9]", "_", tp)
  cat("--- Processing", tp, "---
")

  # Root signaling
  root_tp <- root_data[root_data$time_bin == tp, ]
  if (nrow(root_tp) > 0) {
    tryCatch({
      root_map <- ggPm.At.roottip.longitudinal
      root_merged <- merge_map_values(root_map, root_tp)
      p_root <- ggPlantmap.heatmap(root_merged, value.quant = value) +
        stress_scale +
        labs(title = paste0("Root signaling - ", tp)) +
        theme_bw() + svg_font()
      ggsave(paste0("/mnt/results/zenodo_bundle/results/figures/ggplantmap_root_", tp_safe, ".svg"), p_root,
             width = 5, height = 4, device = svglite::svglite)
      ggsave(paste0("/mnt/results/zenodo_bundle/results/figures/ggplantmap_root_", tp_safe, ".png"), p_root,
             width = 5, height = 4, dpi = 150)
      cat("  root OK
")
    }, error = function(e) cat("  root FAILED:", conditionMessage(e), "
"))
  }

  # Leaf signaling
  leaf_tp <- leaf_data[leaf_data$time_bin == tp, ]
  if (nrow(leaf_tp) > 0) {
    tryCatch({
      leaf_map <- ggPm.At.leaf.crosssection
      leaf_merged <- merge_map_values(leaf_map, leaf_tp)
      p_leaf <- ggPlantmap.heatmap(leaf_merged, value.quant = value) +
        stress_scale +
        labs(title = paste0("Leaf signaling - ", tp)) +
        theme_bw() + svg_font()
      ggsave(paste0("/mnt/results/zenodo_bundle/results/figures/ggplantmap_leaf_", tp_safe, ".svg"), p_leaf,
             width = 5, height = 4, device = svglite::svglite)
      ggsave(paste0("/mnt/results/zenodo_bundle/results/figures/ggplantmap_leaf_", tp_safe, ".png"), p_leaf,
             width = 5, height = 4, dpi = 150)
      cat("  leaf OK
")
    }, error = function(e) cat("  leaf FAILED:", conditionMessage(e), "
"))
  }

  # Seedling RRI (whole-plant RRI overlaid on seedling map)
  seed_tp <- seedling_data[seedling_data$time_bin == tp, ]
  if (nrow(seed_tp) > 0) {
    tryCatch({
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
      ggsave(paste0("/mnt/results/zenodo_bundle/results/figures/ggplantmap_seedling_", tp_safe, ".svg"), p_seed,
             width = 5, height = 4, device = svglite::svglite)
      ggsave(paste0("/mnt/results/zenodo_bundle/results/figures/ggplantmap_seedling_", tp_safe, ".png"), p_seed,
             width = 5, height = 4, dpi = 150)
      cat("  seedling OK
")
    }, error = function(e) cat("  seedling FAILED:", conditionMessage(e), "
"))
  }

  cat("Saved figures for", tp, "
")
}

cat("All ggPlantmap figures saved
")
