> **Historical record, not current results.** This is the plan the Biomni agent wrote on
> 23 Jul 2026 for the ligand-receptor specificity extension, which became
> [`code/12_lr_pair_specificity.py`](../code/12_lr_pair_specificity.py). It is kept unchanged below
> for provenance. Its sample counts and runtimes are as stated at planning time; for the results, see
> [`README_pipeline.md`](../README_pipeline.md) and the manuscript.

# Plan: LR Pair Specificity Across Radiation Qualities + Visualization Strengthening

## Summary

Extend the CellChat signaling analysis to compute **per-radiation-quality, per-timepoint signaling networks** for all 5 qualities in the cohort (gamma, GCR, spaceflight-LEO, HZE-Fe, UV-B), then test which ligand-receptor pairs show quality-specific signaling patterns using Kruskal-Wallis differential signaling tests. Generate 7 new visualizations that strengthen the kinetic narrative. Integrate findings into both existing Section 4 and a new Section 12 in the manuscript, and update the README.

**Key discovery from data inspection:** HZE-Fe (OSD-46/OSD-320) and UV-B (OSD-296) microarray data is available on disk at `/mnt/shared-workspace/shared/raw/counts/`, already log-normalized (range 2.5–14.3), with sample IDs matching metadata exactly, 21,267 TAIR genes (95.7% overlap with expression matrix), and 2,672 LR pairs (73.1% of full database). This means we can run **empirical** CellChat on all 5 qualities — no GP-AE decoder interpolation needed. The "hybrid" approach simplifies to: empirical CellChat for all 5 qualities, with microarray qualities (HZE-Fe, UV-B) using a slightly reduced LR pair set (2,672 vs 3,605).

## Data Constraints (verified)

| Quality | Assay | Samples | Time bins available | LR pairs |
|---------|-------|---------|-------------------|----------|
| gamma | RNA-seq | 108 (68 with time>0) | 0-0.5h(16), 0.5-2h(14), 2-6h(8), 6-12h(8), 12-30h(10) | 3,605 |
| GCR | RNA-seq | 14 (8 with time) | 2-6h(8) only | 3,605 |
| spaceflight-LEO | RNA-seq | 53 (no time) | static only | 3,605 |
| HZE-Fe | microarray | 36 (OSD-46+OSD-320) | 0.5-2h(8), 2-6h(12), 6-12h(8), 12-30h(8) | 2,672 |
| UV-B | microarray | 36 (OSD-296) | 0.5-2h(9), 2-6h(9), 30h+(18) | 2,672 |

**Cross-quality comparison timepoints:**
- 2-6h: gamma, GCR, HZE-Fe, UV-B (4-way) — richest comparison
- 0.5-2h: gamma, HZE-Fe, UV-B (3-way)
- 6-12h: gamma, HZE-Fe (2-way)
- 12-30h: gamma, HZE-Fe (2-way)

## Implementation Steps

### Step 1: Write `12_lr_pair_specificity.py`

A single self-contained Python script that:

**1a. Microarray deconvolution (HZE-Fe + UV-B)**
- Load OSD-46, OSD-320, OSD-296 microarray probeset expression from `/mnt/shared-workspace/shared/raw/counts/`
- Extract TAIR gene column + sample columns; data is already log-normalized
- Run NNLS deconvolution using the same marker gene sets as script 04 (CELL_TYPE_MARKERS + STRESS_MARKERS), with fallback markers for the 3 missing genes (AT1G03220→SCR alt, AT1G73890→WOL alt, AT1G29910→CAB1 alt)
- Output: `results/deconvolution/proportions_microarray.csv` and `results/deconvolution/pseudo_celltype_expr_microarray.csv`
- Merge with metadata to get quality + time bin per sample

**1b. Per-quality, per-timepoint CellChat**
- For each (quality, time_bin) combination with ≥2 non-control samples:
  - Compute mean pseudo-cell-type expression (proportions-weighted scaling of cell-type profiles, same as script 06)
  - Run mass-action scoring (Hill-type, Kh=0.5) with available LR pairs (3,605 for RNA-seq, 2,672 for microarray)
  - 200-permutation p-values
  - Track per-LR-pair signaling strength (not just aggregated 24×24 matrix)
- Also compute static profiles for spaceflight-LEO (no time)
- Output: `results/cellchat/signaling_flow_per_quality_timepoint.csv` (long format: Quality, Time, Source, Target, SignalStrength, nLRpairs, pvalue)
- Output: `results/cellchat/lr_pair_strength_per_quality.csv` (per-LR-pair total signaling strength per quality-timepoint)

**1c. LR pair specificity analysis (differential signaling)**
- At each shared timepoint (2-6h primary, 0.5-2h secondary):
  - For each LR pair, collect signaling strength across qualities
  - Kruskal-Wallis test (non-parametric, handles unequal sample sizes)
  - BH-FDR correction across all LR pairs
  - Classify pairs as:
    - **Quality-specific** (FDR < 0.05, max/min ratio > 2)
    - **Quality-shared** (FDR > 0.1 or max/min ratio < 1.5)
    - **Intermediate** (everything else)
  - For quality-specific pairs, identify which quality drives the signal (post-hoc pairwise Mann-Whitney U)
- Output: `results/cellchat/lr_pair_specificity.csv` (LR_pair, timepoint, KW_stat, pvalue, padj, specificity_class, enriched_quality, max_quality, min_quality, max_min_ratio)

### Step 2: Run the script

Execute `12_lr_pair_specificity.py` on worker-0. Expected runtime ~2-3 hours (14 quality-timepoint CellChat runs + specificity tests). Run in background if >60 min.

### Step 3: Generate 7 visualizations

All saved as SVG + PNG to `results/figures/`:

1. **LR pair specificity heatmap** (`lr_pair_specificity_heatmap.svg/png`)
   - Rows: top 50 quality-specific LR pairs (ranked by KW statistic at 2-6h)
   - Columns: 5 radiation qualities (at 2-6h) + 3 qualities (at 0.5-2h)
   - Color: signaling strength (z-scored per row)
   - Annotated with enriched quality labels

2. **Pathway kinetic line plot** (`pathway_kinetic_lines.svg/png`)
   - X: time (0-30h, continuous)
   - Y: pathway activation z-score
   - Lines: 8 pathways (DNA repair, oxidative stress, SA, JA, ETH, ABA, auxin, cell cycle/DDR)
   - Vertical dashed line at 2-6h nadir
   - Uses existing `pathway_by_timepoint.csv` data

3. **RRI dose-response curve** (`rri_dose_response.svg/png`)
   - X: dose (Gy, log scale)
   - Y: RRI
   - Points colored by radiation quality
   - Separate panels for each timepoint bin with ≥2 dose levels
   - Uses existing `rri_per_sample.csv` data

4. **Chord diagram for cell-type signaling** (`signaling_chord_diagram.svg/png`)
   - Circular layout: 24 cell types as segments
   - Arcs: top 30 signaling flows at 2-6h (crisis nadir)
   - Arc width ∝ signal strength, color by source cell type
   - Generated for the 2-6h timepoint (most disrupted)

5. **Module dynamics streamgraph** (`module_streamgraph.svg/png`)
   - X: time bins (0-0.5h through 12-30h)
   - Y: stacked area showing module eigengene magnitude
   - 3 layers: blue (early-response), turquoise (sustained), grey
   - Shows the temporal handoff from early-response to sustained modules

6. **RRI component radar per quality** (`rri_radar_per_quality.svg/png`)
   - 5 radar charts (one per quality) arranged in a panel
   - Axes: RRI_latent, RRI_pathway, RRI_module, and 3 pathway sub-components
   - Shows the "fingerprint" of each quality's resilience profile

7. **2D dose-time RRI heatmap** (`rri_dose_time_heatmap.svg/png`)
   - X: time (log scale, 0.1-24h)
   - Y: dose (log scale, 0.1-100 Gy)
   - Color: RRI (from GP-AE surface data, already computed)
   - One panel per quality (gamma, GCR, HZE-Fe)
   - Overlaid measured data points as black dots
   - Uses existing `rri_surface_{gamma,gcr,hze_fe}.csv` data

### Step 4: Update manuscript

**Extend Section 4** (Tissue-Resolved Signaling Kinetics):
- Add subsection 4.3 "Radiation-Quality-Specific Signaling"
- Summarize key quality-specific LR pairs found
- Reference the new specificity heatmap and chord diagram
- Note the microarray deconvolution extension (HZE-Fe, UV-B now included)

**Add Section 12** "Ligand-Receptor Pair Specificity Across Radiation Qualities":
- 12.1: Methodology (per-quality CellChat, Kruskal-Wallis differential signaling, BH-FDR)
- 12.2: Quality-specific LR pairs (table of top pairs, enriched quality, function)
- 12.3: Cross-quality comparison at 2-6h (4-way: gamma vs GCR vs HZE-Fe vs UV-B)
- 12.4: Shared vs specific signaling architecture
- 12.5: New visualizations (7 figures with interpretation)

### Step 5: Update README

- Add `12_lr_pair_specificity.py` to code listing
- Add new data files: `proportions_microarray.csv`, `pseudo_celltype_expr_microarray.csv`, `signaling_flow_per_quality_timepoint.csv`, `lr_pair_strength_per_quality.csv`, `lr_pair_specificity.csv`
- Add 7 new figures to figures listing
- Add microarray deconvolution caveat to Data Dictionary (PseudoCellTypeProfile section)

## Compute Estimate

- **Total runtime:** ~2-3 hours
- **Peak memory:** ~4 GB (permutation matrices during CellChat)
- **Disk:** ~200 MB (new CSV outputs + figures)
- **Execution target:** worker-0 (8 CPU, 32 GB max) — sufficient, no HPC needed
- **Strategy:** Run script 12 in background (run_in_background=true) since CellChat permutations across 14 quality-timepoint combos will exceed 60 min. Generate visualizations as a separate foreground step after the background script completes.

## Assumptions

1. Microarray data is directly comparable to RNA-seq for CellChat mass-action scoring after log-normalization (both are log-scale, non-negative). Cross-platform signaling strength comparisons use z-scoring within each platform to normalize scale differences.
2. The 3 missing marker genes in microarray (AT1G03220, AT1G73890, AT1G29910) are replaced with functionally equivalent alternatives from the same cell type.
3. Kruskal-Wallis is appropriate for LR pair specificity testing because signaling strengths are non-normally distributed and sample sizes are unequal across qualities.
4. BH-FDR < 0.05 with max/min ratio > 2 defines "quality-specific" — this is a discovery-level threshold, not a validated cutoff.
5. The 2,672 LR pairs available in microarray (vs 3,605 in RNA-seq) are sufficient for meaningful specificity analysis — the missing 933 pairs are those with genes not on the ATH1 array.
6. spaceflight-LEO is included as a static (no-time) comparison only, not in time-resolved analyses.
