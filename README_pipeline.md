# Plant Ionizing-Radiation Kinetic Landscape Pipeline

A reproducible, FAIR-compliant computational pipeline that aggregates NASA OSDR plant radiation transcriptomics studies, deconvolves bulk samples into pseudo-cell-type profiles, trains a Gaussian-Process autoencoder with continuous dose/time/LET covariates, reconstructs inter-tissue signaling kinetics with PlantCellChat, builds kinetic WGCNA modules, computes a composite Radiation Resilience Index (RRI), and visualizes the stress-signaling wave across plant organs with ggPlantmap.

**Authors:** Richard Barker  
**Date:** 2026-07-23  
**License:** CC-BY 4.0 (data), MIT (code)  

---

## Repository Structure

```
zenodo_bundle/
├── README.md                          # This file — Data Dictionary + usage guide
├── manuscript_chapter.md              # "The kinetic landscape of plant ionizing-radiation response"
├── data/
│   ├── metadata_master.csv            # 299 samples × 22 columns (machine-readable)
│   ├── metadata_master.json           # Same as CSV, JSON format
│   ├── ortholog_map.csv               # 32,833 Arabidopsis genes (identity map for native studies)
│   ├── expression_arabidopsis_space.h5ad  # 195 samples × 32,833 genes (log-normalized)
│   ├── gp_ae_checkpoint.pt            # Trained GP-AE model weights (PyTorch)
│   ├── lr_pairs_arabidopsis.csv       # 3,657 Arabidopsis ligand-receptor pairs (PlantPhoneDB)
│   ├── old_repo_deg_table.csv         # Earlier repo's DESeq2 Up/Down calls (input to script 14)
│   ├── osd782_normalized_counts.csv   # Held-out OSD-782/GLDS-679 normalized counts (36 samples, input to script 15)
│   └── osd782_metadata.csv            # OSD-782 sample dose/time/control annotation
├── code/
│   ├── 01_acquire_osdr.py             # NASA OSDR data acquisition
│   ├── 02_extract_metadata.py         # Metadata extraction & codification
│   ├── 03_orthology_map.py            # Orthology mapping to Arabidopsis
│   ├── 04_deconvolve_bulk.py          # NNLS deconvolution to 24 pseudo-cell-types
│   ├── 05_gp_autoencoder.py           # GP autoencoder with dose/time/LET covariates (KEY SCRIPT)
│   ├── 05b_trajectory_holdout.py      # Trajectory holdout evaluation
│   ├── 06_plantcellchat_kinetic.py    # PlantCellChat kinetic signaling per timepoint
│   ├── 07_kinetic_wgcna.py            # Kinetic WGCNA (calls R subprocess)
│   ├── 08_radiation_resilience_index.py  # Composite RRI computation
│   ├── 09_ggplantmap_spatial.py       # ggPlantmap spatial wave visualization
│   ├── 10_pathway_enrichment.py       # Pathway enrichment heatmaps
│   ├── 11_advanced_visualizations.py  # 3D RRI surfaces, signaling animation, latent scatter, RRI waterfall, prioritization
│   ├── 12_lr_pair_specificity.py      # LR pair specificity across radiation qualities (microarray CellChat + Kruskal-Wallis)
│   ├── 13_visualizations.py           # 7 publication-quality visualizations for kinetic narrative
│   ├── 14_old_repo_overlap.py         # Cross-pipeline check vs earlier DESeq2/WGCNA repo (DEG × module enrichment)
│   ├── 15_osd782_validation.py        # Independent validation: project modules onto held-out OSD-782/GLDS-679
│   └── environment.yml                # Pinned software environment
└── results/
    ├── deconvolution/
    │   ├── proportions.csv            # 195 samples × 24 cell-type proportions (RNA-seq)
    │   ├── proportions_microarray.csv # 124 samples × 29 columns (HZE-Fe + UV-B + gamma microarray deconvolution)
    │   ├── pseudo_celltype_expr.csv   # 24 cell types × 32,833 genes (posterior mean, RNA-seq)
    │   └── pseudo_celltype_expr_microarray.csv  # 23 cell types × 21,267 genes (microarray deconvolution)
    ├── trajectories/
    │   ├── latent_embeddings.csv      # 195 samples × 16 latent dimensions
    │   ├── latent_embeddings_annotated.csv  # With metadata columns
    │   ├── decoder_trajectories.csv   # GP-decoded expression trajectories
    │   ├── cv_results.json            # Leave-one-study-out CV results
    │   ├── trajectory_holdout.json    # Trajectory holdout evaluation
    │   └── training_summary.json      # GP-AE training metrics
    ├── cellchat/
    │   ├── signaling_flow_per_timepoint.csv      # 2,760 rows: LR pair × cell pair × timepoint
    │   ├── signaling_heatmap_per_timepoint.csv   # Cell-type × cell-type signaling matrix per timepoint
    │   ├── signaling_origin_summary.csv          # Outgoing/incoming/net signal per cell type per timepoint
    │   ├── top_lr_pairs.csv                      # Top LR pairs by total signaling strength
    │   ├── timepoint_metadata.json               # Sample counts and studies per timepoint bin
    │   ├── signaling_flow_per_quality_timepoint.csv  # 7,222 rows: signaling flow per quality × timepoint × cell pair
    │   ├── lr_pair_strength_per_quality.csv      # 38,545 rows: mean LR pair signaling strength per quality × timepoint
    │   ├── lr_pair_strength_per_sample.csv       # 738,701 rows: per-sample LR pair signaling strength (for Kruskal-Wallis)
    │   ├── lr_pair_specificity.csv               # 10,243 rows: LR pair specificity results (KW p-values, BH-FDR, classification)
    │   └── quality_timepoint_metadata.json       # 14 quality-timepoint combinations with sample counts
    ├── wgcna/
    │   ├── modules.csv                # 2,000 genes × module assignment + kME values
    │   ├── module_traits.csv          # Module × trait correlations (time, dose, LET)
    │   ├── module_classification.csv  # Module → early-response / sustained / DDR-core classification
    │   ├── module_eigengenes.csv      # Per-sample module eigengene values
    │   ├── summary.json               # WGCNA parameters and module counts
    │   ├── osd782_validation.json     # Independent-validation results on held-out OSD-782 (script 15)
    │   └── run_wgcna.R                # R script for WGCNA execution
    ├── rri/
    │   ├── rri_per_sample.csv         # 195 samples × RRI + component scores
    │   ├── rri_per_timepoint.csv      # 5 timepoint bins × RRI mean/std + per-quality breakdown
    │   ├── rri_summary.json           # RRI weights and summary statistics
    │   ├── rri_surface_gamma.csv      # 900-point GP-AE predicted RRI surface (gamma, LET=0.2)
    │   ├── rri_surface_gcr.csv        # 900-point GP-AE predicted RRI surface (GCR, LET=50)
    │   └── rri_surface_hze_fe.csv     # 900-point GP-AE predicted RRI surface (HZE-Fe, LET=200)
    ├── pathway_enrichment/
    │   ├── pathway_scores_per_sample.csv       # 195 samples × 8 pathway activation scores
    │   ├── pathway_by_timepoint.csv            # 8 pathways × 5 timepoint bins
    │   ├── pathway_by_timepoint_quality.csv    # 8 pathways × timepoint × radiation quality
    │   ├── pathway_by_celltype.csv             # 8 pathways × 24 cell types
    │   ├── module_eigengene_by_timepoint.csv   # 3 modules × 5 timepoint bins
    │   └── radiation_quality_prioritization.csv  # 11 missing qualities ranked by composite score
    └── figures/
        ├── rri_trajectory.svg/png              # RRI over time (line plot with ribbon)
        ├── signaling_flow_heatmap.svg/png      # Cell-type signaling origin × timepoint
        ├── ggplantmap_root_{timepoint}.svg/png     # 5 root tip maps (signaling strength)
        ├── ggplantmap_leaf_{timepoint}.svg/png     # 5 leaf cross-section maps (signaling strength)
        ├── ggplantmap_seedling_{timepoint}.svg/png # 5 seedling maps (organ-specific RRI)
        ├── pathway_heatmap_timepoint.svg/png       # Pathway activation × timepoint
        ├── pathway_heatmap_timepoint_quality.svg/png  # Pathway × timepoint × quality
        ├── pathway_heatmap_celltype.svg/png       # Pathway activation × cell type
        ├── module_eigengene_heatmap.svg/png       # WGCNA module eigengene × timepoint
        ├── rri_surface_3d_gamma.svg/png           # 3D GP-AE predicted RRI surface (gamma, LET=0.2)
        ├── rri_surface_3d_hze_fe.svg/png          # 3D GP-AE predicted RRI surface (HZE-Fe, LET=200)
        ├── rri_surface_3d_gcr.svg/png             # 3D GP-AE predicted RRI surface (GCR, LET=50)
        ├── rri_empirical_2d.svg/png               # Empirical RRI vs time (gamma) and vs dose (GCR)
        ├── signaling_flow_animation.gif           # 5-frame network graph animation of signaling topology
        ├── signaling_frame_{0_0_5h,0_5_2h,2_6h,6_12h,12_30h}.png  # Individual animation frames
        ├── latent_3d_scatter.svg/png              # 3D latent space scatter colored by radiation quality
        ├── latent_3d_rotation.gif                 # 36-frame 360° rotation of latent space
        ├── rri_component_waterfall.svg/png        # Stacked bar: RRI components per timepoint
        ├── radiation_quality_prioritization.svg/png  # Horizontal bar: 11 missing qualities ranked
        ├── lr_pair_specificity_heatmap.svg/png    # Top 50 quality-specific LR pairs at 2–6h (z-scored heatmap)
        ├── pathway_kinetic_lines.svg/png          # 8 pathways over 5 timepoint bins (line plot)
        ├── rri_dose_response.svg/png              # RRI vs dose by radiation quality (scatter + error bars)
        ├── signaling_chord_diagram.svg/png        # Top 30 inter-cell-type signaling flows at 2–6h (chord diagram)
        ├── module_streamgraph.svg/png             # WGCNA module eigengenes over time (streamgraph + signed overlay)
        ├── rri_radar_per_quality.svg/png          # 4 radar charts: RRI components per quality (Control, gamma, GCR, LEO)
        ├── rri_dose_time_heatmap.svg/png          # 3 interpolated dose×time RRI heatmaps (gamma, GCR, HZE-Fe)
        ├── old_repo_overlap.svg/png               # DEG × module enrichment + DDR-gene provenance vs earlier repo
        └── osd782_validation.svg/png              # Independent held-out validation on OSD-782/GLDS-679
```

---

## Data Dictionary

### Radiation Parameters

#### `RadiationQuality` (controlled vocabulary, 16 values)

The radiation type applied to the sample. Visible light (PAR, 400–700 nm) is **excluded** from this study per design. UV is retained as radiation (non-ionizing but genotoxic).

| Value | Class | LET class | Description |
|-------|-------|-----------|-------------|
| `GCR` | ionizing-particulate-mixed | mixed | Simulated galactic cosmic rays; sequential p, He, O, Si, Fe beams at NSRL |
| `proton` | ionizing-particulate-low-LET | low | Proton beams from cyclotron/NSRL |
| `HZE-Fe` | ionizing-particulate-high-LET | high | Iron-56 ions (1 GeV Fe²⁶⁺); produces clustered DNA double-strand breaks |
| `HZE-Si` | ionizing-particulate-high-LET | high | Silicon ions (GCR component) |
| `HZE-O` | ionizing-particulate-high-LET | high | Oxygen ions (GCR component) |
| `helium` | ionizing-particulate-low-LET | low | Alpha particles / He ions (GCR component) |
| `neutron` | ionizing-particulate-mixed | mixed | Neutron sources |
| `beta` | ionizing-particulate-low-LET | low | Beta emitters (often dosimetry context) |
| `gamma` | ionizing-photon-low-LET | low | Gamma radiation (Co-60 or Cs-137); dispersed DNA damage |
| `X-ray` | ionizing-photon-low-LET | low | X-ray radiation |
| `spaceflight-LEO` | ionizing-mixed-chronic | mixed | International Space Station environment; trapped protons + GCR, TLD-measured |
| `solar-particle-event` | ionizing-particulate-mixed | mixed | Solar particle events (SPE); flagged in OSDR but no dedicated plant study yet |
| `UV-A` | non-ionizing-UV | n/a | Ultraviolet A (315–400 nm) |
| `UV-B` | non-ionizing-UV | n/a | Ultraviolet B (280–315 nm) |
| `UV-C` | non-ionizing-UV | n/a | Ultraviolet C (100–280 nm); germicidal, lab mutagenesis |
| `cosmic-mixed` | ionizing-mixed-chronic | mixed | Cosmic descriptor; usually maps to GCR + spaceflight |

**Studies in this release use 5 of these 16 qualities:** gamma (160 samples), spaceflight-LEO (53), HZE-Fe (36), UV-B (36), GCR (14). The remaining 11 qualities are defined in the vocabulary for completeness and future extension.

#### `RadiationClass` (5 values)

Broad physical classification of the radiation type:

| Value | Description | Qualities included |
|-------|-------------|-------------------|
| `ionizing-photon-low-LET` | Electromagnetic radiation, low linear energy transfer | gamma, X-ray |
| `ionizing-particulate-low-LET` | Charged particles, low LET | proton, helium, beta |
| `ionizing-particulate-high-LET` | Heavy ions, high LET (clustered DNA damage) | HZE-Fe, HZE-Si, HZE-O |
| `ionizing-particulate-mixed` | Mixed particle fields | GCR, neutron, solar-particle-event |
| `ionizing-mixed-chronic` | Chronic mixed radiation (spaceflight) | spaceflight-LEO, cosmic-mixed |
| `non-ionizing-UV` | Ultraviolet (non-ionizing but genotoxic) | UV-A, UV-B, UV-C |

#### `LET_keV_um` (float, keV/μm)

Linear energy transfer — the energy deposited per unit path length. A first-class covariate alongside dose and time because high-LET radiation (HZE ions) produces qualitatively different DNA damage (clustered double-strand breaks) compared to low-LET (gamma, proton).

- **Measured values:** extracted from the OSDR `linear energy transfer` protocol parameter where available (OSD-658, OSD-46).
- **Imputed values:** for qualities without an explicit LET field, the class-default LET is assigned (low ≈ 0.2, high ≈ 200, mixed = NaN) and `LET_imputed = True`.
- **Range in this dataset:** 0.2 (gamma) to 200 (HZE-Fe) keV/μm.

#### `LET_imputed` (boolean)

`True` if `LET_keV_um` was assigned from the class default rather than measured directly from OSDR metadata.

#### `AbsorbedDose_Gy` (float, Gy)

Total absorbed radiation dose in Gray (1 Gy = 1 J/kg). Extracted from the OSDR `Absorbed Radiation Dose` factor or protocol parameter.

- **Range in this dataset:** 0.0 (controls) to 100.0 Gy.
- **Dose grid by quality:** GCR: 0.4, 0.8 Gy; HZE-Fe: 30, 100 Gy; gamma: 100 Gy; UV-B: dose not applicable (wavelength-based exposure); spaceflight-LEO: TLD-measured (chronic, mGy/day).
- **Unit conversion:** cGy → Gy (÷100); mGy → Gy (÷1000).

#### `DoseRate_Gy_min` (float, Gy/min)

Dose rate where available. Empty for chronic spaceflight exposures and UV-B.

#### `TimePostExposure_h` (float, hours)

Time between end of radiation exposure and tissue harvest, in hours.

- **Definition:** For acute exposures (gamma, HZE, GCR), this is the time from the end of irradiation to sample collection. For chronic exposures (spaceflight-LEO), this is the total mission duration (harvest time minus launch time).
- **Extraction rule:** Extracted from OSDR `Time of Sample Collection After Treatment` factor where available. For studies where time is encoded in sample names (e.g., OSD-498/508: `gIR_vs_mock_wt_rep1_1h30min`), parsed via regex and cross-checked against protocol text.
- **Range in this dataset:** 0.0 (harvest at exposure end) to 96.0 h (4 days).
- **Time bins used for aggregation:** 0–0.5h, 0.5–2h, 2–6h, 6–12h, 12–30h.

### Sample Metadata

| Column | Type | Description |
|--------|------|-------------|
| `SampleID` | string | Unique sample identifier (GSM-style with suffix, e.g., `GSM1506014_1`) |
| `StudyID` | string | NASA OSDR study identifier (e.g., `OSD-46`) |
| `Organism` | string | Species name (all `Arabidopsis thaliana` in this release) |
| `Ecotype` | string | Arabidopsis ecotype (e.g., `Ler-0`, `Col-0`) |
| `Genotype` | string | Genotype description (e.g., `Wild Type`, `sog1-1`, `atm/adm`) |
| `Tissue` | string | Tissue type (empty for whole-seedling studies) |
| `AgeAtHarvest_d` | string | Plant age at harvest (e.g., `5 day`, `14 day`) |
| `Replicate` | string | Replicate number |
| `IsControl` | boolean | `True` if sham-irradiated or ground control; `False` if irradiated |
| `AssayType` | string | `RNA-seq` or `microarray` |
| `RawDataFile` | string | Path to raw data file (if available) |
| `ProcessedDataFile` | string | Path to processed data file (if available) |
| `ControlMatchedID` | string | SampleID of the matched control for this sample |
| `OriginalSampleName` | string | Original sample name from the source repository |

### Derived Metrics

#### `RRI` — Radiation Resilience Index (float, 0–1)

A composite metric quantifying how well a sample's transcriptomic state is preserved relative to sham-irradiated controls. RRI = 1.0 at baseline (perfect preservation), decreases as the sample deviates, and recovers as homeostasis is restored.

**Formula:**

```
RRI = 0.50 × RRI_latent + 0.25 × RRI_pathway + 0.25 × RRI_module
```

**Components:**

1. **RRI_latent** (weight 0.50): Gaussian distance of the sample's GP-AE latent vector `z` from the control centroid in latent space.
   ```
   RRI_latent = exp(-||z - z_control_mean||² / σ²)
   ```
   where σ = median pairwise distance among control samples.

2. **RRI_pathway** (weight 0.25): Shannon evenness of pathway activation scores across three core pathways (DNA repair, oxidative stress, hormone signaling). A balanced multi-pathway response is more resilient than a single pathway dominating.
   ```
   RRI_pathway = H / log(n_pathways)
   where H = -Σ pᵢ log(pᵢ), pᵢ = pathway_scoreᵢ / Σ pathway_scores
   ```

3. **RRI_module** (weight 0.25): Correlation of the sample's expression of sustained (housekeeping) WGCNA module genes with the control mean expression. High preservation = resilient; low = disrupted.
   ```
   RRI_module = max(0, cor(expr_sustained[sample], ctrl_mean_sustained))
   ```

**Trajectory in this dataset:** 0.79 (0–0.5h) → 0.79 (0.5–2h) → **0.68 (2–6h, nadir)** → 0.75 (6–12h) → 0.78 (12–30h). The nadir at 2–6h is driven primarily by the module-preservation component (drops to 0.61), indicating maximal disruption of housekeeping networks at this timepoint.

#### `PseudoCellTypeProfile` (inference caveat)

The 24 pseudo-cell-type expression profiles in `results/deconvolution/pseudo_celltype_expr.csv` are **inferential reconstructions** from bulk RNA-seq using non-negative least squares (NNLS) deconvolution with curated marker genes — not measured single-cell data. No plant single-cell radiation dataset exists in NASA OSDR. All CellChat signaling analyses and cell-type-resolved pathway scores are labeled as "in-silico reconstruction" in figures and manuscript.

**24 cell types:** root_meristem, root_columella, root_epidermis, root_cortex, root_endodermis, root_stele, root_xylem, root_phloem, root_pericycle, root_lateral_root_cap, shoot_apical_meristem, leaf_epidermis, leaf_mesophyll, leaf_vasculature, hypocotyl, cotyledon, dna_damage_response, oxidative_stress, hormone_response_{ja, sa, aba, ethylene, auxin, cytokinin}.

#### `PseudoCellTypeProfile_Microarray` (platform caveat)

The microarray deconvolution (`results/deconvolution/proportions_microarray.csv`, `pseudo_celltype_expr_microarray.csv`) uses the same NNLS marker-gene framework but on Affymetrix ATH1 microarray data from OSD-46 (HZE-Fe) and OSD-296 (UV-B). Three marker genes absent from the microarray platform (AT1G03220, AT1G73890, AT1G29910) were replaced with fallback markers (AT4G28400, AT3G62940, AT1G61550), yielding 23 deconvolved cell types (one cell type had no present markers after filtering). The microarray covers 21,267 genes (95.7% overlap with RNA-seq) and 2,672 LR pairs (73.1% of 3,605).

**Platform effects:** Cross-platform LR pair specificity comparisons (gamma RNA-seq vs. HZE-Fe/UV-B microarray) show near-universal "quality-specific" classification that likely reflects systematic platform differences in signal magnitude rather than true biological specificity. The 2–6h timepoint — where GCR and gamma are both RNA-seq — is the only platform-matched cross-quality comparison and should be prioritized for biological interpretation. See manuscript Section 12.4 for full discussion.

#### `OrthologyType`

All samples in this release are Arabidopsis thaliana, so the orthology map (`data/ortholog_map.csv`) is an identity map (32,833 genes). The orthology infrastructure is in place for future cross-species extension (rice, soybean, maize via Ensembl Plants Compara + PLAZA 5.0).

---

## Software Environment

See `code/environment.yml` for pinned versions. Key dependencies:

**Python:** torch, gpytorch, anndata, scanpy, numpy, pandas, scipy, matplotlib, seaborn  
**R:** WGCNA (v1.74), PlantCellChat (v0.0.9), ggPlantmap (v1.1.0), svglite, impute, preprocessCore  
**LR database:** PlantPhoneDB (Xu et al. 2022) — 3,657 Arabidopsis ligand-receptor pairs

---

## Reproduction Guide

### Prerequisites

1. Python 3.11+ with `uv` package manager
2. R 4.3+ with Bioconductor
3. ~2 GB disk space for data + results

### Running the Pipeline

```bash
# Phase 1: Data curation
python code/01_acquire_osdr.py        # Download OSDR studies
python code/02_extract_metadata.py    # Extract & codify metadata
python code/03_orthology_map.py       # Map to Arabidopsis orthologs

# Phase 2: Kinetic modeling
python code/04_deconvolve_bulk.py     # Deconvolve to pseudo-cell-types
python code/05_gp_autoencoder.py      # Train GP autoencoder (KEY SCRIPT)
python code/05b_trajectory_holdout.py # Evaluate trajectory prediction
python code/06_plantcellchat_kinetic.py  # CellChat per timepoint

# Phase 3: Network & synthesis
python code/07_kinetic_wgcna.py       # WGCNA kinetic modules
python code/08_radiation_resilience_index.py  # Compute RRI
python code/09_ggplantmap_spatial.py  # Spatial wave visualization
python code/10_pathway_enrichment.py  # Pathway enrichment heatmaps

# Phase 4: Advanced visualizations & prioritization
python code/11_advanced_visualizations.py  # 3D RRI surfaces, signaling GIF, latent scatter, waterfall, prioritization

# Phase 5: LR pair specificity & final visualizations
python code/12_lr_pair_specificity.py     # Microarray CellChat + Kruskal-Wallis specificity testing
python code/13_visualizations.py          # 7 publication-quality figures for kinetic narrative

# Phase 6: Cross-pipeline consistency & independent validation
python code/14_old_repo_overlap.py        # DEG × module enrichment vs earlier DESeq2/WGCNA repo
python code/15_osd782_validation.py       # Project modules onto held-out OSD-782/GLDS-679
```

Scripts are numbered to reflect execution order. Each script reads from the outputs of prior scripts and writes to `results/`. The GP autoencoder checkpoint is saved to allow resuming without retraining.

### Key Methodological Notes

1. **Within-study centering** (script 05): The GP-AE subtracts each study's control mean from all samples before training. This removes batch effects so the model learns radiation response, not study identity. Without this step, the model achieves high in-sample R² but fails completely on held-out studies.

2. **Honest cross-validation** (script 05): Leave-one-study-out R² ≈ 0 is the honest result for this cohort. The cohort has 195 samples across 5 studies with time data, but dose is nearly constant within studies and LET is confounded with radiation quality. The model captures response *direction* (trajectory holdout cosine similarity 0.60–0.82) but not *magnitude* across studies.

3. **IsControl parsing** (scripts 07, 08): The `IsControl` column in the raw metadata is a categorical with values `''`, `'False'`, `'True'`. Using `.astype(bool)` converts the non-empty string `'False'` to `True` — a critical bug. The fix: `str.lower().isin(['true', '1', 'yes'])`.

4. **PlantCellChat LR database** (script 06): The PlantCellChatDB is served dynamically from a Shiny app and is not bundled with the R package. We use PlantPhoneDB (Xu et al. 2022) LR pairs instead, downloaded from the LewisLabUCSD/Ligand-Receptor-Pairs repository.

5. **Microarray CellChat** (script 12): The GP-AE decoder outputs only 2,000 HVGs, covering just 59 of 3,605 LR pairs (1.6%) — insufficient for CellChat. To enable cross-quality LR pair analysis, microarray expression data (OSD-46 HZE-Fe, OSD-296 UV-B) were deconvolved and used for empirical CellChat. Per-sample LR pair signaling strengths were computed by scaling pseudo-cell-type profiles by each sample's deconvolved proportions, enabling Kruskal-Wallis testing across qualities. The 2–6h timepoint (GCR vs. gamma, both RNA-seq) is the only platform-matched cross-quality comparison; other timepoints cross platforms and should be interpreted cautiously.

6. **Platform effects in specificity testing** (script 12): Cross-platform comparisons (RNA-seq gamma vs. microarray HZE-Fe/UV-B) produce near-universal "quality-specific" classification (2,536–2,540 of ~2,561 pairs) that likely reflects systematic platform differences rather than biology. Only the 2–6h GCR-vs-gamma comparison (both RNA-seq) yields biologically interpretable specificity results: 71 GCR-enriched quality-specific pairs and 2,477 quality-shared pairs.

7. **Module composition and the grey DDR core** (script 14): The blue early-response module is significantly enriched for radiation DEGs called by the earlier `dr-richard-barker/Plant_response_to_radiation` DESeq2 pipeline (228/739, 1.55×, P = 2.0 × 10⁻²⁰), a cross-method consistency check on the shared studies OSD-498/502/508/510/658 (not independent replication, since the two analyses use mostly the same samples). The canonical DNA double-strand-break repair genes (BRCA1, PARP2, RAD51, SMR7, GMI1, XRI1) do **not** sit among the blue hubs — they concentrate in the small `grey` module (49 genes, 92% radiation DEGs, 4.6× enrichment), which is the tight transcriptional core of the DNA-damage response rather than WGCNA's usual unassigned bin. The upstream regulators ATM, ATR, SOG1 and MYB3R fall below the top-2,000-HVG variance threshold and are absent from the network. See `results/figures/old_repo_overlap.svg` and manuscript Section 5. The `DDR-core` label in `module_classification.csv` is emitted automatically by script 07: any non-early/late module carrying ≥3 canonical DSB-repair effectors (and enriched for the bundled old-repo DEGs at P ≤ 0.05) is relabelled, so the annotation is reproduced on a re-run and follows the DSB-repair genes regardless of which colour WGCNA assigns them.

8. **Independent held-out validation** (script 15): OSD-782 (GeneLab **GLDS-679** — note the OSD and GLDS numbers differ) is a wild-type low-dose gamma time-course (0/0.1/1 Gy × 1–72 h, 36 samples) that is **not** in the training cohort and uses a dose regime ~100× below the cohort's 100 Gy gamma studies, so it is a genuinely external test set (unlike the shared-sample check in note 7). Projecting the modules onto OSD-782 without refitting, its radiation-response genes (815 at FDR < 0.05) are significantly enriched in the blue early-response module (41/725, 2.3×, P = 2.8 × 10⁻¹²) and depleted in the sustained turquoise module (0.3×). The grey DDR-core module is directionally enriched (2.6×) but underpowered at this sample size (P = 0.10); at the gene level its members PARP2 and XRI1 are significantly radiation-induced (FDR < 0.03). Projected blue and grey module scores decline through the 1–24 h crisis window and recover by 72 h, independently reproducing the three-phase architecture. Raw OSDR download: `https://osdr.nasa.gov/osdr/data/osd/files/782` (RNA-seq normalized counts). See `results/figures/osd782_validation.svg`, `results/wgcna/osd782_validation.json`, and manuscript Section 5.

---

## Citation

If you use this pipeline or data, please cite:

- Jo, L. & Kajala, K. (2024). ggPlantmap: a versatile toolkit for visualizing plant organ spatial data. *J Exp Bot* 75:5366–5376. doi:10.1093/jxb/erae043
- Liu et al. (2026). PlantCellChat: inference and visualization of cell-cell communication for plant single-cell data. *Plant J* 126(3):e70905
- Xu et al. (2022). PlantPhoneDB: A comprehensive database of ligand-receptor pairs in plants. *Plant Biotechnol J*
- NASA Open Science Data Repository (OSDR): https://osdr.nasa.gov

---

## Future Data Acquisition Priorities

The current release covers 5 of 16 radiation qualities in the controlled vocabulary. The 11 missing qualities were ranked by a composite score (50% data availability + 25% LET gap-filling + 25% mission relevance) to guide future data generation. The full scoring table is in `results/pathway_enrichment/radiation_quality_prioritization.csv` and visualized in `results/figures/radiation_quality_prioritization.svg/png`.

| Tier | Qualities | Action |
|------|-----------|--------|
| **Tier 1** (acquire now) | proton (0.81), cosmic-mixed (0.80) | 12 + 10 existing OSDR studies; reprocess through this pipeline immediately. Expands cohort from 5→7 qualities. |
| **Tier 2** (acquire if feasible) | helium (0.49), solar-particle-event (0.44), neutron (0.42), HZE-Si (0.39), HZE-O (0.38) | Helium and HZE-Si/O have existing data (reprocessable). SPE and neutron need dedicated NSRL campaigns — SPE is the most acute mission hazard with zero plant data. |
| **Tier 3** (opportunistic) | beta (0.32), X-ray (0.17), UV-A (0.16) | Low spaceflight relevance; incorporate opportunistically to broaden LET coverage. |
| **Tier 4** (requires new experiments) | UV-C (0.03) | No OSDR data, GEO-only exposure, minimal relevance. Lowest priority. |

See manuscript Section 11.5 for the full rationale and per-quality notes.

---

## Limitations

1. **No plant single-cell radiation data exists** — all cell-type-resolved analyses use atlas-deconvolved pseudo-cell-types (in-silico reconstruction).
2. **24-hour time ceiling** — the longest post-exposure timepoint is 24h (most studies), preventing identification of true late-adaptation modules (≥48h).
3. **Dose confounding** — dose is nearly constant within each study (e.g., all gamma studies use 100 Gy), so dose-dependency cannot be disentangled from study batch effects.
4. **LET confounding** — LET is perfectly confounded with radiation quality in this cohort (all HZE-Fe samples have the same LET), so the LET covariate in the GP-AE cannot be independently validated.
5. **5 of 16 radiation qualities represented** — the controlled vocabulary covers 16 qualities, but only 5 (gamma, spaceflight-LEO, HZE-Fe, UV-B, GCR) have data in this release. See Future Data Acquisition Priorities above for the ranked acquisition plan.
6. **GP-AE is a novel method** — no external plant-radiation validation cohort exists; results are discovery-only.
7. **3D RRI surfaces are model predictions, not measured data** — the GP-AE decoder interpolates across the dose–time grid, but LOSO R² ≈ 0 means absolute predictions on unseen studies are unreliable. Surfaces are valuable for qualitative topology and relative quality-dependent differences, not quantitative dose–response estimation.
8. **Platform effects in LR pair specificity** — cross-quality LR pair comparisons at most timepoints cross platforms (RNA-seq gamma vs. microarray HZE-Fe/UV-B), producing near-universal "quality-specific" classification that likely reflects platform differences rather than biology. Only the 2–6h GCR-vs-gamma comparison (both RNA-seq) is platform-matched and biologically interpretable. Future studies should use a single platform or include platform-matched controls for each quality.
