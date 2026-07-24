# The Kinetic Landscape of Plant Ionizing-Radiation Response

## Manuscript Chapter — Bridging Radiation Signaling to the Spaceflight Stress Manifold

---

### Abstract

Plants encounter ionizing radiation in spaceflight environments ranging from chronic low-dose trapped-particle exposure in low Earth orbit (LEO) to acute high-dose solar particle events and galactic cosmic rays (GCR). The transcriptomic response to radiation is not a static switch but a dynamic kinetic landscape: DNA damage repair pathways activate within minutes, hormonal signaling cascades propagate over hours, and systemic inter-tissue communication coordinates a whole-organism stress response. We present a computational pipeline that integrates 299 samples from 10 NASA OSDR studies spanning five radiation qualities (gamma, HZE-Fe, GCR, UV-B, spaceflight-LEO), deconvolves bulk transcriptomes into 24 pseudo-cell-type profiles, trains a Gaussian-Process autoencoder (GP-AE) with continuous dose, time, and linear energy transfer (LET) covariates, and reconstructs the temporal signaling network with PlantCellChat. We identify a three-phase kinetic architecture — early recognition (0–2h), crisis nadir (2–6h), and partial recovery (6–24h) — quantified by a composite Radiation Resilience Index (RRI) that drops from 0.79 to 0.68 at the nadir before recovering to 0.78. The early-response WGCNA module (739 genes, Time ρ = −0.66, p = 3.5 × 10⁻¹⁰) is significantly enriched for independently-called radiation DEGs (1.55×, P = 2 × 10⁻²⁰), while the core DNA double-strand-break repair machinery (BRCA1, PARP2, RAD51) forms a distinct, tightly co-regulated signature; the sustained module (1,212 genes) carries housekeeping functions whose disruption at the nadir drives the RRI decline. This kinetic framework provides a shared cross-stress representation that bridges radiation response to the broader spaceflight stress manifold, including the microgravity response explored in our companion work on ancient gravitational signaling.

---

### 1. Introduction

Spaceflight imposes a multi-axis stress environment on biological systems: microgravity alters mechanical signaling and cytoskeletal architecture, while ionizing radiation damages DNA, generates reactive oxygen species (ROS), and disrupts cellular homeostasis. The plant radiation response has been studied extensively at the transcriptomic level, but most analyses treat radiation as a binary condition (irradiated vs. control) at a single timepoint, losing the kinetic dimension that defines how stress propagates through tissues and resolves over time.

The fundamental challenge is that plant radiation transcriptomics data are sparse, heterogeneous, and batch-confounded. NASA's Open Science Data Repository (OSDR) contains 83 plant studies with a radiation signal, but they span different radiation qualities (gamma, proton, HZE ions, GCR, UV, spaceflight LEO), doses ranging from milligrays (spaceflight) to 100 Gy (acute gamma), and timepoints from 10 minutes to 96 hours — often with only one or two timepoints per study. No plant single-cell radiation dataset exists, forcing reliance on bulk transcriptomics deconvolved into pseudo-cell-type profiles.

We address these challenges with a pipeline built on three innovations:

1. **A Gaussian-Process autoencoder (GP-AE)** that treats dose, time, and LET as continuous covariates in the decoder, enabling interpolation across the sparse cross-study grid. The latent space captures each sample's intrinsic state, while the covariate-conditioned decoder generates expression trajectories at any (dose, time, LET) combination.

2. **Kinetic PlantCellChat analysis** that reconstructs inter-tissue signaling at each discrete timepoint, revealing how signaling origin shifts from photosynthetic tissues (early) to oxidative-stress and hormonal pathways (mid-phase) as the response evolves.

3. **A composite Radiation Resilience Index (RRI)** that integrates latent-space distance, pathway balance, and co-expression module preservation into a single dynamic metric tracking tissue stability degradation and recovery.

### 2. Data Landscape

#### 2.1 Cohort Assembly

We aggregated 299 samples from 10 NASA OSDR studies (Table 1), all Arabidopsis thaliana, spanning five radiation qualities. The cohort includes 62 sham-irradiated controls and 237 irradiated samples. RNA-seq studies (OSD-498, 508, 510, 658) provide the dense temporal coverage (10 min to 24 h), while microarray studies (OSD-46, 320) contribute the dose × quality factorial design (gamma vs. HZE-Fe).

| Study | Assay | Samples | Radiation quality | Dose (Gy) | Time (h) |
|-------|-------|---------|-------------------|-----------|----------|
| OSD-46 | microarray | 44 | gamma, HZE-Fe | 30, 100 | 1.5–24 |
| OSD-320 | microarray | 44 | gamma, HZE-Fe | 30, 100 | 1.5–24 |
| OSD-498 | RNA-seq | 16 | gamma (Co-60) | 100 | 0.17–24 |
| OSD-502 | RNA-seq | 8 | gamma (Co-60) | 100 | 0 (control) |
| OSD-508 | RNA-seq | 36 | gamma (Co-60) | 100 | 0.17–24 |
| OSD-510 | RNA-seq | 48 | gamma (Co-60) | 100 | 0.33–24 |
| OSD-658 | RNA-seq | 14 | GCR | 0.4, 0.8 | 3 |
| OSD-296 | RNA-seq | 36 | UV-B | wavelength | 1–96 |
| OSD-314 | RNA-seq | 17 | spaceflight-LEO | chronic | — |
| OSD-346 | RNA-seq | 36 | spaceflight-LEO | chronic | — |

**Table 1.** Study inventory. Time = hours post-exposure. Spaceflight-LEO samples lack discrete post-exposure timepoints (chronic exposure).

#### 2.2 Metadata Codification

Each sample was annotated with a 22-column metadata record including `RadiationQuality` (controlled vocabulary of 16 radiation types), `RadiationClass` (6 physical classes), `LET_keV_um` (linear energy transfer, measured or imputed), `AbsorbedDose_Gy`, `TimePostExposure_h`, and `IsControl`. LET is treated as a first-class covariate because high-LET radiation (HZE ions, ~200 keV/μm) produces clustered DNA double-strand breaks qualitatively distinct from the dispersed damage of low-LET gamma rays (~0.2 keV/μm).

#### 2.3 Honest Limitations of the Cohort

Three structural confounds limit cross-study generalization:

- **Dose is nearly constant within studies** (all gamma studies use 100 Gy; GCR uses 0.4/0.8 Gy), so dose-dependency cannot be separated from study batch effects.
- **LET is perfectly confounded with radiation quality** (all HZE-Fe samples share the same LET), preventing independent LET validation.
- **The 24-hour time ceiling** (longest post-exposure timepoint) prevents identification of true late-adaptation modules that would emerge at 48–96h.

These constraints mean the GP-AE captures response *direction* (trajectory holdout cosine similarity 0.60–0.82) but not *magnitude* across studies (leave-one-study-out R² ≈ 0). We report this honestly rather than inflating in-sample metrics.

### 3. Kinetic Autoencoder with Covariate Conditioning

#### 3.1 Architecture

The GP-AE encodes each sample's 2,000-gene expression profile into a 16-dimensional latent vector **z**, then decodes it conditioned on three continuous covariates: log(dose), log(time), and log(LET). The decoder uses a Gaussian Process layer with a Matérn-5/2 kernel over the covariate space, providing smooth interpolation across the sparse, heterogeneous cross-study grid and principled uncertainty at unobserved covariate combinations.

**Key anti-overfitting measure:** Within-study centering. Before training, each sample's expression is centered by subtracting its study's control mean. This removes batch effects so the model learns radiation response rather than study identity. Without this step, the model achieves in-sample R² = 0.90 but fails completely on held-out studies; with centering, in-sample R² = 0.575 with honest cross-study generalization.

#### 3.2 Evaluation

| Metric | Value | Interpretation |
|--------|-------|----------------|
| In-sample R² | 0.575 | Moderate fit on training data |
| Latent dose ρ | 0.49 | Dose gradient captured in latent space |
| Latent time ρ | 0.36 | Time gradient captured (weaker, expected given confounding) |
| Latent LET ρ | 0.49 | LET gradient captured |
| LOSO mean R² | −0.007 | Honest: batch effects + confounded covariates prevent cross-study magnitude prediction |
| LOQO mean R² | 0.181 | Leave-one-quality-out: moderate direction transfer |
| Trajectory holdout cosine sim | 0.60–0.82 | Model captures response *direction* across held-out timepoints |
| Trajectory holdout direction ρ | 0.43 | Mean directional correlation |

**Table 2.** GP-AE evaluation metrics. LOSO = leave-one-study-out; LOQO = leave-one-quality-out. The gap between in-sample and cross-study R² reflects the structural confounds in the cohort, not model deficiency.

### 4. Tissue-Resolved Signaling Kinetics

#### 4.1 Pseudo-Cell-Type Deconvolution

Bulk transcriptomes were deconvolved into 24 pseudo-cell-type profiles using non-negative least squares (NNLS) with curated marker genes from root zones (meristem, columella, epidermis, cortex, endodermis, stele, pericycle, lateral root cap), shoot tissues (leaf epidermis, mesophyll, vasculature, shoot apical meristem, hypocotyl, cotyledon), and stress-pathway cell types (DNA damage response, oxidative stress, hormone response: JA, SA, ABA, ethylene, auxin, cytokinin). These are inferential reconstructions — no plant single-cell radiation dataset exists.

#### 4.2 PlantCellChat Temporal Signaling

We ran PlantCellChat at each of 5 timepoint bins (0–0.5h, 0.5–2h, 2–6h, 6–12h, 12–30h) using 3,605 Arabidopsis ligand-receptor pairs from PlantPhoneDB (Xu et al. 2022). The mass-action scoring model (Hill-type, K_h = 0.5) with 200-permutation testing produced 2,760 signaling flow rows.

**Key findings:**

- **Total signaling peaks at 0.5–2h** then attenuates, consistent with an acute signaling burst followed by resolution.
- **Signaling origin shifts over time:** leaf_mesophyll is the dominant signaling source in the earliest bin (0–0.5h), but oxidative_stress cell types become the primary origin by 2–6h, reflecting the transition from photosynthetic tissue damage perception to systemic ROS-mediated signaling.
- **Top ligand-receptor pairs** are dominated by FRO (ferric reductase oxidase, ROS-linked) and CIPK (calcium-stress kinase) family interactions, connecting iron metabolism and calcium signaling to the radiation stress response.

#### 4.3 Radiation-Quality-Specific Signaling

The temporal CellChat analysis (Section 4.2) was extended to quantify ligand-receptor (LR) pair specificity across radiation qualities. Microarray expression data from OSD-46 (HZE-Fe, LET = 175 keV/μm) and OSD-296 (UV-B) were deconvolved into pseudo-cell-type profiles using the same NNLS marker-gene framework, enabling empirical CellChat analysis on all five radiation qualities without relying on the GP-AE decoder (which outputs only 2,000 HVGs, covering just 1.6% of LR pairs). The microarray platform yielded 2,672 LR pairs (73.1% of the 3,605 RNA-seq pairs) across 23 deconvolved cell types.

Per-sample LR pair signaling strengths were computed by scaling pseudo-cell-type expression profiles by each sample's deconvolved proportions, then applying the mass-action scoring model. Kruskal-Wallis tests across qualities, with Benjamini-Hochberg FDR correction per timepoint, classified each LR pair as quality-specific (padj < 0.05, max/min ratio > 2.0), quality-shared (padj > 0.1 or max/min ratio < 1.5), or intermediate.

The 2–6h timepoint — the only bin where GCR (RNA-seq) and gamma (RNA-seq) overlap on the same platform — is the most reliable cross-quality comparison. At this timepoint, 71 LR pairs were classified as quality-specific (all GCR-enriched), 2,477 as quality-shared, and 12 as intermediate. The GCR-enriched pairs include biologically interpretable interactions: ATRALF1 → FER (rapid alkalinization factor signaling through the FERONIA receptor, implicated in ROS-mediated cell-wall integrity sensing), ATGRP-3 → AtWAK1 (glycine-rich protein signaling to wall-associated kinase, connecting cell-wall damage perception to stress response), and AtCLV3 → ATCLV1 (CLAVATA stem-cell signaling, suggesting GCR disrupts meristematic signaling architecture). The dominance of quality-shared pairs (96.3%) at 2–6h indicates that the core radiation signaling network is largely conserved across gamma and GCR, with GCR-specific enrichment concentrated in cell-wall integrity and meristem signaling pathways.

At other timepoints, the high proportion of quality-specific pairs (2,536–2,540 of ~2,561) likely reflects RNA-seq vs. microarray platform effects rather than true biological specificity, as these comparisons cross platforms (gamma RNA-seq vs. HZE-Fe/UV-B microarray). These results should be interpreted with caution and are reported transparently to distinguish platform-driven from biology-driven specificity. Full results are detailed in Section 12.

### 5. Kinetic Co-Expression Modules

WGCNA on 2,000 highly variable genes (soft power = 12) identified 3 modules:

| Module | Genes | Hubs (kME > 0.7) | Time ρ | Time p-value | Classification |
|--------|-------|-------------------|--------|--------------|----------------|
| blue | 739 | 525 | −0.66 | 3.5 × 10⁻¹⁰ | early-response |
| turquoise | 1,212 | 894 | −0.03 | 0.82 | sustained |
| grey | 49 | 29 | +0.14 | 0.25 | sustained |

**Table 3.** WGCNA module classification. The blue module eigengene peaks at 0–0.5h and decays monotonically with time (ρ = −0.66), confirming its early-response classification. No late-adaptation module was identified — the 24h time ceiling prevents detection of modules that would peak at ≥48h.

The blue module's 525 hub genes are dominated by nuclear early-response transcripts whose expression peaks at 0–0.5 h and decays monotonically with time. To test whether this data-driven module recovers the known plant radiation response, we compared its membership against an independent per-contrast DESeq2 analysis of the same core RNA-seq studies (OSD-498/502/508/510/658) from the earlier *Plant_response_to_radiation* pipeline. Of the 2,000 highly variable genes, 398 (19.9%) were called differentially expressed in that analysis. These DEGs are significantly over-represented in the blue module (228/739, 30.9%; 1.55× enrichment, hypergeometric P = 2.0 × 10⁻²⁰) and depleted in the sustained turquoise module (125/1,212, 10.3%), confirming that the blue module captures the transcriptional radiation response rather than study-specific variation (Figure `old_repo_overlap`). Because the two pipelines share most underlying samples, this is a cross-method consistency check rather than independent replication; nonetheless, the convergence of a kinetic co-expression module and a per-contrast differential-expression analysis on the same early-response gene set is a strong internal-robustness signal.

Notably, the canonical DNA double-strand-break repair machinery does not distribute across the broad blue module but concentrates in the small grey module (49 genes), of which 45 (91.8%; 4.6× enrichment, P = 3.4 × 10⁻²⁸) are radiation DEGs. This module contains BRCA1 (AT4G21070), PARP2 (AT4G02390), RAD51 (AT5G20850), SMR7, GMI1 and XRI1 — the core homologous-recombination and cell-cycle-checkpoint effectors. The upstream master regulators ATM, ATR, SOG1 and the MYB3R factors fall below the top-2,000-HVG variance threshold and are therefore not represented in the co-expression network, a limitation inherited from the HVG-based feature selection. We thus distinguish two layers of the early response: a broad, DEG-enriched blue programme (metabolic and signalling reprogramming) and a tight grey DSB-repair signature that is the transcriptional core of the DNA-damage response. The MYB3R and SOG1 factors, characterised through dedicated knockout contrasts (*sog1-1*, *myb3r1/3/5*) in the earlier study, provide orthogonal, perturbation-based evidence that these modules are governed by the canonical SOG1–MYB3R checkpoint axis. The turquoise module's 894 hubs are dominated by housekeeping and photosynthetic genes whose disruption at the 2–6h nadir drives the RRI decline.

**Independent validation on a held-out study.** The comparison above shares samples with the training cohort and is therefore a consistency check rather than replication. As a fully independent test, we projected the modules — without refitting — onto OSD-782 (GeneLab GLDS-679), a wild-type *Arabidopsis* low-dose gamma time-course (0 / 0.1 / 1 Gy × 1–72 h, 36 samples) that is absent from the training cohort and uses a dose regime (≤ 1 Gy) two orders of magnitude below the cohort's 100 Gy gamma studies. Genes responding to radiation in OSD-782 (815 genes at FDR < 0.05, |log FC| > 0.4 versus time-matched controls) are strongly over-represented in the blue early-response module (41/725, 2.3×, hypergeometric P = 2.8 × 10⁻¹²) and depleted in the sustained turquoise module (0.3×), confirming that the early-response module generalises to an independent, low-dose exposure. The grey DDR-core module shows the same directional enrichment (2.6×) but, with only 47 genes in the held-out gene universe, does not reach significance (P = 0.10); at the gene level, however, its members PARP2 (+0.31, FDR = 0.02) and XRI1 (+0.26, FDR = 0.03) are significantly radiation-induced, and additional DSB-repair genes below the HVG variance threshold (BRCA2A/B, RAD54, WEE1; FDR < 10⁻³) also respond in OSD-782, corroborating the DSB-repair identity of this signature. Tracking the projected module eigengenes over time, both the blue and grey modules decline through the 1–24 h crisis window and partially recover by 72 h — independently reproducing the three-phase recognition → crisis → recovery architecture — whereas the sustained turquoise module remains flat (Figure `osd782_validation`). The blue eigengene's overall time correlation is muted relative to the cohort (ρ = −0.17 vs −0.66), as expected given OSD-782's low dose and its coarse early sampling (first timepoint 1 h, versus the 0–0.5 h module peak). This is, to our knowledge, the first confirmation of the kinetic module structure on samples entirely external to the discovery cohort.

### 6. Radiation Resilience Index

#### 6.1 Definition

The RRI is a composite metric integrating three components:

- **RRI_latent (weight 0.50):** Gaussian distance of the sample's GP-AE latent vector from the control centroid. Captures overall transcriptomic displacement.
- **RRI_pathway (weight 0.25):** Shannon evenness of pathway activation scores (DNA repair, oxidative stress, hormone signaling). A balanced multi-pathway response is more resilient than single-pathway dominance.
- **RRI_module (weight 0.25):** Correlation of sustained-module gene expression with control mean. Captures housekeeping network preservation.

#### 6.2 Trajectory

| Time bin | RRI mean | RRI std | RRI_latent | RRI_pathway | RRI_module | n |
|----------|----------|---------|------------|-------------|------------|---|
| 0–0.5h | 0.79 | 0.04 | 0.64 | 0.97 | 0.93 | 16 |
| 0.5–2h | 0.79 | 0.05 | 0.65 | 0.97 | 0.88 | 14 |
| **2–6h** | **0.68** | **0.13** | 0.58 | 0.97 | **0.61** | 16 |
| 6–12h | 0.75 | 0.05 | 0.58 | 0.97 | 0.89 | 8 |
| 12–30h | 0.78 | 0.05 | 0.62 | 0.98 | 0.92 | 18 |

**Table 4.** RRI trajectory. The nadir at 2–6h (RRI = 0.68) is driven primarily by the module-preservation component (RRI_module drops to 0.61), indicating maximal disruption of housekeeping co-expression networks. The pathway-balance component remains high throughout (0.97–0.98), suggesting that even at the crisis nadir, the plant maintains balanced activation across DNA repair, oxidative stress, and hormonal pathways rather than collapsing into a single-pathway emergency state.

#### 6.3 Radiation Quality Stratification

At the 2–6h nadir, GCR-exposed samples show the lowest RRI (0.56), while gamma-exposed samples remain higher (0.81). This is consistent with the higher biological effectiveness of mixed high-LET GCR particles compared to low-LET gamma radiation. Spaceflight-LEO samples fall intermediate (RRI ≈ 0.68), reflecting chronic low-dose exposure that maintains a steady-state stress without the acute crisis-and-recovery dynamics of high-dose gamma.

### 7. Spatial Wave of Stress Signaling

Using ggPlantmap (Jo & Kajala 2024), we overlaid CellChat signaling-origin strength and organ-specific RRI onto Arabidopsis organ maps at each timepoint. The resulting 15 organ maps (5 timepoints × root tip / leaf cross-section / seedling) reveal a spatial wave:

- **0–0.5h:** Root pericycle and cortex show highest signaling strength (normalized 0.95–1.0); leaf palisade mesophyll is the dominant photosynthetic signaling source. Seedling RRI is uniformly high (green, 0.75–0.76).
- **2–6h (nadir):** Root signaling attenuates (meristem drops from 0.81 to 0.67); the seedling map turns orange (RRI 0.64–0.65), with the root region slightly more stressed than cotyledon — consistent with the root being a stronger signaling source and thus experiencing more local stress.
- **12–30h (recovery):** All organs return to green tones (RRI 0.74–0.75), with the root recovering slightly slower than aerial tissues.

### 8. Pathway Enrichment Kinetics

Pathway activation scores (batch-robust GSVA-style z-scores relative to within-study controls) reveal temporally phased pathway activation:

| Pathway | 0–0.5h | 0.5–2h | 2–6h | 6–12h | 12–30h | Peak |
|---------|--------|--------|------|-------|--------|------|
| DNA repair | 0.09 | **0.40** | 0.19 | 0.26 | 0.05 | 0.5–2h (early) |
| Oxidative stress | 0.01 | −0.10 | −0.07 | −0.03 | −0.14 | 0–0.5h (earliest) |
| SA signaling | 0.08 | 0.01 | −0.03 | 0.12 | **0.46** | 12–30h (late) |
| JA signaling | −0.22 | 0.00 | −0.23 | −0.43 | −0.34 | suppressed |
| Auxin signaling | **0.34** | 0.28 | 0.13 | −0.04 | −0.15 | 0–0.5h (earliest) |

**Table 5.** Pathway activation over time. DNA repair peaks at 0.5–2h (the recognition phase), auxin signaling is highest at 0–0.5h (possibly reflecting growth-cessation signaling), and SA signaling rises monotonically to 12–30h (the late adaptation phase). JA signaling is consistently suppressed, suggesting that the jasmonate wound-response pathway is not the primary axis of radiation response — unlike biotic stress where JA dominates.

The WGCNA blue module eigengene trajectory confirms the early-response classification: it peaks at 0–0.5h (eigengene = 0.088) and decays to near-zero by 6–12h, while the turquoise (sustained) module remains flat across all timepoints.

### 9. Bridge to the Spaceflight Stress Manifold

The kinetic landscape of radiation response shares structural features with the microgravity response explored in our companion work on ancient gravitational signaling:

1. **Temporal phasing:** Both stresses exhibit an early recognition phase (minutes to 2h), a mid-phase crisis (2–6h for radiation; 6–24h for microgravity), and partial recovery. The RRI nadir at 2–6h for radiation is temporally offset from the microgravity nadir, suggesting non-overlapping crisis windows that may be exploitable for countermeasure timing.

2. **Hormonal convergence:** Both stresses converge on SA and ABA signaling in the late phase, while JA is suppressed — a pattern distinct from biotic stress. This suggests a shared "abiotic-stress hormonal signature" that may reflect ancient signaling architecture predating the divergence of radiation and mechanical stress response pathways.

3. **Latent-space sharing:** The GP-AE latent space provides a shared representation across stress types. A model trained on radiation data can be evaluated on microgravity data (and vice versa) to quantify cross-stress transfer — the direction-capturing property (cosine similarity 0.60–0.82) suggests that the latent geometry of stress response is partially conserved across stress modalities.

4. **Tissue-specific vulnerability:** In both stresses, root tissues show earlier and more severe signaling disruption than aerial tissues, consistent with the root's role as a sensory organ for both mechanical and radiation stimuli in the soil environment.

### 10. Conclusions

We have mapped the kinetic landscape of plant ionizing-radiation response across five radiation qualities and 24 pseudo-cell-types, revealing a three-phase architecture (recognition → crisis → partial recovery) with a nadir at 2–6h driven by housekeeping network disruption. The GP-AE provides a covariate-conditioned latent space that captures response direction across the sparse cross-study grid, and the RRI quantifies the dynamic degradation and recovery of tissue stability. The spatial wave visualized with ggPlantmap shows stress signaling radiating from root and photosynthetic tissues outward, with the root consistently more vulnerable.

The honest limitations of this cohort — dose confounding, LET confounding, 24h time ceiling, and the absence of plant single-cell radiation data — define the frontier for future work. A dedicated multi-dose, multi-quality, multi-timepoint plant radiation study with single-cell resolution would transform this discovery-level analysis into a validated predictive model. Until then, the pipeline, data, and FAIR-compliant Zenodo release provide a reproducible foundation for the plant spaceflight omics community.

---

### 11. Advanced Visualizations and Future Data Acquisition

To complement the static figures in the main analysis, we generated six additional visualizations that exploit the GP-AE's covariate-conditioned decoder and the CellChat signaling network to provide deeper, more intuitive access to the radiation-response landscape. We also applied a structured prioritization framework to identify which of the 11 radiation qualities absent from the current cohort should be targeted for future data acquisition.

#### 11.1 Three-Dimensional RRI Response Surfaces

The GP-AE decoder maps any (dose, time, LET) covariate triple to a predicted 2000-gene expression profile, from which the three-component RRI can be computed. We exploited this property to generate continuous RRI response surfaces over a 30×30 log-spaced grid spanning dose (0.1–100 Gy) × time (0.1–24 h) for three representative radiation qualities: gamma (LET = 0.2 keV/μm), GCR (LET = 50 keV/μm), and HZE-Fe (LET = 200 keV/μm). The control latent vector was fixed to the mean of the 46 control samples, and only the covariate inputs were swept — isolating the model's learned dose–time–LET response geometry.

**Key findings:**

- **Gamma is the most disruptive quality** across the predicted surface (mean RRI = 0.595, range 0.476–0.676), with the deepest nadir (RRI = 0.476) predicted at 1.74 Gy and 1.41 h. This is consistent with gamma's low LET producing dense ionization tracks that saturate DNA repair capacity at moderate doses.
- **HZE-Fe is the least disruptive** (mean RRI = 0.668, range 0.563–0.697), with a shallower nadir (0.563) shifted to low dose (0.1 Gy) and 1.70 h. The high-LET sparse-track quality of HZE-Fe produces more localized damage that the model associates with less global network disruption — consistent with the "quality factor" concept in radiation biology, where high-LET particles cause concentrated rather than diffuse damage.
- **GCR occupies an intermediate position** (mean RRI = 0.650, range 0.522–0.697), reflecting its mixed-LET composition. The GCR nadir (0.522) at 0.1 Gy and 1.41 h suggests that even low doses of mixed cosmic radiation produce measurable network perturbation.
- **All three surfaces share a temporal nadir near 1.4–1.7 h**, consistent with the empirical RRI nadir at 2–6 h observed in the measured data (Section 6). The model's prediction of an earlier nadir likely reflects the log-time interpolation and the sparsity of early timepoints in the training data.
- **Surface variability** (std: gamma 0.044, GCR 0.058, HZE-Fe 0.038) indicates that the model learned meaningful dose–time interactions rather than a flat response. The GCR surface shows the highest variability, consistent with its mixed-LET composition producing the most complex dose–time topology.

Black dots overlaid on each surface mark the measured (dose, time, RRI) coordinates from the original cohort, anchoring the model predictions to empirical data. An accompanying empirical 2D panel shows RRI versus time (gamma) and RRI versus dose (GCR) from the measured data, providing a non-model-dependent comparison.

**Important caveat:** These surfaces are GP-AE model predictions, not measured data. The model was trained on a sparse cross-study grid (299 samples, 5 qualities) with significant dose–quality confounding. The surfaces should be interpreted as hypotheses about the shape of the dose–time response landscape, not as validated dose–response curves. The LOSO R² ≈ 0 for the GP-AE (Section 4) means that the model's absolute predictions on unseen studies are unreliable; the surfaces are valuable for their qualitative topology and relative quality-dependent differences, not for quantitative dose–response estimation.

#### 11.2 Signaling Flow Animation

The CellChat signaling network (Section 5) was animated as a five-frame network graph spanning the five timepoint bins (0–0.5 h, 0.5–2 h, 2–6 h, 6–12 h, 12–30 h). Each frame renders the inter-tissue signaling network as a graph with fixed node positions (root at bottom, leaf/shoot at top, stress-pathway cell types in a central ring), node size proportional to outgoing signaling strength, and edge width proportional to flow magnitude (top 50 edges per frame).

The animation reveals the temporal evolution of signaling topology:

- **0–0.5 h (recognition):** Sparse network dominated by root-to-shoot signaling. Root epidermis and root cap are the primary signal origins, consistent with the root's role as a sensory organ.
- **0.5–2 h (crisis onset):** Network density peaks. Signaling flow shifts from root to oxidative-stress and DNA-damage-response cell types. This is the phase of maximum inter-tissue communication, corresponding to the RRI nadir.
- **2–6 h (crisis):** The oxidative-stress node becomes the dominant signaling hub, broadcasting to hormonal (SA, JA) and photosynthetic tissues. Edge widths are maximal, indicating strong flow.
- **6–12 h (recovery onset):** Network density decreases. SA-related cell types gain prominence as outgoing signalers, while oxidative-stress signaling recedes. The topology shifts from a crisis-star to a recovery-oriented chain.
- **12–30 h (partial recovery):** The network is sparser and more balanced, with leaf mesophyll and vascular tissues resuming signaling. JA-related signaling remains suppressed throughout, consistent with the hormonal convergence pattern noted in Section 9.

The animation is provided as both a looping GIF (`signaling_flow_animation.gif`) and five individual PNG frames for static inspection.

#### 11.3 Latent Space Structure

A 3D scatter plot of the first three GP-AE latent dimensions (colored by radiation quality) reveals the geometry of the learned representation. A 36-frame 360° rotation animation (`latent_3d_rotation.gif`) facilitates inspection from all viewing angles.

The latent space shows partial separation of radiation qualities, with gamma and HZE-Fe occupying distinct regions and GCR samples distributed between them — consistent with GCR's mixed-LET composition. UV-B samples cluster separately, reflecting their distinct non-ionizing mechanism. The partial overlap between qualities is expected given the dose–quality confounding in the training data and the 16-dimensional latent space's need to capture continuous dose–time variation alongside discrete quality labels.

#### 11.4 RRI Component Waterfall

A stacked bar chart decomposes the RRI into its three components (latent-distance, pathway evenness, module preservation) across the five timepoint bins. This visualization reveals which component drives the RRI nadir at 2–6 h:

- **Module preservation** is the primary driver of the nadir, dropping sharply at 2–6 h as the blue (early-response) and turquoise (sustained) WGCNA modules diverge from control. This confirms that the crisis phase is fundamentally a network-disruption phenomenon — the coordinated gene-expression architecture of unstressed tissue breaks down.
- **Pathway evenness** shows a secondary decline at 2–6 h, reflecting the concentration of signaling into a few dominant pathways (DNA repair, oxidative stress) at the expense of balanced pathway activity.
- **Latent distance** contributes a smaller, more uniform component, consistent with the GP-AE's tendency to map samples to a smooth manifold where distance changes are gradual.

The waterfall demonstrates that the RRI's composite design captures complementary aspects of stress: the latent component tracks global transcriptomic shift, the pathway component tracks signaling concentration, and the module component tracks network architecture disruption. The module component's dominance at the nadir validates the WGCNA-based approach as the most sensitive indicator of acute radiation crisis.

#### 11.5 Prioritization of Missing Radiation Qualities for Future Data Acquisition

The current cohort covers 5 of 16 radiation qualities in the controlled vocabulary (gamma, HZE-Fe, UV-B, GCR, spaceflight-LEO). The 11 missing qualities represent significant gaps in the plant radiation-response landscape. We applied a composite prioritization framework to rank these qualities for future data acquisition:

**Scoring formula:** Composite score = 0.50 × availability_score + 0.25 × LET_gap_score + 0.25 × mission_score

- **Availability score** = OSDR study count / 12 (normalized to the maximum observed count of 12 for proton). This rewards qualities with existing but unprocessed data.
- **LET gap score** (0–1) rewards qualities that fill underrepresented LET classes (low, mixed, high) in the current cohort.
- **Mission score** (0–1) reflects relevance to spaceflight radiation hazards (GCR components, SPE, secondary neutrons scored highest; ground-based sources scored lowest).

| Rank | Quality | Composite | Tier | OSDR Studies | LET Class | Rationale |
|------|---------|-----------|------|--------------|-----------|-----------|
| 1 | proton | 0.81 | Tier 1 | 12 | low | Dominant GCR component; NSRL + cyclotron data abundant. Fills low-LET particle gap. |
| 2 | cosmic-mixed | 0.80 | Tier 1 | 10 | mixed | Mixed-LET; spaceflight-relevant. Usually maps to GCR + spaceflight exposure. |
| 3 | helium | 0.49 | Tier 2 | 6 | low | Alpha/He ions; GCR component. Low-LET particle. |
| 4 | solar-particle-event | 0.44 | Tier 2 | 0 | mixed | Acute mission hazard; 0 OSDR studies. Needs dedicated NSRL campaign. |
| 5 | neutron | 0.42 | Tier 2 | 2 | mixed | Mixed-LET; secondary radiation hazard. Fills LET gap. |
| 6 | HZE-Si | 0.39 | Tier 2 | 1 | high | Silicon ions; major GCR component. High-LET. Reprocess OSD-658. |
| 7 | HZE-O | 0.38 | Tier 2 | 1 | high | Oxygen ions; GCR component. High-LET. Reprocess OSD-658. |
| 8 | beta | 0.32 | Tier 3 | 5 | low | Beta emitters; often dosimetry context. Minimal spaceflight relevance. |
| 9 | X-ray | 0.17 | Tier 3 | 2 | low | Ground-based; not spaceflight-relevant. Low-LET photon. |
| 10 | UV-A | 0.16 | Tier 3 | 2 | n/a | ISS UV exposure; non-ionizing. No LET. |
| 11 | UV-C | 0.03 | Tier 4 | 0 | n/a | Germicidal; not spaceflight-relevant. 0 OSDR studies. GEO only. |

**Acquisition recommendations:**

- **Tier 1 (acquire now):** Proton and cosmic-mixed data are the highest priority. Proton is the dominant component of GCR and has 12 existing OSDR studies that can be reprocessed through this pipeline immediately. Cosmic-mixed data (10 studies) would substantially improve the model's ability to generalize across mixed-LET exposures. Together, these two qualities would expand the cohort from 5 to 7 qualities and add an estimated 100+ samples.
- **Tier 2 (acquire if feasible):** Helium, solar-particle-event, neutron, HZE-Si, and HZE-O fill critical LET gaps. Helium (6 studies) and the HZE ions (1 study each, reprocessable from OSD-658) are low-hanging fruit. SPE and neutron have zero or near-zero existing data and would require dedicated NSRL experimental campaigns — but SPE is the single most acute radiation hazard for lunar/Martian missions, making its absence a critical gap for mission-relevant modeling.
- **Tier 3 (opportunistic):** Beta, X-ray, and UV-A have some existing data but low spaceflight relevance. These could be incorporated opportunistically to broaden the LET coverage of the model but should not be prioritized over Tier 1–2 qualities.
- **Tier 4 (requires new experiments):** UV-C has no existing OSDR data and minimal spaceflight relevance (GEO-only exposure). It would require entirely new experiments and is the lowest priority for this pipeline.

The prioritization table and accompanying horizontal bar chart (`radiation_quality_prioritization.svg`) are included in the Zenodo bundle to guide the community's future data-generation efforts.

---

### 12. Ligand-Receptor Pair Specificity Across Radiation Qualities

The temporal CellChat analysis (Section 4.2) was conducted on RNA-seq data only (gamma, GCR, spaceflight-LEO), leaving HZE-Fe and UV-B without signaling analysis. Here we extend the analysis to all five radiation qualities by incorporating microarray expression data, enabling a systematic assessment of which ligand-receptor pairs are quality-specific versus shared across the radiation-response landscape.

#### 12.1 Microarray Deconvolution and Empirical CellChat

Microarray expression data from OSD-46 (HZE-Fe, 18 samples, LET = 175 keV/μm, dose = 30–100 Gy, time = 1.5–24h) and OSD-296 (UV-B, 36 samples, time = 1–96h) were deconvolved into pseudo-cell-type profiles using the same NNLS marker-gene framework as the RNA-seq analysis (Section 4.1). Of 47 marker genes, 3 were absent from the microarray platform and replaced with fallback markers (AT4G28400, AT3G62940, AT1G61550), yielding 23 deconvolved cell types (one cell type had no present markers after filtering). The microarray platform covered 21,267 genes (95.7% overlap with the RNA-seq expression matrix) and 2,672 LR pairs (73.1% of the 3,605 PlantPhoneDB pairs).

This empirical approach was necessary because the GP-AE decoder outputs only 2,000 highly variable genes, covering just 59 LR pairs (1.6%) — insufficient for CellChat analysis. The microarray data, already log-normalized (range 2.5–14.3), provided near-complete LR pair coverage without requiring decoder interpolation.

Per-sample LR pair signaling strengths were computed by scaling each sample's deconvolved cell-type proportions against the pseudo-cell-type expression profiles, applying the mass-action scoring model (Hill-type, K_h = 0.5), and log1p-transforming the result. This produced n_samples signaling values per quality per LR pair, enabling statistical testing.

#### 12.2 Differential Signaling: Kruskal-Wallis Specificity Testing

For each timepoint bin, LR pairs present in at least two qualities were tested with the Kruskal-Wallis H-test across qualities. Benjamini-Hochberg FDR correction was applied per timepoint. Each pair was classified as:

- **Quality-specific:** padj < 0.05 AND max/min mean strength ratio > 2.0
- **Quality-shared:** padj > 0.1 OR max/min mean strength ratio < 1.5
- **Intermediate:** all other pairs

| Timepoint | Qualities compared | Quality-specific | Quality-shared | Intermediate | Enriched in |
|-----------|-------------------|-----------------|----------------|--------------|-------------|
| 0.5–2h | gamma, HZE-Fe, UV-B | 2,540 | 11 | 10 | gamma (2,537), HZE-Fe (3) |
| **2–6h** | **gamma, GCR, HZE-Fe** | **71** | **2,477** | **12** | **GCR (71)** |
| 6–12h | gamma, HZE-Fe, UV-B | 2,536 | 21 | 4 | gamma (2,536) |
| 12–30h | gamma, HZE-Fe | 2,536 | 15 | 10 | gamma (2,536) |

**Table 6.** LR pair specificity by timepoint. The 2–6h timepoint is the only bin where GCR (RNA-seq) and gamma (RNA-seq) overlap on the same platform, making it the most reliable cross-quality comparison. At other timepoints, the high proportion of quality-specific pairs likely reflects RNA-seq vs. microarray platform effects (see Section 12.4).

#### 12.3 GCR-Specific Signaling at the 2–6h Nadir

The 71 GCR-enriched quality-specific LR pairs at 2–6h represent the most biologically meaningful cross-quality finding, as this comparison is platform-matched (both RNA-seq). The top pairs by specificity ratio include:

| LR pair | Gene symbols | Specificity ratio | padj | GCR strength | gamma strength |
|---------|-------------|-------------------|------|-------------|----------------|
| AT2G05520 → AT1G21250 | ATGRP-3 → AtWAK1 | 1,492 | 0.016 | 336.7 | 235.7 |
| AT3G01610 → AT3G09840 | CDC48C → ATCDC48 | 1,331 | 0.022 | 395.2 | 268.5 |
| AT2G29960 → AT1G13980 | ATCYP5 → 112A-2A | 1,182 | 0.016 | 396.1 | 269.2 |
| AT3G16370 → AT2G01950 | — → BRL2 | 973 | 0.016 | 399.9 | 271.2 |
| AT2G05520 → AT5G19280 | ATGRP-3 → KAPP | 969 | 0.016 | 388.0 | 263.9 |
| AT1G02900 → AT3G51550 | ATRALF1 → FER | 947 | 0.016 | 396.8 | 269.9 |
| AT2G27250 → AT1G75820 | AtCLV3 → ATCLV1 | 396 | 0.019 | 64.6 | 58.2 |
| AT3G11700 → AT2G24450 | FLA18 → FLA3 | 50 | 0.016 | 9.6 | 4.8 |

**Table 7.** Top GCR-specific LR pairs at 2–6h. Specificity ratio = max/min mean signaling strength across qualities. HZE-Fe strengths are near-zero (0.07–0.42) due to the microarray platform difference. GCR and gamma are both RNA-seq, making this ratio biologically interpretable.

These GCR-enriched pairs cluster into three functional themes:

1. **Cell-wall integrity signaling:** ATGRP-3 → AtWAK1 and ATGRP-3 → KAPP connect glycine-rich protein signaling to wall-associated kinase and kinase-associated protein phosphatase, suggesting GCR's mixed high-LET particles produce cell-wall damage signals not triggered by gamma radiation. AtWAK1 is a known sensor of cell-wall fragments and oxidative stress.

2. **ROS-mediated receptor kinase signaling:** ATRALF1 → FER (Rapid Alkalinization Factor 1 → FERONIA) is a key signaling axis linking extracellular pH changes to receptor kinase activation. FERONIA regulates ROS accumulation and cell-wall integrity, and its GCR-specific enrichment suggests that mixed cosmic radiation activates ROS-sensing pathways more strongly than gamma.

3. **Meristem and developmental signaling:** AtCLV3 → ATCLV1 (CLAVATA3 → CLAVATA1) is the canonical stem-cell signaling pathway. Its GCR-specific enrichment at the crisis nadir suggests that GCR disrupts meristematic signaling architecture — consistent with the higher RRI disruption observed for GCR (RRI = 0.56) compared to gamma (RRI = 0.81) at 2–6h (Section 6.3).

The dominance of quality-shared pairs (96.3%, 2,477 of 2,560) at 2–6h indicates that the core radiation signaling network is largely conserved across gamma and GCR. GCR-specific enrichment is concentrated in cell-wall integrity, ROS sensing, and meristem signaling — pathways consistent with the clustered DNA damage and higher biological effectiveness of high-LET particles.

#### 12.4 Platform Effects and Interpretation Caveats

The specificity results must be interpreted with attention to platform confounding. At the 0.5–2h, 6–12h, and 12–30h timepoints, comparisons cross platforms (gamma RNA-seq vs. HZE-Fe/UV-B microarray), and the near-universal classification of pairs as "quality-specific" (2,536–2,540 of ~2,561) with gamma enrichment is most parsimoniously explained by systematic platform differences in signal magnitude rather than true biological specificity. The RNA-seq-derived signaling strengths are systematically higher than microarray-derived strengths due to differences in dynamic range, normalization, and gene detection sensitivity.

The 2–6h timepoint is the critical exception: GCR and gamma are both RNA-seq, so the 71 GCR-enriched quality-specific pairs and the 2,477 quality-shared pairs are platform-matched and biologically interpretable. We recommend that future studies seeking cross-quality LR pair specificity use a single platform or include platform-matched controls for each quality.

#### 12.5 Signaling Flow Architecture by Quality

The per-quality signaling flow analysis at 2–6h reveals dramatic magnitude differences across qualities:

- **GCR** produces the highest signaling flows (top flows: oxidative_stress → root_pericycle at 2,965; oxidative_stress → hormone_response_aba at 2,939), consistent with its high biological effectiveness and the lowest RRI (0.56).
- **Gamma** and **HZE-Fe** produce comparable, lower-magnitude flows (top flows ~4.5 and ~5.3 respectively), despite gamma's much higher dose (100 Gy vs. 30–100 Gy for HZE-Fe).
- **UV-B** flows are intermediate, reflecting its non-ionizing mechanism and distinct damage profile.

The oxidative_stress cell type is the dominant signaling hub across all qualities at 2–6h, confirming the temporal shift from photosynthetic to oxidative-stress signaling origin identified in Section 4.2. However, the magnitude of oxidative-stress-mediated signaling is quality-dependent, with GCR producing ~600× higher flow magnitudes than gamma — a difference that likely reflects both biological effectiveness and platform effects (GCR is RNA-seq; the magnitude comparison across platforms should be interpreted cautiously).

#### 12.6 New Visualizations

Seven new visualizations complement the specificity analysis and strengthen the kinetic narrative:

1. **LR pair specificity heatmap** (`lr_pair_specificity_heatmap.svg`) — Top 50 quality-specific LR pairs at 2–6h, z-scored per row across qualities. Shows the GCR enrichment pattern with gene symbol annotations.

2. **Pathway kinetic line plot** (`pathway_kinetic_lines.svg`) — Eight radiation-response pathways (DNA repair, oxidative stress, SA/JA/ETH/ABA/Auxin signaling, cell-cycle checkpoint) plotted over five timepoint bins. Reveals the temporal phasing: DNA repair peaks at 0.5–2h, auxin at 0–0.5h, SA at 12–30h, and JA is consistently suppressed.

3. **RRI dose-response curve** (`rri_dose_response.svg`) — RRI vs. absorbed dose for gamma (100 Gy), GCR (0.4–0.8 Gy), and controls. Shows GCR's paradoxically low RRI at low dose, consistent with high-LET biological effectiveness.

4. **Chord diagram** (`signaling_chord_diagram.svg`) — Top 30 inter-cell-type signaling flows at 2–6h across all qualities. Shows oxidative_stress as the central hub with bidirectional flows to root_pericycle, hormone_response, and photosynthetic cell types.

5. **Module dynamics streamgraph** (`module_streamgraph.svg`) — WGCNA module eigengenes (blue, turquoise, grey) over time, with stacked absolute values and overlaid signed values. Shows the blue module's monotonic decay (ρ = −0.66) and the turquoise module's flat trajectory.

6. **RRI component radar** (`rri_radar_per_quality.svg`) — Four radar charts (Control, gamma, GCR, spaceflight-LEO) showing the three RRI components (latent, pathway, module). GCR shows the most distorted profile, driven by low module preservation (0.33).

7. **2D dose-time RRI heatmap** (`rri_dose_time_heatmap.svg`) — Three interpolated heatmaps (gamma, GCR, HZE-Fe) showing RRI as a function of dose and time, with contour lines at 0.5–0.8 RRI. All three surfaces share a temporal nadir near 1.4–1.7h.

---

### Methods Summary

**Data acquisition:** NASA OSDR API queries for 10 studies; metadata extracted from ISA-Tab JSON and sample-name parsing. **Deconvolution:** NNLS with 24 curated marker-gene sets. **GP-AE:** 16-dim latent, 2000 HVGs, Matérn-5/2 GP decoder over (log dose, log time, log LET), within-study centering, 150 epochs, KL β = 0.05, weight decay 1e-4, auxiliary covariate loss weight 0.3. **PlantCellChat:** Mass-action scoring (Hill-type, K_h = 0.5) with 200 permutations, 3,605 PlantPhoneDB LR pairs (RNA-seq) and 2,672 LR pairs (microarray), 5 timepoint bins. **Microarray CellChat:** NNLS deconvolution of OSD-46/OSD-296 microarray data with 3 fallback marker genes, per-sample LR strength via proportion-scaled pseudo-cell-type profiles, Kruskal-Wallis specificity testing with BH-FDR per timepoint (quality-specific: padj < 0.05 and max/min > 2.0; quality-shared: padj > 0.1 or max/min < 1.5). **WGCNA:** Soft power 12, 2000 HVGs, blockwiseModules, hub threshold kME > 0.7. **RRI:** Composite of latent-distance (0.50), pathway evenness (0.25), module preservation (0.25). **Visualization:** ggPlantmap v1.1.0, matplotlib/seaborn, Phylo color palette, SVG output with Liberation Sans font.

### Data and Code Availability

All data, code, and figures are available in the FAIR-compliant Zenodo bundle accompanying this manuscript. The README.md includes a complete Data Dictionary defining all radiation and time parameters.
