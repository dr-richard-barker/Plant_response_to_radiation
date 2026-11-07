# Technical Benchmark Report: Cross-Species Spaceflight Cosmic Radiation & Redox Decoders

> **Status: draft under audit. Not for citation.** An AI agent (Google Antigravity) generated this report on 2–4 Nov 2026. On 2026-11-07 an ABAI check compared it with the scripts and tables in this repository, and it was corrected to match them. Every number below is read from the named output file. Items that could not be verified are marked **TODO (Richard)**. See §5 of the companion manuscript for the audit record.

**Project:** Cross-species radiation decoding (*Brassica rapa* GCR dataset × *Arabidopsis thaliana* tools)  
**Companion manuscript:** `docs/manuscript_brappa_radiation_decoding.md`  
**Date:** November 2026 (audited 2026-11-07)  
**Pipeline scripts:** `code/16_brappa_llgcss_ingest.py` through `code/20_generate_figures_cross_species.py`

---

## 1. Executive Summary

This report documents scripts 16–20. They score *B. rapa* GCR contrasts (40 vs 0 cGy; WT vs *anthocyaninless*; DRS vs RNAlater fixative) with tools built on *Arabidopsis thaliana*. The dataset's 39 libraries, 31 of them in the processed matrix used here, carry a weak dose signal. Seven genes reach padj < 0.05 in WT, none in `anl`, one for dose × genotype and none for dose × fixative.

### Benchmark matrix
| Tool | What it actually uses | Target in *B. rapa* | Result (from the output file) | Output file |
|---|---|---|---|---|
| **DSRS** (script 17) | 27 *Arabidopsis* reference contrasts × 123 seed panels (`decoder_nes_matrix_v7.csv`, local to `deepspace-seed-stress-decoder`) | WT 40 cGy, `anl` 40 cGy, dose × genotype, dose × fixative | WT: top match `ug_root_light` (ρ = 0.259); radiation families rank 9th–11th of 13. Dose × genotype: top `rad_100Gy_co60_90m_a` (ρ = 0.254), `radiation_GCR` family first (0.150). No significance test. | `results/dsrs_brappa_stressor_rankings.csv`, `results/dsrs_brappa_family_rankings.csv` |
| **Redox heuristic** (script 18) | A hand-assigned weight table for 18 sensor genes; **not** the Redox_decoder CVAE | Same four contrasts | WT: H₂O₂ 51.1%, 3.84 h. `anl`: paraquat 24.9% (H₂O₂ 23.7%), 2.23 h. No sensor gene is DE. | `results/redox_decoder_brappa_predictions.csv` |
| **RRI** (script 19) | Expression-space distance, pathway evenness and housekeeping correlation, all referenced to WT 0 cGy; **not** the *Arabidopsis* GP-AE | 31 libraries | WT 0 cGy 0.762; WT 40 cGy 0.636; `anl` 0.637 / 0.637. Lower in every non-reference group by construction. | `results/brappa_rri_sample_scores.csv` |
| **Fixative comparison** (script 19) | WT dose effect computed within DRS vs within RL samples | 30,756 expressed genes | Pearson r = 0.297 (P < 10⁻¹⁵); 0 genes with a significant dose × fixative interaction | console output of script 19; RRI by fixative from `results/brappa_rri_sample_scores.csv` |
| **DNA-repair and flavonoid orthologs** (script 19) | Six orthologs looked up in the DESeq2 table | WT and `anl` 40 cGy | BRCA1 +2.03, PARP2 +0.86, F3H +1.51 log₂FC in WT; **none significant** (padj 0.72–1.0) | `results/brappa_preservative_and_ddr_summary.csv` |
| **Cell-type shares** (script 18) | Hand-assigned per-cell-type scores for the same 18 sensors | Same four contrasts | Shares lie in 4.6–13.2%; highest are stele/vascular (11.3–13.2%). Not data-derived. | `results/redox_decoder_organ_cell_scores.csv` |

---

## 2. Data Dictionary

### 2.1 Ingested datasets
- **`data/brappa_master_radiation_contrasts.csv`** (42,278 genes; 18,770 with an ortholog):
  - `BRA_id`: *B. rapa* locus in the radiation tables (Ensembl Plants Brapa_1.0 / Chiifu `Bra######`, upper-cased).
  - `arabidopsis_agi`: TAIR10 locus from the 18,770 one-to-one reciprocal best hits.
  - `WT_40cGy_vs_0cGy_log2FC`, `_padj`: DESeq2 (iDEP) WT dose contrast; padj as reported by iDEP.
  - `antho_less_40cGy_vs_0cGy_log2FC`, `_padj`: `anl` dose contrast.
  - `inter_dose_x_genotype_log2FC`, `_padj`: interaction term `condition_40_CGY × genotype_ANTHO_LESS`.
  - `inter_dose_x_preservative_log2FC`, `_padj`: interaction term `condition_40_CGY × preservative_RL`.
- **`data/brappa_sig_*.csv`**: AGI-keyed query signatures (16,653 genes for the WT and `anl` contrasts; 16,796 for the interactions).

### 2.2 Output results
- **`results/dsrs_brappa_stressor_rankings.csv`**: `query`, `stressor` (reference contrast), `stressor_class` (family), `spearman_rho` and `pearson_r` across the 123 panel NES values.
- **`results/brappa_rri_sample_scores.csv`**: `sample_id`, `condition`, `genotype`, `preservative`, `is_control` (WT 0 cGy = the reference), `RRI`, `RRI_latent` (an expression-space distance score, despite the name), `RRI_pathway`, `RRI_module`.
- **`results/brappa_preservative_and_ddr_summary.csv`**: `symbol`, `agi`, `description`, `bra_id`, `WT_40cGy_log2FC`, `AnthoLess_40cGy_log2FC`, `Inter_Genotype_log2FC`, `Inter_Preservative_log2FC`. This file holds only the DNA-repair/flavonoid genes; the fixative correlation is printed by script 19, not saved.

---

## 3. Algorithms (as implemented)

### 3.1 DSRS projection (script 17)
Genes are sorted by signed log₂FC r_i. For panel S_k:

P_hit(i) = Σ_{j ≤ i, g_j ∈ S_k} |r_j| / Σ_{g_j ∈ S_k} |r_j|,  P_miss(i) = #{j ≤ i : g_j ∉ S_k} / (N − |S_k|)

ES(S_k) is the value of P_hit − P_miss at its largest absolute deviation (signed), not an argmax. NES divides ES by the mean same-signed ES from 60 permutations. Each permutation shuffles the fold-change values across fixed gene positions, so gene labels stay in place. The query's 123-panel NES vector is then correlated with each column of the reference matrix M ∈ ℝ^{123×27}. No p-values are computed.

### 3.2 Redox heuristic (script 18)
For each sensor i with an ortholog fold change (16 of 18 genes), the script sets a_i = max(0, LFC_i) + 0.1·|LFC_i| and computes

score_k = Σ_i a_i · W_{i,k},  P(stim_k) = softmax(score)_k

W is a hand-assigned table (six stimuli). Timing is Σ a_i·t_i / Σ a_i over hand-assigned response times t_i. Cell-type shares are Σ a_i·c_{i,type}, normalised across 11 types. The Redox_decoder CVAE (its 20,869-gene input and its latent space) is not loaded.

### 3.3 Radiation Resilience Index (script 19)
RRI = 0.50 · exp(−‖x_i − x̄_ref‖² / 2σ²) + 0.25 · H(p)/log K + 0.25 · max(0, corr(x_{i,hk}, x̄_{ref,hk}))

- x_i is the sample's processed expression over ortholog-mapped genes.
- "ref" is the WT 0 cGy group (n = 8), and σ is the median pairwise distance within it.
- p holds the normalised mean expression of three hand-picked gene lists (K = 3).
- "hk" is the 20% of genes with the lowest variance within the reference group.

---

## 4. Results Summary (copied from the output files)

### 4.1 DSRS family ranking: dose × genotype interaction (`results/dsrs_brappa_family_rankings.csv`)
```
Stressor family      Mean Spearman ρ    Mean Pearson r
------------------------------------------------------
radiation_GCR             0.1504            0.1522
radiation_lowdose         0.1321            0.1697
hypergravity              0.1011            0.0671
radiation_acute           0.0914            0.0884
ethylene                  0.0556            0.1688
desiccation               0.0457            0.1875
osmotic                   0.0411            0.0790
temperature               0.0070            0.0577
tropism_photo            -0.0273           -0.0405
microgravity             -0.0374            0.0403
tropism_gravi            -0.1508           -0.1405
low_oxygen               -0.1572           -0.1219
uv                       -0.2368           -0.1774
```
For the **WT 40 cGy** query the same table ranks uv (0.197), tropism_gravi (0.119) and tropism_photo (0.112) first. radiation_lowdose (−0.031), radiation_GCR (−0.051) and radiation_acute (−0.081) rank 9th–11th of 13.

### 4.2 Redox heuristic scores (%) and weighted response time (`results/redox_decoder_brappa_predictions.csv`)
```
Contrast                         Top score        Score (%)  Runner-up (%)       Weighted time
-----------------------------------------------------------------------------------------------
WT_GCR_40cGy                     H2O2             51.1       paraquat 16.2       3.84 h
AnthoLess_GCR_40cGy              paraquat         24.9       H2O2 23.7           2.23 h
Interaction_Dose_x_Genotype      paraquat         21.8       high_light 21.7     5.09 h
Interaction_Dose_x_Preservative  H2O2             23.0       high_light 18.6     3.90 h
```
These scores come from hand-assigned weights, not from a trained model.

### 4.3 RRI by genotype and dose (`results/brappa_rri_sample_scores.csv`)
```
Genotype     Dose     n   RRI (mean ± SD)    RRI_latent   RRI_pathway   RRI_module
-----------------------------------------------------------------------------------
WT           0_cGy    8   0.7621 ± 0.0468    0.7876       0.9922        0.4809
WT           40_cGy   7   0.6361 ± 0.0408    0.7135       0.9912        0.1260
antho_less   0_cGy    8   0.6371 ± 0.0234    0.7392       0.9907        0.0794
antho_less   40_cGy   8   0.6372 ± 0.0384    0.7292       0.9871        0.1031
```
WT 0 cGy is the reference group for every component.

### 4.4 RRI by fixative and dose (both genotypes)
```
Fixative   Dose     n   Mean RRI   RRI SD   RRI_latent   RRI_module
--------------------------------------------------------------------
DRS        0_cGy    8   0.7060     0.078    0.7748       0.2830
DRS        40_cGy   8   0.6551     0.017    0.7541       0.1228
RL         0_cGy    8   0.6932     0.074    0.7521       0.2772
RL         40_cGy   7   0.6156     0.045    0.6851       0.1035
```
The WT dose effect in DRS vs in RL correlates at r = 0.2968 across 30,756 expressed genes (P < 10⁻¹⁵; printed by script 19). No gene reaches padj < 0.05 for dose × fixative. **TODO (Richard):** run an equivalence test before claiming fixative parity.

---

## 5. Reproduction

The scripts read local files outside this repository: `~/Documents/B_rappa_LLGCSS` (DESeq2 tables, git-ignored there; ortholog map), `~/Documents/deepspace-seed-stress-decoder` (panels, `decoder_nes_matrix_v7.csv`) and `~/Documents/Redox_decoder` (named in script 18, which does not read from it). They will not run elsewhere until those inputs are public. **TODO (Richard):** give the deposit or accession for each.

```bash
python3 code/16_brappa_llgcss_ingest.py         # orthology join, query signatures
python3 code/17_run_dsrs_classifier.py          # DSRS projection (123 panels × 60 permutations)
python3 code/18_run_redox_cvae_decoder.py       # hand-weighted redox heuristic
python3 code/19_brappa_rri_and_antioxidant_shield.py   # RRI, fixative correlation, DNA-repair lookup
python3 code/20_generate_figures_cross_species.py      # Figures 1-6 (PNG + PDF)
```

On 2026-11-07, scripts 16–19 were re-run in a scratch copy and reproduced every committed table in `data/brappa_*.csv` and `results/*.csv` exactly.
