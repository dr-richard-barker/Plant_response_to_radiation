#!/usr/bin/env python3
"""
20_generate_figures_cross_species.py
Phase 4 — Generate publication-quality figures (Figures 1–6) for the
pure Cross-Species Radiation, Redox & Fixative Decoding manuscript.

Figures produced (100% floral scent-free):
- Fig 1: Conceptual architecture of the cross-species radiation diagnostic stack
- Fig 2: DSRS Stress Recognition & Cross-Species Spaceflight Concordance
- Fig 3: Redox CVAE Decoding & The "Antioxidant Shield Deficit" (WT vs antho_less)
- Fig 4: Radiation Resilience Index (RRI) Degradation & Housekeeping Collapse
- Fig 5: Preservative Fixative Parity (DRS vs RNAlater) & DNA Damage Repair / Flavonoid Activation
- Fig 6: Spatial Anatomical Vulnerability Across Leaf, Root, and Floral Tissues

Author: Antigravity Assistant & Dr. Richard Barker
Date: November 2026
"""

import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import textwrap

WORKSPACE_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = WORKSPACE_DIR / "results"
FIG_DIR = RESULTS_DIR / "figures"
FIG_DIR.mkdir(exist_ok=True, parents=True)

# Styling configuration
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.titlesize": 13,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight"
})

print("=" * 70)
print("20_generate_figures_cross_species: Generating Pure Radiation Publication Figures")
print("=" * 70)

# Load pure radiation results
dsrs_ranks = pd.read_csv(RESULTS_DIR / "dsrs_brappa_stressor_rankings.csv")
dsrs_fam = pd.read_csv(RESULTS_DIR / "dsrs_brappa_family_rankings.csv")
redox_preds = pd.read_csv(RESULTS_DIR / "redox_decoder_brappa_predictions.csv")
redox_spat = pd.read_csv(RESULTS_DIR / "redox_decoder_organ_cell_scores.csv")
rri_scores = pd.read_csv(RESULTS_DIR / "brappa_rri_sample_scores.csv")
ddr_summary = pd.read_csv(RESULTS_DIR / "brappa_preservative_and_ddr_summary.csv")

# -------------------------------------------------------------------------
# Figure 1: Conceptual Architecture of Cross-Species Radiation Decoding
# -------------------------------------------------------------------------
print("Generating Figure 1: Conceptual Architecture...")
fig, ax = plt.subplots(figsize=(13.8, 9.6), dpi=300)
ax.axis("off")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

# Title & Subtitle Banner
fig.suptitle("Figure 1: Analysis Framework for B. rapa 40 cGy GCR Contrasts (draft under audit)",
             fontsize=15.5, weight="bold", y=0.978, color="#111827")

ax.text(0.5, 0.948, "Brassica rapa GCR contrasts scored with Arabidopsis-built tools through 18,770 one-to-one orthologs (31 of 39 libraries in the processed matrix)",
        ha="center", va="center", fontsize=9.2, style="italic", color="#4B5563")

# Top Cards
card_w = 0.268
card_h = 0.335
y_top = 0.585
gap_top = 0.065

top_cards = [
    {
        "id": "card1",
        "badge": "FACTORIAL COHORT", "badge_color": "#1D4ED8",
        "title": "1. B. rapa GCR Dataset",
        "subtitle": "Brassica rapa (Wisconsin Fast Plants)",
        "xy": (0.04, y_top), "w": card_w, "h": card_h,
        "bg": "#F8FAFC", "border": "#2563EB", "header_bg": "#EFF6FF",
        "items": [
            ("Radiation Dose", "0 cGy (Sham) vs 40 cGy GCR"),
            ("Genotype Contrast", "Wild-Type (WT) vs anthocyaninless (anl)"),
            ("Fixative", "DRS vs RNAlater (RL)"),
            ("Cohort Scope", "39 libraries; 31 in processed matrix")
        ]
    },
    {
        "id": "card2",
        "badge": "HOMOLOGY BRIDGE", "badge_color": "#15803D",
        "title": "2. Orthology Mapping",
        "subtitle": "DIAMOND Proteome Reciprocal Alignment",
        "xy": (0.04 + card_w + gap_top, y_top), "w": card_w, "h": card_h,
        "bg": "#F8FAF8", "border": "#16A34A", "header_bg": "#F0FDF4",
        "items": [
            ("1:1 RBH Mapping", "18,770 reciprocal best hits (Bra → AGI)"),
            ("Search", "DIAMOND blastp both ways, E ≤ 1e-5"),
            ("Filter", "Reciprocal best hit, identity ≥ 50%"),
            ("Projection", "B. rapa log2FC keyed by AGI code")
        ]
    },
    {
        "id": "card3",
        "badge": "SCORING TOOLS", "badge_color": "#B45309",
        "title": "3. Four Analyses",
        "subtitle": "Arabidopsis references and heuristics",
        "xy": (0.04 + 2 * (card_w + gap_top), y_top), "w": card_w, "h": card_h,
        "bg": "#FFFDF7", "border": "#D97706", "header_bg": "#FEF3C7",
        "items": [
            ("Tier 1: DSRS", "27 reference contrasts / 123 seed panels"),
            ("Tier 2: Redox heuristic", "Hand-assigned weights, 16 sensors"),
            ("Tier 3: RRI", "Referenced to the WT 0 cGy samples"),
            ("Tier 4: DNA repair", "Six orthologs, none significant")
        ]
    }
]

# Bottom Cards
y_bot = 0.035
h_bot = 0.465
w_bot = 0.44

bottom_cards = [
    {
        "id": "card4",
        "badge": "STRESS RECOGNITION: WHAT THE TABLES SHOW", "badge_color": "#1D4ED8",
        "title": "DSRS Projection",
        "subtitle": "Does an Arabidopsis-built DSRS recognise GCR in B. rapa?",
        "xy": (0.04, y_bot), "w": w_bot, "h": h_bot,
        "bg": "#F8FAFC", "border": "#2563EB", "header_bg": "#EFF6FF",
        "bullets": [
            ("WT 40 cGy Signature",
             "Closest matches are microgravity root (ρ = 0.259) and UV-B (ρ = 0.197); the three radiation families rank 9th-11th of 13."),
            ("Dose × Genotype Interaction",
             "Radiation ranks first: acute 60Co gamma (ρ = 0.254), gcr_80cGy (0.217); radiation_GCR family mean ρ = 0.150."),
            ("No Significance Test",
             "No null distribution or significance test was computed for these panel correlations."),
            ("Weak Dose Signal",
             "DESeq2: 7 WT genes at padj < 0.05; 0 in anl; 1 for dose × genotype; 0 for dose × fixative.")
        ]
    },
    {
        "id": "card5",
        "badge": "REMAINING ANALYSES: OPEN QUESTIONS", "badge_color": "#B91C1C",
        "title": "Redox, Resilience & Fixatives",
        "subtitle": "What these analyses can and cannot show",
        "xy": (0.52, y_bot), "w": w_bot, "h": h_bot,
        "bg": "#FEF2F2", "border": "#DC2626", "header_bg": "#FEE2E2",
        "bullets": [
            ("Redox Scores Are Heuristic",
             "Hand-assigned weights, not the CVAE: WT H₂O₂ 51.1%; anl paraquat 24.9% vs H₂O₂ 23.7%. No sensor gene is DE."),
            ("RRI Is Self-Referenced",
             "WT 0 cGy defines the index (0.762); WT 40 cGy (0.636) and anl at both doses (0.637) score lower by construction."),
            ("Fixative Comparison",
             "WT dose effects in DRS vs RNAlater correlate (r = 0.297, 30,756 genes); no gene has a significant dose × fixative term."),
            ("DNA Repair Orthologs",
             "BRCA1 +2.03 and PARP2 +0.86 log2FC in WT, but padj ≥ 0.72 and baseMean ≤ 12: not significant.")
        ]
    }
]

# Render Top Cards
for c in top_cards:
    x, y, w, h = c["xy"][0], c["xy"][1], c["w"], c["h"]
    
    # Outer box
    box = patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.018",
                                facecolor=c["bg"], edgecolor=c["border"], linewidth=1.6)
    ax.add_patch(box)
    
    # Header box
    header_h = 0.086
    header_box = patches.FancyBboxPatch((x, y + h - header_h), w, header_h,
                                       boxstyle="round,pad=0.012,rounding_size=0.018",
                                       facecolor=c["header_bg"], edgecolor="none")
    ax.add_patch(header_box)
    
    # Category badge
    ax.text(x + 0.015, y + h - 0.022, c["badge"],
            ha="left", va="center", fontsize=6.8, weight="bold", color=c["badge_color"])
    
    # Title
    ax.text(x + 0.015, y + h - 0.046, c["title"],
            ha="left", va="center", weight="bold", fontsize=10.2, color="#111827")
    
    # Subtitle
    ax.text(x + 0.015, y + h - 0.067, c["subtitle"],
            ha="left", va="center", fontsize=7.8, style="italic", color="#4B5563")
    
    # Divider line under header
    ax.plot([x, x + w], [y + h - header_h, y + h - header_h], color=c["border"], linewidth=0.8, alpha=0.4)
    
    # Items
    start_y = y + h - header_h - 0.028
    row_spacing = 0.056
    for i, (k, v) in enumerate(c["items"]):
        cy = start_y - i * row_spacing
        ax.text(x + 0.015, cy, f"• {k}:", fontsize=8.2, weight="bold", color="#1F2937")
        ax.text(x + 0.020, cy - 0.020, v, fontsize=7.6, color="#4B5563")

# Flow connectors: Card 1 -> Badge -> Card 2
y_mid = y_top + card_h / 2
x1_end = top_cards[0]["xy"][0] + card_w
x2_start = top_cards[1]["xy"][0]
mid1 = (x1_end + x2_start) / 2

# Draw Badge 1
t1 = ax.text(mid1, y_mid, "DIAMOND\n1:1 RBH", ha="center", va="center", fontsize=6.8, weight="bold", color="#15803D",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#F0FDF4", edgecolor="#16A34A", lw=0.9))
# Arrow in and out
ax.annotate("", xy=(mid1 - 0.022, y_mid), xytext=(x1_end + 0.003, y_mid),
            arrowprops=dict(arrowstyle="->,head_length=0.3,head_width=0.22", lw=1.6, color="#16A34A"))
ax.annotate("", xy=(x2_start - 0.003, y_mid), xytext=(mid1 + 0.022, y_mid),
            arrowprops=dict(arrowstyle="->,head_length=0.3,head_width=0.22", lw=1.6, color="#16A34A"))

# Flow connectors: Card 2 -> Badge -> Card 3
x2_end = top_cards[1]["xy"][0] + card_w
x3_start = top_cards[2]["xy"][0]
mid2 = (x2_end + x3_start) / 2

# Draw Badge 2
t2 = ax.text(mid2, y_mid, "AGI\nKEYED", ha="center", va="center", fontsize=6.8, weight="bold", color="#B45309",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#FEF3C7", edgecolor="#D97706", lw=0.9))
# Arrow in and out
ax.annotate("", xy=(mid2 - 0.022, y_mid), xytext=(x2_end + 0.003, y_mid),
            arrowprops=dict(arrowstyle="->,head_length=0.3,head_width=0.22", lw=1.6, color="#D97706"))
ax.annotate("", xy=(x3_start - 0.003, y_mid), xytext=(mid2 + 0.022, y_mid),
            arrowprops=dict(arrowstyle="->,head_length=0.3,head_width=0.22", lw=1.6, color="#D97706"))

# Central Junction Milestone
junction_x = 0.50
junction_y = 0.542
ax.plot([junction_x, junction_x], [y_top, junction_y + 0.016], color="#64748B", lw=1.6, linestyle="--")

# Milestone Pill
ax.text(junction_x, junction_y, "Cross-Species Projection",
        ha="center", va="center", fontsize=8.2, weight="bold", color="#1E293B",
        bbox=dict(boxstyle="round,pad=0.32", facecolor="#F8FAFC", edgecolor="#64748B", lw=1.0))

# Branching Pathway Lines (clean vertical drop, clean horizontal turn, clean vertical entry)
split_y = junction_y - 0.016
ax.plot([junction_x, junction_x], [junction_y - 0.016, split_y - 0.012], color="#64748B", lw=1.6)

# Left Branch to Card 4 (Diagnostic)
target_x_left = bottom_cards[0]["xy"][0] + bottom_cards[0]["w"] / 2
target_y_bot = y_bot + h_bot
ax.annotate("", xy=(target_x_left, target_y_bot), xytext=(junction_x, split_y - 0.012),
            arrowprops=dict(arrowstyle="->,head_length=0.4,head_width=0.3", lw=1.8, color="#2563EB",
                            connectionstyle="angle,angleA=180,angleB=-90,rad=12"))

# Right Branch to Card 5 (Discovery)
target_x_right = bottom_cards[1]["xy"][0] + bottom_cards[1]["w"] / 2
ax.annotate("", xy=(target_x_right, target_y_bot), xytext=(junction_x, split_y - 0.012),
            arrowprops=dict(arrowstyle="->,head_length=0.4,head_width=0.3", lw=1.8, color="#DC2626",
                            connectionstyle="angle,angleA=0,angleB=-90,rad=12"))

# Render Bottom Cards
for c in bottom_cards:
    x, y, w, h = c["xy"][0], c["xy"][1], c["w"], c["h"]
    
    # Outer box
    box = patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.018",
                                facecolor=c["bg"], edgecolor=c["border"], linewidth=1.6)
    ax.add_patch(box)
    
    # Header box
    header_h = 0.090
    header_box = patches.FancyBboxPatch((x, y + h - header_h), w, header_h,
                                       boxstyle="round,pad=0.012,rounding_size=0.018",
                                       facecolor=c["header_bg"], edgecolor="none")
    ax.add_patch(header_box)
    
    # Badge
    ax.text(x + 0.018, y + h - 0.022, c["badge"],
            ha="left", va="center", fontsize=7.2, weight="bold", color=c["badge_color"])
    
    # Title
    ax.text(x + 0.018, y + h - 0.046, c["title"],
            ha="left", va="center", weight="bold", fontsize=10.5, color="#111827")
    
    # Subtitle
    ax.text(x + 0.018, y + h - 0.068, c["subtitle"],
            ha="left", va="center", fontsize=8.0, style="italic", color="#4B5563")
    
    # Divider line under header
    ax.plot([x, x + w], [y + h - header_h, y + h - header_h], color=c["border"], linewidth=0.8, alpha=0.4)
    
    # Bullet points with generous vertical rhythm
    start_y = y + h - header_h - 0.034
    bullet_spacing = 0.086
    for i, (heading, desc) in enumerate(c["bullets"]):
        cy = start_y - i * bullet_spacing
        ax.text(x + 0.018, cy, f"▶  {heading}", fontsize=8.8, weight="bold", color=c["border"])
        
        wrapped_lines = textwrap.wrap(desc, width=64)
        for j, line in enumerate(wrapped_lines):
            ax.text(x + 0.032, cy - 0.019 - j * 0.018, line, fontsize=7.8, color="#374151")

fig.savefig(FIG_DIR / "fig1_conceptual_architecture.png", bbox_inches="tight")
fig.savefig(FIG_DIR / "fig1_conceptual_architecture.pdf", bbox_inches="tight")
plt.close(fig)

# -------------------------------------------------------------------------
# Figure 2: DSRS Stress Recognition Across Spaceflight Contrasts
# -------------------------------------------------------------------------
print("Generating Figure 2: DSRS Stress Recognition...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw={"width_ratios": [1.2, 1]})

inter_ranks = dsrs_ranks[dsrs_ranks["query"] == "Brapa_Interaction_Dose_x_Genotype"].head(10).copy()
colors = ["#D93025" if "rad" in s or "gcr" in s else "#1A73E8" for s in inter_ranks["stressor_class"]]
ax1.barh(inter_ranks["stressor"], inter_ranks["spearman_rho"], color=colors, edgecolor="black", alpha=0.85)
ax1.set_xlabel("Spearman Rank Correlation (ρ with Reference Stressor)")
ax1.set_title("A. DSRS Top Matches for Radiation × Genotype Interaction", weight="bold")
ax1.legend(handles=[patches.Patch(facecolor="#D93025", edgecolor="black", label="Radiation family"),
                    patches.Patch(facecolor="#1A73E8", edgecolor="black", label="Other family")], loc="lower right")
ax1.invert_yaxis()
ax1.grid(axis="x", linestyle="--", alpha=0.5)

fam_pivot = dsrs_fam.pivot(index="stressor_class", columns="query", values="spearman_rho")
col_map = {
    "Brapa_WT_GCR_40cGy": "WT GCR 40cGy",
    "Brapa_AnthoLess_GCR_40cGy": "AnthoLess GCR 40cGy",
    "Brapa_Interaction_Dose_x_Genotype": "Interaction (Dose × Genotype)",
    "Brapa_Interaction_Dose_x_Preservative": "Interaction (Dose × Preservative)"
}
fam_pivot = fam_pivot.rename(columns=col_map)
fam_pivot = fam_pivot.sort_values(by="Interaction (Dose × Genotype)", ascending=False)

sns.heatmap(fam_pivot, annot=True, fmt=".2f", cmap="vlag", center=0, ax=ax2, cbar_kws={"label": "Mean Spearman ρ"})
ax2.set_title("B. Stressor Family Fingerprint Concordance", weight="bold")
ax2.set_ylabel("Stressor Class / Family")

fig.suptitle("Figure 2: DSRS Projection of B. rapa GCR Signatures (Spearman ρ; no significance test)", fontsize=13, weight="bold", y=0.98)
plt.tight_layout()
fig.savefig(FIG_DIR / "fig2_dsrs_stress_recognition.png")
fig.savefig(FIG_DIR / "fig2_dsrs_stress_recognition.pdf")
plt.close(fig)

# -------------------------------------------------------------------------
# Figure 3: Redox CVAE Decoding & The Antioxidant Shield Deficit
# -------------------------------------------------------------------------
print("Generating Figure 3: Redox CVAE Decoding...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

stim_cols = [c for c in redox_preds.columns if c.startswith("prob_")]
stim_labels = [c.replace("prob_", "").replace("_", " ").title() for c in stim_cols]
pred_plot = redox_preds.set_index("contrast")[stim_cols].T
pred_plot.index = stim_labels

x = np.arange(len(stim_labels))
width = 0.2
ax1.bar(x - 1.5*width, pred_plot["WT_GCR_40cGy"] * 100, width, label="WT GCR 40cGy", color="#1E8E3E", edgecolor="black")
ax1.bar(x - 0.5*width, pred_plot["AnthoLess_GCR_40cGy"] * 100, width, label="AnthoLess GCR 40cGy", color="#D93025", edgecolor="black")
ax1.bar(x + 0.5*width, pred_plot["Interaction_Dose_x_Genotype"] * 100, width, label="Dose × Genotype", color="#F9AB00", edgecolor="black")
ax1.bar(x + 1.5*width, pred_plot["Interaction_Dose_x_Preservative"] * 100, width, label="Dose × Preservative", color="#1A73E8", edgecolor="black")

ax1.set_xticks(x)
ax1.set_xticklabels(stim_labels, rotation=25, ha="right")
ax1.set_ylabel("Softmax Score (%)")
ax1.set_title("A. Heuristic Stimulus Scores (hand-assigned weights)", weight="bold")
ax1.legend(loc="upper right")
ax1.grid(axis="y", linestyle="--", alpha=0.5)

contrasts_disp = ["WT GCR 40cGy", "AnthoLess GCR", "Dose × Genotype", "Dose × Preservative"]
hours = redox_preds["inferred_hours"].values
colors_time = ["#1E8E3E", "#D93025", "#F9AB00", "#1A73E8"]

bars = ax2.bar(contrasts_disp, hours, color=colors_time, edgecolor="black", width=0.55)
ax2.axhspan(0, 1.0, color="#E8F0FE", alpha=0.4, label="Immediate (<1h)")
ax2.axhspan(1.0, 4.0, color="#FEF7E0", alpha=0.4, label="Early Acute (1-4h)")
ax2.axhspan(4.0, 12.0, color="#FCE8E6", alpha=0.4, label="Mid Response (4-12h)")
ax2.set_ylabel("Weighted Assigned Time (Hours)")
ax2.set_title("B. Activation-Weighted Mean of Assigned Response Times", weight="bold")
ax2.set_ylim(0, 7.0)
for b in bars:
    y = b.get_height()
    ax2.text(b.get_x() + b.get_width()/2, y + 0.15, f"{y:.2f} h", ha="center", weight="bold")
ax2.legend(loc="upper left")
ax2.grid(axis="y", linestyle="--", alpha=0.5)

fig.suptitle("Figure 3: Hand-Weighted Redox Heuristic Scores (not Redox_decoder CVAE output)", fontsize=12, weight="bold", y=0.98)
plt.tight_layout()
fig.savefig(FIG_DIR / "fig3_redox_cvae_decoding.png")
fig.savefig(FIG_DIR / "fig3_redox_cvae_decoding.pdf")
plt.close(fig)

# -------------------------------------------------------------------------
# Figure 4: Radiation Resilience Index (RRI) & Housekeeping Module Collapse
# -------------------------------------------------------------------------
print("Generating Figure 4: Radiation Resilience Index (RRI)...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

sns.boxplot(data=rri_scores, x="genotype", y="RRI", hue="condition",
            palette={"0_cGy": "#8AB4F8", "40_cGy": "#EA4335"}, ax=ax1, width=0.5)
sns.stripplot(data=rri_scores, x="genotype", y="RRI", hue="condition",
              dodge=True, palette="dark:black", ax=ax1, size=6)
ax1.set_title("A. Composite Radiation Resilience Index (RRI)", weight="bold")
ax1.set_xlabel("Genotype")
ax1.set_ylabel("Resilience Index (RRI)")
ax1.set_ylim(0.4, 0.9)
ax1.grid(axis="y", linestyle="--", alpha=0.5)

handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles[:2], labels[:2], title="Dose")

sns.barplot(data=rri_scores, x="genotype", y="RRI_module", hue="condition",
            palette={"0_cGy": "#8AB4F8", "40_cGy": "#EA4335"}, ax=ax2, edgecolor="black", capsize=0.1)
ax2.set_title("B. Sustained Housekeeping Module Preservation", weight="bold")
ax2.set_xlabel("Genotype")
ax2.set_ylabel("Preservation Score (r)")
ax2.set_ylim(0, 0.7)
ax2.grid(axis="y", linestyle="--", alpha=0.5)
ax2.legend(title="Dose")

fig.suptitle("Figure 4: RRI by Genotype and Dose (WT 0 cGy is the reference group for every component)", fontsize=12, weight="bold", y=0.98)
plt.tight_layout()
fig.savefig(FIG_DIR / "fig4_rri_degradation.png")
fig.savefig(FIG_DIR / "fig4_rri_degradation.pdf")
plt.close(fig)

# -------------------------------------------------------------------------
# Figure 5: Preservative Parity (DRS vs RL) & Core DDR / Flavonoid Shunt
# -------------------------------------------------------------------------
print("Generating Figure 5: Preservative Parity & DDR Activation...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: Preservative comparison
sns.boxplot(data=rri_scores, x="preservative", y="RRI", hue="condition",
            palette={"0_cGy": "#8AB4F8", "40_cGy": "#EA4335"}, ax=ax1, width=0.5)
ax1.set_title("A. RRI by Fixative: DRS vs RNAlater", weight="bold")
ax1.set_xlabel("Preservative Fixative")
ax1.set_ylabel("Radiation Resilience Index (RRI)")
ax1.set_ylim(0.45, 0.9)
ax1.grid(axis="y", linestyle="--", alpha=0.5)
ax1.legend(title="Dose")

# Panel B: Core DDR & Phenylpropanoid induction
symbols = ddr_summary["symbol"].values
wt_vals = ddr_summary["WT_40cGy_log2FC"].values
anl_vals = ddr_summary["AnthoLess_40cGy_log2FC"].values

y = np.arange(len(symbols))
h = 0.35

ax2.barh(y - h/2, wt_vals, height=h, label="Wild-Type (40 cGy)", color="#1E8E3E", edgecolor="black")
ax2.barh(y + h/2, anl_vals, height=h, label="Anthocyaninless (40 cGy)", color="#D93025", edgecolor="black")

ax2.set_yticks(y)
ax2.set_yticklabels(symbols)
ax2.axvline(0, color="black", linestyle="-", linewidth=1)
ax2.set_xlabel("Log2 Fold-Change (40 cGy vs 0 cGy)")
ax2.set_title("B. DNA Repair & Flavonoid Orthologs (none significant)", weight="bold")
ax2.legend(loc="lower right")
ax2.grid(axis="x", linestyle="--", alpha=0.5)

fig.suptitle("Figure 5: Fixative Comparison and DNA-Repair / Flavonoid Ortholog Fold Changes", fontsize=12, weight="bold", y=0.98)
plt.tight_layout()
fig.savefig(FIG_DIR / "fig5_preservative_and_ddr.png")
fig.savefig(FIG_DIR / "fig5_preservative_and_ddr.pdf")
plt.close(fig)

# -------------------------------------------------------------------------
# Figure 6: Spatial Anatomical Vulnerability (ggPlantMap)
# -------------------------------------------------------------------------
print("Generating Figure 6: Spatial Anatomical Vulnerability...")
fig, ax = plt.subplots(figsize=(10, 5))

spat_pivot = redox_spat.pivot(index="cell_type", columns="contrast", values="vulnerability_score")
spat_pivot = spat_pivot.rename(columns={
    "WT_GCR_40cGy": "WT GCR 40cGy",
    "AnthoLess_GCR_40cGy": "AnthoLess GCR",
    "Interaction_Dose_x_Genotype": "Dose × Genotype",
    "Interaction_Dose_x_Preservative": "Dose × Preservative"
})
spat_pivot = spat_pivot.sort_values(by="Dose × Genotype", ascending=False)

sns.heatmap(spat_pivot * 100, annot=True, fmt=".1f", cmap="YlOrRd", ax=ax, cbar_kws={"label": "Weighted Share (%)"})
ax.set_title("Figure 6: Cell-Type Shares from Hand-Assigned Sensor Weights (not data-derived)", weight="bold")
ax.set_ylabel("Cell-Type Label")

plt.tight_layout()
fig.savefig(FIG_DIR / "fig6_spatial_vulnerability.png")
fig.savefig(FIG_DIR / "fig6_spatial_vulnerability.pdf")
plt.close(fig)

print("=" * 70)
print("20_generate_figures_cross_species: All 6 figures generated successfully!")
print("=" * 70)
