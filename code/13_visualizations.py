#!/usr/bin/env python3
"""
13_visualizations.py
====================
Generate 7 publication-quality visualizations for the radiation-quality
signaling analysis and kinetic narrative.

Figures:
  1. LR pair specificity heatmap (top 50 quality-specific pairs at 2-6h)
  2. Pathway kinetic line plot (8 pathways over time)
  3. RRI dose-response curve (dose vs RRI by quality)
  4. Chord diagram for cell-type signaling (top 30 flows at 2-6h)
  5. Module dynamics streamgraph (WGCNA module eigengenes over time)
  6. RRI component radar per quality (5 radar charts)
  7. 2D dose-time RRI heatmap (gamma, GCR, HZE-Fe surfaces)

All figures saved as both .svg and .png to results/figures/.
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
rc_params = matplotlib.rcParams
from matplotlib.patches import FancyArrowPatch, Arc
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.gridspec as gridspec
from scipy.interpolate import griddata

# ── Global style ──────────────────────────────────────────────────────────────
rc_params['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
rc_params['svg.fonttype'] = 'none'  # Keep SVG text editable
rc_params['pdf.fonttype'] = 42
rc_params['axes.unicode_minus'] = False
rc_params['figure.dpi'] = 150
rc_params['savefig.dpi'] = 300
rc_params['axes.spines.top'] = False
rc_params['axes.spines.right'] = False

# Phylo color palette
PHYLO_COLORS = {
    'black': '#000000',
    'cream': '#ECE9E2',
    'white': '#FAF9F3',
    'yellow': '#E9ED4C',
    'orange': '#FF9400',
    'green': '#75A025',
    'pink': '#FD9BED',
    'blue': '#0279EE',
}

# Quality colors
QUALITY_COLORS = {
    'gamma': '#FF9400',
    'GCR': '#0279EE',
    'HZE-Fe': '#FD9BED',
    'UV-B': '#75A025',
    'spaceflight-LEO': '#E9ED4C',
}

# Pathway colors
PATHWAY_COLORS = {
    'DNA repair': '#0279EE',
    'Oxidative stress': '#FF9400',
    'SA signaling': '#FD9BED',
    'JA signaling': '#75A025',
    'ETH signaling': '#E9ED4C',
    'ABA signaling': '#8B4513',
    'Auxin signaling': '#9370DB',
    'Cell cycle / DDR checkpoint': '#DC143C',
}

BASE = '/mnt/results/zenodo_bundle'
FIG_DIR = os.path.join(BASE, 'results', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

# Time bin ordering
TIME_ORDER = ['0-0.5h', '0.5-2h', '2-6h', '6-12h', '12-30h', '30h+']
TIME_LABELS = ['0–0.5h', '0.5–2h', '2–6h', '6–12h', '12–30h', '30h+']


def save_fig(fig, name):
    """Save figure as both SVG and PNG."""
    svg_path = os.path.join(FIG_DIR, f'{name}.svg')
    png_path = os.path.join(FIG_DIR, f'{name}.png')
    fig.savefig(svg_path, format='svg', bbox_inches='tight', facecolor='white')
    fig.savefig(png_path, format='png', bbox_inches='tight', facecolor='white', dpi=300)
    plt.close(fig)
    print(f'  Saved: {name}.svg, {name}.png')


def load_tair_symbol_map():
    """Load TAIR -> gene symbol mapping from microarray annotation."""
    ma_path = '/mnt/shared-workspace/shared/raw/counts/OSD-46/GLDS-46_array_normalized_expression_probeset_GLmicroarray.csv'
    try:
        ma = pd.read_csv(ma_path, usecols=['TAIR', 'SYMBOL'])
        ma['sym'] = ma['SYMBOL'].str.split('|').str[0]
        return ma.drop_duplicates('TAIR').set_index('TAIR')['sym'].to_dict()
    except Exception:
        return {}


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 1: LR Pair Specificity Heatmap (top 50 quality-specific at 2-6h)
# ═══════════════════════════════════════════════════════════════════════════════
def fig1_lr_specificity_heatmap():
    print('Figure 1: LR pair specificity heatmap')
    spec = pd.read_csv(os.path.join(BASE, 'results/cellchat/lr_pair_specificity.csv'))
    sym_map = load_tair_symbol_map()

    # Focus on 2-6h — the most reliable cross-quality comparison (GCR vs gamma, both RNA-seq)
    spec26 = spec[spec['timepoint'] == '2-6h'].copy()
    qs = spec26[spec26['specificity'] == 'quality-specific'].sort_values('max_min_ratio', ascending=False).head(50)

    # Parse ligand/receptor and get symbols
    qs['ligand'] = qs['LR_pair'].str.split('->').str[0]
    qs['receptor'] = qs['LR_pair'].str.split('->').str[1]
    qs['ligand_sym'] = qs['ligand'].map(sym_map).fillna(qs['ligand'])
    qs['receptor_sym'] = qs['receptor'].map(sym_map).fillna(qs['receptor'])
    qs['label'] = qs['ligand_sym'] + ' → ' + qs['receptor_sym']

    # Build matrix: rows = LR pairs, cols = qualities
    quality_cols = ['mean_GCR', 'mean_gamma', 'mean_HZE-Fe', 'mean_UV-B']
    quality_labels = ['GCR', 'gamma', 'HZE-Fe', 'UV-B']

    mat = qs[quality_cols].values.astype(float)
    # Replace NaN with 0 for visualization
    mat = np.nan_to_num(mat, nan=0.0)

    # Z-score per row
    row_mean = mat.mean(axis=1, keepdims=True)
    row_std = mat.std(axis=1, keepdims=True)
    row_std[row_std == 0] = 1
    mat_z = (mat - row_mean) / row_std

    fig, ax = plt.subplots(figsize=(8, 12))
    cmap = LinearSegmentedColormap.from_list('phylo_div',
        ['#0279EE', '#FAF9F3', '#FF9400'], N=256)
    im = ax.imshow(mat_z, aspect='auto', cmap=cmap, vmin=-2, vmax=2,
                   interpolation='nearest')

    ax.set_xticks(range(len(quality_labels)))
    ax.set_xticklabels(quality_labels, fontsize=11, fontweight='bold')
    ax.set_yticks(range(len(qs)))
    ax.set_yticklabels(qs['label'].values, fontsize=7, fontfamily='monospace')
    ax.tick_params(axis='x', top=True, bottom=False, labeltop=True, labelbottom=False)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=0.4, pad=0.02)
    cbar.set_label('Z-scored signaling strength', fontsize=10)

    ax.set_title('Top 50 Quality-Specific LR Pairs (2–6h)\nGCR-enriched, same-platform RNA-seq comparison',
                 fontsize=13, fontweight='bold', pad=15)
    ax.set_xlabel('')
    ax.set_ylabel('')

    # Add grid
    ax.set_xticks(np.arange(-0.5, len(quality_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(qs), 0.5), minor=True)
    ax.grid(which='minor', color='white', linewidth=0.5)
    ax.tick_params(which='minor', size=0)

    save_fig(fig, 'lr_pair_specificity_heatmap')


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 2: Pathway Kinetic Line Plot (8 pathways over time)
# ═══════════════════════════════════════════════════════════════════════════════
def fig2_pathway_kinetics():
    print('Figure 2: Pathway kinetic line plot')
    pt = pd.read_csv(os.path.join(BASE, 'results/pathway_enrichment/pathway_by_timepoint.csv'))

    # Order time bins
    pt['time_bin'] = pd.Categorical(pt['time_bin'], categories=TIME_ORDER[:5], ordered=True)
    pt = pt.sort_values('time_bin')

    pathways = ['DNA repair', 'Oxidative stress', 'SA signaling', 'JA signaling',
                'ETH signaling', 'ABA signaling', 'Auxin signaling', 'Cell cycle / DDR checkpoint']
    x = np.arange(len(pt))

    fig, ax = plt.subplots(figsize=(10, 7))

    for pw in pathways:
        ax.plot(x, pt[pw].values, 'o-', linewidth=2.2, markersize=7,
                color=PATHWAY_COLORS[pw], label=pw, alpha=0.85)

    ax.axhline(0, color='#CCCCCC', linewidth=0.8, linestyle='--')
    ax.set_xticks(x)
    ax.set_xticklabels([TIME_LABELS[i] for i in range(len(pt))], fontsize=11)
    ax.set_xlabel('Time post-exposure', fontsize=12)
    ax.set_ylabel('Pathway enrichment score (z-scored GSVA)', fontsize=12)
    ax.set_title('Radiation-Response Pathway Kinetics Over Time\n(gamma RNA-seq, n=72 irradiated samples)',
                 fontsize=13, fontweight='bold', pad=12)
    ax.legend(loc='upper right', fontsize=8.5, framealpha=0.9, ncol=1)
    ax.set_ylim(-0.55, 0.55)
    ax.set_xlim(-0.3, len(pt) - 0.7)

    # Highlight nadir region
    ax.axvspan(1.8, 2.2, alpha=0.08, color='#FF9400')
    ax.annotate('RRI nadir\n(2–6h)', xy=(2, -0.5), fontsize=8, ha='center',
                color='#FF9400', fontweight='bold')

    save_fig(fig, 'pathway_kinetic_lines')


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 3: RRI Dose-Response Curve
# ═══════════════════════════════════════════════════════════════════════════════
def fig3_rri_dose_response():
    print('Figure 3: RRI dose-response curve')
    r = pd.read_csv(os.path.join(BASE, 'results/rri/rri_per_sample.csv'))

    # Non-control samples with dose and quality
    rc = r[(r['IsControl'] == False) & r['RadiationQuality'].notna() & r['AbsorbedDose_Gy'].notna()].copy()

    fig, ax = plt.subplots(figsize=(9, 7))

    # Plot per quality
    for q in ['gamma', 'GCR', 'HZE-Fe', 'UV-B']:
        sub = rc[rc['RadiationQuality'] == q]
        if len(sub) == 0:
            continue
        # Jitter dose slightly for visibility
        jitter = np.random.normal(0, 0.3, len(sub)) if q == 'gamma' else np.random.normal(0, 0.01, len(sub))
        ax.scatter(sub['AbsorbedDose_Gy'] + jitter, sub['RRI'],
                   color=QUALITY_COLORS[q], s=50, alpha=0.7, edgecolors='white',
                   linewidth=0.5, label=f'{q} (n={len(sub)})', zorder=3)

        # Mean ± SEM
        if len(sub) > 1:
            mean_rri = sub['RRI'].mean()
            sem_rri = sub['RRI'].std() / np.sqrt(len(sub))
            ax.errorbar(sub['AbsorbedDose_Gy'].mean(), mean_rri,
                        yerr=sem_rri, fmt='D', color=QUALITY_COLORS[q],
                        markersize=10, markeredgecolor='black', markeredgewidth=1,
                        capsize=5, capthick=1.5, zorder=4)

    # Controls
    ctrl = r[r['IsControl'] == True]
    if len(ctrl) > 0:
        ax.axhline(ctrl['RRI'].mean(), color='#888888', linestyle=':', linewidth=1.5,
                   label=f'Control mean ({ctrl["RRI"].mean():.3f})', zorder=2)
        ax.fill_between([0, 105], ctrl['RRI'].mean() - ctrl['RRI'].std(),
                        ctrl['RRI'].mean() + ctrl['RRI'].std(), alpha=0.1, color='#888888')

    ax.set_xlabel('Absorbed Dose (Gy)', fontsize=12)
    ax.set_ylabel('Radiation Response Index (RRI)', fontsize=12)
    ax.set_title('RRI Dose-Response by Radiation Quality', fontsize=13, fontweight='bold', pad=12)
    ax.legend(loc='lower right', fontsize=9, framealpha=0.9)
    ax.set_xlim(-2, 105)
    ax.set_ylim(0.3, 1.0)

    # Annotate GCR low-dose cluster
    gcr = rc[rc['RadiationQuality'] == 'GCR']
    if len(gcr) > 0:
        ax.annotate(f'GCR: low dose,\nhigh LET\n(RRI={gcr["RRI"].mean():.3f})',
                    xy=(gcr['AbsorbedDose_Gy'].mean(), gcr['RRI'].mean()),
                    xytext=(20, 0.45), fontsize=8, color=QUALITY_COLORS['GCR'],
                    arrowprops=dict(arrowstyle='->', color=QUALITY_COLORS['GCR'], lw=1.2))

    save_fig(fig, 'rri_dose_response')


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 4: Chord Diagram for Cell-Type Signaling (top 30 flows at 2-6h)
# ═══════════════════════════════════════════════════════════════════════════════
def fig4_chord_diagram():
    print('Figure 4: Chord diagram for cell-type signaling')
    sf = pd.read_csv(os.path.join(BASE, 'results/cellchat/signaling_flow_per_quality_timepoint.csv'))

    # Top 30 flows at 2-6h across all qualities
    sf26 = sf[sf['Time'] == '2-6h'].copy()
    top30 = sf26.nlargest(30, 'SignalStrength').copy()

    # Get unique cell types
    all_cts = sorted(set(top30['Source'].unique()) | set(top30['Target'].unique()))
    n_ct = len(all_cts)
    ct_idx = {ct: i for i, ct in enumerate(all_cts)}

    # Aggregate flows between same source-target pairs
    flow_mat = np.zeros((n_ct, n_ct))
    for _, row in top30.iterrows():
        i, j = ct_idx[row['Source']], ct_idx[row['Target']]
        flow_mat[i, j] += row['SignalStrength']

    # Layout: circular
    angles = np.linspace(0, 2 * np.pi, n_ct, endpoint=False) + np.pi / 2
    radius = 1.0

    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': None})
    ax.set_aspect('equal')
    ax.axis('off')

    # Draw arcs for each cell type
    arc_width = 2 * np.pi / n_ct * 0.8
    max_flow = flow_mat.sum()

    # Compute arc lengths proportional to total flow
    out_flow = flow_mat.sum(axis=1)
    in_flow = flow_mat.sum(axis=0)
    total_flow = out_flow + in_flow
    total_flow[total_flow == 0] = 1

    for i, ct in enumerate(all_cts):
        a = angles[i]
        # Arc segment
        seg_width = 2 * np.pi / n_ct * 0.7
        theta_start = a - seg_width / 2
        theta_end = a + seg_width / 2

        # Outer arc
        theta = np.linspace(theta_start, theta_end, 30)
        x_outer = radius * np.cos(theta)
        y_outer = radius * np.sin(theta)
        x_inner = 0.92 * radius * np.cos(theta[::-1])
        y_inner = 0.92 * radius * np.sin(theta[::-1])

        color = plt.cm.tab20(i / n_ct)
        ax.fill(np.concatenate([x_outer, x_inner]),
                np.concatenate([y_outer, y_inner]), color=color, alpha=0.8)

        # Label
        label_r = 1.15 * radius
        ha = 'left' if np.cos(a) > 0.1 else 'right' if np.cos(a) < -0.1 else 'center'
        va = 'bottom' if np.sin(a) > 0.1 else 'top' if np.sin(a) < -0.1 else 'center'
        label = ct.replace('_', ' ').replace('hormone response ', 'HR-').replace('dna damage response', 'DDR')
        ax.text(label_r * np.cos(a), label_r * np.sin(a), label,
                fontsize=7, ha=ha, va=va, fontweight='bold', color=color)

    # Draw chords (bezier curves)
    for i in range(n_ct):
        for j in range(n_ct):
            if flow_mat[i, j] > 0:
                a_i, a_j = angles[i], angles[j]
                # Control points through center
                p0 = np.array([radius * np.cos(a_i), radius * np.sin(a_i)])
                p1 = np.array([0.0, 0.0])
                p2 = np.array([radius * np.cos(a_j), radius * np.sin(a_j)])

                # Quadratic bezier
                t = np.linspace(0, 1, 50)
                curve_x = (1 - t)**2 * p0[0] + 2 * (1 - t) * t * p1[0] + t**2 * p2[0]
                curve_y = (1 - t)**2 * p0[1] + 2 * (1 - t) * t * p1[1] + t**2 * p2[1]

                alpha = min(0.6, flow_mat[i, j] / max_flow * 3)
                lw = max(0.5, np.log10(flow_mat[i, j] + 1) * 0.8)
                color_i = plt.cm.tab20(i / n_ct)
                ax.plot(curve_x, curve_y, color=color_i, alpha=alpha, linewidth=lw, zorder=1)

    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)
    ax.set_title('Cell-Type Signaling Network at 2–6h\n(Top 30 flows across all radiation qualities)',
                 fontsize=13, fontweight='bold', pad=20)

    save_fig(fig, 'signaling_chord_diagram')


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 5: Module Dynamics Streamgraph
# ═══════════════════════════════════════════════════════════════════════════════
def fig5_module_streamgraph():
    print('Figure 5: Module dynamics streamgraph')
    mebt = pd.read_csv(os.path.join(BASE, 'results/pathway_enrichment/module_eigengene_by_timepoint.csv'))

    # Order time bins
    mebt['time_bin'] = pd.Categorical(mebt['time_bin'], categories=TIME_ORDER[:5], ordered=True)
    mebt = mebt.sort_values('time_bin')

    modules = ['MEblue', 'MEturquoise', 'MEgrey']
    module_colors = {'MEblue': '#0279EE', 'MEturquoise': '#FF9400', 'MEgrey': '#888888'}
    module_labels = {'MEblue': 'Blue module', 'MEturquoise': 'Turquoise module', 'MEgrey': 'Grey module'}

    x = np.arange(len(mebt))
    # For streamgraph, use absolute values stacked
    vals = np.array([mebt[m].abs().values for m in modules]).T

    fig, ax = plt.subplots(figsize=(10, 6))

    # Stackplot
    ax.stackplot(x, vals.T,
                 colors=[module_colors[m] for m in modules],
                 labels=[module_labels[m] for m in modules],
                 alpha=0.7, edgecolor='white', linewidth=1.5)

    # Overlay actual eigengene values as lines
    ax2 = ax.twinx()
    ax2.spines['top'].set_visible(False)
    for m in modules:
        ax2.plot(x, mebt[m].values, 'o-', color=module_colors[m],
                 linewidth=2, markersize=6, alpha=0.9, label=f'{module_labels[m]} (signed)')

    ax.set_xticks(x)
    ax.set_xticklabels([TIME_LABELS[i] for i in range(len(mebt))], fontsize=11)
    ax.set_xlabel('Time post-exposure', fontsize=12)
    ax.set_ylabel('|Module eigengene| (stacked)', fontsize=12)
    ax2.set_ylabel('Signed eigengene value', fontsize=12, color='#555555')
    ax.set_title('WGCNA Module Dynamics Over Time\n(Blue module: time-correlated, ρ = −0.66, p = 3.5×10⁻¹⁰)',
                 fontsize=13, fontweight='bold', pad=12)
    ax.legend(loc='upper left', fontsize=9, framealpha=0.9)
    ax2.legend(loc='upper right', fontsize=8, framealpha=0.9)
    ax.set_xlim(-0.3, len(mebt) - 0.7)

    save_fig(fig, 'module_streamgraph')


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 6: RRI Component Radar Per Quality
# ═══════════════════════════════════════════════════════════════════════════════
def fig6_rri_radar():
    print('Figure 6: RRI component radar per quality')
    r = pd.read_csv(os.path.join(BASE, 'results/rri/rri_per_sample.csv'))
    rc = r[(r['IsControl'] == False) & r['RadiationQuality'].notna()].copy()

    components = ['RRI_latent', 'RRI_pathway', 'RRI_module']
    comp_labels = ['Latent\n(GP-AE)', 'Pathway\n(GSVA)', 'Module\n(WGCNA)']

    # Compute mean per quality
    qualities = ['gamma', 'GCR', 'spaceflight-LEO']
    quality_data = {}
    for q in qualities:
        sub = rc[rc['RadiationQuality'] == q]
        if len(sub) > 0:
            quality_data[q] = [sub[c].mean() for c in components]

    # Also add control
    ctrl = r[r['IsControl'] == True]
    if len(ctrl) > 0:
        quality_data['Control'] = [ctrl[c].mean() for c in components]

    n_comp = len(components)
    angles = np.linspace(0, 2 * np.pi, n_comp, endpoint=False).tolist()
    angles += angles[:1]  # Close the polygon

    fig, axes = plt.subplots(1, 4, figsize=(16, 5), subplot_kw=dict(polar=True))

    plot_order = ['Control', 'gamma', 'GCR', 'spaceflight-LEO']
    for idx, q in enumerate(plot_order):
        ax = axes[idx]
        if q not in quality_data:
            ax.axis('off')
            continue

        values = quality_data[q]
        values_closed = values + values[:1]

        color = QUALITY_COLORS.get(q, '#888888')
        ax.fill(angles, values_closed, alpha=0.25, color=color)
        ax.plot(angles, values_closed, 'o-', linewidth=2, color=color, markersize=6)

        # Add reference: control outline
        if q != 'Control' and 'Control' in quality_data:
            ctrl_vals = quality_data['Control'] + quality_data['Control'][:1]
            ax.plot(angles, ctrl_vals, '--', color='#888888', linewidth=1, alpha=0.5)

        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(comp_labels, fontsize=9)
        ax.set_ylim(0, 1.05)
        ax.set_title(f'{q}\n(n={len(rc[rc["RadiationQuality"]==q]) if q != "Control" else len(ctrl)})',
                     fontsize=11, fontweight='bold', pad=15)
        ax.grid(True, alpha=0.3)

    fig.suptitle('RRI Component Decomposition by Radiation Quality',
                 fontsize=14, fontweight='bold', y=1.02)
    fig.tight_layout()
    save_fig(fig, 'rri_radar_per_quality')


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 7: 2D Dose-Time RRI Heatmap
# ═══════════════════════════════════════════════════════════════════════════════
def fig7_dose_time_heatmap():
    print('Figure 7: 2D dose-time RRI heatmap')
    surfaces = {}
    for q in ['gamma', 'gcr', 'hze_fe']:
        path = os.path.join(BASE, 'results/rri', f'rri_surface_{q}.csv')
        if os.path.exists(path):
            surfaces[q] = pd.read_csv(path)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    titles = {'gamma': 'Gamma (LET ~0.2 keV/μm)', 'gcr': 'GCR (LET ~1.5 keV/μm)', 'hze_fe': 'HZE-Fe (LET 175 keV/μm)'}
    cmaps = {'gamma': 'Oranges', 'gcr': 'Blues', 'hze_fe': 'RdPu'}

    for idx, (q, df) in enumerate(surfaces.items()):
        ax = axes[idx]
        # Pivot to dose x time grid
        pivot = df.pivot_table(index='dose_Gy', columns='time_h', values='RRI', aggfunc='mean')

        # Interpolate to regular grid for smooth heatmap
        dose_unique = np.linspace(df['dose_Gy'].min(), df['dose_Gy'].max(), 50)
        time_unique = np.linspace(df['time_h'].min(), df['time_h'].max(), 50)
        D, T = np.meshgrid(dose_unique, time_unique)

        points = df[['dose_Gy', 'time_h']].values
        values = df['RRI'].values
        Z = griddata(points, values, (D, T), method='linear')

        im = ax.pcolormesh(D, T, Z.T, cmap=cmaps[q], shading='auto',
                           vmin=0.3, vmax=1.0)
        ax.set_xlabel('Dose (Gy)', fontsize=11)
        ax.set_ylabel('Time post-exposure (h)', fontsize=11)
        ax.set_title(titles[q], fontsize=12, fontweight='bold')

        # Add contour lines
        contour_levels = [0.5, 0.6, 0.7, 0.8]
        cs = ax.contour(D, T, Z.T, levels=contour_levels, colors='black',
                        linewidths=0.8, alpha=0.5)
        ax.clabel(cs, fontsize=7, fmt='%.2f')

    # Shared colorbar
    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.94, 0.15, 0.015, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label('RRI', fontsize=11)

    fig.suptitle('Radiation Response Index: Dose–Time Response Surfaces',
                 fontsize=14, fontweight='bold', y=1.02)
    fig.tight_layout(rect=[0, 0, 0.93, 0.98])
    save_fig(fig, 'rri_dose_time_heatmap')


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print('=' * 70)
    print('Generating 7 publication-quality visualizations')
    print('=' * 70)

    fig1_lr_specificity_heatmap()
    fig2_pathway_kinetics()
    fig3_rri_dose_response()
    fig4_chord_diagram()
    fig5_module_streamgraph()
    fig6_rri_radar()
    fig7_dose_time_heatmap()

    print('=' * 70)
    print('All 7 figures generated successfully.')
    print(f'Output directory: {FIG_DIR}')
    print('=' * 70)
