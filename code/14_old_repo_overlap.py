#!/usr/bin/env python3
"""
14_old_repo_overlap.py

Cross-pipeline consistency check between this kinetic pipeline and the earlier
per-contrast DESeq2 / WGCNA analysis in the dr-richard-barker/Plant_response_to_radiation
GitHub repo (which was run on largely the same OSDR RNA-seq studies:
OSD-498/502/508/510/658).

Produces a two-panel figure:
  A) Enrichment of old-repo radiation DEGs within each new WGCNA module.
  B) Provenance of the canonical DNA-damage-response (DDR) genes: which new
     module each lands in (correcting the manuscript's blue-hub gene list).

Inputs (paths via CLI or defaults):
  OLD_DEG   : old repo data/processed/rnaseq/Up_and_down_genes_DESeq2.csv
  NEW_MOD   : this bundle results/wgcna/modules.csv
Outputs:
  results/figures/old_repo_overlap.{svg,png}
"""
import csv, sys, os
from math import comb

# ---------------------------------------------------------------- paths
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # zenodo_bundle
NEW_MOD = sys.argv[1] if len(sys.argv) > 1 else f"{HERE}/results/wgcna/modules.csv"
# Bundled copy of the earlier repo's DESeq2 Up/Down calls (dr-richard-barker/
# Plant_response_to_radiation, data/processed/rnaseq/Up_and_down_genes_DESeq2.csv),
# included here so this check is reproducible from the Zenodo bundle alone.
OLD_DEG = sys.argv[2] if len(sys.argv) > 2 else f"{HERE}/data/old_repo_deg_table.csv"
OUTDIR = f"{HERE}/results/figures"
os.makedirs(OUTDIR, exist_ok=True)

# ---------------------------------------------------------------- load new modules
rows = list(csv.DictReader(open(NEW_MOD)))
new_mod = {r["Gene"].upper(): r["Module"] for r in rows}
universe = set(new_mod)                                  # 2000 HVGs = background
mod_genes = {m: {g for g, mm in new_mod.items() if mm == m}
             for m in ("blue", "turquoise", "grey")}

# ---------------------------------------------------------------- load old DEGs
deg_rows = list(csv.DictReader(open(OLD_DEG)))
contrasts = [c for c in deg_rows[0] if c != "gene_id"]
old_deg = {r["gene_id"].upper() for r in deg_rows
           if any(r[c] in ("Up", "Down") for c in contrasts)}
deg_in_univ = old_deg & universe                         # 398 genes

# ---------------------------------------------------------------- stats
def hyper_ge(k, K, n, N):
    return sum(comb(K, i) * comb(N - K, n - i)
               for i in range(k, min(K, n) + 1)) / comb(N, N and n)

N = len(universe)
n = len(deg_in_univ)
bg = n / N
stats = {}
for m, gs in mod_genes.items():
    k = len(deg_in_univ & gs)
    K = len(gs)
    exp = K * n / N
    p = hyper_ge(k, K, n, N)
    stats[m] = dict(k=k, K=K, frac=k / K, exp_frac=exp / K, fold=(k / exp) if exp else 0, p=p)

# ---------------------------------------------------------------- DDR gene provenance
DDR = [("ATM", "AT3G48190"), ("ATR", "AT5G40820"), ("SOG1", "AT1G25580"),
       ("MYB3R1", "AT4G32730"), ("MYB3R3", "AT3G09370"),
       ("BRCA1", "AT4G21070"), ("BRCA2A", "AT4G00020"),
       ("PARP2", "AT4G02390"), ("RAD51", "AT5G20850"), ("RAD54", "AT3G19210"),
       ("SMR7", "AT3G27630"), ("GMI1", "AT5G24280"), ("XRI1", "AT5G48720"),
       ("KU70", "AT1G16970"), ("WEE1", "AT1G02970")]

# ---------------------------------------------------------------- plot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none",
                     "font.size": 10, "axes.edgecolor": "#444"})
COL = {"blue": "#3B6EA5", "turquoise": "#3FB6A8", "grey": "#8A8A8A",
       "absent": "#D9534F"}

fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 5.2),
                               gridspec_kw={"width_ratios": [1.05, 1.0]})

# --- Panel A: DEG enrichment per module ---
order = ["blue", "turquoise", "grey"]
labels = ["blue\n(early-response)", "turquoise\n(sustained)", "grey\n(DDR core)"]
vals = [stats[m]["frac"] * 100 for m in order]
bars = axA.bar(range(3), vals, color=[COL[m] for m in order],
               edgecolor="#222", width=0.62, zorder=3)
axA.axhline(bg * 100, ls="--", lw=1.4, color="#333", zorder=2)
axA.text(2.46, bg * 100 + 1.5, f"background\n{bg*100:.1f}%", ha="right",
         va="bottom", fontsize=8.5, color="#333", style="italic")
for i, m in enumerate(order):
    s = stats[m]
    axA.text(i, vals[i] + 1.6, f"{s['k']}/{s['K']}", ha="center",
             va="bottom", fontsize=9, fontweight="bold")
# significance annotation for blue
sb = stats["blue"]
axA.annotate(f"{sb['fold']:.2f}× enrichment\nP = {sb['p']:.1e}",
             xy=(0, vals[0]), xytext=(0.15, 66),
             fontsize=9, color=COL["blue"], fontweight="bold",
             ha="left", va="top")
axA.set_xticks(range(3))
axA.set_xticklabels(labels)
axA.set_ylabel("% of module genes that are\nold-repo radiation DEGs")
axA.set_ylim(0, 100)
axA.set_title("A  Old DESeq2 DEGs are enriched in the\nnew early-response & DDR modules",
              fontsize=10.5, loc="left", fontweight="bold")
axA.grid(axis="y", ls=":", color="#ccc", zorder=0)
for s in ("top", "right"):
    axA.spines[s].set_visible(False)

# --- Panel B: DDR gene provenance ---
def provenance(agi):
    a = agi.upper()
    if a not in universe:
        return "absent"
    return new_mod[a]  # blue / turquoise / grey

yB = list(range(len(DDR)))[::-1]
for y, (sym, agi) in zip(yB, DDR):
    prov = provenance(agi)
    is_deg = agi.upper() in old_deg
    axB.scatter([0.5], [y], s=340, color=COL.get(prov, "#999"),
                edgecolor="#222", zorder=3)
    tag = {"absent": "not in HVG panel"}.get(prov, f"{prov} module")
    axB.text(0.62, y, f"{sym}", va="center", ha="left",
             fontsize=9.5, fontweight="bold")
    axB.text(1.35, y, tag, va="center", ha="left", fontsize=8.6,
             color=COL.get(prov, "#666"))
    axB.text(2.75, y, "DEG" if is_deg else "–", va="center", ha="center",
             fontsize=8.6, color="#222" if is_deg else "#bbb",
             fontweight="bold" if is_deg else "normal")
axB.text(2.75, len(DDR) - 0.35, "old DEG?", ha="center", fontsize=8.2,
         style="italic", color="#555")
axB.set_xlim(0.2, 3.1)
axB.set_ylim(-0.8, len(DDR) - 0.2)
axB.axis("off")
axB.set_title("B  Canonical DDR genes land in grey (or are\nabsent) — not among blue hubs",
              fontsize=10.5, loc="left", fontweight="bold")
leg = [Patch(fc=COL["grey"], ec="#222", label="grey module (DDR core)"),
       Patch(fc=COL["blue"], ec="#222", label="blue module"),
       Patch(fc=COL["absent"], ec="#222", label="not in 2000-HVG panel")]
fig.legend(handles=leg, loc="lower center", ncol=3, fontsize=8.6,
           frameon=False, bbox_to_anchor=(0.5, 0.005))

fig.suptitle("Cross-pipeline consistency with the earlier Plant_response_to_radiation DESeq2/WGCNA analysis "
             "(shared studies OSD-498/502/508/510/658)",
             fontsize=9.2, y=0.995, color="#555")
fig.tight_layout(rect=[0, 0.05, 1, 0.96])
fig.savefig(f"{OUTDIR}/old_repo_overlap.svg")
fig.savefig(f"{OUTDIR}/old_repo_overlap.png", dpi=200)
print("wrote", f"{OUTDIR}/old_repo_overlap.svg/.png")

# ---------------------------------------------------------------- console summary
print(f"\nbackground DEG rate in 2000-HVG universe: {n}/{N} = {bg*100:.1f}%")
for m in order:
    s = stats[m]
    print(f"  {m:10s} {s['k']:3d}/{s['K']:4d} = {s['frac']*100:5.1f}%  "
          f"fold={s['fold']:.2f}  P={s['p']:.2e}")
print("\nDDR provenance:")
for sym, agi in DDR:
    print(f"  {sym:8s} {agi} -> {provenance(agi):10s} "
          f"{'DEG' if agi.upper() in old_deg else '.'}")
