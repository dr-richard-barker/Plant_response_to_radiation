#!/usr/bin/env python3
"""
15_osd782_validation.py

INDEPENDENT validation of the kinetic WGCNA modules on a held-out study.

OSD-782 (GeneLab GLDS-679) is a wild-type Arabidopsis gamma time-course
(0 / 0.1 / 1 Gy x 60 / 180 / 1440 / 4320 min, 3 reps = 36 samples) that is NOT
part of the 10-study cohort used to build the GP-AE / WGCNA modules, and uses a
low-dose regime (<=1 Gy) distinct from the cohort's 100 Gy gamma studies. It is
therefore a genuinely independent test set (unlike the cross-pipeline check in
script 14, which shares samples with the training cohort).

Two projection tests, neither of which refits the model:
  1. Enrichment - are OSD-782 radiation-responsive genes over-represented in the
     new blue (early-response) and grey (DDR-core) modules?
  2. Eigengene trajectory - projecting each new module onto OSD-782 and
     correlating its module score with time post-exposure, does the blue module
     independently reproduce its early-response signature (module score declining
     with time, matching the cohort's Time rho = -0.66)?

Inputs (bundled):
  data/osd782_normalized_counts.csv   GLDS-679 normalized counts (genes x 36 GSM)
  data/osd782_metadata.csv            SampleID, Dose_Gy, Time_min, Time_h, IsControl
  results/wgcna/modules.csv           new module assignments (2000 HVGs)
Outputs:
  results/wgcna/osd782_validation.json
  results/figures/osd782_validation.{svg,png}
"""
import csv, json, os, sys
from math import comb
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COUNTS = sys.argv[1] if len(sys.argv) > 1 else f"{HERE}/data/osd782_normalized_counts.csv"
META = sys.argv[2] if len(sys.argv) > 2 else f"{HERE}/data/osd782_metadata.csv"
MODS = sys.argv[3] if len(sys.argv) > 3 else f"{HERE}/results/wgcna/modules.csv"
OUTJSON = f"{HERE}/results/wgcna/osd782_validation.json"
FIG = f"{HERE}/results/figures/osd782_validation"

CANONICAL_DDR = {"AT4G21070": "BRCA1", "AT4G00020": "BRCA2A", "AT5G01630": "BRCA2B",
                 "AT5G20850": "RAD51", "AT3G19210": "RAD54", "AT2G31320": "PARP1",
                 "AT4G02390": "PARP2", "AT1G07500": "SMR5", "AT3G27630": "SMR7",
                 "AT5G24280": "GMI1", "AT5G48720": "XRI1", "AT1G16970": "KU70",
                 "AT1G48050": "KU80", "AT5G54260": "MRE11", "AT5G66130": "RAD17",
                 "AT1G02970": "WEE1"}


def hyper_ge(k, K, n, N):
    if n == 0 or N == 0:
        return 1.0
    return sum(comb(K, i) * comb(N - K, n - i)
               for i in range(k, min(K, n) + 1)) / comb(N, n)


def spearman(x, y):
    x = pd.Series(x).rank().values
    y = pd.Series(y).rank().values
    x = x - x.mean(); y = y - y.mean()
    d = np.sqrt((x**2).sum() * (y**2).sum())
    return float((x * y).sum() / d) if d else 0.0


# ---------------------------------------------------------------- load
expr = pd.read_csv(COUNTS, index_col=0)
expr.index = expr.index.str.upper()
meta = pd.read_csv(META).set_index("SampleID")
meta = meta.loc[[s for s in expr.columns if s in meta.index]]
expr = expr[meta.index]                       # align columns to metadata order
log = np.log1p(expr)                           # normalized counts -> log space
mods = pd.read_csv(MODS)
mods["Gene"] = mods["Gene"].str.upper()
new_mod = dict(zip(mods["Gene"], mods["Module"]))

ctrl = meta.index[meta["IsControl"]]
irr = meta.index[~meta["IsControl"]]
times_min = sorted(meta["Time_min"].unique())

# ---- within-study centering vs TIME-MATCHED controls (removes circadian/time baseline)
centered = pd.DataFrame(index=log.index, columns=log.columns, dtype=float)
for t in times_min:
    cc = [s for s in ctrl if meta.loc[s, "Time_min"] == t]
    ss = [s for s in meta.index if meta.loc[s, "Time_min"] == t]
    ref = log[cc].mean(axis=1)
    centered[ss] = log[ss].sub(ref, axis=0)

# ---------------------------------------------------------------- (1) response genes
# radiation response = irradiated samples deviate from time-matched control baseline.
# one-sample t across the 24 irradiated centered values; effect = mean centered log2-ish.
from scipy import stats
irr_vals = centered[irr]
tstat, pval = stats.ttest_1samp(irr_vals.T.values, 0.0, axis=0)
resp = pd.DataFrame({"gene": centered.index, "effect": irr_vals.mean(axis=1).values,
                     "t": tstat, "p": pval}).set_index("gene")
resp["padj"] = stats.false_discovery_control(np.nan_to_num(resp["p"], nan=1.0))
resp_genes = set(resp.index[(resp["padj"] < 0.05) & (resp["effect"].abs() > np.log(1.5))])

# ---------------------------------------------------------------- (2) enrichment vs modules
universe = set(new_mod) & set(expr.index)      # new-module genes measured in OSD-782
deg_u = resp_genes & universe
N, n = len(universe), len(deg_u)
bg = n / N if N else 0
enrich = {}
for m in ("blue", "turquoise", "grey"):
    mg = {g for g, mm in new_mod.items() if mm == m} & universe
    k = len(deg_u & mg)
    exp = len(mg) * n / N if N else 0
    enrich[m] = dict(k=k, module_n=len(mg), frac=k / len(mg) if mg else 0,
                     fold=(k / exp) if exp else 0, p=hyper_ge(k, len(mg), n, N))

# ---------------------------------------------------------------- (3) eigengene trajectory
# project each module onto OSD-782: module score = mean centered expr over module genes.
traj = {}
irr_time = meta.loc[irr, "Time_h"].values
for m in ("blue", "turquoise", "grey"):
    mg = [g for g, mm in new_mod.items() if mm == m and g in centered.index]
    score = centered.loc[mg, irr].mean(axis=0).values
    traj[m] = dict(n_genes=len(mg), time_rho=spearman(irr_time, score),
                   score_by_time={float(t): float(centered.loc[mg,
                       [s for s in irr if meta.loc[s, "Time_h"] == t]].mean().mean())
                       for t in sorted(set(irr_time))})

# ---------------------------------------------------------------- DDR gene behaviour
ddr = {}
for agi, sym in CANONICAL_DDR.items():
    if agi in resp.index:
        ddr[sym] = dict(agi=agi, effect=round(float(resp.loc[agi, "effect"]), 3),
                        padj=float(resp.loc[agi, "padj"]),
                        module=new_mod.get(agi, "not-in-HVG-panel"),
                        responsive=bool(agi in resp_genes))

result = dict(
    study="OSD-782 / GLDS-679", n_samples=int(expr.shape[1]),
    n_control=int(len(ctrl)), n_irradiated=int(len(irr)),
    doses_Gy=sorted(meta["Dose_Gy"].unique()), times_h=sorted(meta["Time_h"].unique()),
    n_genes_measured=int(expr.shape[0]),
    n_response_genes=int(len(resp_genes)),
    universe_n=N, response_in_universe=n, background_frac=round(bg, 4),
    enrichment=enrich, module_trajectory=traj, ddr_genes=ddr,
    cohort_blue_time_rho=-0.66)
os.makedirs(os.path.dirname(OUTJSON), exist_ok=True)
json.dump(result, open(OUTJSON, "w"), indent=2)

# ---------------------------------------------------------------- console
print(f"OSD-782 (GLDS-679): {expr.shape[1]} samples, {expr.shape[0]} genes, "
      f"{len(ctrl)} control / {len(irr)} irradiated")
print(f"radiation-response genes (FDR<0.05, |FC|>1.5): {len(resp_genes)}  "
      f"({n}/{N} in new-module universe, bg {bg*100:.1f}%)")
print("\nEnrichment in new modules (independent held-out data):")
for m in ("blue", "turquoise", "grey"):
    e = enrich[m]
    print(f"  {m:10s} {e['k']:3d}/{e['module_n']:4d} = {e['frac']*100:5.1f}%  "
          f"fold={e['fold']:.2f}  P={e['p']:.2e}")
print("\nProjected module-score vs time among irradiated (early-response = negative):")
for m in ("blue", "turquoise", "grey"):
    print(f"  {m:10s} Time rho = {traj[m]['time_rho']:+.3f}  (n_genes={traj[m]['n_genes']})")
print(f"  [cohort blue Time rho was -0.66]")

# ---------------------------------------------------------------- figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none", "font.size": 10})
COL = {"blue": "#3B6EA5", "turquoise": "#3FB6A8", "grey": "#8A8A8A"}
fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.8))

order = ["blue", "turquoise", "grey"]
labels = ["blue\n(early-response)", "turquoise\n(sustained)", "grey\n(DDR-core)"]
vals = [enrich[m]["frac"] * 100 for m in order]
axA.bar(range(3), vals, color=[COL[m] for m in order], edgecolor="#222", width=0.62, zorder=3)
axA.axhline(bg * 100, ls="--", lw=1.4, color="#333", zorder=2)
axA.text(2.46, bg * 100 + 0.6, f"background {bg*100:.1f}%", ha="right", va="bottom",
         fontsize=8.5, style="italic", color="#333")
for i, m in enumerate(order):
    e = enrich[m]
    star = " *" if e["p"] < 0.05 else ""
    axA.text(i, vals[i] + 0.6, f"{e['k']}/{e['module_n']}\n{e['fold']:.2f}x{star}",
             ha="center", va="bottom", fontsize=8.6,
             fontweight="bold" if e["p"] < 0.05 else "normal")
axA.set_xticks(range(3)); axA.set_xticklabels(labels)
axA.set_ylabel("% of module genes responsive in OSD-782")
axA.set_ylim(0, max(vals) * 1.35 + 2)
axA.set_title("A  Held-out OSD-782 response genes recover the\nearly-response & DDR modules",
              fontsize=10.3, loc="left", fontweight="bold")
axA.grid(axis="y", ls=":", color="#ccc", zorder=0)
for s in ("top", "right"): axA.spines[s].set_visible(False)

for m in order:
    sd = traj[m]["score_by_time"]
    ts = sorted(sd); ys = [sd[t] for t in ts]
    axB.plot(ts, ys, "-o", color=COL[m], lw=2, ms=6,
             label=f"{m} (rho={traj[m]['time_rho']:+.2f})", zorder=3)
axB.axhline(0, color="#999", lw=1, ls=":")
axB.set_xlabel("time post-exposure (h)")
axB.set_ylabel("projected module score\n(centered vs time-matched control)")
axB.set_title("B  Blue module score declines with time in OSD-782\n(independent early-response signature)",
              fontsize=10.3, loc="left", fontweight="bold")
axB.legend(fontsize=8.5, frameon=False, loc="best")
for s in ("top", "right"): axB.spines[s].set_visible(False)
axB.grid(ls=":", color="#eee", zorder=0)

fig.suptitle("Independent validation on held-out study OSD-782 / GLDS-679 "
             "(low-dose gamma, not in the training cohort)", fontsize=9.2, y=0.99, color="#555")
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(f"{FIG}.svg"); fig.savefig(f"{FIG}.png", dpi=200)
print(f"\nwrote {OUTJSON}\nwrote {FIG}.svg/.png")
