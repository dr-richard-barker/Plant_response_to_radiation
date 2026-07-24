#!/usr/bin/env python3
"""
07_kinetic_wgcna.py
Phase 3.1 — Kinetic WGCNA: gene co-expression network with time-resolved
module identification.

Goal: Build a gene co-expression network on the radiation-response expression
matrix, identify modules, and classify them as "early-response" vs "late-
adaptation" based on their correlation with time post-exposure. This directly
addresses the user's request for kinetic modules that distinguish the immediate
DNA-damage/oxidative response from the slower hormonal adaptation phase.

METHOD
  1. Load expression (within-study centered, top HVGs) for non-control samples
     with valid time post-exposure.
  2. Run WGCNA (R) via rpy2 or subprocess:
     - Soft-thresholding power selection (scale-free topology fit)
     - Block-wise module detection (dynamic tree cut)
     - Module eigengenes (MEs) per sample
  3. Correlate MEs with time post-exposure (Spearman) to classify:
     - Early-response modules: ME negatively correlated with time
       (high expression early, decays later) — immediate stress response
     - Late-adaptation modules: ME positively correlated with time
       (low early, rises later) — recovery/adaptation
     - Sustained modules: no significant time correlation
  4. Also correlate MEs with dose and LET for additional trait associations.
  5. Identify hub genes per module (module membership > 0.8).

OUTPUTS
  - wgcna_modules.csv: gene -> module assignment, module membership, hub flag
  - wgcna_module_traits.csv: module -> correlation with time/dose/LET, p-value
  - wgcna_module_classification.csv: module -> early/late/sustained classification
  - wgcna_module_eigengenes.csv: sample -> ME values per module
  - wgcna_soft_threshold.json: selected soft-thresholding power
  - wgcna_summary.json: summary statistics
"""
import json, os, subprocess, sys
from pathlib import Path
import numpy as np
import pandas as pd

DATA_DIR = Path("/mnt/results/zenodo_bundle/data")
RESULTS_DIR = Path("/mnt/results/zenodo_bundle/results")
OUT_DIR = RESULTS_DIR / "wgcna"
EXPR_PATH = Path("/workspace/expression_raw.h5ad")
N_HVG = 2000  # reduced for n=72 samples (WGCNA needs n >> genes/blocks)
MIN_MODULE_SIZE = 20

# --- DDR-core reclassification (see script 14 and manuscript Section 5) ---
# A tight module carrying the canonical DNA double-strand-break repair /
# homologous-recombination / cell-cycle-checkpoint effectors is the
# transcriptional core of the DNA-damage response, not WGCNA's usual grey bin.
# It is relabelled "DDR-core" even when its time correlation is non-significant
# (the sharp, transient DDR genes do not track the smooth early-response decay).
OLD_DEG_PATH = DATA_DIR / "old_repo_deg_table.csv"  # bundled; input to script 14
MIN_DDR_MARKERS = 3          # min canonical effectors present in a module
DDR_DEG_P_MAX = 0.05         # module must also be enriched for old-repo DEGs
CANONICAL_DDR = {            # AGI id -> symbol (DSB repair / HR / checkpoint effectors)
    "AT4G21070": "BRCA1", "AT4G00020": "BRCA2A", "AT5G01630": "BRCA2B",
    "AT5G20850": "RAD51", "AT3G19210": "RAD54", "AT2G31320": "PARP1",
    "AT4G02390": "PARP2", "AT1G07500": "SMR5", "AT3G27630": "SMR7",
    "AT5G24280": "GMI1", "AT5G48720": "XRI1", "AT1G16970": "KU70",
    "AT1G48050": "KU80", "AT5G54260": "MRE11", "AT5G66130": "RAD17",
    "AT1G02970": "WEE1",
}


def _hyper_ge(k, K, n, N):
    """P(X >= k): hypergeometric upper tail. K successes in pop N, n draws."""
    from math import comb
    if n == 0 or N == 0:
        return 1.0
    return sum(comb(K, i) * comb(N - K, n - i)
               for i in range(k, min(K, n) + 1)) / comb(N, n)


def prepare_expression():
    """Load expression, within-study center, select HVGs, filter to non-control
    samples with valid time. Save as CSV for the R script."""
    import anndata as ad
    adata = ad.read_h5ad(EXPR_PATH)
    obs = adata.obs.copy()
    expr = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)

    # Filter to non-control samples with valid time
    # IsControl is a categorical with values '', 'False', 'True'
    ic_str = obs["IsControl"].astype(str)
    is_ctrl = ic_str.str.lower().isin(["true", "1", "yes"])
    obs["IsControl"] = is_ctrl
    obs["TimePostExposure_h"] = pd.to_numeric(obs["TimePostExposure_h"], errors="coerce")
    mask = (~obs["IsControl"]) & (~obs["TimePostExposure_h"].isna())
    expr = expr.loc[mask]
    obs = obs.loc[mask]
    print(f"  Non-control samples with time: {len(expr)}")

    # Select top HVGs
    expr = expr.loc[:, expr.var(axis=0).sort_values(ascending=False).head(N_HVG).index]
    print(f"  HVGs selected: {expr.shape[1]}")

    # Within-study centering (response relative to control)
    expr_log = np.log1p(expr.values.astype(float))
    study = obs["StudyID"].values
    is_ctrl_arr = is_ctrl.loc[mask].values  # use the correctly-computed is_ctrl
    study_mean = np.zeros_like(expr_log)
    for sid in np.unique(study):
        m = study == sid
        cm = m & is_ctrl_arr
        ref = expr_log[cm].mean(axis=0) if cm.sum() > 0 else expr_log[m].mean(axis=0)
        study_mean[m] = ref
    expr_centered = expr_log - study_mean

    # Build expression DataFrame with sample IDs
    expr_df = pd.DataFrame(expr_centered, index=expr.index, columns=expr.columns)
    expr_df.index.name = "Sample"

    # Trait matrix: time, dose, LET (for module-trait correlation)
    traits = pd.DataFrame({
        "Time_h": obs["TimePostExposure_h"].values,
        "Dose_Gy": pd.to_numeric(obs["AbsorbedDose_Gy"], errors="coerce").fillna(0).values,
        "LET_keV_um": pd.to_numeric(obs["LET_keV_um"], errors="coerce").fillna(0).values,
    }, index=expr.index)
    traits.index.name = "Sample"

    return expr_df, traits, obs


def run_wgcna_r(expr_path, trait_path, out_dir):
    """Run WGCNA via R subprocess. Returns paths to output files."""
    r_script = f"""
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# Load data (rows = samples, cols = genes)
expr <- read.csv("{expr_path}", row.names=1, check.names=FALSE)
traits <- read.csv("{trait_path}", row.names=1, check.names=FALSE)
cat("Expression (samples x genes):", dim(expr), "\\n")
cat("Traits:", dim(traits), "\\n")

# WGCNA datExpr: rows = samples, cols = genes (standard WGCNA convention)
datExpr <- as.matrix(expr)
datTraits <- as.data.frame(traits)
# Ensure sample order matches
datTraits <- datTraits[rownames(datExpr), ]
cat("datExpr (samples x genes):", dim(datExpr), "\\n")

gsg <- goodSamplesGenes(datExpr, verbose=3)
if (!gsg$allOK) {{
  if (sum(!gsg$goodGenes) > 0) datExpr <- datExpr[gsg$goodGenes, ]
  if (sum(!gsg$goodSamples) > 0) {{ datExpr <- datExpr[, gsg$goodSamples]
    datTraits <- datTraits[gsg$goodSamples, ] }}
}}
cat("After QC (genes x samples):", dim(datExpr), "\\n")

# Soft-thresholding power selection
powers <- c(1:10, seq(12, 20, by=2))
sft <- pickSoftThreshold(datExpr, powerVector=powers, networkType="signed", verbose=2)
power <- sft$powerEst
if (is.na(power)) power <- 12
cat("Selected soft power:", power, "\\n")
writeLines(as.character(power), "{out_dir}/soft_power.txt")

# Block-wise module detection
# With n=72 samples, use a lower merge height and smaller min module to
# get more granular modules
net <- blockwiseModules(datExpr, power=power, networkType="signed",
  TOMType="signed", minModuleSize={MIN_MODULE_SIZE},
  mergeCutHeight=0.35, numericLabels=TRUE,
  saveTOMs=FALSE, verbose=3)

moduleColors <- labels2colors(net$colors)
moduleLabels <- net$colors
nModules <- length(unique(moduleColors))
cat("Modules found (incl grey):", nModules, "\\n")
cat("Module sizes:\\n")
print(table(moduleColors))
cat("length(moduleColors):", length(moduleColors), "\\n")
cat("nrow(datExpr):", nrow(datExpr), "ncol(datExpr):", ncol(datExpr), "\\n")

# If only grey module found, try with unsigned network (more permissive)
if (nModules <= 1) {{
  cat("Only grey module found with signed network. Trying unsigned...\\n")
  net <- blockwiseModules(datExpr, power=power, networkType="unsigned",
    TOMType="unsigned", minModuleSize={MIN_MODULE_SIZE},
    mergeCutHeight=0.35, numericLabels=TRUE,
    saveTOMs=FALSE, verbose=3)
  moduleColors <- labels2colors(net$colors)
  moduleLabels <- net$colors
  nModules <- length(unique(moduleColors))
  cat("Modules found (unsigned, incl grey):", nModules, "\\n")
  print(table(moduleColors))
}}

# Module eigengenes (one per sample per module)
# datExpr: rows=samples, cols=genes. moduleEigengenes expects this orientation.
# Returns: eigengenes with rows=samples, cols=modules
ME_list <- moduleEigengenes(datExpr, moduleColors)
MEs <- ME_list$eigengenes  # rows = samples, cols = MEs
cat("MEs dim (samples x modules):", dim(MEs), "\\n")
cat("MEs rownames (first 5):", head(rownames(MEs)), "\\n")
MEs_ordered <- orderMEs(MEs)

# Save module eigengenes (samples x modules)
MEs_df <- data.frame(Sample=rownames(MEs_ordered), MEs_ordered, check.names=FALSE)
write.csv(MEs_df, "{out_dir}/module_eigengenes.csv", row.names=FALSE)

# Module-trait correlations (Spearman) — MEs rows = samples, traits rows = samples
moduleTraitCor <- cor(MEs_ordered, datTraits, use="p", method="spearman")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples=nrow(datTraits))
mt_df <- data.frame(Module=colnames(MEs_ordered))
for (i in 1:ncol(moduleTraitCor)) {{
  mt_df[[paste0(colnames(moduleTraitCor)[i], "_rho")]] <- moduleTraitCor[,i]
  mt_df[[paste0(colnames(moduleTraitCor)[i], "_pvalue")]] <- moduleTraitPvalue[,i]
}}
write.csv(mt_df, "{out_dir}/module_traits.csv", row.names=FALSE)

# Gene-module membership (kME = correlation of each gene with each ME)
# datExpr: rows=samples, cols=genes. cor(datExpr, MEs) correlates columns:
#   result rows = genes (cols of datExpr), cols = MEs
geneModuleMembership <- cor(datExpr, MEs_ordered, use="p", method="spearman")
colnames(geneModuleMembership) <- paste0("kME_", colnames(MEs_ordered))
genes_df <- data.frame(Gene=colnames(datExpr), Module=moduleColors, stringsAsFactors=FALSE)
genes_df <- cbind(genes_df, geneModuleMembership)
# Hub gene flag: |kME in own module| > 0.7
# Column names are kME_ME<color>, module is <color>, so kME col = paste0("kME_ME", color)
genes_df$IsHub <- FALSE
for (i in 1:nrow(genes_df)) {{
  kME_col <- paste0("kME_ME", genes_df$Module[i])
  if (kME_col %in% colnames(genes_df)) {{
    genes_df$IsHub[i] <- abs(genes_df[[kME_col]][i]) > 0.7
  }}
}}
write.csv(genes_df, "{out_dir}/modules.csv", row.names=FALSE)

cat("Done. Modules:", nModules, "\\n")
cat("Module sizes:\\n")
print(table(moduleColors))
"""
    r_path = out_dir / "run_wgcna.R"
    with open(r_path, "w") as f:
        f.write(r_script)
    print(f"  R script -> {r_path}")
    result = subprocess.run(["Rscript", str(r_path)], capture_output=True, text=True, timeout=1800)
    print(result.stdout[-2000:] if len(result.stdout) > 2000 else result.stdout)
    if result.returncode != 0:
        print("STDERR:", result.stderr[-2000:])
        raise RuntimeError(f"WGCNA R script failed (exit {result.returncode})")
    return out_dir


def classify_modules(out_dir, obs):
    """Classify modules as early-response, late-adaptation, or sustained
    based on module-trait correlations with time."""
    traits_df = pd.read_csv(out_dir / "module_traits.csv")
    modules_df = pd.read_csv(out_dir / "modules.csv")

    classifications = []
    for _, row in traits_df.iterrows():
        module_me = row["Module"]  # e.g. "MEblue"
        module = module_me.replace("ME", "", 1)  # e.g. "blue"
        time_rho = row.get("Time_h_rho", 0)
        time_p = row.get("Time_h_pvalue", 1)
        dose_rho = row.get("Dose_Gy_rho", 0)
        let_rho = row.get("LET_keV_um_rho", 0)
        n_genes = (modules_df["Module"] == module).sum()
        n_hubs = ((modules_df["Module"] == module) & (modules_df["IsHub"])).sum()
        if time_p < 0.05 and time_rho < -0.3:
            cls = "early-response"
        elif time_p < 0.05 and time_rho > 0.3:
            cls = "late-adaptation"
        elif time_p < 0.05:
            cls = "time-associated"
        else:
            cls = "sustained"
        classifications.append({"Module": module, "ME_name": module_me,
                                "Classification": cls,
                                "n_genes": int(n_genes), "n_hubs": int(n_hubs),
                                "Time_rho": time_rho, "Time_pvalue": time_p,
                                "Dose_rho": dose_rho, "LET_rho": let_rho})
    cls_df = pd.DataFrame(classifications)
    cls_df = _reclassify_ddr_core(cls_df, modules_df)
    cls_df.to_csv(out_dir / "module_classification.csv", index=False)
    return cls_df


def _reclassify_ddr_core(cls_df, modules_df):
    """Relabel a non-early/late module as 'DDR-core' if it carries the canonical
    DSB-repair effectors and is enriched for old-repo radiation DEGs. Keeps the
    label stable across re-runs (deterministic, marker-driven)."""
    universe = set(modules_df["Gene"].str.upper())
    # old-repo radiation DEGs (any Up/Down contrast), if the bundled table exists
    deg_in_univ = set()
    if OLD_DEG_PATH.exists():
        deg = pd.read_csv(OLD_DEG_PATH)
        contrasts = [c for c in deg.columns if c != "gene_id"]
        is_deg = deg[contrasts].isin(["Up", "Down"]).any(axis=1)
        deg_in_univ = set(deg.loc[is_deg, "gene_id"].str.upper()) & universe
    N, n = len(universe), len(deg_in_univ)

    for i, row in cls_df.iterrows():
        if row["Classification"] in ("early-response", "late-adaptation"):
            continue  # never override a significant time signature
        genes = set(modules_df.loc[modules_df["Module"] == row["Module"], "Gene"].str.upper())
        n_markers = len(genes & set(CANONICAL_DDR))
        if n_markers < MIN_DDR_MARKERS:
            continue
        # optional DEG-enrichment gate (only when the DEG table is bundled)
        deg_p = _hyper_ge(len(genes & deg_in_univ), len(genes), n, N) if n else 0.0
        if deg_p <= DDR_DEG_P_MAX:
            markers = sorted(CANONICAL_DDR[g] for g in genes & set(CANONICAL_DDR))
            cls_df.at[i, "Classification"] = "DDR-core"
            print(f"  {row['Module']}: reclassified sustained -> DDR-core "
                  f"({n_markers} DSB-repair effectors: {', '.join(markers)}; "
                  f"DEG-enrichment P={deg_p:.1e})")
    return cls_df


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("Preparing expression for WGCNA...")
    expr_df, traits, obs = prepare_expression()
    expr_path = OUT_DIR / "expression_input.csv"
    trait_path = OUT_DIR / "traits_input.csv"
    expr_df.to_csv(expr_path)
    traits.to_csv(trait_path)
    print(f"  Expression -> {expr_path}  ({expr_df.shape})")
    print(f"  Traits -> {trait_path}  ({traits.shape})")

    print("\nRunning WGCNA in R...")
    run_wgcna_r(str(expr_path), str(trait_path), OUT_DIR)

    print("\nClassifying modules...")
    cls_df = classify_modules(OUT_DIR, obs)
    print("\n=== Module classification ===")
    print(cls_df.groupby("Classification").size().to_string())
    print()
    print("=== Early-response modules ===")
    early = cls_df[cls_df["Classification"] == "early-response"].sort_values("Time_rho")
    print(early[["Module","n_genes","n_hubs","Time_rho","Time_pvalue"]].to_string(index=False))
    print()
    print("=== Late-adaptation modules ===")
    late = cls_df[cls_df["Classification"] == "late-adaptation"].sort_values("Time_rho", ascending=False)
    print(late[["Module","n_genes","n_hubs","Time_rho","Time_pvalue"]].to_string(index=False))

    # Summary
    soft_power = (OUT_DIR / "soft_power.txt").read_text().strip() if (OUT_DIR / "soft_power.txt").exists() else "NA"
    summary = {
        "n_samples": int(expr_df.shape[0]),
        "n_genes": int(expr_df.shape[1]),
        "soft_power": soft_power,
        "n_modules": int(len(cls_df)),
        "n_early_response": int((cls_df["Classification"] == "early-response").sum()),
        "n_late_adaptation": int((cls_df["Classification"] == "late-adaptation").sum()),
        "n_sustained": int((cls_df["Classification"] == "sustained").sum()),
        "n_ddr_core": int((cls_df["Classification"] == "DDR-core").sum()),
    }
    with open(OUT_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSummary -> {OUT_DIR / 'summary.json'}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
