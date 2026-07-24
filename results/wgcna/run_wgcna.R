
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# Load data (rows = samples, cols = genes)
expr <- read.csv("/mnt/results/zenodo_bundle/results/wgcna/expression_input.csv", row.names=1, check.names=FALSE)
traits <- read.csv("/mnt/results/zenodo_bundle/results/wgcna/traits_input.csv", row.names=1, check.names=FALSE)
cat("Expression (samples x genes):", dim(expr), "\n")
cat("Traits:", dim(traits), "\n")

# WGCNA datExpr: rows = samples, cols = genes (standard WGCNA convention)
datExpr <- as.matrix(expr)
datTraits <- as.data.frame(traits)
# Ensure sample order matches
datTraits <- datTraits[rownames(datExpr), ]
cat("datExpr (samples x genes):", dim(datExpr), "\n")

gsg <- goodSamplesGenes(datExpr, verbose=3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) datExpr <- datExpr[gsg$goodGenes, ]
  if (sum(!gsg$goodSamples) > 0) { datExpr <- datExpr[, gsg$goodSamples]
    datTraits <- datTraits[gsg$goodSamples, ] }
}
cat("After QC (genes x samples):", dim(datExpr), "\n")

# Soft-thresholding power selection
powers <- c(1:10, seq(12, 20, by=2))
sft <- pickSoftThreshold(datExpr, powerVector=powers, networkType="signed", verbose=2)
power <- sft$powerEst
if (is.na(power)) power <- 12
cat("Selected soft power:", power, "\n")
writeLines(as.character(power), "/mnt/results/zenodo_bundle/results/wgcna/soft_power.txt")

# Block-wise module detection
# With n=72 samples, use a lower merge height and smaller min module to
# get more granular modules
net <- blockwiseModules(datExpr, power=power, networkType="signed",
  TOMType="signed", minModuleSize=20,
  mergeCutHeight=0.35, numericLabels=TRUE,
  saveTOMs=FALSE, verbose=3)

moduleColors <- labels2colors(net$colors)
moduleLabels <- net$colors
nModules <- length(unique(moduleColors))
cat("Modules found (incl grey):", nModules, "\n")
cat("Module sizes:\n")
print(table(moduleColors))
cat("length(moduleColors):", length(moduleColors), "\n")
cat("nrow(datExpr):", nrow(datExpr), "ncol(datExpr):", ncol(datExpr), "\n")

# If only grey module found, try with unsigned network (more permissive)
if (nModules <= 1) {
  cat("Only grey module found with signed network. Trying unsigned...\n")
  net <- blockwiseModules(datExpr, power=power, networkType="unsigned",
    TOMType="unsigned", minModuleSize=20,
    mergeCutHeight=0.35, numericLabels=TRUE,
    saveTOMs=FALSE, verbose=3)
  moduleColors <- labels2colors(net$colors)
  moduleLabels <- net$colors
  nModules <- length(unique(moduleColors))
  cat("Modules found (unsigned, incl grey):", nModules, "\n")
  print(table(moduleColors))
}

# Module eigengenes (one per sample per module)
# datExpr: rows=samples, cols=genes. moduleEigengenes expects this orientation.
# Returns: eigengenes with rows=samples, cols=modules
ME_list <- moduleEigengenes(datExpr, moduleColors)
MEs <- ME_list$eigengenes  # rows = samples, cols = MEs
cat("MEs dim (samples x modules):", dim(MEs), "\n")
cat("MEs rownames (first 5):", head(rownames(MEs)), "\n")
MEs_ordered <- orderMEs(MEs)

# Save module eigengenes (samples x modules)
MEs_df <- data.frame(Sample=rownames(MEs_ordered), MEs_ordered, check.names=FALSE)
write.csv(MEs_df, "/mnt/results/zenodo_bundle/results/wgcna/module_eigengenes.csv", row.names=FALSE)

# Module-trait correlations (Spearman) — MEs rows = samples, traits rows = samples
moduleTraitCor <- cor(MEs_ordered, datTraits, use="p", method="spearman")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples=nrow(datTraits))
mt_df <- data.frame(Module=colnames(MEs_ordered))
for (i in 1:ncol(moduleTraitCor)) {
  mt_df[[paste0(colnames(moduleTraitCor)[i], "_rho")]] <- moduleTraitCor[,i]
  mt_df[[paste0(colnames(moduleTraitCor)[i], "_pvalue")]] <- moduleTraitPvalue[,i]
}
write.csv(mt_df, "/mnt/results/zenodo_bundle/results/wgcna/module_traits.csv", row.names=FALSE)

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
for (i in 1:nrow(genes_df)) {
  kME_col <- paste0("kME_ME", genes_df$Module[i])
  if (kME_col %in% colnames(genes_df)) {
    genes_df$IsHub[i] <- abs(genes_df[[kME_col]][i]) > 0.7
  }
}
write.csv(genes_df, "/mnt/results/zenodo_bundle/results/wgcna/modules.csv", row.names=FALSE)

cat("Done. Modules:", nModules, "\n")
cat("Module sizes:\n")
print(table(moduleColors))
