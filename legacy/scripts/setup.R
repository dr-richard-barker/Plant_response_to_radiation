#!/usr/bin/env Rscript

# setup.R: Install and verify dependencies for SBGNview and pathview analysis

# Function to install a package if not already installed
install_if_missing <- function(pkg, method = "cran") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste0("Installing package: ", pkg, "\n"))
    if (method == "cran") {
      install.packages(pkg, repos = "https://cran.r-project.org")
    } else if (method == "bioc") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cran.r-project.org")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
  } else {
    cat(paste0("Package already installed: ", pkg, "\n"))
  }
}

# Essential CRAN packages
cran_packages <- c("BiocManager", "openxlsx", "dplyr", "ggplot2")
for (pkg in cran_packages) {
  install_if_missing(pkg, method = "cran")
}

# Essential Bioconductor packages
bioc_packages <- c("SBGNview", "pathview", "KEGGREST")
for (pkg in bioc_packages) {
  install_if_missing(pkg, method = "bioc")
}

cat("\nVerification of essential packages:\n")
all_loaded <- TRUE
for (pkg in c(cran_packages, bioc_packages)) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(paste0("[OK] ", pkg, "\n"))
  } else {
    cat(paste0("[ERROR] Failed to load ", pkg, "\n"))
    all_loaded <- FALSE
  }
}

if (all_loaded) {
  cat("\nAll dependencies are successfully installed.\n")
} else {
  cat("\nSome dependencies are missing. Please check the error messages above.\n")
  quit(status = 1)
}
