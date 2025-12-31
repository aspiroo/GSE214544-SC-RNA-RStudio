# ==============================================================================
# Complete Analysis Pipeline for Exercise-Induced Skeletal Muscle scRNA-seq
# Based on:  "Single-cell sequencing deconvolutes cellular responses to exercise 
# in human skeletal muscle" (Comm Biol 2022)
# Dataset: GSE214544
# ==============================================================================

# Clear environment
rm(list = ls())
gc()

# ==============================================================================
# Install and Load Required Packages
# ==============================================================================

cat("Checking and installing required packages...\n")

# CRAN packages
required_packages <- c(
  "Seurat", "dplyr", "ggplot2", "patchwork", "viridis",
  "pheatmap", "RColorBrewer", "gridExtra", "tidyr",
  "cowplot", "reshape2", "scales", "data.table"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste0("Installing ", pkg, "...\n"))
    install.packages(pkg)
  }
}

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c("limma", "ComplexHeatmap", "enrichplot", "clusterProfiler")
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste0("Installing ", pkg, " from Bioconductor...\n"))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

# Load core libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

cat("✓ All packages loaded successfully!\n\n")

# ==============================================================================
# Set Working Directory and Create Project Structure
# ==============================================================================

# Set to your data directory
setwd("E:/Compilers and Libraries/R Library/GSE214544_RAW")

# Create project directories
dirs <- c("results", "figures", "tables", "scripts")
for (dir in dirs) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}

cat("✓ Project directories created:\n")
cat(paste0("  - ", paste(dirs, collapse = "\n  - "), "\n\n"))

# ==============================================================================
# Sample Metadata Based on GSE214544
# ==============================================================================

# Based on your folder names
sample_metadata <- data.frame(
  sample_id = c(
    "GSM6611295_P15306_5001",
    "GSM6611296_P15306_5002",
    "GSM6611297_P14601_4004",
    "GSM6611298_P14601_4005",
    "GSM6611299_P15306_5003",
    "GSM6611300_P15306_5004",
    "GSM6611301_10X_20_003",
    "GSM6611302_10X_20_004"
  ),
  participant = c(
    "P15306", "P15306",
    "P14601", "P14601",
    "P15306", "P15306",
    "10X_20", "10X_20"
  ),
  timepoint = c(
    "Pre", "Post",
    "Pre", "Post",
    "Pre", "Post",
    "Pre", "Post"
  ),
  stringsAsFactors = FALSE
)

# Save metadata
write.csv(sample_metadata, "tables/sample_metadata.csv", row.names = FALSE)

cat("✓ Sample metadata created:\n\n")
print(sample_metadata)

cat("\n========================================\n")
cat("SETUP COMPLETE!\n")
cat("========================================\n")
cat("Your data structure:\n")
cat("E:/Compilers and Libraries/R Library/GSE214544_RAW/\n")
cat("  ├── data/\n")
cat("  │   ├── GSM6611295_P15306_5001/ (matrix, barcodes, features)\n")
cat("  │   ├── GSM6611296_P15306_5002/\n")
cat("  │   └── ...  (6 more samples)\n")
cat("  ├── results/    (analysis outputs will be saved here)\n")
cat("  ├── figures/    (plots will be saved here)\n")
cat("  ├── tables/     (tables will be saved here)\n")
cat("  └── scripts/    (you can save R scripts here)\n")
cat("\n")
cat("Next steps:\n")
cat("  1. Run data_loading_and_qc_UPDATED.R\n")
cat("  2. Run normalization_integration.R\n")
cat("  3. Run clustering_annotation. R\n")
cat("========================================\n")