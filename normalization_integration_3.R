# ==============================================================================
# Part 3: Normalization, Integration, and Dimensionality Reduction
# ==============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Load the merged Seurat object
merged_seurat <- readRDS("../results/merged_seurat_raw.rds")

cat("Starting with:", ncol(merged_seurat), "cells\n")

# ==============================================================================
Increase memory limit for parallel processing
# ==============================================================================

# Option 1: Increase the future globals size (recommended)
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB

# Option 2: Disable parallelization (slower but uses less memory)
# plan("sequential")

cat("Memory limit increased for large dataset processing\n")

# Load the merged Seurat object
merged_seurat <- readRDS("../results/merged_seurat_raw.rds")

cat("Starting with:", ncol(merged_seurat), "cells\n")

# ==============================================================================
# STEP 1: Normalize and Find Variable Features
# ==============================================================================

cat("\n1.Normalizing data...\n")

# Split by sample for integration
seurat_list <- SplitObject(merged_seurat, split.by = "sample_id")

# Normalize and find variable features for each sample
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  return(x)
})

cat("✓ Normalization complete\n")

# ==============================================================================
# STEP 2: Integration (Critical for removing batch effects)
# ==============================================================================

cat("\n2.Integrating samples to remove batch effects...\n")

# Select integration features
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 2000)

# Find integration anchors
anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  anchor.features = features,
  dims = 1:30,
  verbose = FALSE
)

# Integrate data
integrated_seurat <- IntegrateData(anchorset = anchors, dims = 1:30, verbose = FALSE)

cat("✓ Integration complete\n")

# ==============================================================================
# STEP 3: Scale Data and Run PCA
# ==============================================================================

cat("\n3.Scaling data and running PCA...\n")

DefaultAssay(integrated_seurat) <- "integrated"

# Scale the data
integrated_seurat <- ScaleData(integrated_seurat, verbose = FALSE)

# Run PCA
integrated_seurat <- RunPCA(integrated_seurat, npcs = 50, verbose = FALSE)

# Visualize PCA
pdf("../figures/03_PCA_analysis.pdf", width = 14, height = 10)

# Elbow plot
print(ElbowPlot(integrated_seurat, ndims = 50) + 
        ggtitle("Elbow Plot - Determine Number of PCs"))

# PCA by timepoint
print(DimPlot(integrated_seurat, reduction = "pca", group.by = "timepoint", 
              cols = c("Pre" = "#4DBBD5FF", "Post" = "#E64B35FF")) +
        ggtitle("PCA - Pre vs Post Exercise"))

# PCA by sample
print(DimPlot(integrated_seurat, reduction = "pca", group.by = "sample_id") +
        ggtitle("PCA - By Sample"))

# PCA by participant
print(DimPlot(integrated_seurat, reduction = "pca", group.by = "participant") +
        ggtitle("PCA - By Participant"))

dev.off()

cat("✓ PCA complete\n")

# ==============================================================================
# STEP 4: UMAP and tSNE
# ==============================================================================

cat("\n4.Running UMAP and tSNE...\n")

# Run UMAP (using first 30 PCs based on elbow plot)
integrated_seurat <- RunUMAP(integrated_seurat, reduction = "pca", dims = 1:30, verbose = FALSE)

# Run tSNE
integrated_seurat <- RunTSNE(integrated_seurat, reduction = "pca", dims = 1:30, verbose = FALSE)

cat("✓ Dimensionality reduction complete\n")

# ==============================================================================
# STEP 5: Clustering
# ==============================================================================

cat("\n5.Performing clustering...\n")

# Find neighbors
integrated_seurat <- FindNeighbors(integrated_seurat, reduction = "pca", dims = 1:30, verbose = FALSE)

# Find clusters at multiple resolutions
resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)
for (res in resolutions) {
  integrated_seurat <- FindClusters(integrated_seurat, resolution = res, verbose = FALSE)
}

cat("✓ Clustering complete\n")

# ==============================================================================
# STEP 6: Visualize UMAP with Different Groupings
# ==============================================================================

cat("\n6.Creating UMAP visualizations...\n")

# Use resolution 0.5 as default
Idents(integrated_seurat) <- "integrated_snn_res.0.5"

pdf("../figures/04_UMAP_clustering.pdf", width = 16, height = 12)

# UMAP by clusters
p1 <- DimPlot(integrated_seurat, reduction = "umap", label = TRUE, label.size = 6) +
  ggtitle("UMAP - Clusters (Resolution 0.5)")

# UMAP by timepoint
p2 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "timepoint",
              cols = c("Pre" = "#4DBBD5FF", "Post" = "#E64B35FF")) +
  ggtitle("UMAP - Pre vs Post Exercise")

# UMAP by sample
p3 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "sample_id") +
  ggtitle("UMAP - By Sample (Integration Quality Check)")

# UMAP by participant
p4 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "participant") +
  ggtitle("UMAP - By Participant")

# Combined plot
print((p1 | p2) / (p3 | p4))

# Split by timepoint
print(DimPlot(integrated_seurat, reduction = "umap", split.by = "timepoint", label = TRUE) +
        ggtitle("UMAP - Clusters Split by Timepoint"))

dev.off()

cat("✓ UMAP plots saved\n")

# ==============================================================================
# STEP 7: Compare Different Clustering Resolutions
# ==============================================================================

pdf("../figures/05_clustering_resolutions.pdf", width = 18, height = 12)

plot_list <- list()
for (i in 1:length(resolutions)) {
  res <- resolutions[i]
  plot_list[[i]] <- DimPlot(integrated_seurat, reduction = "umap", 
                            group.by = paste0("integrated_snn_res.", res),
                            label = TRUE, label.size = 4) +
    ggtitle(paste("Resolution", res)) +
    NoLegend()
}

print(wrap_plots(plot_list, ncol = 3))

dev.off()

# ==============================================================================
# STEP 8: Save Processed Object
# ==============================================================================

saveRDS(integrated_seurat, "../results/integrated_seurat_clustered.rds")

cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("INTEGRATION & CLUSTERING COMPLETE\n")
cat(paste(rep("=", 50), collapse = ""), "\n")
cat(paste0("Total cells:  ", ncol(integrated_seurat), "\n"))
cat(paste0("Number of clusters (res 0.5): ", length(unique(integrated_seurat$integrated_snn_res.0.5)), "\n"))
cat(paste(rep("=", 50), collapse = ""), "\n\n")

cat("✓ Analysis complete!\n")
cat("Next step: Run 04_cell_type_annotation.R\n")