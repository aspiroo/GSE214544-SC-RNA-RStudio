# ==============================================================================
# Part 2: Data Loading and Quality Control 
# ==============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Set working directory to where your data folders are
setwd("E:/Compilers and Libraries/R Library/GSE214544_RAW/data")

# Load the sample metadata created in part 1
sample_metadata <- read.csv("../tables/sample_metadata.csv")

# ==============================================================================
# STEP 1: Load Data from Each Folder
# ==============================================================================

# Initialize list to store Seurat objects
seurat_objects <- list()

# Loop through samples and create Seurat objects
for (i in 1:nrow(sample_metadata)) {
  sample_id <- sample_metadata$sample_id[i]
  
  cat(paste0("\n", paste(rep("=", 50), collapse = ""), "\n"))
  cat(paste0("Loading sample ", i, "/", nrow(sample_metadata), ": ", sample_id, "\n"))
  cat(paste0(paste(rep("=", 50), collapse = ""), "\n"))
  
  # Construct the folder path
  data_dir <- paste0("./", sample_id)
  
  # Check if directory exists
  if (! dir.exists(data_dir)) {
    cat(paste0("✗ Directory not found: ", data_dir, "\n"))
    next
  }
  
  # List files in directory to verify structure
  files <- list.files(data_dir)
  cat("Files found:\n")
  print(files)
  
  tryCatch({
    # Read 10X data from the folder
    # This will automatically look for matrix.mtx(.gz), barcodes.tsv(.gz), features.tsv(.gz) or genes.tsv(.gz)
    counts <- Read10X(data.dir = data_dir)
    
    # Create Seurat object
    seurat_obj <- CreateSeuratObject(
      counts = counts,
      project = sample_metadata$participant[i],
      min.cells = 3,        # Include features detected in at least 3 cells
      min.features = 200    # Include cells with at least 200 features
    )
    
    # Add metadata
    seurat_obj$sample_id <- sample_id
    seurat_obj$participant <- sample_metadata$participant[i]
    seurat_obj$timepoint <- sample_metadata$timepoint[i]
    
    # Store in list
    seurat_objects[[sample_id]] <- seurat_obj
    
    cat(paste0("✓ Successfully loaded ", sample_id, "\n"))
    cat(paste0("  - Cells: ", ncol(seurat_obj), "\n"))
    cat(paste0("  - Features:  ", nrow(seurat_obj), "\n"))
    
  }, error = function(e) {
    cat(paste0("✗ Error loading ", sample_id, ": ", e$message, "\n"))
  })
}

cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat(paste0("Successfully loaded ", length(seurat_objects), "/", nrow(sample_metadata), " samples\n"))
cat(paste(rep("=", 50), collapse = ""), "\n\n")

# ==============================================================================
# STEP 2: Calculate QC Metrics
# ==============================================================================

cat("Calculating QC metrics...\n")

for (sample_id in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[sample_id]]
  
  # Calculate mitochondrial percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Calculate ribosomal percentage
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
  
  # Calculate hemoglobin percentage (common in muscle biopsies)
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^HB[^(P)]")
  
  # Store back
  seurat_objects[[sample_id]] <- seurat_obj
}

# ==============================================================================
# STEP 3: QC Visualization (Before Filtering)
# ==============================================================================

cat("Generating QC plots...\n")

# Combine all samples for visualization
combined_qc <- data.frame()

for (sample_id in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[sample_id]]
  
  qc_data <- data.frame(
    sample_id = sample_id,
    participant = seurat_obj$participant,
    timepoint = seurat_obj$timepoint,
    nCount_RNA = seurat_obj$nCount_RNA,
    nFeature_RNA = seurat_obj$nFeature_RNA,
    percent.mt = seurat_obj$percent.mt,
    percent.ribo = seurat_obj$percent.ribo
  )
  
  combined_qc <- rbind(combined_qc, qc_data)
}

# Create output directories if they don't exist
dir.create("../results", showWarnings = FALSE)
dir.create("../figures", showWarnings = FALSE)
dir.create("../tables", showWarnings = FALSE)

# QC Violin plots
p1 <- ggplot(combined_qc, aes(x = sample_id, y = nFeature_RNA, fill = timepoint)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "Number of Genes Detected", y = "Genes per cell", x = "") +
  geom_hline(yintercept = c(200, 6000), linetype = "dashed", color = "red", alpha = 0.5) +
  scale_fill_manual(values = c("Pre" = "#4DBBD5FF", "Post" = "#E64B35FF"))

p2 <- ggplot(combined_qc, aes(x = sample_id, y = nCount_RNA, fill = timepoint)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "Total UMI Counts", y = "UMIs per cell", x = "") +
  geom_hline(yintercept = 500, linetype = "dashed", color = "red", alpha = 0.5) +
  scale_fill_manual(values = c("Pre" = "#4DBBD5FF", "Post" = "#E64B35FF"))

p3 <- ggplot(combined_qc, aes(x = sample_id, y = percent.mt, fill = timepoint)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "Mitochondrial Gene %", y = "% MT genes", x = "") +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red", alpha = 0.5) +
  scale_fill_manual(values = c("Pre" = "#4DBBD5FF", "Post" = "#E64B35FF"))

p4 <- ggplot(combined_qc, aes(x = sample_id, y = percent.ribo, fill = timepoint)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "Ribosomal Gene %", y = "% Ribosomal genes", x = "") +
  scale_fill_manual(values = c("Pre" = "#4DBBD5FF", "Post" = "#E64B35FF"))

# Save QC plots
pdf("../figures/01_QC_before_filtering.pdf", width = 16, height = 12)
print((p1 | p2) / (p3 | p4))
dev.off()

cat("✓ Saved:  ../figures/01_QC_before_filtering.pdf\n")

# Scatter plots showing relationships
p5 <- ggplot(combined_qc, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
  geom_point(alpha = 0.3, size = 0.5) +
  scale_color_viridis_c() +
  theme_classic() +
  facet_wrap(~sample_id, ncol = 4) +
  labs(title = "QC Metrics:  UMI vs Genes (colored by MT%)",
       x = "Total UMI counts",
       y = "Number of genes detected") +
  geom_hline(yintercept = c(200, 6000), linetype = "dashed", color = "red", alpha = 0.3) +
  geom_vline(xintercept = 500, linetype = "dashed", color = "red", alpha = 0.3)

pdf("../figures/02_QC_scatter.pdf", width = 14, height = 10)
print(p5)
dev.off()

cat("✓ Saved: ../figures/02_QC_scatter.pdf\n")

# QC statistics by timepoint
qc_stats <- combined_qc %>%
  group_by(timepoint) %>%
  summarise(
    n_cells = n(),
    mean_genes = mean(nFeature_RNA),
    median_genes = median(nFeature_RNA),
    mean_UMIs = mean(nCount_RNA),
    median_UMIs = median(nCount_RNA),
    mean_mt = mean(percent.mt),
    median_mt = median(percent.mt)
  )

cat("\nQC Statistics by Timepoint:\n")
print(qc_stats)

# ==============================================================================
# STEP 4: Apply QC Filters
# ==============================================================================

cat("\nApplying QC filters...\n")

# QC thresholds for skeletal muscle scRNA-seq
min_features <- 200
max_features <- 6000
max_mt_percent <- 20
min_counts <- 500

filter_summary <- data.frame()

for (sample_id in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[sample_id]]
  
  # Store pre-filter cell count
  cells_before <- ncol(seurat_obj)
  
  # Apply filters
  seurat_obj <- subset(seurat_obj,
                       subset = nFeature_RNA > min_features &
                         nFeature_RNA < max_features &
                         percent.mt < max_mt_percent &
                         nCount_RNA > min_counts)
  
  # Store post-filter cell count
  cells_after <- ncol(seurat_obj)
  percent_retained <- round(100 * cells_after / cells_before, 1)
  
  cat(paste0(sample_id, ":\n"))
  cat(paste0("  Before: ", cells_before, " cells\n"))
  cat(paste0("  After:  ", cells_after, " cells\n"))
  cat(paste0("  Retained: ", percent_retained, "%\n\n"))
  
  # Store for summary table
  filter_summary <- rbind(filter_summary, data.frame(
    sample_id = sample_id,
    cells_before = cells_before,
    cells_after = cells_after,
    percent_retained = percent_retained
  ))
  
  seurat_objects[[sample_id]] <- seurat_obj
}

write.csv(filter_summary, "../tables/filtering_summary.csv", row.names = FALSE)
cat("✓ Saved: ../tables/filtering_summary.csv\n")

# ==============================================================================
# STEP 5:  Merge All Samples
# ==============================================================================

cat("\nMerging all samples...\n")

# Merge all samples into one Seurat object
merged_seurat <- merge(x = seurat_objects[[1]],
                       y = seurat_objects[2:length(seurat_objects)],
                       add.cell.ids = names(seurat_objects),
                       project = "GSE214544_Exercise")

# Save the merged object
saveRDS(merged_seurat, "../results/merged_seurat_raw.rds")
cat("✓ Saved:  ../results/merged_seurat_raw.rds\n")

# ==============================================================================
# STEP 6: Generate Summary Statistics
# ==============================================================================

# Summary table
qc_summary <- data.frame(
  Sample = names(seurat_objects),
  Participant = sapply(seurat_objects, function(x) unique(x$participant)),
  Timepoint = sapply(seurat_objects, function(x) unique(x$timepoint)),
  Cells = sapply(seurat_objects, ncol),
  Mean_Genes = sapply(seurat_objects, function(x) round(mean(x$nFeature_RNA))),
  Median_Genes = sapply(seurat_objects, function(x) round(median(x$nFeature_RNA))),
  Mean_UMIs = sapply(seurat_objects, function(x) round(mean(x$nCount_RNA))),
  Median_UMIs = sapply(seurat_objects, function(x) round(median(x$nCount_RNA))),
  Mean_MT_Percent = sapply(seurat_objects, function(x) round(mean(x$percent.mt), 2)),
  Median_MT_Percent = sapply(seurat_objects, function(x) round(median(x$percent.mt), 2))
)

write.csv(qc_summary, "../tables/qc_summary.csv", row.names = FALSE)
cat("✓ Saved:  ../tables/qc_summary.csv\n")

# Print final summary
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("DATA LOADING COMPLETE - SUMMARY\n")
cat(paste(rep("=", 50), collapse = ""), "\n")
cat(paste0("Total samples loaded: ", length(seurat_objects), "\n"))
cat(paste0("Total cells after QC: ", ncol(merged_seurat), "\n"))
cat(paste0("Total features:  ", nrow(merged_seurat), "\n"))
cat(paste0("Pre-exercise cells: ", sum(merged_seurat$timepoint == "Pre"), "\n"))
cat(paste0("Post-exercise cells: ", sum(merged_seurat$timepoint == "Post"), "\n"))
cat(paste(rep("=", 50), collapse = ""), "\n\n")

cat("Sample Summary:\n")