# ==============================================================================
# Script: 01-load-and-qc.R
# Description: Load 10x data and perform quality control for scRNA-seq
# Author: Satish Patnaik
# Date: 2024
# ==============================================================================

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Set options
options(stringsAsFactors = FALSE)
set.seed(42)

# Define paths
data_path <- "data/filtered_feature_bc_matrix/"  # 10x output directory
output_dir <- "results"
figures_dir <- "figures"

# Create output directories
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

# Load data
# ------------------------------------------------------------------------------
cat("Loading 10x Genomics data...\n")

# Load 10x data
counts <- Read10X(data.dir = data_path)

# Create Seurat object
seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "scRNA_analysis",
  min.cells = 3,      # Include genes detected in at least 3 cells
  min.features = 200  # Include cells with at least 200 detected genes
)

cat("Initial dataset:\n")
cat("  Number of cells:", ncol(seurat_obj), "\n")
cat("  Number of genes:", nrow(seurat_obj), "\n")

# Calculate QC metrics
# ------------------------------------------------------------------------------
cat("\nCalculating QC metrics...\n")

# Calculate mitochondrial percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Calculate ribosomal percentage
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

# View QC metrics
cat("\nQC metrics summary:\n")
summary(seurat_obj@meta.data$nFeature_RNA)
summary(seurat_obj@meta.data$nCount_RNA)
summary(seurat_obj@meta.data$percent.mt)

# Visualize QC metrics
# ------------------------------------------------------------------------------
cat("\nGenerating QC plots...\n")

# Violin plots for QC metrics
p1 <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

ggsave(
  filename = file.path(figures_dir, "qc_violin_plots.pdf"),
  plot = p1,
  width = 12,
  height = 4
)

# Scatter plots to visualize relationships
p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
p3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

combined_scatter <- p2 + p3

ggsave(
  filename = file.path(figures_dir, "qc_scatter_plots.pdf"),
  plot = combined_scatter,
  width = 12,
  height = 5
)

# Filter cells based on QC metrics
# ------------------------------------------------------------------------------
cat("\nFiltering cells...\n")

# Define filtering thresholds
min_features <- 200
max_features <- 6000
max_mt_percent <- 15

# Apply filters
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > min_features & 
           nFeature_RNA < max_features & 
           percent.mt < max_mt_percent
)

cat("After filtering:\n")
cat("  Number of cells:", ncol(seurat_obj), "\n")
cat("  Number of genes:", nrow(seurat_obj), "\n")

# Save filtered data
# ------------------------------------------------------------------------------
cat("\nSaving filtered Seurat object...\n")
saveRDS(seurat_obj, file = file.path(output_dir, "seurat_filtered.rds"))

# Export metadata
write.csv(
  seurat_obj@meta.data,
  file = file.path(output_dir, "cell_metadata_filtered.csv"),
  row.names = TRUE
)

# Generate QC report
# ------------------------------------------------------------------------------
qc_report <- data.frame(
  Metric = c(
    "Initial cells",
    "Cells after QC",
    "Cells removed",
    "Initial genes",
    "Genes after filtering",
    "Mean features per cell",
    "Median features per cell",
    "Mean UMI per cell",
    "Median UMI per cell",
    "Mean MT%",
    "Median MT%"
  ),
  Value = c(
    ncol(CreateSeuratObject(counts, min.cells = 3, min.features = 200)),
    ncol(seurat_obj),
    ncol(CreateSeuratObject(counts, min.cells = 3, min.features = 200)) - ncol(seurat_obj),
    nrow(seurat_obj),
    nrow(seurat_obj),
    round(mean(seurat_obj$nFeature_RNA), 2),
    median(seurat_obj$nFeature_RNA),
    round(mean(seurat_obj$nCount_RNA), 2),
    median(seurat_obj$nCount_RNA),
    round(mean(seurat_obj$percent.mt), 2),
    round(median(seurat_obj$percent.mt), 2)
  )
)

write.csv(
  qc_report,
  file = file.path(output_dir, "qc_summary.csv"),
  row.names = FALSE
)

# Print session info
# ------------------------------------------------------------------------------
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo_qc.txt")
)

cat("\n=== QC complete! ===\n")
cat("Filtered Seurat object saved to:", file.path(output_dir, "seurat_filtered.rds"), "\n")
cat("QC plots saved to:", figures_dir, "\n")
