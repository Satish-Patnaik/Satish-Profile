# ==============================================================================
# Script: 03-clustering.R
# Description: PCA, UMAP, and clustering for scRNA-seq data
# Author: Satish Patnaik
# Date: 2024
# ==============================================================================

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Load normalized data
seurat_obj <- readRDS("results/seurat_normalized.rds")

# Define paths
output_dir <- "results"
figures_dir <- "figures"

# Dimensionality reduction
# ------------------------------------------------------------------------------
cat("Performing PCA...\n")

# Run PCA on scaled data
seurat_obj <- RunPCA(
  seurat_obj,
  features = VariableFeatures(seurat_obj),
  npcs = 50,
  verbose = FALSE
)

# Visualize PCA results
p1 <- DimPlot(seurat_obj, reduction = "pca")
p2 <- ElbowPlot(seurat_obj, ndims = 50)

combined_pca <- p1 + p2

ggsave(
  filename = file.path(figures_dir, "pca_plots.pdf"),
  plot = combined_pca,
  width = 12,
  height = 5
)

# Print top genes associated with PCs
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)

# Determine dimensionality
# ------------------------------------------------------------------------------
# Use elbow plot and JackStraw (optional, computationally intensive)

# Uncomment to run JackStraw (takes time)
# seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
# seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:20)
# JackStrawPlot(seurat_obj, dims = 1:20)

# Based on elbow plot, choose number of PCs
n_pcs <- 30  # Adjust based on your data

cat("Using", n_pcs, "principal components for downstream analysis\n")

# Clustering
# ------------------------------------------------------------------------------
cat("\nPerforming clustering...\n")

# Find neighbors
seurat_obj <- FindNeighbors(
  seurat_obj,
  dims = 1:n_pcs,
  verbose = FALSE
)

# Find clusters at multiple resolutions
resolutions <- c(0.4, 0.6, 0.8, 1.0, 1.2)

for (res in resolutions) {
  seurat_obj <- FindClusters(
    seurat_obj,
    resolution = res,
    verbose = FALSE
  )
}

# Use default resolution (0.8) for main analysis
Idents(seurat_obj) <- "RNA_snn_res.0.8"

cat("Number of clusters (res=0.8):", length(unique(Idents(seurat_obj))), "\n")

# UMAP
# ------------------------------------------------------------------------------
cat("\nRunning UMAP...\n")

seurat_obj <- RunUMAP(
  seurat_obj,
  dims = 1:n_pcs,
  verbose = FALSE
)

# Visualize UMAP
p_umap <- DimPlot(
  seurat_obj,
  reduction = "umap",
  label = TRUE,
  label.size = 6
) + NoLegend()

ggsave(
  filename = file.path(figures_dir, "umap_clusters.pdf"),
  plot = p_umap,
  width = 8,
  height = 6
)

# Compare different resolutions
# ------------------------------------------------------------------------------
cat("\nComparing clustering resolutions...\n")

resolution_plots <- lapply(resolutions, function(res) {
  DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by = paste0("RNA_snn_res.", res),
    label = TRUE,
    label.size = 4
  ) + ggtitle(paste("Resolution:", res)) + NoLegend()
})

combined_res <- wrap_plots(resolution_plots, ncol = 3)

ggsave(
  filename = file.path(figures_dir, "umap_resolutions.pdf"),
  plot = combined_res,
  width = 15,
  height = 10
)

# Feature plots for QC
# ------------------------------------------------------------------------------
cat("\nGenerating feature plots...\n")

qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

p_features <- FeaturePlot(
  seurat_obj,
  features = qc_features,
  ncol = 3
)

ggsave(
  filename = file.path(figures_dir, "umap_qc_features.pdf"),
  plot = p_features,
  width = 15,
  height = 5
)

# Save clustered object
# ------------------------------------------------------------------------------
cat("\nSaving clustered Seurat object...\n")
saveRDS(seurat_obj, file = file.path(output_dir, "seurat_clustered.rds"))

# Export cluster information
cluster_info <- data.frame(
  cell_id = colnames(seurat_obj),
  cluster = Idents(seurat_obj),
  seurat_obj@meta.data
)

write.csv(
  cluster_info,
  file = file.path(output_dir, "cluster_assignments.csv"),
  row.names = FALSE
)

# Cluster statistics
# ------------------------------------------------------------------------------
cluster_stats <- table(Idents(seurat_obj))
cluster_df <- data.frame(
  Cluster = names(cluster_stats),
  Count = as.numeric(cluster_stats),
  Percentage = round(100 * as.numeric(cluster_stats) / sum(cluster_stats), 2)
)

write.csv(
  cluster_df,
  file = file.path(output_dir, "cluster_statistics.csv"),
  row.names = FALSE
)

cat("\nCluster distribution:\n")
print(cluster_df)

# Print session info
# ------------------------------------------------------------------------------
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo_clustering.txt")
)

cat("\n=== Clustering complete! ===\n")
cat("Clustered object saved to:", file.path(output_dir, "seurat_clustered.rds"), "\n")
cat("UMAP plots saved to:", figures_dir, "\n")
