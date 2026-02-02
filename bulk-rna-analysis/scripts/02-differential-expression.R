# ==============================================================================
# Script: 02-differential-expression.R
# Description: Differential expression analysis using DESeq2
# Author: Satish Patnaik
# Date: 2024
# ==============================================================================

# Load setup and data
source("scripts/00-setup.R")

# Load normalized data (from previous step)
# Alternatively, load raw counts and metadata directly
# counts <- load_counts("data/counts.txt")
# metadata <- load_metadata("data/metadata.txt")

# Create DESeq2 object
# ------------------------------------------------------------------------------
cat("Creating DESeq2 object...\n")

# Example metadata structure:
# sample_id  condition  batch
# sample1    control    1
# sample2    control    1
# sample3    treatment  1
# sample4    treatment  1

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition  # Adjust design formula as needed
)

# Pre-filtering: remove genes with very low counts
keep <- rowSums(counts(dds)) >= params$min_counts
dds <- dds[keep, ]

cat("Number of genes after filtering:", nrow(dds), "\n")

# Run DESeq2 analysis
# ------------------------------------------------------------------------------
cat("Running differential expression analysis...\n")
dds <- DESeq(dds)

# Get results for specific contrast
res <- results(
  dds,
  contrast = c("condition", "treatment", "control"),
  alpha = params$padj_threshold
)

# Order by adjusted p-value
res <- res[order(res$padj), ]

# Summarize results
cat("\n=== DESeq2 Results Summary ===\n")
summary(res)

# Extract significant genes
sig_genes <- subset(
  res,
  padj < params$padj_threshold & abs(log2FoldChange) > params$lfc_threshold
)

cat("\nNumber of significantly differentially expressed genes:", nrow(sig_genes), "\n")
cat("Upregulated:", sum(sig_genes$log2FoldChange > 0), "\n")
cat("Downregulated:", sum(sig_genes$log2FoldChange < 0), "\n")

# Save results
# ------------------------------------------------------------------------------

# Save all results
write.csv(
  as.data.frame(res),
  file = file.path(results_dir, "deseq2_results_all.csv"),
  row.names = TRUE
)

# Save significant genes only
write.csv(
  as.data.frame(sig_genes),
  file = file.path(results_dir, "deseq2_results_significant.csv"),
  row.names = TRUE
)

# Save normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(
  normalized_counts,
  file = file.path(results_dir, "normalized_counts.csv"),
  row.names = TRUE
)

# Create visualizations
# ------------------------------------------------------------------------------

# MA plot
cat("\nGenerating MA plot...\n")
pdf(file.path(figures_dir, "ma_plot.pdf"), width = 8, height = 6)
plotMA(res, ylim = c(-5, 5))
dev.off()

# Volcano plot
cat("Generating volcano plot...\n")
volcano_plot <- EnhancedVolcano(
  res,
  lab = rownames(res),
  x = 'log2FoldChange',
  y = 'padj',
  pCutoff = params$padj_threshold,
  FCcutoff = params$lfc_threshold,
  title = 'Differential Expression',
  subtitle = paste0('padj < ', params$padj_threshold, ', |log2FC| > ', params$lfc_threshold),
  legendPosition = 'right',
  legendLabSize = 10,
  legendIconSize = 4.0
)

ggsave(
  file.path(figures_dir, "volcano_plot.pdf"),
  plot = volcano_plot,
  width = 10,
  height = 8
)

# PCA plot
cat("Generating PCA plot...\n")
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  theme(legend.position = "top")

ggsave(
  file.path(figures_dir, "pca_plot.pdf"),
  plot = pca_plot,
  width = 8,
  height = 6
)

# Heatmap of top genes
cat("Generating heatmap...\n")
top_genes <- head(rownames(sig_genes), 50)
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)  # Center the data

pdf(file.path(figures_dir, "heatmap_top50.pdf"), width = 10, height = 12)
pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  annotation_col = as.data.frame(colData(dds)[, "condition", drop = FALSE]),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Top 50 Differentially Expressed Genes"
)
dev.off()

# Save R objects for future use
# ------------------------------------------------------------------------------
save(dds, res, sig_genes, file = file.path(results_dir, "deseq2_objects.RData"))

# Print session info
# ------------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), session_file)

cat("\n=== Analysis complete! ===\n")
cat("Results saved to:", results_dir, "\n")
cat("Figures saved to:", figures_dir, "\n")
