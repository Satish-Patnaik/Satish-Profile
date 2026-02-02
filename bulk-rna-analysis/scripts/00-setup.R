# ==============================================================================
# Script: 00-setup.R
# Description: Setup script for bulk RNA-seq analysis
# Author: Satish Patnaik
# Date: 2024
# ==============================================================================

# Load required libraries
# ------------------------------------------------------------------------------

# Core analysis packages
library(DESeq2)
library(edgeR)

# Annotation and databases
library(AnnotationDbi)
library(org.Hs.eg.db)  # Human annotation, change if working with other species

# Visualization
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)

# Data manipulation
library(tidyverse)
library(here)

# Set global options
# ------------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
set.seed(123)  # For reproducibility

# Define paths (relative paths using here package)
# ------------------------------------------------------------------------------
data_dir <- here("data")
results_dir <- here("results")
figures_dir <- here("figures")

# Create directories if they don't exist
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

# Helper functions
# ------------------------------------------------------------------------------

#' Load count matrix from file
#' @param file_path Path to count matrix file
#' @return Count matrix with genes as rows and samples as columns
load_counts <- function(file_path) {
  counts <- read.table(file_path, header = TRUE, row.names = 1, sep = "\t")
  return(as.matrix(counts))
}

#' Load sample metadata
#' @param file_path Path to metadata file
#' @return Data frame with sample information
load_metadata <- function(file_path) {
  metadata <- read.table(file_path, header = TRUE, sep = "\t")
  return(metadata)
}

#' Save plot to file
#' @param plot_obj ggplot object
#' @param filename Output filename
#' @param width Plot width in inches
#' @param height Plot height in inches
save_plot <- function(plot_obj, filename, width = 8, height = 6) {
  ggsave(
    filename = file.path(figures_dir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300
  )
}

# Analysis parameters
# ------------------------------------------------------------------------------
params <- list(
  padj_threshold = 0.05,
  lfc_threshold = 1,
  min_counts = 10,
  norm_method = "DESeq2"
)

# Print session info
# ------------------------------------------------------------------------------
cat("=== Setup Complete ===\n")
cat("R version:", R.version.string, "\n")
cat("Working directory:", getwd(), "\n")
cat("Results will be saved to:", results_dir, "\n")
cat("Figures will be saved to:", figures_dir, "\n")

# Session info will be saved at the end of analysis
session_file <- file.path(results_dir, "sessionInfo.txt")
