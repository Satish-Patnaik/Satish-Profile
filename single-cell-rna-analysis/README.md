# Single Cell RNA-seq Analysis

## Overview
This directory contains workflows and scripts for single cell/nucleus RNA sequencing analysis, including clustering, cell type identification, and trajectory analysis.

## Directory Structure
```
single-cell-rna-analysis/
├── scripts/           # R scripts for analysis pipeline
├── data/             # Raw and processed data (not tracked by git)
├── results/          # Analysis outputs (not tracked by git)
├── figures/          # Generated plots (not tracked by git)
└── README.md         # This file
```

## Analysis Pipeline

### 1. Data Loading and QC
- Load 10x Genomics data or other formats
- Quality control metrics
- Filtering low-quality cells
- Doublet detection

### 2. Normalization and Feature Selection
- Log normalization or SCTransform
- Highly variable gene selection
- Scaling and centering

### 3. Dimensionality Reduction
- PCA
- UMAP or t-SNE
- Batch effect correction (if needed)

### 4. Clustering
- Graph-based clustering
- Cluster annotation
- Marker gene identification

### 5. Cell Type Identification
- Automated annotation tools
- Manual annotation using known markers
- Validation with reference datasets

### 6. Downstream Analysis
- Differential expression between clusters
- Trajectory analysis (pseudotime)
- Cell-cell communication
- Gene regulatory networks

## Requirements

### R Packages (Seurat Workflow)
```r
# Core scRNA-seq analysis
install.packages("Seurat")
install.packages("BiocManager")
BiocManager::install(c(
  "SingleR",
  "celldex",
  "scater",
  "scran",
  "SingleCellExperiment"
))

# Visualization
install.packages(c(
  "ggplot2",
  "patchwork",
  "RColorBrewer",
  "viridis"
))

# Utilities
install.packages(c(
  "tidyverse",
  "here",
  "Matrix"
))
```

### Python Alternative (Scanpy)
```bash
# Create conda environment
conda create -n scanpy python=3.9
conda activate scanpy

# Install packages
pip install scanpy python-igraph leidenalg
pip install numpy pandas matplotlib seaborn
```

## Usage

### Seurat Workflow (R)
```r
# Set working directory
setwd("single-cell-rna-analysis/")

# Load Seurat
library(Seurat)

# Run standard workflow
source("scripts/01-load-and-qc.R")
source("scripts/02-normalization.R")
source("scripts/03-clustering.R")
source("scripts/04-cell-type-annotation.R")
source("scripts/05-differential-expression.R")
```

### Scanpy Workflow (Python)
```python
import scanpy as sc
import pandas as pd

# Load data
adata = sc.read_10x_mtx('data/filtered_feature_bc_matrix/')

# Preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# Analysis
sc.pp.highly_variable_genes(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
```

## Input Data Formats

### 10x Genomics Output
```
data/
├── filtered_feature_bc_matrix/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
```

### Other Formats
- H5AD (AnnData/Scanpy)
- H5 (10x HDF5)
- RDS (Seurat object)
- Loom files
- CSV/TSV count matrices

## Template Scripts

See `scripts/` directory for template analysis scripts:
- `01-load-and-qc.R`: Data loading and quality control
- `02-normalization.R`: Normalization and scaling
- `03-clustering.R`: PCA, UMAP, and clustering
- `04-cell-type-annotation.R`: Identify cell types
- `05-differential-expression.R`: Find marker genes
- `06-advanced-analysis.R`: Trajectory, velocity, etc.

## QC Metrics

### Key Metrics to Monitor
- Number of genes per cell (nFeature_RNA)
- Total counts per cell (nCount_RNA)
- Mitochondrial percentage (percent.mt)
- Ribosomal percentage (percent.ribo)

### Typical Filtering Thresholds
```r
# Filter cells
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 200 & 
           nFeature_RNA < 6000 & 
           percent.mt < 15
)
```

## Visualization Examples

### UMAP with Cell Types
```r
DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE
)
```

### Feature Plots
```r
FeaturePlot(
  seurat_obj,
  features = c("PAX6", "RBPMS", "RLBP1", "RHO"),
  ncol = 2
)
```

### Dot Plot
```r
DotPlot(
  seurat_obj,
  features = marker_genes,
  group.by = "cell_type"
) + RotatedAxis()
```

## Cell Type Annotation

### Retinal Cell Markers (Example)
```r
retinal_markers <- list(
  "Photoreceptors" = c("RHO", "OPN1SW", "OPN1MW"),
  "Bipolar_cells" = c("VSX2", "GRIK1", "SCGN"),
  "Retinal_ganglion" = c("RBPMS", "POU4F2", "SNCG"),
  "Müller_glia" = c("RLBP1", "SLC1A3", "GLUL"),
  "Horizontal_cells" = c("ONECUT1", "LHX1", "CALB1"),
  "Amacrine_cells" = c("GAD1", "PAX6", "TFAP2A")
)
```

### Automated Annotation
```r
library(SingleR)
library(celldex)

# Use reference dataset
ref <- celldex::HumanPrimaryCellAtlasData()

# Run SingleR
predictions <- SingleR(
  test = GetAssayData(seurat_obj),
  ref = ref,
  labels = ref$label.main
)
```

## Best Practices

### Computational Efficiency
- Use sparse matrices for large datasets
- Consider downsampling for visualization
- Use batch processing for multiple samples
- Store intermediate objects

### Batch Effect Correction
```r
# Harmony integration
library(harmony)
seurat_obj <- RunHarmony(
  seurat_obj,
  group.by.vars = "batch"
)

# Or use Seurat integration
seurat_list <- SplitObject(seurat_obj, split.by = "batch")
features <- SelectIntegrationFeatures(seurat_list)
anchors <- FindIntegrationAnchors(seurat_list, anchor.features = features)
integrated <- IntegrateData(anchors)
```

### Reproducibility
```r
# Set random seed
set.seed(42)

# Document analysis parameters
params <- list(
  resolution = 0.8,
  n_pcs = 30,
  min_genes = 200,
  max_genes = 6000,
  max_mt_percent = 15
)

# Save session info
sessionInfo()
```

## Common Issues and Solutions

### High Mitochondrial Content
- May indicate dying cells
- Typical threshold: < 15-20%
- Adjust based on tissue type

### Doublets
```r
# Use DoubletFinder
library(DoubletFinder)

# Or Scrublet (Python)
import scrublet as scr
scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
```

### Low Cell Numbers
- Check QC thresholds
- Verify sequencing depth
- Consider cell viability issues

## Example Complete Analysis

See `scripts/seurat-complete-workflow.R` for a full example analyzing retinal tissue single-cell data.

## Resources

### Tutorials
- Seurat: https://satijalab.org/seurat/
- Scanpy: https://scanpy.readthedocs.io/
- Bioconductor: https://bioconductor.org/books/release/OSCA/

### Datasets
- Human Cell Atlas: https://www.humancellatlas.org/
- Retinal Atlas: https://singlecell.broadinstitute.org/

## Citation
If you use this code, please cite:
```
[Your publication details will go here]
```

## Contact
For questions about this analysis:
- Email: satishbiochem1@gmail.com
- GitHub: @Satish-Patnaik

## License
[Specify license after publication - e.g., MIT, GPL-3]

## Changelog
- 2024-XX-XX: Initial version
