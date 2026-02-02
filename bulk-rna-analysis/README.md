# Bulk RNA-seq Analysis

## Overview
This directory contains workflows and scripts for bulk RNA sequencing analysis, focusing on differential gene expression and pathway analysis in ocular tissues.

## Directory Structure
```
bulk-rna-analysis/
├── scripts/           # R scripts for analysis pipeline
├── data/             # Raw and processed data (not tracked by git)
├── results/          # Analysis outputs (not tracked by git)
├── figures/          # Generated plots (not tracked by git)
└── README.md         # This file
```

## Analysis Pipeline

### 1. Data Preprocessing
- Quality control of raw reads
- Alignment and counting
- Normalization

### 2. Differential Expression
- DESeq2 or edgeR analysis
- Multiple testing correction
- Result filtering and annotation

### 3. Visualization
- PCA plots
- Volcano plots
- Heatmaps
- MA plots

### 4. Pathway Analysis
- GSEA
- GO enrichment
- KEGG pathway analysis

## Requirements

### R Packages
```r
# Core analysis
install.packages("BiocManager")
BiocManager::install(c(
  "DESeq2",
  "edgeR",
  "limma",
  "clusterProfiler",
  "org.Hs.eg.db",
  "EnhancedVolcano"
))

# Visualization
install.packages(c(
  "ggplot2",
  "pheatmap",
  "RColorBrewer",
  "ggrepel"
))

# Utilities
install.packages(c(
  "tidyverse",
  "here",
  "data.table"
))
```

### System Requirements
- R >= 4.0.0
- RStudio (recommended)
- Sufficient RAM (8GB minimum, 16GB+ recommended)

## Usage

### Quick Start
```r
# Set working directory
setwd("bulk-rna-analysis/")

# Load required libraries
source("scripts/00-setup.R")

# Run preprocessing
source("scripts/01-preprocessing.R")

# Differential expression analysis
source("scripts/02-differential-expression.R")

# Generate visualizations
source("scripts/03-visualization.R")

# Pathway analysis
source("scripts/04-pathway-analysis.R")
```

### Input Data Format
Expected input is a count matrix with:
- Rows: Genes
- Columns: Samples
- Format: Tab-delimited or CSV

```
GeneID    Sample1    Sample2    Sample3    Sample4
ENSG001   150        200        180        220
ENSG002   50         60         55         65
...
```

## Template Scripts

See `scripts/` directory for template R scripts:
- `00-setup.R`: Package loading and helper functions
- `01-preprocessing.R`: QC and normalization
- `02-differential-expression.R`: DESeq2 analysis
- `03-visualization.R`: Plotting functions
- `04-pathway-analysis.R`: GSEA and enrichment

## Best Practices

### Data Management
- Keep raw data separate from processed data
- Use relative paths (not absolute paths)
- Document data sources and versions
- Use consistent naming conventions

### Code Organization
- One analysis step per script
- Comment your code thoroughly
- Include session info at end of analysis
- Use version control (git)

### Reproducibility
```r
# Always include session info
sessionInfo()

# Use set.seed() for reproducible results
set.seed(123)

# Document parameters
params <- list(
  padj_cutoff = 0.05,
  lfc_threshold = 1,
  norm_method = "DESeq2"
)
```

## Example Analysis

### DESeq2 Differential Expression
```r
# Load DESeq2
library(DESeq2)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition
)

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(dds, contrast = c("condition", "treatment", "control"))

# Filter significant genes
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
```

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
