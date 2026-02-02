# Guide: Managing Private RNA Analysis Code Until Publication

## Overview
This guide explains how to upload and manage your Bulk RNA-seq and Single Cell RNA-seq analysis codes from R Studio while keeping them private until you're ready to publish.

## Quick Start: Steps to Keep Code Private

### Option 1: Using Private Repository (RECOMMENDED)
1. **Create a separate private repository for unpublished work:**
   ```bash
   # On GitHub:
   # - Go to https://github.com/new
   # - Name it (e.g., "RNA-Analysis-Private")
   # - Select "Private" visibility
   # - Create repository
   ```

2. **Clone and work in the private repo:**
   ```bash
   git clone https://github.com/YOUR-USERNAME/RNA-Analysis-Private.git
   cd RNA-Analysis-Private
   # Add your R scripts and analysis files
   git add .
   git commit -m "Add RNA analysis scripts"
   git push
   ```

3. **When ready to publish, transfer to public repo:**
   ```bash
   # Copy specific files you want to make public
   cp -r bulk-rna-analysis/ /path/to/Satish-Profile/
   cd /path/to/Satish-Profile/
   git add bulk-rna-analysis/
   git commit -m "Add published bulk RNA-seq analysis"
   git push
   ```

### Option 2: Using .gitignore in This Repository
1. **Create private directories (already in .gitignore):**
   ```bash
   mkdir -p bulk-rna-analysis/private
   mkdir -p single-cell-rna-analysis/private
   ```

2. **Store sensitive code in private folders:**
   - Put unpublished scripts in `/private/` subdirectories
   - The `.gitignore` file will prevent them from being committed

3. **Move to public when ready:**
   ```bash
   # Move files from private to public directory
   mv bulk-rna-analysis/private/my-analysis.R bulk-rna-analysis/
   git add bulk-rna-analysis/my-analysis.R
   git commit -m "Publish bulk RNA-seq analysis"
   git push
   ```

### Option 3: Using Git Branches (For This Repository)
1. **Create a private branch:**
   ```bash
   git checkout -b private-analysis
   # Add all your analysis code
   git add .
   git commit -m "Add private RNA analysis code"
   git push origin private-analysis
   ```

2. **Keep the branch private:**
   - Do NOT create a Pull Request or merge to main
   - Only you can see this branch (unless others have repository access)

3. **When ready to publish:**
   ```bash
   # Cherry-pick specific commits to main
   git checkout main
   git cherry-pick <commit-hash>
   git push origin main
   ```

## What to Keep Private

### Always Keep Private:
- ❌ Raw sequencing data (.fastq, .bam files)
- ❌ Patient/sample identifiers
- ❌ Unpublished results and figures
- ❌ Preliminary analyses
- ❌ Count matrices from unpublished experiments
- ❌ API keys, credentials, database connections

### Can Make Public (After Publication):
- ✅ Analysis scripts and workflows
- ✅ Documentation and methods
- ✅ Visualization code
- ✅ Package requirements and environment files
- ✅ Published figures and supplementary materials
- ✅ Processed, de-identified data (with appropriate permissions)

## Best Practices for R Studio Users

### 1. Project Organization
```
your-analysis-project/
├── data/              # Raw data (keep private)
│   ├── bulk/
│   └── single-cell/
├── scripts/           # Analysis scripts
│   ├── 01-preprocessing.R
│   ├── 02-quality-control.R
│   ├── 03-differential-expression.R
│   └── 04-visualization.R
├── results/           # Output (keep private until published)
├── figures/           # Plots (keep private until published)
└── README.md          # Documentation
```

### 2. Remove Sensitive Paths from Code
Before publishing, replace absolute paths:
```r
# DON'T publish this:
data <- read.csv("/Users/satish/Documents/private-project/data.csv")

# DO publish this:
data <- read.csv("data/expression_data.csv")
# Or use relative paths with here package
library(here)
data <- read.csv(here("data", "expression_data.csv"))
```

### 3. Document Your Methods
Create clear documentation for reproducibility:
```r
# Analysis: Bulk RNA-seq differential expression
# Author: Your Name
# Date: 2024-XX-XX
# Description: DESeq2 analysis of retinal tissue samples
#
# Requirements:
# - R >= 4.0
# - DESeq2, ggplot2, pheatmap
#
# Input: Raw count matrix
# Output: Differential expression results, QC plots
```

### 4. Version Control Before Publishing
```bash
# Before publishing, review what you're about to share:
git status
git diff

# Only add files you want to publish:
git add scripts/01-preprocessing.R
git add scripts/02-analysis.R
git add README.md

# Don't accidentally commit everything:
# git add .  # Be careful with this!
```

## GitHub Repository Settings

### Making Repository Private
1. Go to your repository on GitHub
2. Click "Settings"
3. Scroll to "Danger Zone"
4. Click "Change repository visibility"
5. Select "Make private"

### Making Repository Public (When Ready to Publish)
1. Go to repository Settings
2. Scroll to "Danger Zone"
3. Click "Change repository visibility"
4. Select "Make public"
5. Type repository name to confirm

## Pre-Publication Checklist

Before making code public, verify:
- [ ] Remove all patient/sample identifiers
- [ ] Remove absolute file paths
- [ ] Remove API keys and credentials
- [ ] Update README with clear documentation
- [ ] Test code runs on fresh environment
- [ ] Remove large data files (use data repositories like GEO, SRA)
- [ ] Add license (e.g., MIT, GPL-3)
- [ ] Get co-author approval
- [ ] Check journal/funding agency data sharing policies
- [ ] Add citation to published paper (if available)

## Recommended Tools

### For Large Data Management
- **Git LFS**: Store large files (not recommended for truly private data)
- **GEO/SRA**: Public repositories for published omics data
- **Zenodo/Figshare**: DOI-citable data repositories

### For Reproducible Environments
- **renv**: R package version management
- **conda**: Environment management
- **Docker**: Containerized analysis environments

### For Documentation
- **R Markdown**: Literate programming and reports
- **Jupyter Notebooks**: Interactive analysis documentation
- **pkgdown**: R package documentation websites

## Example Workflow

### Upload to Private Repository
```bash
# 1. Create private repo on GitHub (RNA-Analysis-Private)
# 2. Clone it locally
git clone https://github.com/YOUR-USERNAME/RNA-Analysis-Private.git
cd RNA-Analysis-Private

# 3. Add your R Studio project
cp -r ~/RStudio-Projects/BulkRNA-Analysis/ .
cp -r ~/RStudio-Projects/scRNA-Analysis/ .

# 4. Commit and push
git add .
git commit -m "Add RNA-seq analysis pipelines"
git push
```

### Publish to Public Repository (After Paper Acceptance)
```bash
# 1. Review and clean code
cd RNA-Analysis-Private
# Remove sensitive data, update paths, add documentation

# 2. Clone public profile repo
cd ..
git clone https://github.com/Satish-Patnaik/Satish-Profile.git
cd Satish-Profile

# 3. Copy published analysis
mkdir -p rna-seq-analyses
cp -r ../RNA-Analysis-Private/BulkRNA-Analysis/ rna-seq-analyses/
cp -r ../RNA-Analysis-Private/scRNA-Analysis/ rna-seq-analyses/

# 4. Commit to public repo
git add rna-seq-analyses/
git commit -m "Add published RNA-seq analysis pipelines

Published in: Journal Name, Year
DOI: 10.xxxx/xxxxx"
git push
```

## Getting Help

If you need assistance:
1. Check GitHub documentation: https://docs.github.com
2. R-specific version control: https://happygitwithr.com/
3. Ask in bioinformatics communities (Biostars, SEQanswers)

## Additional Resources

- **GitHub Private Repositories**: https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/managing-repository-settings/setting-repository-visibility
- **Git Ignore Files**: https://git-scm.com/docs/gitignore
- **R Studio and Git**: https://support.rstudio.com/hc/en-us/articles/200532077-Version-Control-with-Git-and-SVN
- **Reproducible Research**: https://ropensci.github.io/reproducibility-guide/
