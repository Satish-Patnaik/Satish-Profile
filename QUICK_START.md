# Quick Start Guide: Uploading Your RNA Analysis Code

## What Was Set Up

This repository now has a complete structure for managing your Bulk RNA and Single Cell RNA analysis codes while keeping them private until publication.

## Directory Structure Created

```
Satish-Profile/
├── .gitignore                          # Protects sensitive data from being committed
├── PRIVATE_CODE_GUIDE.md              # Comprehensive privacy management guide
├── README.md                          # Updated with links to RNA analysis
├── bulk-rna-analysis/
│   ├── README.md                      # Documentation for bulk RNA-seq
│   └── scripts/
│       ├── 00-setup.R                 # Setup and helper functions
│       └── 02-differential-expression.R  # DESeq2 analysis template
└── single-cell-rna-analysis/
    ├── README.md                      # Documentation for scRNA-seq
    └── scripts/
        ├── 01-load-and-qc.R          # Quality control template
        └── 03-clustering.R            # Clustering workflow template
```

## Next Steps: How to Upload Your R Studio Code

### Option 1: Keep Everything Private (RECOMMENDED)

1. **Create a separate private repository for your unpublished work:**
   - Go to https://github.com/new
   - Name it something like "RNA-Analysis-Private"
   - Select **"Private"** visibility
   - Create the repository

2. **Clone it and add your R Studio projects:**
   ```bash
   git clone https://github.com/YOUR-USERNAME/RNA-Analysis-Private.git
   cd RNA-Analysis-Private
   
   # Copy your R Studio projects
   cp -r ~/path/to/your/BulkRNA-Project/ .
   cp -r ~/path/to/your/scRNA-Project/ .
   
   # Commit and push
   git add .
   git commit -m "Add RNA analysis code"
   git push
   ```

3. **When ready to publish:**
   - Review and clean your code
   - Copy public-ready scripts to this Satish-Profile repository
   - Push to make them public

### Option 2: Use This Repository with .gitignore

1. **Create private folders:**
   ```bash
   cd bulk-rna-analysis
   mkdir private
   cd ../single-cell-rna-analysis
   mkdir private
   ```

2. **Add your R scripts to the private folders:**
   - Put unpublished code in `bulk-rna-analysis/private/`
   - Put scRNA code in `single-cell-rna-analysis/private/`
   - These folders are ignored by git (won't be committed)

3. **Move to public when ready:**
   ```bash
   # Move from private to public directory
   mv bulk-rna-analysis/private/my-analysis.R bulk-rna-analysis/scripts/
   git add bulk-rna-analysis/scripts/my-analysis.R
   git commit -m "Publish analysis code"
   git push
   ```

### Option 3: Use Git Branches

1. **Create a private branch:**
   ```bash
   git checkout -b private-analysis
   # Add all your code here
   git push origin private-analysis
   ```

2. **Keep branch private** - don't merge to main until ready

3. **When ready to publish:**
   ```bash
   git checkout main
   git cherry-pick <specific-commits>
   git push
   ```

## Important Reminders

### ❌ NEVER Commit These:
- Raw sequencing data (.fastq, .bam files)
- Count matrices from unpublished experiments  
- Patient/sample identifiers
- API keys or credentials
- Large result files

The `.gitignore` file protects you from accidentally committing these.

### ✅ Safe to Make Public (after publication):
- Analysis R scripts
- Documentation and methods
- Published figures
- Package requirements
- De-identified summary data

## Using the Template Scripts

The template R scripts in the `scripts/` directories are starting points:

1. **Copy and modify them** for your specific analysis
2. **Update file paths** to match your data
3. **Adjust parameters** for your experimental design
4. **Add your specific visualizations**

## Getting Help

- Read the comprehensive guide: `PRIVATE_CODE_GUIDE.md`
- Check directory READMEs for detailed instructions
- GitHub docs on private repos: https://docs.github.com/en/repositories

## Example Workflow

Here's how you might work with your code:

```bash
# In your private repository or private branch:
# 1. Develop and test your analysis
# 2. Keep all data and preliminary results private

# When paper is accepted:
# 3. Review and clean your code
# 4. Remove absolute paths and sensitive information
# 5. Copy to this public repository:

cd /path/to/Satish-Profile
cp ~/private-repo/bulk-rna/my-deseq2-analysis.R bulk-rna-analysis/scripts/
git add bulk-rna-analysis/scripts/my-deseq2-analysis.R
git commit -m "Add published DESeq2 analysis

Published in: Journal Name, 2024
DOI: 10.xxxx/xxxxx"
git push
```

## Questions?

If you have questions about:
- How to structure your specific analysis
- Which files to keep private
- How to organize multiple projects
- Data sharing requirements

Refer to `PRIVATE_CODE_GUIDE.md` for detailed answers.

---

**Remember:** The safest approach is Option 1 - keep everything in a private repository until publication, then selectively move code to this public profile repository.
