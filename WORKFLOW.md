# Analysis Workflow Guide

## Overview

This guide provides detailed instructions for running the oilcane_T3 analysis pipeline. The pipeline follows best practices from the germs-lab repositories and provides a structured approach to microbial community analysis.

## Prerequisites

### Required R Packages

The analysis requires the following R packages:

**Core packages:**
- `phyloseq` - Microbial community analysis framework
- `vegan` - Community ecology package
- `tidyverse` - Data manipulation and visualization
- `Biostrings` - DNA sequence manipulation

**Analysis packages:**
- `iNEXT` - Interpolation and extrapolation for diversity
- `microbiome` - Microbiome analytics
- `metagMisc` - Metagenomic utilities
- `ggpubr` - Publication-ready plots

**Utilities:**
- `conflicted` - Manage namespace conflicts
- `janitor` - Data cleaning
- `here` - Project-relative paths
- `future` - Parallel processing framework
- `future.apply` - Parallel apply functions

### Installation

If using `renv` (recommended):
```r
renv::restore()
```

Manual installation:
```r
# Install from CRAN
install.packages(c("tidyverse", "vegan", "janitor", "here", 
                   "conflicted", "ggpubr", "readr", "readxl", 
                   "stringr", "future", "future.apply"))

# Install from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("phyloseq", "Biostrings", "microbiome"))

# Install iNEXT
install.packages("iNEXT")

# Install metagMisc (if available)
# devtools::install_github("vmikk/metagMisc")
```

## Running the Analysis

### Step 0: Setup Environment

Open the R project:
```r
# In RStudio, open oilcane_T3.Rproj
# Or from the command line:
setwd("/path/to/oilcane_T3")
```

### Step 1: Data Import

Run the import script to load and process raw data:

```r
source("R/001_import.R")
```

**What this does:**
1. Loads metadata from `data/input/16S_three_timepoints/psssu_metadata.csv`
2. Loads ASV table and taxonomy for Timepoint 3
3. Loads phyloseq object for Timepoints 1 & 2
4. Creates individual phyloseq objects for each timepoint
5. Merges all timepoints into a single object
6. Saves processed objects to `data/output/rdata/`

**Expected outputs:**
- `data/output/rdata/oilcane_t3_physeq.rda`
- `data/output/rdata/oilcane_t1_t2_physeq.rda`
- `data/output/rdata/oilcane_merged_physeq.rda`
- `data/output/rdata/oilcane_physeq.rda` (main object)

**Console output should include:**
- Summary of each phyloseq object
- Number of samples and taxa
- Sample variables and taxonomic ranks

### Step 2: Exploratory Data Analysis

Run the EDA script:

```r
source("R/04_eda.R")
```

**What this does:**
1. Loads the merged phyloseq object
2. Calculates basic summary statistics
3. Analyzes read counts and Good's coverage
4. Generates visualization plots
5. Runs parallel iNEXT for rarefaction curves
6. Provides taxonomic composition overview

**Expected outputs:**
- `data/output/rdata/eda_results.rda`
- Multiple plots displayed in R graphics device:
  - Read count density plot
  - Read count histogram
  - Read count boxplot
  - Ranked read count distribution
  - Good's coverage vs sequencing depth
  - Rarefaction curves (iNEXT)
  - Phylum-level composition

**Notes:**
- iNEXT computation may take several minutes depending on sample size
- Default is set to 4 cores; adjust `nCores` parameter based on your system
- Large datasets may require increasing `future.globals.maxSize`

### Step 3: Alpha and Beta Diversity Analysis

Run the diversity analysis script:

```r
source("R/05_alpha_beta_div.R")
```

**What this does:**
1. Calculates alpha diversity metrics (Observed, Shannon, Simpson, InvSimpson)
2. Creates alpha diversity visualizations
3. Calculates beta diversity ordination (PCoA with Bray-Curtis)
4. Creates beta diversity plots
5. Provides templates for PERMANOVA analysis

**Expected outputs:**
- `data/output/rdata/alpha_beta_results.rda`
- Multiple plots displayed:
  - Individual alpha diversity metrics (4 plots)
  - Combined alpha diversity panel
  - PCoA ordination plot(s)

**Customization needed:**
- PERMANOVA analysis requires specifying grouping variables from your metadata
- Edit the script to uncomment and adjust the PERMANOVA section
- Grouping variables for colored ordination plots may need adjustment

### Step 4: Beta Diversity Partitioning (betapart)

Run the betapart analysis script:

```r
source("R/06_betapart_analysis.R")
```

**What this does:**
1. Decomposes Bray-Curtis dissimilarity into two components:
   - **Balanced variation**: Species replacement (turnover)
   - **Abundance gradient**: Differences due to abundance changes
2. Compares oilcane lineages vs wild type
3. Analyzes temporal patterns across T1, T2, T3
4. Explores genotype and sample type effects
5. Generates comprehensive visualizations

**Key Research Question:**
Do oilcane lineages have more balanced (turnover) Bray-Curtis dissimilarity than they have abundance gradients? This helps understand whether community differences are driven by species replacement or abundance changes.

**Expected outputs:**
- `data/output/rdata/betapart_results.rda`
- Multiple plots displayed:
  - Overall beta diversity partitioning
  - Oilcane vs Wild type comparison
  - Temporal patterns (T1, T2, T3)
  - Genotype-specific patterns
  - Sample type patterns
  - Combined oilcane x timepoint analysis

**Interpretation guide:**
- High balanced proportion: Communities differ mainly due to species replacement
- High gradient proportion: Communities differ mainly due to abundance changes
- Compare oilcane vs wild type to determine if transgenic modifications affect community assembly patterns

## Working with the Results

### Loading Saved Results

To load and work with saved results in a new R session:

```r
# Load setup (sources functions and loads packages)
source("R/utils/000_setup.R")

# Load specific results
load("data/output/rdata/oilcane_physeq.rda")
load("data/output/rdata/eda_results.rda")
load("data/output/rdata/alpha_beta_results.rda")
load("data/output/rdata/betapart_results.rda")

# Access results
print(oilcane_physeq)
names(eda_results)
names(alpha_beta_results)
names(betapart_results)
```

### Accessing Individual Components

```r
# Alpha diversity data frame
alpha_div <- alpha_beta_results$alpha_diversity
head(alpha_div)

# Alpha diversity plots
print(alpha_beta_results$alpha_plots$shannon)
print(alpha_beta_results$alpha_combined)

# Beta diversity ordination
beta_div <- alpha_beta_results$beta_diversity
print(beta_div$ordination)

# iNEXT results from EDA
inext_result <- eda_results$inext_result
print(eda_results$inext_plot)

# Betapart results
print(betapart_results$summary)
print(betapart_results$by_oilcane$group_comparison)
print(betapart_results$plots$combined)
```

## Troubleshooting

### Common Issues

1. **Package loading errors**
   - Ensure all required packages are installed
   - Try `renv::restore()` to sync package versions
   - Check for package conflicts with `conflicts()`

2. **Memory issues with iNEXT**
   - Reduce number of bootstrap replicates: `nboot = 50`
   - Increase memory limit: `options(future.globals.maxSize = 2000 * 1024^2)`
   - Reduce number of cores: `nCores = 2`

3. **Path issues**
   - Ensure working directory is set to project root
   - Use `here::here()` for all file paths
   - Check that input data exists: `list.files("data/input/16S_three_timepoints")`

4. **Missing metadata columns**
   - Check available columns: `sample_variables(oilcane_physeq)`
   - Modify grouping variables in scripts as needed
   - Some visualizations may need adjustment based on your metadata structure

### Validation Steps

To verify the pipeline is working correctly:

```r
# Check phyloseq object structure
source("R/utils/000_setup.R")
load("data/output/rdata/oilcane_physeq.rda")

# Basic checks
ntaxa(oilcane_physeq)  # Should return number of taxa
nsamples(oilcane_physeq)  # Should return number of samples
sample_variables(oilcane_physeq)  # List available metadata
rank_names(oilcane_physeq)  # List taxonomic ranks

# Check data quality
min(sample_sums(oilcane_physeq))  # Minimum reads per sample
max(sample_sums(oilcane_physeq))  # Maximum reads per sample
mean(sample_sums(oilcane_physeq))  # Mean reads per sample
```

## Next Steps

After running the basic pipeline, consider:

1. **Data transformation**
   - Rarefaction for equal sampling depth
   - Relative abundance transformation
   - Log transformation for visualization

2. **Additional analyses**
   - Differential abundance testing
   - Core microbiome analysis
   - Network analysis
   - Machine learning models

3. **Advanced visualization**
   - Heatmaps of abundant taxa
   - Stacked bar plots by sample groups
   - Network plots
   - Custom ordination plots with metadata

4. **Statistical testing**
   - PERMANOVA for beta diversity
   - Kruskal-Wallis for alpha diversity
   - Indicator species analysis
   - LEfSe analysis

5. **Beta diversity partitioning**
   - Run `R/06_betapart_analysis.R` for betapart analysis
   - Compare balanced (turnover) vs gradient (abundance) components

## References

- McMurdie, P.J. & Holmes, S. (2013) phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLoS ONE 8(4): e61217.
- Chao, A. et al. (2014) Rarefaction and extrapolation with Hill numbers: a framework for sampling and estimation in species diversity studies. Ecological Monographs 84(1): 45-67.
- Baselga, A. (2010) Partitioning the turnover and nestedness components of beta diversity. Global Ecology and Biogeography 19: 134-143.

## Support

For issues specific to this pipeline:
- Check the README.md for general information
- Review function documentation in R/functions/
- Consult reference repositories: germs_miscanthus and lightSABR

For package-specific issues:
- phyloseq: https://joey711.github.io/phyloseq/
- iNEXT: https://github.com/JohnsonHsieh/iNEXT
- vegan: https://github.com/vegandevs/vegan
- betapart: https://cran.r-project.org/package=betapart
