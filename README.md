# Oilcane T3 - Microbial Community Analysis

This repository contains the analysis pipeline for 16S rRNA sequencing data from the Oilcane T3 project, covering three timepoints of microbial community data.

## Repository Structure

```
oilcane_T3/
├── R/
│   ├── functions/           # Reusable R functions
│   │   ├── p_iNEXT.R       # Parallel iNEXT computation
│   │   ├── explore_phyloseq.R
│   │   ├── calculate_alpha_diversity.R
│   │   ├── calculate_beta_diversity.R
│   │   ├── calculate_betapart.R  # Beta diversity partitioning
│   │   └── analyze_read_counts.R
│   ├── utils/              # Utility scripts
│   │   └── 000_setup.R     # Common setup and library loading
│   ├── 001_import.R        # Data import and phyloseq creation
│   ├── 04_eda.R            # Exploratory data analysis
│   ├── 05_alpha_beta_div.R # Alpha and beta diversity analysis
│   └── 06_betapart_analysis.R  # Beta diversity partitioning
├── data/
│   ├── input/              # Input data files
│   │   └── 16S_three_timepoints/
│   └── output/             # Analysis outputs (gitignored)
│       └── processed/
│           ├── rdata/      # R data objects
│           └── metadata/   # Processed metadata
└── README.md
```

## Analysis Pipeline

### 1. Setup (`R/utils/000_setup.R`)

Common setup script that:
- Loads required packages (phyloseq, vegan, tidyverse, iNEXT, etc.)
- Sources all function files
- Loads processed data objects
- Resolves package conflicts

### 2. Data Import (`R/001_import.R`)

Imports and processes raw data:
- Loads metadata, ASV tables, and taxonomy from three timepoints
- Creates individual phyloseq objects for each timepoint
- Merges all timepoints into a single phyloseq object
- Saves processed objects to `data/output/rdata/`

**Output files:**
- `oilcane_t3_physeq.rda` - Timepoint 3 only
- `oilcane_t1_t2_physeq.rda` - Timepoints 1 & 2
- `oilcane_merged_physeq.rda` - All timepoints combined
- `oilcane_physeq.rda` - Main analysis object (alias for merged)

### 3. Exploratory Data Analysis (`R/04_eda.R`)

Performs initial data exploration:
- Basic phyloseq object exploration
- Read count analysis and visualization
- Good's coverage estimation
- Rarefaction curves using parallel iNEXT
- Taxonomic composition overview

**Output:** `data/output/rdata/eda_results.rda`

### 4. Alpha and Beta Diversity (`R/05_alpha_beta_div.R`)

Calculates and visualizes diversity metrics:
- Alpha diversity (Observed, Shannon, Simpson, InvSimpson)
- Alpha diversity plots and statistical comparisons
- Beta diversity ordination (PCoA with Bray-Curtis)
- PERMANOVA analysis (template provided)

**Output:** `data/output/rdata/alpha_beta_results.rda`

### 5. Beta Diversity Partitioning (`R/06_betapart_analysis.R`)

Uses the `betapart` package to decompose Bray-Curtis dissimilarity into:
- **Balanced variation (turnover)**: Differences due to species replacement
- **Abundance gradient (nestedness-like)**: Differences due to abundance changes

Key analyses include:
- Oilcane vs Wild type comparison
- Temporal patterns across T1, T2, T3
- Genotype-specific patterns (Lineages 1-5 vs WT)
- Sample type patterns (Roots, Soil, Leaves, Stalks)

**Output:** `data/output/rdata/betapart_results.rda`

**Key Question:** Do oilcane lineages have more balanced (turnover) Bray-Curtis dissimilarity compared to abundance gradients?

## Usage

### Running the Analysis

1. Clone the repository
2. Open the R project file `oilcane_T3.Rproj`
3. Restore the R environment (if using renv):
   ```r
   renv::restore()
   ```
4. Run the analysis scripts in order:
   ```r
   source("R/001_import.R")    # Import and process data
   source("R/04_eda.R")        # Exploratory analysis
   source("R/05_alpha_beta_div.R")  # Diversity analysis
   source("R/06_betapart_analysis.R")  # Beta diversity partitioning
   ```

### Key Functions

#### Exploration
- `explore_phyloseq(physeq, name)` - Print summary statistics
- `analyze_read_counts(physeq)` - Extract and analyze read counts

#### Diversity Analysis
- `calculate_alpha_diversity(physeq, measures)` - Calculate alpha diversity
- `calculate_beta_diversity(physeq, method, distance)` - Calculate beta diversity ordination
- `calculate_betapart(physeq, index.family)` - Partition beta diversity into balanced and gradient components
- `calculate_betapart_by_group(physeq, group_var)` - Compare betapart across groups
- `p_iNEXT(x, q, nCores, ...)` - Parallel iNEXT rarefaction curves

## Data

Input data is located in `data/input/16S_three_timepoints/`:
- `psssu_metadata.csv` - Sample metadata
- `seqtab.nochim.rds` - ASV table (Timepoint 3)
- `taxa.rds` - Taxonomy assignments (Timepoint 3)
- `ps-ssu.rds` - Phyloseq object (Timepoints 1 & 2)
- `merged.rds` - Merged phyloseq object
- `asv_merge.fasta` - Merged ASV sequences

## Dependencies

Main R packages:
- phyloseq - Microbial community analysis
- vegan - Community ecology analysis
- tidyverse - Data manipulation and visualization
- iNEXT - Rarefaction and extrapolation
- betapart - Beta diversity partitioning
- microbiome - Microbiome analysis tools
- metagMisc - Metagenomic analysis utilities
- ggpubr - Publication-ready plots
- Biostrings - DNA sequence manipulation

## Credits

Analysis pipeline adapted from:
- [germs-lab/germs_miscanthus](https://github.com/germs-lab/germs_miscanthus)
- [germs-lab/lightSABR](https://github.com/germs-lab/lightSABR)

Author: Bolívar Aponte Rolón  
Date: 2025-11-17
