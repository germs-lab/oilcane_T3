###################################################################
# Import of Data for Oilcane T3 Analysis
#
# This script imports phyloseq objects and metadata files from
# three timepoints of 16S rRNA sequencing data.
#
# Author: Adapted from Lorie's original script
# Date: 2025-11-17
###################################################################

# Load setup
source("R/utils/000_setup.R")

# Load required libraries
library(Biostrings) # for renaming ASVs

#--------------------------------------------------------
# Import Data Files
#--------------------------------------------------------

# Set paths relative to project root
input_dir <- here::here("data/input/16S_three_timepoints")

# Import metadata
samdf <- read.csv(file.path(input_dir, "psssu_metadata.csv"))

# Import Timepoint 3 data
seqtab.nochim <- readRDS(file.path(input_dir, "seqtab.nochim.rds"))
taxa <- readRDS(file.path(input_dir, "taxa.rds"))

# Import Jihoon's timepoints 1 and 2
psssu <- readRDS(file.path(input_dir, "ps-ssu.rds"))

#--------------------------------------------------------
# Create Phyloseq Objects
#--------------------------------------------------------

# Timepoint 3
row.names(samdf) <- samdf$Sequence_ID
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(taxa),
  sample_data(samdf)
)

# Add DNA sequences as refseq
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)

# Rename taxa to ASV IDs
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

cat("\n### Timepoint 3 Phyloseq Object ###\n")
print(ps)

# Add sample data to Jihoon's phyloseq object
# Note: This may need adjustment based on sample ID matching
psssu_with_metadata <- merge_phyloseq(psssu, sample_data(samdf))

cat("\n### Timepoints 1 & 2 Phyloseq Object ###\n")
print(psssu_with_metadata)

#--------------------------------------------------------
# Merge All Studies
#--------------------------------------------------------

ps_merged <- merge_phyloseq(psssu_with_metadata, ps)

cat("\n### Merged Phyloseq Object (All 3 Timepoints) ###\n")
print(ps_merged)

#--------------------------------------------------------
# Explore Merged Data
#--------------------------------------------------------

# Basic exploration
explore_phyloseq(ps_merged, name = "Merged All Timepoints")

# Check sample variables
cat("\nSample variables:\n")
print(sample_variables(ps_merged))

# Check taxonomic ranks
cat("\nTaxonomic ranks:\n")
print(rank_names(ps_merged))

#--------------------------------------------------------
# Save Processed Objects
#--------------------------------------------------------

# Create output directories if they don't exist
dir.create(
  here::here("data/output/processed/rdata"),
  recursive = TRUE,
  showWarnings = FALSE
)

# Save individual phyloseq objects
save(
  ps,
  file = here::here("data/output/processed/rdata/oilcane_t3_physeq.rda")
)

save(
  psssu_with_metadata,
  file = here::here("data/output/processed/rdata/oilcane_t1_t2_physeq.rda")
)

# Save merged object
save(
  ps_merged,
  file = here::here("data/output/processed/rdata/oilcane_merged_physeq.rda")
)

# Also save the merged object with a simpler name for easy loading
oilcane_physeq <- ps_merged
save(
  oilcane_physeq,
  file = here::here("data/output/processed/rdata/oilcane_physeq.rda")
)

cat("\n### Data Import Complete ###\n")
cat("Saved objects:\n")
cat("  - oilcane_t3_physeq.rda (Timepoint 3)\n")
cat("  - oilcane_t1_t2_physeq.rda (Timepoints 1 & 2)\n")
cat("  - oilcane_merged_physeq.rda (All timepoints)\n")
cat("  - oilcane_physeq.rda (Main analysis object)\n")
