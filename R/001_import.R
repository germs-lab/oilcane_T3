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
samdf <- read.csv(file.path(input_dir, "psssu_metadata.csv")) %>%
  janitor::clean_names()
row.names(samdf) <- samdf$sequence_id

# Import Timepoint 3 data
seqtab.nochim <- readRDS(file.path(input_dir, "seqtab.nochim.rds"))
taxa <- readRDS(file.path(input_dir, "taxa.rds"))

# Import Jihoon's timepoints 1 and 2
psssu <- readRDS(file.path(input_dir, "ps-ssu.rds"))

#--------------------------------------------------------
# Create Phyloseq Objects
#--------------------------------------------------------

# Timepoint 3
oilcane_t3_physeq <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(taxa),
  sample_data(samdf)
)

# Add DNA sequences as refseq
dna <- Biostrings::DNAStringSet(taxa_names(oilcane_t3_physeq))
names(dna) <- taxa_names(oilcane_t3_physeq)
oilcane_t3_physeq <- merge_phyloseq(oilcane_t3_physeq, dna)

# Rename taxa to ASV IDs
taxa_names(oilcane_t3_physeq) <- paste0("ASV", seq(ntaxa(oilcane_t3_physeq)))

cat("\n### Timepoint 3 Phyloseq Object ###\n")
print(oilcane_t3_physeq)

# Add sample data to Jihoon's phyloseq object
# Note: This may need adjustment based on sample ID matching
oilcane_t1_t2_physeq <- merge_phyloseq(psssu, sample_data(samdf))

cat("\n### Timepoints 1 & 2 Phyloseq Object ###\n")
print(oilcane_t1_t2_physeq)

#--------------------------------------------------------
# Merge All Studies
#--------------------------------------------------------

main_oilcane_physeq <- merge_phyloseq(oilcane_t1_t2_physeq, oilcane_t3_physeq)

cat("\n### Merged Phyloseq Object (All 3 Timepoints) ###\n")
print(main_oilcane_physeq)

#--------------------------------------------------------
# Explore Merged Data
#--------------------------------------------------------

# Basic exploration
explore_phyloseq(main_oilcane_physeq, name = "Merged All Timepoints")

# Check sample variables
cat("\nSample variables:\n")
print(sample_variables(main_oilcane_physeq))

# Check taxonomic ranks
cat("\nTaxonomic ranks:\n")
print(rank_names(main_oilcane_physeq))

#--------------------------------------------------------
# Save Processed Objects
#--------------------------------------------------------

# Create output directories if they don't exist
dir.create(
  here::here("data/output/rdata"),
  recursive = TRUE,
  showWarnings = FALSE
)

# Save individual phyloseq objects

save(
  oilcane_t3_physeq,
  file = here::here("data/output/rdata/phyloseq/oilcane_t3_physeq.rda")
)

save(
  oilcane_t1_t2_physeq,
  file = here::here("data/output/rdata/phyloseq/oilcane_t1_t2_physeq.rda")
)

# Save merged object
save(
  main_oilcane_physeq,
  file = here::here("data/output/rdata/phyloseq/main_oilcane_physeq.rda")
)

cat("\n### Data Import Complete ###\n")
cat("Saved objects:\n")
cat("  - oilcane_t3_physeq.rda (Timepoint 3)\n")
cat("  - oilcane_t1_t2_physeq.rda (Timepoints 1 & 2)\n")
cat("  - main_oilcane_physeq.rda (All timepoints & Main analysis object)\n")
