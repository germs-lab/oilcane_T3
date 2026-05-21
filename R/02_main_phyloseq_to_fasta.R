###################################################################
# From Phyloseq to FASTA
# Author: Bolívar Aponte Rolón
# Date: 2025-10-23
###################################################################

source("R/utils/00_setup.R")


# Create output dir
fs::dir_create(
  here::here("data/output/sequences"),
  recurse = TRUE,
  showWarnings = FALSE
)

# Phyloseq to FASTA

refseq2fasta(
  main_oilcane_physeq,
  phyloseq_list = FALSE,
  out_dir = "data/output/sequences"
)
