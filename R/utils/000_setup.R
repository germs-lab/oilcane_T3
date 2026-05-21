#####################################################################
# Common Functions and Utilities for Oilcane T3 Microbial Analysis
#
# This file contains shared functions, constants, and utilities that
# are used across multiple analysis scripts.
#
# Author: Bolívar Aponte Rolón
# Date: 2025-11-17
#####################################################################

# Package and Environment setup

invisible(
  lapply(
    c(
      "conflicted",
      "phyloseq",
      "vegan",
      "tidyverse",
      "data.table",
      "janitor",
      "microbiome",
      "metagMisc",
      "ggtext",
      "readr",
      "readxl",
      "stringr",
      "iNEXT",
      "ggpubr",
      "here"
    ),
    library,
    character.only = TRUE
  )
)


# List files and source each
list.files(here::here("R/functions"), pattern = "\\.R$", full.names = TRUE) %>%
  purrr::map(source)


# Load processed data objects if they exist
if (dir.exists(here::here("data/output/rdata"))) {
  list.files(
    here::here("data/output/rdata/phyloseq"),
    full.names = TRUE,
    recursive = FALSE,
    pattern = "\\.rda$"
  ) %>%
    purrr::walk(~ load(.x, envir = .GlobalEnv))
}


# Solve known conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("right_join", "dplyr")
conflict_prefer("intersect", "base")
conflict_prefer("setdiff", "Biostrings")
