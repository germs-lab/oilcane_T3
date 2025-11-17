# Testing Checklist for Oilcane T3 Analysis Pipeline

## Purpose
This document provides a testing checklist to validate the analysis pipeline after setup.

## Environment Setup Tests

### 1. Package Installation
```r
# Check if all required packages are available
required_packages <- c(
  "conflicted", "phyloseq", "vegan", "tidyverse", "data.table",
  "janitor", "microbiome", "metagMisc", "ggtext", "readr", "readxl",
  "stringr", "iNEXT", "ggpubr", "here", "Biostrings", "future",
  "future.apply", "parallelly"
)

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("Missing packages:\n")
  print(missing_packages)
  cat("\nInstall with:\n")
  cat('install.packages(c("', paste(missing_packages, collapse = '", "'), '"))\n', sep = '')
} else {
  cat("All required packages are installed!\n")
}
```

### 2. Directory Structure Test
```r
# Check if required directories exist
required_dirs <- c(
  "R/functions",
  "R/utils",
  "data/input/16S_three_timepoints",
  "data/output/processed/rdata",
  "data/output/processed/metadata"
)

for (dir in required_dirs) {
  if (dir.exists(here::here(dir))) {
    cat("✓", dir, "exists\n")
  } else {
    cat("✗", dir, "MISSING\n")
  }
}
```

### 3. Input Data Test
```r
# Check if input data files exist
input_files <- c(
  "data/input/16S_three_timepoints/psssu_metadata.csv",
  "data/input/16S_three_timepoints/seqtab.nochim.rds",
  "data/input/16S_three_timepoints/taxa.rds",
  "data/input/16S_three_timepoints/ps-ssu.rds"
)

for (file in input_files) {
  if (file.exists(here::here(file))) {
    cat("✓", file, "exists\n")
  } else {
    cat("✗", file, "MISSING\n")
  }
}
```

## Function Tests

### 4. Setup Script Test
```r
# Test if setup script loads without errors
tryCatch({
  source("R/utils/000_setup.R")
  cat("✓ Setup script loaded successfully\n")
}, error = function(e) {
  cat("✗ Setup script failed:\n")
  print(e)
})
```

### 5. Function Loading Test
```r
# Test if all functions are available
function_files <- list.files("R/functions", pattern = "\\.R$", full.names = TRUE)

for (file in function_files) {
  tryCatch({
    source(file)
    cat("✓", basename(file), "loaded\n")
  }, error = function(e) {
    cat("✗", basename(file), "failed:\n")
    print(e)
  })
}
```

### 6. Helper Function Tests
```r
# Test explore_phyloseq function
test_explore_phyloseq <- function() {
  tryCatch({
    # Load a test phyloseq object
    data(GlobalPatterns, package = "phyloseq")
    result <- explore_phyloseq(GlobalPatterns, name = "Test")
    cat("✓ explore_phyloseq works\n")
    return(TRUE)
  }, error = function(e) {
    cat("✗ explore_phyloseq failed:\n")
    print(e)
    return(FALSE)
  })
}

# Test calculate_alpha_diversity function
test_calculate_alpha <- function() {
  tryCatch({
    data(GlobalPatterns, package = "phyloseq")
    alpha <- calculate_alpha_diversity(GlobalPatterns)
    cat("✓ calculate_alpha_diversity works\n")
    return(TRUE)
  }, error = function(e) {
    cat("✗ calculate_alpha_diversity failed:\n")
    print(e)
    return(FALSE)
  })
}

# Test calculate_beta_diversity function
test_calculate_beta <- function() {
  tryCatch({
    data(GlobalPatterns, package = "phyloseq")
    beta <- calculate_beta_diversity(GlobalPatterns)
    cat("✓ calculate_beta_diversity works\n")
    return(TRUE)
  }, error = function(e) {
    cat("✗ calculate_beta_diversity failed:\n")
    print(e)
    return(FALSE)
  })
}

# Run all helper function tests
cat("\n### Testing Helper Functions ###\n")
test_explore_phyloseq()
test_calculate_alpha()
test_calculate_beta()
```

## Pipeline Tests

### 7. Import Script Test
```r
# Test data import (dry run check)
test_import <- function() {
  cat("\n### Testing Import Script ###\n")
  
  # Check if input files exist
  input_dir <- here::here("data/input/16S_three_timepoints")
  
  required_files <- c(
    "psssu_metadata.csv",
    "seqtab.nochim.rds",
    "taxa.rds",
    "ps-ssu.rds"
  )
  
  all_exist <- TRUE
  for (file in required_files) {
    if (!file.exists(file.path(input_dir, file))) {
      cat("✗ Missing:", file, "\n")
      all_exist <- FALSE
    }
  }
  
  if (all_exist) {
    cat("✓ All input files present\n")
    cat("  Ready to run: source('R/001_import.R')\n")
  } else {
    cat("✗ Cannot run import - missing input files\n")
  }
  
  return(all_exist)
}

test_import()
```

### 8. Full Pipeline Test (if data available)
```r
# Run full pipeline test
test_full_pipeline <- function() {
  cat("\n### Running Full Pipeline Test ###\n")
  
  # Step 1: Import
  cat("\nStep 1: Running import...\n")
  tryCatch({
    source("R/001_import.R")
    cat("✓ Import completed\n")
  }, error = function(e) {
    cat("✗ Import failed:\n")
    print(e)
    return(FALSE)
  })
  
  # Step 2: EDA
  cat("\nStep 2: Running EDA...\n")
  tryCatch({
    source("R/04_eda.R")
    cat("✓ EDA completed\n")
  }, error = function(e) {
    cat("✗ EDA failed:\n")
    print(e)
    return(FALSE)
  })
  
  # Step 3: Diversity
  cat("\nStep 3: Running diversity analysis...\n")
  tryCatch({
    source("R/05_alpha_beta_div.R")
    cat("✓ Diversity analysis completed\n")
  }, error = function(e) {
    cat("✗ Diversity analysis failed:\n")
    print(e)
    return(FALSE)
  })
  
  cat("\n✓ Full pipeline completed successfully!\n")
  return(TRUE)
}

# Uncomment to run full test (requires data and all packages):
# test_full_pipeline()
```

## Output Validation

### 9. Check Output Files
```r
# After running the pipeline, check if outputs were created
check_outputs <- function() {
  cat("\n### Checking Output Files ###\n")
  
  expected_outputs <- c(
    "data/output/processed/rdata/oilcane_physeq.rda",
    "data/output/processed/rdata/oilcane_t3_physeq.rda",
    "data/output/processed/rdata/oilcane_t1_t2_physeq.rda",
    "data/output/processed/rdata/oilcane_merged_physeq.rda",
    "data/output/processed/rdata/eda_results.rda",
    "data/output/processed/rdata/alpha_beta_results.rda"
  )
  
  for (file in expected_outputs) {
    if (file.exists(here::here(file))) {
      cat("✓", file, "exists\n")
    } else {
      cat("✗", file, "NOT FOUND\n")
    }
  }
}

# Run after pipeline execution:
# check_outputs()
```

### 10. Validate Output Objects
```r
# Validate the structure of output objects
validate_outputs <- function() {
  cat("\n### Validating Output Objects ###\n")
  
  # Load main phyloseq object
  load(here::here("data/output/processed/rdata/oilcane_physeq.rda"))
  
  # Basic validation
  cat("Phyloseq object validation:\n")
  cat("  Taxa:", ntaxa(oilcane_physeq), "\n")
  cat("  Samples:", nsamples(oilcane_physeq), "\n")
  cat("  Total reads:", sum(sample_sums(oilcane_physeq)), "\n")
  
  # Check components
  if (!is.null(otu_table(oilcane_physeq, errorIfNULL = FALSE))) {
    cat("  ✓ OTU table present\n")
  }
  if (!is.null(tax_table(oilcane_physeq, errorIfNULL = FALSE))) {
    cat("  ✓ Taxonomy table present\n")
  }
  if (!is.null(sample_data(oilcane_physeq, errorIfNULL = FALSE))) {
    cat("  ✓ Sample data present\n")
  }
  
  # Load and validate EDA results
  load(here::here("data/output/processed/rdata/eda_results.rda"))
  cat("\nEDA results validation:\n")
  cat("  Components:", paste(names(eda_results), collapse = ", "), "\n")
  
  # Load and validate diversity results
  load(here::here("data/output/processed/rdata/alpha_beta_results.rda"))
  cat("\nDiversity results validation:\n")
  cat("  Components:", paste(names(alpha_beta_results), collapse = ", "), "\n")
  
  cat("\n✓ All outputs validated successfully\n")
}

# Run after pipeline execution:
# validate_outputs()
```

## Summary Test Runner

### Run All Tests
```r
# Master test function
run_all_tests <- function() {
  cat("\n" , rep("=", 60), "\n")
  cat("OILCANE T3 PIPELINE TEST SUITE\n")
  cat(rep("=", 60), "\n\n")
  
  # Environment tests
  cat("PART 1: Environment Setup\n")
  cat(rep("-", 40), "\n")
  
  # Run tests...
  # (Include all test functions above)
  
  cat("\n", rep("=", 60), "\n")
  cat("TEST SUITE COMPLETE\n")
  cat(rep("=", 60), "\n")
}

# Uncomment to run all tests:
# run_all_tests()
```

## Manual Verification Steps

After running the automated tests, manually verify:

1. **Visual Inspection**
   - Check that plots are generated correctly
   - Verify rarefaction curves look reasonable
   - Examine alpha/beta diversity plots

2. **Data Quality**
   - Review read count distributions
   - Check Good's coverage values (should be >0.95)
   - Verify sample sizes match expectations

3. **Statistical Results**
   - Review alpha diversity summary statistics
   - Check ordination axis variance explained
   - Verify PERMANOVA results (if applicable)

4. **File Organization**
   - Confirm output files are properly named
   - Check that .gitignore prevents tracking output files
   - Verify documentation is complete

## Troubleshooting Guide

If tests fail, check:

1. **Package issues**: Run `renv::restore()` or reinstall packages
2. **Path issues**: Ensure working directory is project root
3. **Data issues**: Verify input files exist and are not corrupted
4. **Memory issues**: Reduce iNEXT parameters or increase memory limits
5. **Function issues**: Check for typos or missing dependencies

## Notes

- Some tests require the full dataset and all packages installed
- Tests can be run individually or as a suite
- Modify tests as needed for your specific use case
- Document any custom modifications or additional tests
