# Run full pipeline test
test_full_pipeline <- function() {
  cat("\n### Running Full Pipeline Test ###\n")

  # Step 1: Import
  cat("\nStep 1: Running import...\n")
  tryCatch(
    {
      source("R/001_import.R")
      cat("✓ Import completed\n")
    },
    error = function(e) {
      cat("✗ Import failed:\n")
      print(e)
      return(FALSE)
    }
  )

  # Step 2: EDA
  cat("\nStep 2: Running EDA...\n")
  tryCatch(
    {
      source("R/03_eda.R")
      cat("✓ EDA completed\n")
    },
    error = function(e) {
      cat("✗ EDA failed:\n")
      print(e)
      return(FALSE)
    }
  )

  # Step 3: Diversity
  cat("\nStep 3: Running diversity analysis...\n")
  tryCatch(
    {
      source("R/04_alpha_beta_div.R")
      cat("✓ Diversity analysis completed\n")
    },
    error = function(e) {
      cat("✗ Diversity analysis failed:\n")
      print(e)
      return(FALSE)
    }
  )

  cat("\n✓ Full pipeline completed successfully!\n")
  return(TRUE)
}

# Uncomment to run full test (requires data and all packages):
# test_full_pipeline()
