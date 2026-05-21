#' Explore Phyloseq Object
#'
#' @description
#' Prints basic information about a phyloseq object including
#' number of taxa, samples, and total reads.
#'
#' @param physeq A phyloseq object
#' @param name Optional name for the phyloseq object (for printing)
#'
#' @return Invisibly returns a list with summary statistics
#'
#' @export
explore_phyloseq <- function(physeq, name = NULL) {
  if (!is.null(name)) {
    cat("\n### Exploring:", name, "###\n")
  }
  
  cat("Number of taxa:", ntaxa(physeq), "\n")
  cat("Number of samples:", nsamples(physeq), "\n")
  cat("Total reads:", sum(sample_sums(physeq)), "\n")
  cat("Min reads per sample:", min(sample_sums(physeq)), "\n")
  cat("Max reads per sample:", max(sample_sums(physeq)), "\n")
  cat("Mean reads per sample:", round(mean(sample_sums(physeq)), 2), "\n")
  cat("Median reads per sample:", median(sample_sums(physeq)), "\n")
  
  if (!is.null(tax_table(physeq, errorIfNULL = FALSE))) {
    cat("Taxonomic ranks:", paste(rank_names(physeq), collapse = ", "), "\n")
  }
  
  if (!is.null(sample_data(physeq, errorIfNULL = FALSE))) {
    cat("Sample variables:", paste(sample_variables(physeq), collapse = ", "), "\n")
  }
  
  invisible(list(
    n_taxa = ntaxa(physeq),
    n_samples = nsamples(physeq),
    total_reads = sum(sample_sums(physeq)),
    min_reads = min(sample_sums(physeq)),
    max_reads = max(sample_sums(physeq)),
    mean_reads = mean(sample_sums(physeq)),
    median_reads = median(sample_sums(physeq))
  ))
}


#' Explore Nested Phyloseq List
#'
#' @description
#' Explores a nested list structure containing phyloseq objects
#'
#' @param physeq_list A nested list of phyloseq objects
#'
#' @return A data frame with summary statistics for all phyloseq objects
#'
#' @export
explore_nested_phyloseq <- function(physeq_list) {
  purrr::imap_dfr(
    physeq_list,
    function(project_list, project_name) {
      purrr::imap_dfr(
        project_list,
        function(physeq_obj, physeq_name) {
          data.frame(
            project = project_name,
            physeq_name = physeq_name,
            n_taxa = ntaxa(physeq_obj),
            n_samples = nsamples(physeq_obj),
            total_reads = sum(sample_sums(physeq_obj)),
            min_reads = min(sample_sums(physeq_obj)),
            max_reads = max(sample_sums(physeq_obj)),
            mean_reads = mean(sample_sums(physeq_obj)),
            median_reads = median(sample_sums(physeq_obj))
          )
        }
      )
    }
  )
}
