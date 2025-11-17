#' Analyze Read Counts
#'
#' @description
#' Extracts and analyzes read counts from phyloseq object(s)
#'
#' @param physeq A phyloseq object or nested list of phyloseq objects
#'
#' @return A data frame with read count statistics
#'
#' @export
analyze_read_counts <- function(physeq) {
  if ("phyloseq" %in% class(physeq)) {
    # Single phyloseq object
    sample_sums_df <- data.frame(
      sample_id = sample_names(physeq),
      n_seqs = sample_sums(physeq)
    )
    
    # Add sample data if available
    if (!is.null(sample_data(physeq, errorIfNULL = FALSE))) {
      sample_df <- data.frame(sample_data(physeq))
      sample_sums_df <- cbind(sample_sums_df, sample_df)
    }
    
    # Calculate Good's coverage
    otu_mat <- as(otu_table(physeq), "matrix")
    if (!taxa_are_rows(physeq)) {
      otu_mat <- t(otu_mat)
    }
    
    singletons <- colSums(otu_mat == 1)
    goods_coverage <- 1 - (singletons / sample_sums_df$n_seqs)
    sample_sums_df$goods <- goods_coverage
    
    return(sample_sums_df)
    
  } else if (is.list(physeq)) {
    # Nested list structure
    purrr::imap_dfr(
      physeq,
      function(project_list, project_name) {
        purrr::imap_dfr(
          project_list,
          function(physeq_obj, physeq_name) {
            read_counts <- analyze_read_counts(physeq_obj)
            read_counts$project <- project_name
            read_counts$physeq_name <- physeq_name
            return(read_counts)
          }
        )
      }
    )
  } else {
    stop("Input must be a phyloseq object or nested list of phyloseq objects")
  }
}
