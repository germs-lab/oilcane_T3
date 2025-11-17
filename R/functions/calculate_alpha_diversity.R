#' Calculate Alpha Diversity
#'
#' @description
#' Calculates alpha diversity metrics for a phyloseq object
#'
#' @param physeq A phyloseq object
#' @param measures Character vector of diversity measures to calculate.
#'   Options: "Observed", "Shannon", "Simpson", "InvSimpson", "Chao1", "ACE", "Fisher"
#'
#' @return A data frame with alpha diversity metrics
#'
#' @export
calculate_alpha_diversity <- function(physeq, measures = c("Observed", "Shannon", "Simpson", "InvSimpson")) {
  # Calculate alpha diversity
  alpha <- phyloseq::estimate_richness(physeq, measures = measures)
  
  # Add sample data if available
  if (!is.null(sample_data(physeq, errorIfNULL = FALSE))) {
    sample_df <- data.frame(sample_data(physeq))
    alpha <- cbind(alpha, sample_df)
  }
  
  # Clean column names
  alpha <- alpha %>%
    janitor::clean_names() %>%
    tibble::rownames_to_column(var = "sample_id")
  
  return(alpha)
}


#' Calculate Alpha Diversity for Nested Phyloseq Lists
#'
#' @description
#' Calculates alpha diversity metrics for nested phyloseq list structures
#'
#' @param physeq_list A nested list of phyloseq objects
#' @param measures Character vector of diversity measures to calculate
#'
#' @return A nested list with alpha diversity data frames
#'
#' @export
calculate_alpha_diversity_nested <- function(physeq_list, measures = c("Observed", "Shannon", "Simpson", "InvSimpson")) {
  purrr::imap(
    physeq_list,
    function(project_list, project_name) {
      purrr::imap(
        project_list,
        function(physeq_obj, physeq_name) {
          alpha_df <- calculate_alpha_diversity(physeq_obj, measures = measures)
          alpha_df$project <- project_name
          alpha_df$physeq_name <- physeq_name
          return(alpha_df)
        }
      )
    }
  )
}
