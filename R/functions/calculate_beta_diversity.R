#' Calculate Beta Diversity
#'
#' @description
#' Calculates beta diversity ordination for a phyloseq object
#'
#' @param physeq A phyloseq object
#' @param method Ordination method: "PCoA", "NMDS", "DCA", "CCA", "RDA", "DPCoA"
#' @param distance Distance/dissimilarity method (for PCoA/NMDS): "bray", "jaccard", "unifrac", etc.
#'
#' @return A list containing the phyloseq object, ordination result, and method info
#'
#' @export
calculate_beta_diversity <- function(physeq, method = "PCoA", distance = "bray") {
  # Calculate ordination
  ord <- phyloseq::ordinate(physeq, method = method, distance = distance)
  
  # Return results
  list(
    physeq = physeq,
    ordination = ord,
    method = method,
    distance = distance
  )
}


#' Calculate Beta Diversity for Nested Phyloseq Lists
#'
#' @description
#' Calculates beta diversity ordinations for nested phyloseq list structures
#'
#' @param physeq_list A nested list of phyloseq objects
#' @param method Ordination method
#' @param distance Distance/dissimilarity method
#'
#' @return A nested list with beta diversity ordination results
#'
#' @export
calculate_beta_diversity_nested <- function(physeq_list, method = "PCoA", distance = "bray") {
  purrr::imap(
    physeq_list,
    function(project_list, project_name) {
      purrr::imap(
        project_list,
        function(physeq_obj, physeq_name) {
          beta_result <- calculate_beta_diversity(physeq_obj, method = method, distance = distance)
          beta_result$project <- project_name
          beta_result$physeq_name <- physeq_name
          return(beta_result)
        }
      )
    }
  )
}
