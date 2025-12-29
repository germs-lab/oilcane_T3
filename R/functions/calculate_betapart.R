#' Calculate Beta Diversity Partitioning using betapart
#'
#' @description
#' Partitions Bray-Curtis dissimilarity into balanced variation (turnover)
#' and abundance gradient (nestedness-resultant) components using the betapart package.
#' This allows understanding whether community differences are driven by species
#' replacement or abundance changes.
#'
#' @param physeq A phyloseq object
#' @param index.family Either "bray" (default) or "ruzicka" for abundance data
#'
#' @return A list containing:
#'   - bray.bal: Balanced variation component (turnover-like)
#'   - bray.gra: Abundance gradient component (nestedness-like)
#'   - bray: Total Bray-Curtis dissimilarity
#'   - summary: Data frame summarizing mean dissimilarity components
#'
#' @details
#' For abundance data, betapart decomposes Bray-Curtis dissimilarity into:
#' - Balanced variation (bray.bal): Differences due to species replacement
#' - Abundance gradient (bray.gra): Differences due to abundance changes
#' 
#' The sum of bray.bal and bray.gra equals total Bray-Curtis dissimilarity.
#'
#' @export
calculate_betapart <- function(physeq, index.family = "bray") {
  # Require betapart package
  if (!requireNamespace("betapart", quietly = TRUE)) {
    stop("Package 'betapart' is required. Install with install.packages('betapart')")
  }

  # Extract OTU table as data frame (samples as rows, taxa as columns)
  otu_mat <- physeq %>%
    otu_table() %>%
    as.data.frame()

  # Transpose if taxa are rows

  if (taxa_are_rows(physeq)) {
    otu_mat <- t(otu_mat)
  }

  # Convert to data frame for betapart
  otu_mat <- as.data.frame(otu_mat)

  # Calculate beta diversity components
  beta_core <- betapart::beta.pair.abund(otu_mat, index.family = index.family)

  # Extract distance matrices
  bray_bal <- as.matrix(beta_core$beta.bray.bal)
  bray_gra <- as.matrix(beta_core$beta.bray.gra)
  bray_total <- as.matrix(beta_core$beta.bray)

  # Calculate summary statistics
  # Get upper triangle values (excluding diagonal)
  get_upper_tri <- function(mat) {
    mat[upper.tri(mat)]
  }

  summary_df <- data.frame(
    component = c("balanced", "gradient", "total"),
    mean = c(
      mean(get_upper_tri(bray_bal), na.rm = TRUE),
      mean(get_upper_tri(bray_gra), na.rm = TRUE),
      mean(get_upper_tri(bray_total), na.rm = TRUE)
    ),
    sd = c(
      sd(get_upper_tri(bray_bal), na.rm = TRUE),
      sd(get_upper_tri(bray_gra), na.rm = TRUE),
      sd(get_upper_tri(bray_total), na.rm = TRUE)
    ),
    min = c(
      min(get_upper_tri(bray_bal), na.rm = TRUE),
      min(get_upper_tri(bray_gra), na.rm = TRUE),
      min(get_upper_tri(bray_total), na.rm = TRUE)
    ),
    max = c(
      max(get_upper_tri(bray_bal), na.rm = TRUE),
      max(get_upper_tri(bray_gra), na.rm = TRUE),
      max(get_upper_tri(bray_total), na.rm = TRUE)
    )
  )

  # Calculate proportion of balanced vs gradient
  summary_df$proportion_of_total <- summary_df$mean / summary_df$mean[3]

  # Return results
  list(
    bray.bal = beta_core$beta.bray.bal,
    bray.gra = beta_core$beta.bray.gra,
    bray = beta_core$beta.bray,
    summary = summary_df
  )
}


#' Calculate Beta Diversity Partitioning by Groups
#'
#' @description
#' Calculates betapart dissimilarity components within and between groups
#' defined by a sample variable.
#'
#' @param physeq A phyloseq object
#' @param group_var Character string specifying the grouping variable from sample_data
#' @param index.family Either "bray" (default) or "ruzicka"
#'
#' @return A list containing:
#'   - within_group: Betapart results for samples within each group
#'   - between_group: Mean dissimilarity components between groups
#'   - group_comparison: Data frame comparing groups
#'
#' @export
calculate_betapart_by_group <- function(physeq, group_var, index.family = "bray") {
  # Require betapart package
  if (!requireNamespace("betapart", quietly = TRUE)) {
    stop("Package 'betapart' is required. Install with install.packages('betapart')")
  }

  # Get sample data
  sample_df <- data.frame(sample_data(physeq))

  # Check if group variable exists

  if (!group_var %in% colnames(sample_df)) {
    stop(paste("Variable", group_var, "not found in sample data"))
  }

  # Get unique groups
  groups <- unique(sample_df[[group_var]])

  # Calculate within-group betapart
  within_group <- lapply(groups, function(grp) {
    # Subset phyloseq to this group
    samples_in_group <- rownames(sample_df)[sample_df[[group_var]] == grp]
    physeq_subset <- prune_samples(samples_in_group, physeq)

    # Need at least 2 samples
    if (nsamples(physeq_subset) < 2) {
      return(NULL)
    }

    # Calculate betapart
    result <- calculate_betapart(physeq_subset, index.family = index.family)
    result$group <- grp
    result$n_samples <- nsamples(physeq_subset)
    result
  })
  names(within_group) <- groups

  # Remove NULL entries (groups with < 2 samples)
  within_group <- within_group[!sapply(within_group, is.null)]

  # Create summary comparison data frame
  group_comparison <- do.call(rbind, lapply(within_group, function(x) {
    data.frame(
      group = x$group,
      n_samples = x$n_samples,
      mean_balanced = x$summary$mean[1],
      mean_gradient = x$summary$mean[2],
      mean_total = x$summary$mean[3],
      prop_balanced = x$summary$proportion_of_total[1],
      prop_gradient = x$summary$proportion_of_total[2]
    )
  }))
  rownames(group_comparison) <- NULL

  list(
    within_group = within_group,
    group_comparison = group_comparison
  )
}


#' Plot Betapart Components
#'
#' @description
#' Creates visualization of betapart dissimilarity components
#'
#' @param betapart_result Result from calculate_betapart()
#' @param title Plot title
#'
#' @return A ggplot object showing balanced vs gradient contributions
#'
#' @export
plot_betapart_summary <- function(betapart_result, title = "Beta Diversity Partitioning") {
  summary_df <- betapart_result$summary

  # Create stacked bar for proportions
  plot_df <- summary_df %>%
    filter(component != "total") %>%
    mutate(
      component = factor(component, levels = c("gradient", "balanced")),
      component_label = case_when(
        component == "balanced" ~ "Balanced (Turnover)",
        component == "gradient" ~ "Gradient (Abundance)"
      )
    )

  p <- ggplot(plot_df, aes(x = 1, y = mean, fill = component_label)) +
    geom_bar(stat = "identity", width = 0.5) +
    geom_text(
      aes(label = sprintf("%.2f (%.0f%%)", mean, proportion_of_total * 100)),
      position = position_stack(vjust = 0.5),
      color = "white",
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = c("Balanced (Turnover)" = "#2166AC", "Gradient (Abundance)" = "#B2182B"),
      name = "Component"
    ) +
    labs(
      title = title,
      subtitle = sprintf("Total Bray-Curtis: %.3f", summary_df$mean[3]),
      y = "Mean Dissimilarity",
      x = ""
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom"
    ) +
    coord_flip()

  return(p)
}


#' Plot Betapart Group Comparison
#'
#' @description
#' Creates visualization comparing betapart components across groups
#'
#' @param betapart_by_group Result from calculate_betapart_by_group()
#' @param title Plot title
#'
#' @return A ggplot object comparing groups
#'
#' @export
plot_betapart_comparison <- function(betapart_by_group, title = "Beta Diversity by Group") {
  comparison_df <- betapart_by_group$group_comparison

  # Reshape for plotting
  plot_df <- comparison_df %>%
    select(group, mean_balanced, mean_gradient) %>%
    tidyr::pivot_longer(
      cols = c(mean_balanced, mean_gradient),
      names_to = "component",
      values_to = "value"
    ) %>%
    mutate(
      component_label = case_when(
        component == "mean_balanced" ~ "Balanced (Turnover)",
        component == "mean_gradient" ~ "Gradient (Abundance)"
      )
    )

  p <- ggplot(plot_df, aes(x = group, y = value, fill = component_label)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(
      values = c("Balanced (Turnover)" = "#2166AC", "Gradient (Abundance)" = "#B2182B"),
      name = "Component"
    ) +
    labs(
      title = title,
      y = "Mean Dissimilarity",
      x = ""
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )

  return(p)
}
