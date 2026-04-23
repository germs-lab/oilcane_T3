###########################################################################
# Rank Abundance Curves for Oilcane 16S rRNA Data
#
# Produces rank abundance curves for each of three sequencing datasets:
#   - T1: Greenhouse
#   - T2: Field (timepoint 1)
#   - T3: Field + Greenhouse (timepoint 2)
#
# Accessions (genotypes) are used as the treatment variable:
#   Wild type, Lineage 1-5
#
# Author: Bolivar Aponte Rolon
# Date: 2025-04-23
###########################################################################

source("R/utils/000_setup.R")

#--------------------------------------------------------
# SECTION 1: Helper Functions ----
#--------------------------------------------------------

#' Build rank-abundance data frame from a phyloseq object
#'
#' For each sample, taxa are ranked by abundance (rank 1 = most abundant).
#' Returns a long data frame with columns: sample, rank, rel_abundance, genotype.
#'
#' @param physeq A phyloseq object.
#' @param group_var Character. Name of the sample metadata column to use as
#'   the grouping / colour variable.
#' @return A tibble.
build_rank_abundance_df <- function(physeq, group_var = "genotype") {
  otu <- otu_table(physeq)
  if (!taxa_are_rows(physeq)) {
    otu <- t(otu)
  }
  # Relative abundance, columns = samples
  rel_otu <- sweep(otu, 2, colSums(otu), "/")
  rel_otu <- t(rel_otu) # samples as rows

  meta <- data.frame(sample_data(physeq))

  purrr::map_dfr(rownames(rel_otu), function(smp) {
    abundances <- sort(rel_otu[smp, ], decreasing = TRUE)
    abundances <- abundances[abundances > 0] # drop absent taxa
    tibble::tibble(
      sample = smp,
      rank = seq_along(abundances),
      rel_abundance = abundances,
      !!group_var := meta[smp, group_var]
    )
  })
}


#' Plot rank abundance curves coloured by accession
#'
#' Lines represent the mean relative abundance per accession across samples.
#'
#' @param ra_df Data frame produced by \code{build_rank_abundance_df}.
#' @param title Plot title string.
#' @param group_var Column name used for colour / group aesthetic.
#' @param log_y Logical. Log10-transform the y-axis? Default TRUE.
#' @return A ggplot object.
plot_rank_abundance <- function(
  ra_df,
  title = NULL,
  group_var = "genotype",
  log_y = TRUE
) {
  ra_mean <- ra_df |>
    dplyr::summarise(
      mean_rel_abundance = mean(rel_abundance, na.rm = TRUE),
      .by = c(all_of(group_var), rank)
    )

  p <- ggplot(
    ra_mean,
    aes(
      x = rank,
      y = mean_rel_abundance,
      colour = .data[[group_var]],
      group = .data[[group_var]]
    )
  ) +
    geom_line(linewidth = 0.8, alpha = 0.9) +
    scale_colour_brewer(palette = "Dark2", name = "Accession") +
    labs(
      title = title,
      x = "Rank",
      y = if (log_y) {
        "Mean Relative Abundance (log10)"
      } else {
        "Mean Relative Abundance"
      }
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank()
    )

  if (log_y) {
    p <- p + scale_y_log10()
  }

  p
}

#--------------------------------------------------------
# SECTION 2: Subset Phyloseq Objects ----
#--------------------------------------------------------

# T1 -- Greenhouse
physeq_t1_gh <- subset_samples(
  main_oilcane_physeq,
  sampling_time == "T1" & growing_condition == "Greenhouse"
)
physeq_t1_gh <- prune_taxa(taxa_sums(physeq_t1_gh) > 0, physeq_t1_gh)

cat("T1 Greenhouse:\n")
print(physeq_t1_gh)

# T2 -- Field (first field timepoint)
physeq_t2_field <- subset_samples(
  main_oilcane_physeq,
  sampling_time == "T2" & growing_condition == "Field"
)
physeq_t2_field <- prune_taxa(taxa_sums(physeq_t2_field) > 0, physeq_t2_field)

cat("\nT2 Field:\n")
print(physeq_t2_field)

# T3 -- includes both Field and Greenhouse samples (second field timepoint)
physeq_t3 <- subset_samples(
  main_oilcane_physeq,
  sampling_time == "T3"
)
physeq_t3 <- prune_taxa(taxa_sums(physeq_t3) > 0, physeq_t3)

cat("\nT3 (Field + Greenhouse):\n")
print(physeq_t3)

#--------------------------------------------------------
# SECTION 3: Build Rank-Abundance Data Frames ----
#--------------------------------------------------------

ra_t1 <- build_rank_abundance_df(physeq_t1_gh, group_var = "genotype")
ra_t2 <- build_rank_abundance_df(physeq_t2_field, group_var = "genotype")
ra_t3 <- build_rank_abundance_df(physeq_t3, group_var = "genotype")

#--------------------------------------------------------
# SECTION 4: Plot -- Separate Panels per Dataset ----
#--------------------------------------------------------

p_t1 <- plot_rank_abundance(
  ra_t1,
  title = "Rank Abundance - \nT1 Greenhouse",
  group_var = "genotype"
)

p_t2 <- plot_rank_abundance(
  ra_t2,
  title = "Rank Abundance - \nT2 Field",
  group_var = "genotype"
)

p_t3 <- plot_rank_abundance(
  ra_t3,
  title = "Rank Abundance - \nT3 Field + Greenhouse",
  group_var = "genotype"
)

p_t1
p_t2
p_t3

# Combined panel for convenience
ra_combined_plot <- patchwork::wrap_plots(
  p_t1 + theme(legend.position = "none"),
  p_t2 + theme(legend.position = "none"),
  p_t3,
  ncol = 3
) +
  patchwork::plot_annotation(
    title = "Rank Abundance Curves by Accession",
    subtitle = "Mean relative abundance per accession (log10 y-axis)",
    tag_levels = "A"
  )

ra_combined_plot


#--------------------------------------------------------
# SECTION 5: Save Results Object ----
#--------------------------------------------------------

rank_abundance_results <- list(
  data = list(
    t1_greenhouse = ra_t1,
    t2_field = ra_t2,
    t3 = ra_t3
  ),
  plots = list(
    t1_greenhouse = p_t1,
    t2_field = p_t2,
    t3 = p_t3,
    combined = ra_combined_plot
  )
)

save(
  rank_abundance_results,
  file = here::here("data/output/rdata/rank_abundance_results.rda")
)

cat("Results saved to: data/output/rdata/rank_abundance_results.rda\n")

#--------------------------------------------------------
# SECTION 6: Save Plots ----
#--------------------------------------------------------

output_dir <- here::here("data/output/figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

purrr::iwalk(rank_abundance_results$plots, function(plot_obj, name) {
  filename <- file.path(output_dir, paste0("rank_abundance_", name, ".png"))
  ggsave(
    filename = filename,
    plot = plot_obj,
    width = if (name == "combined") 14 else 9,
    height = if (name == "combined") 9 else 5,
    dpi = 300
  )
})

cat("\nPlots saved to:", output_dir, "\n")


cat("\n### Rank Abundance Analysis Complete ###\n")
