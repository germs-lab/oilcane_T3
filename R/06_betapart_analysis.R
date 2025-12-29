###########################################################################
# Exploratory Beta Diversity Partitioning Analysis using betapart
#
# This script uses the betapart package to decompose Bray-Curtis dissimilarity
# into balanced variation (turnover) and abundance gradient components.
#
# Key Questions Addressed:
# 1. Do oilcane lineages have more balanced (turnover) vs abundance gradient
#    dissimilarity compared to wild type?
# 2. How do beta diversity components change across timepoints (T1, T2, T3)?
# 3. Are there differences in community assembly between genotypes?
#
# Background:
# - Balanced variation (bray.bal): Species replacement - communities differ
#   because species are replaced by others with similar abundances
# - Gradient variation (bray.gra): Abundance changes - communities differ
#   due to abundance changes rather than species replacement
#
# Author: Bolívar Aponte Rolón
# Date: 2025-12-29
##########################################################################

source("R/utils/000_setup.R")

# Load betapart package
if (!requireNamespace("betapart", quietly = TRUE)) {
  cat("Installing betapart package...\n")
  install.packages("betapart")
}
library(betapart)

#--------------------------------------------------------
# SECTION 1: Overview and Data Preparation
#--------------------------------------------------------

cat("### Betapart Analysis: Beta Diversity Partitioning ###\n\n")

# Display study design overview
cat("Study Design Overview:\n")
cat("- Genotypes: Lineage 1-5 (oilcane) vs Wild type (control)\n")
cat("- Timepoints: T1, T2, T3 (maturity stages)\n")
cat("- Sample types: Roots, Bulk soil, Root-associated soils, Leaves, Stalks\n\n")

# Check phyloseq object
cat("Phyloseq object summary:\n")
print(main_oilcane_physeq)

# Get sample metadata
sample_df <- data.frame(sample_data(main_oilcane_physeq))

# Create oilcane vs wild type classification
sample_df$is_oilcane <- ifelse(
  sample_df$genotype == "Wild type",
  "Wild type",
  "Oilcane"
)

# Update sample data in phyloseq
sample_data(main_oilcane_physeq)$is_oilcane <- sample_df$is_oilcane

cat("\nSample distribution:\n")
print(table(sample_df$genotype, sample_df$sampling_time))

#--------------------------------------------------------
# SECTION 2: Overall Beta Diversity Partitioning
#--------------------------------------------------------

cat("\n### SECTION 2: Overall Beta Diversity Partitioning ###\n")

# Calculate betapart for entire dataset
overall_betapart <- calculate_betapart(main_oilcane_physeq, index.family = "bray")

cat("\nOverall Bray-Curtis dissimilarity decomposition:\n")
print(overall_betapart$summary)

# Visualize overall decomposition
overall_plot <- plot_betapart_summary(
  overall_betapart,
  title = "Overall Beta Diversity Partitioning"
)
print(overall_plot)

#--------------------------------------------------------
# SECTION 3: Compare Oilcane vs Wild Type
#--------------------------------------------------------

cat("\n### SECTION 3: Oilcane vs Wild Type Comparison ###\n")

# Calculate betapart by oilcane classification
betapart_oilcane <- calculate_betapart_by_group(
  main_oilcane_physeq,
  group_var = "is_oilcane",
  index.family = "bray"
)

cat("\nComparison: Oilcane vs Wild Type\n")
print(betapart_oilcane$group_comparison)

# Plot comparison
oilcane_comparison_plot <- plot_betapart_comparison(
  betapart_oilcane,
  title = "Beta Diversity: Oilcane vs Wild Type"
)
print(oilcane_comparison_plot)

# Interpretation
cat("\n--- Interpretation ---\n")
oilcane_row <- betapart_oilcane$group_comparison[
  betapart_oilcane$group_comparison$group == "Oilcane",
]
wt_row <- betapart_oilcane$group_comparison[
  betapart_oilcane$group_comparison$group == "Wild type",
]

if (nrow(oilcane_row) > 0 && nrow(wt_row) > 0) {
  cat(sprintf(
    "Oilcane: %.1f%% balanced, %.1f%% gradient\n",
    oilcane_row$prop_balanced * 100,
    oilcane_row$prop_gradient * 100
  ))
  cat(sprintf(
    "Wild type: %.1f%% balanced, %.1f%% gradient\n",
    wt_row$prop_balanced * 100,
    wt_row$prop_gradient * 100
  ))

  if (oilcane_row$prop_balanced > wt_row$prop_balanced) {
    cat("\n>> Oilcane communities show MORE balanced (turnover) variation.\n")
    cat("   This suggests species replacement is more common in oilcane.\n")
  } else {
    cat("\n>> Oilcane communities show LESS balanced variation than wild type.\n")
    cat("   This suggests abundance gradients drive more of the dissimilarity.\n")
  }
}

#--------------------------------------------------------
# SECTION 4: Compare by Timepoint
#--------------------------------------------------------

cat("\n### SECTION 4: Temporal Patterns (T1, T2, T3) ###\n")

# Calculate betapart by timepoint
betapart_timepoint <- calculate_betapart_by_group(
  main_oilcane_physeq,
  group_var = "sampling_time",
  index.family = "bray"
)

cat("\nBeta diversity partitioning by timepoint:\n")
print(betapart_timepoint$group_comparison)

# Plot temporal patterns
timepoint_plot <- plot_betapart_comparison(
  betapart_timepoint,
  title = "Beta Diversity Partitioning Across Timepoints"
)
print(timepoint_plot)

#--------------------------------------------------------
# SECTION 5: Compare by Genotype
#--------------------------------------------------------

cat("\n### SECTION 5: Genotype-Specific Patterns ###\n")

# Calculate betapart by genotype
betapart_genotype <- calculate_betapart_by_group(
  main_oilcane_physeq,
  group_var = "genotype",
  index.family = "bray"
)

cat("\nBeta diversity partitioning by genotype:\n")
print(betapart_genotype$group_comparison)

# Plot genotype patterns
genotype_plot <- plot_betapart_comparison(
  betapart_genotype,
  title = "Beta Diversity Partitioning by Genotype"
)
print(genotype_plot)

#--------------------------------------------------------
# SECTION 6: Sample Type Analysis
#--------------------------------------------------------

cat("\n### SECTION 6: Sample Type Patterns ###\n")

# Calculate betapart by sample type
betapart_material <- calculate_betapart_by_group(
  main_oilcane_physeq,
  group_var = "original_materials",
  index.family = "bray"
)

cat("\nBeta diversity partitioning by sample type:\n")
print(betapart_material$group_comparison)

# Plot material patterns
material_plot <- plot_betapart_comparison(
  betapart_material,
  title = "Beta Diversity Partitioning by Sample Type"
)
print(material_plot)

#--------------------------------------------------------
# SECTION 7: Combined Analysis (Oilcane x Timepoint)
#--------------------------------------------------------

cat("\n### SECTION 7: Oilcane vs WT Across Timepoints ###\n")

# Create combined grouping variable
sample_data(main_oilcane_physeq)$oilcane_time <- paste(
  sample_data(main_oilcane_physeq)$is_oilcane,
  sample_data(main_oilcane_physeq)$sampling_time,
  sep = "_"
)

# Calculate betapart for combined groups
betapart_combined <- calculate_betapart_by_group(
  main_oilcane_physeq,
  group_var = "oilcane_time",
  index.family = "bray"
)

cat("\nBeta diversity: Oilcane status x Timepoint\n")
print(betapart_combined$group_comparison)

# Create enhanced visualization
combined_df <- betapart_combined$group_comparison %>%
  tidyr::separate(group, into = c("type", "timepoint"), sep = "_") %>%
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

combined_plot <- ggplot(
  combined_df,
  aes(x = timepoint, y = value, fill = component_label)
) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~type) +
  scale_fill_manual(
    values = c("Balanced (Turnover)" = "#2166AC", "Gradient (Abundance)" = "#B2182B"),
    name = "Component"
  ) +
  labs(
    title = "Beta Diversity Partitioning: Oilcane vs Wild Type Over Time",
    subtitle = "Comparing balanced (species replacement) vs gradient (abundance) variation",
    y = "Mean Dissimilarity",
    x = "Timepoint"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(combined_plot)

#--------------------------------------------------------
# SECTION 8: Statistical Summary & Conclusions
#--------------------------------------------------------

cat("\n### SECTION 8: Summary & Key Findings ###\n")

# Calculate ratio of balanced to gradient for key comparisons
summary_results <- list()

# Overall
summary_results$overall <- data.frame(
  analysis = "Overall",
  n_samples = nsamples(main_oilcane_physeq),
  mean_balanced = overall_betapart$summary$mean[1],
  mean_gradient = overall_betapart$summary$mean[2],
  mean_total = overall_betapart$summary$mean[3],
  balanced_to_gradient_ratio = overall_betapart$summary$mean[1] /
    overall_betapart$summary$mean[2]
)

# Oilcane vs WT
if (nrow(betapart_oilcane$group_comparison) > 0) {
  summary_results$oilcane_vs_wt <- betapart_oilcane$group_comparison %>%
    mutate(
      analysis = paste("Within", group),
      balanced_to_gradient_ratio = mean_balanced / mean_gradient
    ) %>%
    select(
      analysis, n_samples, mean_balanced, mean_gradient,
      mean_total, balanced_to_gradient_ratio
    )
}

# Combine summaries
final_summary <- bind_rows(summary_results)

cat("\n=== Final Summary Table ===\n")
print(final_summary)

# Key findings narrative
cat("\n=== Key Findings ===\n\n")

cat("1. OVERALL COMMUNITY STRUCTURE:\n")
cat(sprintf(
  "   Total Bray-Curtis dissimilarity: %.3f\n",
  overall_betapart$summary$mean[3]
))
cat(sprintf(
  "   Balanced component: %.3f (%.1f%%)\n",
  overall_betapart$summary$mean[1],
  overall_betapart$summary$proportion_of_total[1] * 100
))
cat(sprintf(
  "   Gradient component: %.3f (%.1f%%)\n",
  overall_betapart$summary$mean[2],
  overall_betapart$summary$proportion_of_total[2] * 100
))

cat("\n2. OILCANE VS WILD TYPE:\n")
if (nrow(oilcane_row) > 0 && nrow(wt_row) > 0) {
  bal_diff <- oilcane_row$prop_balanced - wt_row$prop_balanced
  cat(sprintf(
    "   Oilcane balanced proportion: %.1f%%\n",
    oilcane_row$prop_balanced * 100
  ))
  cat(sprintf(
    "   Wild type balanced proportion: %.1f%%\n",
    wt_row$prop_balanced * 100
  ))
  cat(sprintf(
    "   Difference: %.1f percentage points\n",
    abs(bal_diff) * 100
  ))

  if (bal_diff > 0.05) {
    cat("   >> FINDING: Oilcane shows substantially MORE balanced variation.\n")
    cat("      Interpretation: Oilcane microbial communities show more species\n")
    cat("      replacement patterns, suggesting distinct selection pressures.\n")
  } else if (bal_diff < -0.05) {
    cat("   >> FINDING: Oilcane shows substantially LESS balanced variation.\n")
    cat("      Interpretation: Oilcane communities are more driven by\n")
    cat("      abundance gradients rather than species replacement.\n")
  } else {
    cat("   >> FINDING: Similar balance between turnover and abundance.\n")
    cat("      Interpretation: Both oilcane and wild type show similar\n")
    cat("      community assembly patterns.\n")
  }
}

cat("\n3. TEMPORAL PATTERNS:\n")
if (nrow(betapart_timepoint$group_comparison) > 0) {
  for (i in seq_len(nrow(betapart_timepoint$group_comparison))) {
    row <- betapart_timepoint$group_comparison[i, ]
    cat(sprintf(
      "   %s: %.1f%% balanced, %.1f%% gradient\n",
      row$group,
      row$prop_balanced * 100,
      row$prop_gradient * 100
    ))
  }
}

#--------------------------------------------------------
# SECTION 9: Save Results
#--------------------------------------------------------

cat("\n### Saving Results ###\n")

# Create output directory
dir.create(
  here::here("data/output/rdata"),
  recursive = TRUE,
  showWarnings = FALSE
)

# Compile all results
betapart_results <- list(
  overall = overall_betapart,
  by_oilcane = betapart_oilcane,
  by_timepoint = betapart_timepoint,
  by_genotype = betapart_genotype,
  by_material = betapart_material,
  by_oilcane_time = betapart_combined,
  summary = final_summary,
  plots = list(
    overall = overall_plot,
    oilcane_comparison = oilcane_comparison_plot,
    timepoint = timepoint_plot,
    genotype = genotype_plot,
    material = material_plot,
    combined = combined_plot
  )
)

# Save results
save(
  betapart_results,
  file = here::here("data/output/rdata/betapart_results.rda")
)

cat("Results saved to: data/output/rdata/betapart_results.rda\n")

#--------------------------------------------------------
# SECTION 10: Proposed Next Steps
#--------------------------------------------------------

cat("\n### Proposed Next Steps for Further Analysis ###\n\n")

cat("1. STATISTICAL TESTING:\
")
cat("   - PERMANOVA on balanced vs gradient components\n")
cat("   - Test if oilcane vs WT difference is significant\n")
cat("   - Multi-factor analysis: genotype x timepoint x material\n\n")

cat("2. PHYLOGENETIC CONSIDERATIONS:\n")
cat("   - Use UniFrac-based betapart (requires phylogenetic tree)\n")
cat("   - Compare taxonomic vs phylogenetic turnover\n\n")

cat("3. INDICATOR SPECIES:\n")
cat("   - Identify ASVs contributing to balanced variation\n")
cat("   - Identify ASVs contributing to gradient variation\n\n")

cat("4. TEMPORAL DYNAMICS:\n")
cat("   - Track specific lineages across T1 -> T2 -> T3\n")
cat("   - Identify if turnover increases or decreases with time\n\n")

cat("5. FUNCTIONAL IMPLICATIONS:\n")
cat("   - Link betapart results to predicted functions\n")
cat("   - Test if functional diversity follows taxonomic patterns\n\n")

cat("### Betapart Analysis Complete ###\n")
