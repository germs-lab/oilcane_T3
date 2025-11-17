###########################################################################
# Alpha and Beta Diversity Analysis for Oilcane T3
#
# This script performs alpha and beta diversity analyses on the oilcane
# phyloseq object including:
# - Alpha diversity calculation (Observed, Shannon, Simpson, InvSimpson)
# - Alpha diversity visualization and statistical comparisons
# - Beta diversity ordination (PCoA with Bray-Curtis)
# - Beta diversity visualization
# - PERMANOVA analysis
#
# Author: Bolívar Aponte Rolón
# Date: 2025-11-17
##########################################################################

source("R/utils/000_setup.R")

#--------------------------------------------------------
# SECTION 1: Load Data
#--------------------------------------------------------

# Load the main phyloseq object
if (!exists("oilcane_physeq")) {
  load(here::here("data/output/processed/rdata/oilcane_physeq.rda"))
}

#--------------------------------------------------------
# SECTION 2: Alpha Diversity Analysis
#--------------------------------------------------------

cat("### Calculating Alpha Diversity ###\n")

# Calculate alpha diversity
alpha_diversity <- calculate_alpha_diversity(
  oilcane_physeq,
  measures = c("Observed", "Shannon", "Simpson", "InvSimpson")
)

# View summary
summary(alpha_diversity[, c("observed", "shannon", "simpson", "inv_simpson")])

# Print first few rows
head(alpha_diversity)

#--------------------------------------------------------
# SECTION 3: Alpha Diversity Plots
#--------------------------------------------------------

cat("### Creating Alpha Diversity Plots ###\n")

# Check which grouping variables are available
available_vars <- colnames(alpha_diversity)
cat("Available variables for grouping:\n")
print(available_vars)

# Create plots for different alpha diversity metrics
# Note: Adjust grouping variable based on your metadata

# Function to create alpha diversity plot
create_alpha_plot <- function(data, metric, group_var = NULL, title = NULL) {
  if (is.null(group_var) || !group_var %in% colnames(data)) {
    # No grouping - single boxplot
    p <- ggplot(data, aes(x = 1, y = .data[[metric]])) +
      geom_boxplot(alpha = 0.7, fill = "steelblue") +
      geom_jitter(width = 0.2, alpha = 0.5) +
      labs(
        title = ifelse(is.null(title), paste(str_to_title(metric), "Diversity"), title),
        x = "",
        y = str_to_title(metric)
      ) +
      theme_minimal() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    # With grouping
    p <- ggplot(data, aes(x = .data[[group_var]], y = .data[[metric]], fill = .data[[group_var]])) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      labs(
        title = ifelse(is.null(title), paste(str_to_title(metric), "Diversity"), title),
        x = str_to_title(gsub("_", " ", group_var)),
        y = str_to_title(metric)
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
  }
  
  return(p)
}

# Create plots for each diversity metric
alpha_plots <- list(
  observed = create_alpha_plot(alpha_diversity, "observed", title = "Observed Richness"),
  shannon = create_alpha_plot(alpha_diversity, "shannon", title = "Shannon Diversity"),
  simpson = create_alpha_plot(alpha_diversity, "simpson", title = "Simpson Diversity"),
  inv_simpson = create_alpha_plot(alpha_diversity, "inv_simpson", title = "Inverse Simpson Diversity")
)

# Display plots
print(alpha_plots$observed)
print(alpha_plots$shannon)
print(alpha_plots$simpson)
print(alpha_plots$inv_simpson)

# Combined plot
alpha_combined <- ggpubr::ggarrange(
  alpha_plots$observed,
  alpha_plots$shannon,
  alpha_plots$simpson,
  alpha_plots$inv_simpson,
  labels = c("A", "B", "C", "D"),
  nrow = 2,
  ncol = 2
)

print(alpha_combined)

#--------------------------------------------------------
# SECTION 4: Beta Diversity Analysis (PCoA)
#--------------------------------------------------------

cat("\n### Calculating Beta Diversity (PCoA) ###\n")

# Calculate beta diversity using Bray-Curtis dissimilarity
beta_diversity <- calculate_beta_diversity(
  oilcane_physeq,
  method = "PCoA",
  distance = "bray"
)

cat("Beta diversity ordination complete\n")

#--------------------------------------------------------
# SECTION 5: Beta Diversity Plots
#--------------------------------------------------------

cat("### Creating Beta Diversity Plots ###\n")

# Basic PCoA plot
pcoa_plot_basic <- plot_ordination(
  beta_diversity$physeq,
  beta_diversity$ordination,
  type = "samples"
) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    title = "PCoA (Bray-Curtis Dissimilarity)",
    subtitle = "All Samples"
  ) +
  theme_bw()

print(pcoa_plot_basic)

# Check if there are grouping variables available
if (!is.null(sample_data(oilcane_physeq, errorIfNULL = FALSE))) {
  sample_vars <- sample_variables(oilcane_physeq)
  
  # Try to create colored plots if suitable variables exist
  # You may need to adjust these variable names based on your metadata
  potential_color_vars <- intersect(
    sample_vars,
    c("Timepoint", "timepoint", "Treatment", "treatment", "Site", "site", "Crop", "crop")
  )
  
  if (length(potential_color_vars) > 0) {
    color_var <- potential_color_vars[1]
    
    pcoa_plot_colored <- plot_ordination(
      beta_diversity$physeq,
      beta_diversity$ordination,
      type = "samples",
      color = color_var
    ) +
      geom_point(size = 3, alpha = 0.7) +
      labs(
        title = "PCoA (Bray-Curtis Dissimilarity)",
        subtitle = paste("Colored by", str_to_title(color_var))
      ) +
      theme_bw() +
      theme(legend.position = "right")
    
    print(pcoa_plot_colored)
  }
}

#--------------------------------------------------------
# SECTION 6: PERMANOVA Analysis
#--------------------------------------------------------

cat("\n### Running PERMANOVA ###\n")

# Get distance matrix
dist_matrix <- phyloseq::distance(oilcane_physeq, method = "bray")

# Get sample data
sample_df <- data.frame(sample_data(oilcane_physeq))

# Run PERMANOVA if suitable grouping variable exists
# Note: This is a template - adjust the formula based on your metadata
if (ncol(sample_df) > 0) {
  cat("Sample data variables available:\n")
  print(colnames(sample_df))
  cat("\nPERMANOVA analysis requires a grouping variable.\n")
  cat("Please edit this script to add PERMANOVA analysis with appropriate variables.\n")
  
  # Example PERMANOVA (uncomment and adjust based on your data):
  # set.seed(123)
  # permanova_result <- vegan::adonis2(
  #   dist_matrix ~ Timepoint,  # Adjust formula
  #   data = sample_df,
  #   permutations = 999
  # )
  # print(permanova_result)
} else {
  cat("No sample metadata available for PERMANOVA analysis.\n")
}

#--------------------------------------------------------
# SECTION 7: Summary Statistics
#--------------------------------------------------------

cat("\n### Alpha Diversity Summary Statistics ###\n")

# Overall summary
alpha_summary_overall <- alpha_diversity %>%
  summarise(
    n = n(),
    mean_observed = mean(observed, na.rm = TRUE),
    sd_observed = sd(observed, na.rm = TRUE),
    mean_shannon = mean(shannon, na.rm = TRUE),
    sd_shannon = sd(shannon, na.rm = TRUE),
    mean_simpson = mean(simpson, na.rm = TRUE),
    sd_simpson = sd(simpson, na.rm = TRUE),
    mean_inv_simpson = mean(inv_simpson, na.rm = TRUE),
    sd_inv_simpson = sd(inv_simpson, na.rm = TRUE)
  )

print(alpha_summary_overall)

#--------------------------------------------------------
# SECTION 8: Save Results
#--------------------------------------------------------

# Save alpha and beta diversity results
alpha_beta_results <- list(
  alpha_diversity = alpha_diversity,
  alpha_plots = alpha_plots,
  alpha_combined = alpha_combined,
  alpha_summary = alpha_summary_overall,
  beta_diversity = beta_diversity,
  pcoa_plot = pcoa_plot_basic
)

save(
  alpha_beta_results,
  file = here::here("data/output/processed/rdata/alpha_beta_results.rda")
)

cat("\n### Alpha and Beta Diversity Analysis Complete ###\n")
cat("Results saved to: data/output/processed/rdata/alpha_beta_results.rda\n")
