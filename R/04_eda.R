###########################################################################
# Exploratory Data Analysis for Oilcane T3 Phyloseq Objects
#
# This script performs EDA on the oilcane phyloseq object including:
# - Basic exploration of the phyloseq structure
# - Read count analysis and visualization
# - Good's coverage estimation
# - Rarefaction curves using iNEXT
#
# Author: Bolívar Aponte Rolón
# Date: 2025-11-17
##########################################################################

source("R/utils/000_setup.R")

#--------------------------------------------------------
# SECTION 1: Basic Exploration of Phyloseq
#--------------------------------------------------------

# Load the main phyloseq object
if (!exists("oilcane_physeq")) {
  load(here::here("data/output/processed/rdata/oilcane_physeq.rda"))
}

# Examine the structure
cat("### OILCANE PHYLOSEQ OBJECT ###\n")
explore_phyloseq(oilcane_physeq, name = "Oilcane All Timepoints")

# Summary table
physeq_summary <- data.frame(
  n_taxa = ntaxa(oilcane_physeq),
  n_samples = nsamples(oilcane_physeq),
  total_reads = sum(sample_sums(oilcane_physeq)),
  min_reads_per_sample = min(sample_sums(oilcane_physeq)),
  max_reads_per_sample = max(sample_sums(oilcane_physeq)),
  mean_reads_per_sample = mean(sample_sums(oilcane_physeq)),
  median_reads_per_sample = median(sample_sums(oilcane_physeq))
)

print(physeq_summary)

#--------------------------------------------------------
# SECTION 2: Read Count Analysis
#--------------------------------------------------------

# Get read count data
read_counts <- analyze_read_counts(oilcane_physeq)

# Visualize read counts
read_count_plots <- list(
  # Density plot
  density = ggplot(read_counts, aes(x = n_seqs)) +
    geom_density(alpha = 0.6, fill = "steelblue") +
    labs(
      title = "Read Count Distribution",
      x = "Number of Sequences",
      y = "Density"
    ) +
    theme_minimal(),
  
  # Histogram
  histogram = ggplot(read_counts, aes(x = n_seqs)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    labs(
      title = "Read Count Distribution",
      x = "Number of Sequences",
      y = "Count"
    ) +
    theme_minimal(),
  
  # Box plot
  boxplot = ggplot(read_counts, aes(x = 1, y = n_seqs)) +
    geom_boxplot(alpha = 0.7, fill = "steelblue") +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(
      title = "Read Count Distribution",
      x = "",
      y = "Number of Sequences"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()),
  
  # Ranked line plot
  ranked = read_counts %>%
    arrange(n_seqs) %>%
    mutate(sample_rank = row_number()) %>%
    ggplot(aes(x = sample_rank, y = n_seqs)) +
    geom_line(linewidth = 1, color = "steelblue") +
    geom_point(alpha = 0.5) +
    labs(
      title = "Read Count Distribution (Ranked)",
      x = "Sample Rank",
      y = "Number of Sequences"
    ) +
    theme_minimal(),
  
  # Good's coverage
  goods_coverage = read_counts %>%
    ggplot(aes(x = n_seqs, y = goods)) +
    geom_point(alpha = 0.6, color = "steelblue") +
    geom_smooth(se = TRUE, color = "darkblue") +
    labs(
      title = "Good's Coverage vs Sequencing Depth",
      x = "Number of Sequences",
      y = "Good's Coverage"
    ) +
    theme_minimal()
)

# Display plots
print(read_count_plots$density)
print(read_count_plots$histogram)
print(read_count_plots$boxplot)
print(read_count_plots$ranked)
print(read_count_plots$goods_coverage)

#--------------------------------------------------------
# SECTION 3: Rarefaction Curves with iNEXT
#--------------------------------------------------------

cat("\n### Running iNEXT for Rarefaction Curves ###\n")

# Prepare OTU matrix for iNEXT
otu_mat <- oilcane_physeq %>%
  otu_table() %>%
  data.frame() %>%
  as.matrix()

# Transpose if taxa are not in rows
if (!taxa_are_rows(oilcane_physeq)) {
  otu_mat <- t(otu_mat)
}

# Calculate endpoint
max_lib_size <- max(rowSums(t(otu_mat)))
endpoint <- max_lib_size * 2

cat("Max library size:", max_lib_size, "\n")
cat("Endpoint:", endpoint, "\n")

# Run parallel iNEXT
# Note: Adjust nCores based on available resources
options(future.globals.maxSize = 1500 * 1024^2)

inext_result <- p_iNEXT(
  x = t(otu_mat),  # samples in rows, taxa in columns
  q = c(0, 1, 2),
  endpoint = endpoint,
  nboot = 100,
  nCores = 4,
  combine = TRUE,
  verbose = TRUE
)

# Create ggiNEXT plot
library(iNEXT)  # Ensure iNEXT is loaded for ggiNEXT

inext_plot <- ggiNEXT(
  inext_result,
  type = 1,
  facet.var = "Order.q",
  color.var = "Assemblage"
) +
  theme_bw() +
  labs(
    title = "Rarefaction Curves - Oilcane All Timepoints",
    x = "Number of Sequences"
  ) +
  guides(color = "none", shape = "none", fill = "none")

print(inext_plot)

#--------------------------------------------------------
# SECTION 4: Taxonomic Composition Overview
#--------------------------------------------------------

# Phylum-level composition
if (!is.null(tax_table(oilcane_physeq, errorIfNULL = FALSE))) {
  cat("\n### Taxonomic Composition ###\n")
  
  # Get top phyla
  top_phyla <- oilcane_physeq %>%
    tax_glom(taxrank = "Phylum") %>%
    transform_sample_counts(function(x) x / sum(x)) %>%
    psmelt() %>%
    group_by(Phylum) %>%
    summarise(mean_abundance = mean(Abundance)) %>%
    arrange(desc(mean_abundance)) %>%
    head(10)
  
  print(top_phyla)
  
  # Plot top phyla
  phylum_plot <- oilcane_physeq %>%
    tax_glom(taxrank = "Phylum") %>%
    transform_sample_counts(function(x) x / sum(x)) %>%
    psmelt() %>%
    group_by(Phylum) %>%
    mutate(mean_abundance = mean(Abundance)) %>%
    ungroup() %>%
    mutate(Phylum = fct_reorder(Phylum, mean_abundance)) %>%
    ggplot(aes(x = Phylum, y = Abundance, fill = Phylum)) +
    geom_boxplot(alpha = 0.7) +
    coord_flip() +
    labs(
      title = "Phylum-level Relative Abundance",
      x = "Phylum",
      y = "Relative Abundance"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  print(phylum_plot)
}

#--------------------------------------------------------
# SECTION 5: Save Results
#--------------------------------------------------------

# Save EDA results
eda_results <- list(
  physeq_summary = physeq_summary,
  read_counts = read_counts,
  read_count_plots = read_count_plots,
  inext_result = inext_result,
  inext_plot = inext_plot
)

save(
  eda_results,
  file = here::here("data/output/processed/rdata/eda_results.rda")
)

cat("\n### EDA Complete ###\n")
cat("Results saved to: data/output/processed/rdata/eda_results.rda\n")
