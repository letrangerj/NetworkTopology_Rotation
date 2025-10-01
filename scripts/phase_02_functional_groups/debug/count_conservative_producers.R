# count_conservative_producers.R
# Debug script: Count strains by conservative producer categories
# Run from project root: Rscript scripts/phase_02_functional_groups/debug/count_conservative_producers.R

suppressPackageStartupMessages({
  library(dplyr)
})

cat("========================================\n")
cat("Debug: Conservative Producer Counts\n")
cat("========================================\n\n")

# Load updated strain repertoires
strain_repertoires <- readRDS("data/interim/strain_repertoires.rds")

# Compute counts
total_strains <- nrow(strain_repertoires)

# 1. Strains with multiple syn groups where at least two are validated (multi-producer conservative)
multi_conservative <- strain_repertoires %>%
  filter(n_production_groups > 1 & n_validated_production_groups >= 2) %>%
  nrow()

# 2. Strains with multiple syn groups but only one or zero validated (single or non-validated producer)
single_or_non_validated <- strain_repertoires %>%
  filter(n_production_groups > 1 & n_validated_production_groups <= 1) %>%
  nrow()

# 3. Strains with only unvalidated syn groups (not included as producer nodes in conservative network)
only_unvalidated <- strain_repertoires %>%
  filter(n_production_groups > 0 & n_validated_production_groups == 0) %>%
  nrow()

# Additional summary
total_producers_raw <- sum(strain_repertoires$n_production_groups > 0)
total_producers_validated <- sum(strain_repertoires$n_validated_production_groups > 0)

# Generate report
report <- c(
  "Conservative Producer Category Counts",
  paste("Total strains:", total_strains),
  paste("Total producers (raw syn groups):", total_producers_raw),
  paste("Total producers (validated syn groups):", total_producers_validated),
  "",
  "Category 1: Strains with multiple syn groups, at least two validated (multi-producer conservative)",
  paste("  Count:", multi_conservative),
  paste("  Percentage of total strains:", sprintf("%.2f%%", 100 * multi_conservative / total_strains)),
  "",
  "Category 2: Strains with multiple syn groups, but only one or zero validated (single or non-validated producer)",
  paste("  Count:", single_or_non_validated),
  paste("  Percentage of total strains:", sprintf("%.2f%%", 100 * single_or_non_validated / total_strains)),
  "",
  "Category 3: Strains with only unvalidated syn groups (not included as producer nodes in conservative network)",
  paste("  Count:", only_unvalidated),
  paste("  Percentage of total strains:", sprintf("%.2f%%", 100 * only_unvalidated / total_strains)),
  "",
  "Notes:",
  "- Category 1 and 2 are subsets of strains with n_production_groups > 1",
  "- Category 3 includes strains with n_production_groups > 0 but no validated groups",
  "- Conservative network will include producers from Category 1 and validated parts of Category 2"
)

# Print to console
cat(paste(report, collapse = "\n"))
cat("\n")

# Save to file
writeLines(report, "scripts/phase_02_functional_groups/debug/conservative_producer_counts.txt")
cat("Saved report to scripts/phase_02_functional_groups/debug/conservative_producer_counts.txt\n")

cat("========================================\n")
cat("Debug complete!\n")
cat("========================================\n")
