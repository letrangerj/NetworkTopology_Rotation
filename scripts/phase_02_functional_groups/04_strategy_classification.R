# 04_strategy_classification.R
# Phase 2: Functional Group Classification - Step 4
# Classify functional groups by social strategy using conservative criteria
# Run from project root: Rscript scripts/phase_02_functional_groups/04_strategy_classification.R

suppressPackageStartupMessages({
  library(dplyr)
})

# Set random seed for reproducibility
set.seed(42)

cat("========================================\n")
cat("Phase 2 - Step 4: Strategy Classification (Conservative)\n")
cat("========================================\n\n")

# Load conservative functional groups
cat("Loading functional groups...\n")
functional_groups <- readRDS("data/interim/functional_groups_conservative.rds")

cat("Functional groups loaded:\n")
cat("  - Total functional groups:", nrow(functional_groups), "\n")
cat("  - Total strains represented:", sum(functional_groups$strain_count), "\n\n")

# ---- Define strategy classification rules (Conservative) ----
cat("Applying conservative strategy classification rules...\n")

functional_groups_classified <- functional_groups %>%
  mutate(
    # Primary strategy classification
    strategy = case_when(
      # Multi-producer: Validated production of >1 pyoverdine group
      is_validated_multi_producer ~ "Multi-producer",

      # Single-receptor producer: 1 validated production group, 1 usable pyoverdine
      n_validated_production_groups == 1 & n_usable_pyoverdines == 1 ~ "Single-receptor producer",

      # Multi-receptor producer: 1 validated production group, >1 usable pyoverdine
      n_validated_production_groups == 1 & n_usable_pyoverdines > 1 ~ "Multi-receptor producer",

      # Nonproducer: No validated production, but has usable pyoverdines
      n_validated_production_groups == 0 & n_usable_pyoverdines > 0 ~ "Nonproducer",

      # Edge case: Producer with no usable pyoverdines (no receptors that match their production)
      n_validated_production_groups > 0 & n_usable_pyoverdines == 0 ~ "Producer-only (no utilization)",

      # Edge case: Neither producer nor consumer in conservative network
      n_validated_production_groups == 0 & n_usable_pyoverdines == 0 ~ "Unconnected",

      # Fallback for unexpected combinations
      TRUE ~ "Unclassified"
    ),

    # Strategy category for network analysis (simplified)
    strategy_category = case_when(
      strategy %in% c("Multi-producer", "Single-receptor producer", "Multi-receptor producer") ~ "Producer",
      strategy == "Nonproducer" ~ "Nonproducer",
      TRUE ~ "Other"
    ),

    # Detailed description
    strategy_description = paste(
      "Production:", ifelse(n_validated_production_groups > 0,
                           paste(n_validated_production_groups, "validated groups"),
                           "none"),
      "| Utilization:", ifelse(n_usable_pyoverdines > 0,
                              paste(n_usable_pyoverdines, "pyoverdines"),
                              "none")
    )
  )

# ---- Validate classification ----
cat("Validating strategy assignments...\n")

# Count strategies
strategy_counts <- functional_groups_classified %>%
  count(strategy, name = "group_count") %>%
  mutate(strain_count = sapply(strategy, function(s) {
    sum(functional_groups_classified$strain_count[functional_groups_classified$strategy == s])
  }))

cat("Strategy distribution:\n")
for (i in 1:nrow(strategy_counts)) {
  cat(sprintf("  - %-25s: %3d groups (%4d strains, %.1f%% of groups)\n",
              strategy_counts$strategy[i],
              strategy_counts$group_count[i],
              strategy_counts$strain_count[i],
              100 * strategy_counts$group_count[i] / nrow(functional_groups_classified)))
}

# Check for unclassified groups
unclassified <- functional_groups_classified %>% filter(strategy == "Unclassified")
if (nrow(unclassified) > 0) {
  cat(sprintf("\nWARNING: %d groups could not be classified!\n", nrow(unclassified)))
  cat("Sample unclassified signatures:\n")
  print(head(select(unclassified, conservative_canonical_signature, strategy_description), 10))
}

# ---- Save outputs ----
cat("\nSaving outputs...\n")

# Ensure directories exist
dir.create("data/interim", recursive = TRUE, showWarnings = FALSE)
dir.create("docs/phase_02", recursive = TRUE, showWarnings = FALSE)

# Save classified functional groups
saveRDS(functional_groups_classified, "data/interim/functional_groups_classified_conservative.rds")
cat("Saved: data/interim/functional_groups_classified_conservative.rds\n")

# Save strategy summary as CSV for easy viewing
strategy_summary <- functional_groups_classified %>%
  select(functional_group_id, conservative_canonical_signature, strain_count,
         n_validated_production_groups, n_usable_pyoverdines,
         strategy, strategy_category, strategy_description)

write.csv(strategy_summary, "data/interim/strategy_classification_conservative.csv", row.names = FALSE)
cat("Saved: data/interim/strategy_classification_conservative.csv\n")

# ---- Generate comprehensive summary report ----
cat("\nGenerating summary report...\n")

report <- c(
  "Phase 2 - Step 4: Strategy Classification Summary (Conservative)",
  paste("Run timestamp:", Sys.time()),
  "",
  "Classification Rules Applied:",
  "  1. Multi-producer: >1 validated production group",
  "  2. Single-receptor producer: 1 production group, 1 usable pyoverdine",
  "  3. Multi-receptor producer: 1 production group, >1 usable pyoverdine",
  "  4. Nonproducer: 0 production groups, >0 usable pyoverdines",
  "  5. Producer-only: >0 production groups, 0 usable pyoverdines",
  "  6. Unconnected: 0 production groups, 0 usable pyoverdines",
  "",
  "Overall Statistics:",
  paste("  - Total functional groups classified:", nrow(functional_groups_classified)),
  paste("  - Total strains represented:", sum(functional_groups_classified$strain_count)),
  paste("  - Average group size:", sprintf("%.2f", mean(functional_groups_classified$strain_count))),
  "",
  "Strategy Distribution (Groups):"
)

# Add strategy counts to report
for (i in 1:nrow(strategy_counts)) {
  report <- c(report,
    paste(sprintf("  - %-25s: %3d groups (%5.1f%%)",
                  strategy_counts$strategy[i],
                  strategy_counts$group_count[i],
                  100 * strategy_counts$group_count[i] / nrow(functional_groups_classified)))
  )
}

report <- c(report,
  "",
  "Strategy Distribution (Strains):"
)

for (i in 1:nrow(strategy_counts)) {
  report <- c(report,
    paste(sprintf("  - %-25s: %4d strains (%5.1f%%)",
                  strategy_counts$strategy[i],
                  strategy_counts$strain_count[i],
                  100 * strategy_counts$strain_count[i] / sum(functional_groups_classified$strain_count)))
  )
}

report <- c(report,
  "",
  "Network-Ready Categories:",
  paste("  - Producer groups (all types):", sum(functional_groups_classified$strategy_category == "Producer")),
  paste("  - Nonproducer groups:", sum(functional_groups_classified$strategy_category == "Nonproducer")),
  paste("  - Other/Unconnected groups:", sum(functional_groups_classified$strategy_category == "Other")),
  "",
  "Next Steps:",
  "  - Proceed to Phase 3: Network Construction",
  "  - Use functional_groups_classified_conservative.rds for bipartite network building",
  "  - Validate network topology against biological expectations",
  ""
)

writeLines(report, "docs/phase_02/strategy_classification_summary.txt")
cat("Saved summary report to docs/phase_02/strategy_classification_summary.txt\n")

# Print summary to console
cat("\n")
cat(paste(report, collapse = "\n"))
cat("\n")

cat("========================================\n")
cat("Strategy classification complete!\n")
cat("========================================\n")
