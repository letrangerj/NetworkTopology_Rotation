# 03_functional_group_collapse.R
# Phase 2: Functional Group Classification - Step 3
# Collapse strains into functional groups using conservative (validated) criteria
# Run from project root: Rscript scripts/phase_02_functional_groups/03_functional_group_collapse.R

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

# Set random seed for reproducibility
set.seed(42)

cat("========================================\n")
cat("Phase 2 - Step 3: Functional Group Collapse (Conservative)\n")
cat("========================================\n\n")

# Load strain repertoires with conservative validation fields
cat("Loading strain repertoires...\n")
strain_repertoires <- readRDS("data/interim/strain_repertoires.rds")

cat("Strain repertoires loaded:\n")
cat("  - Total strains:", nrow(strain_repertoires), "\n")
cat("  - Strains with validated production:", sum(strain_repertoires$is_validated_producer), "\n")
cat("  - Strains with usable pyoverdines:", sum(strain_repertoires$n_usable_pyoverdines > 0), "\n\n")

# ---- Define conservative network criteria ----
# For conservative network, include strains that are either:
# - validated producers (will have production edges), OR
# - have usable pyoverdines (will have utilization edges), OR
# - both (producer-consumers)

conservative_network_strains <- strain_repertoires %>%
  filter(is_validated_producer | n_usable_pyoverdines > 0)

cat("Strains included in conservative network analysis:\n")
cat("  - Total:", nrow(conservative_network_strains), "\n")
cat("  - Validated producers:", sum(conservative_network_strains$is_validated_producer), "\n")
cat("  - With usable pyoverdines:", sum(conservative_network_strains$n_usable_pyoverdines > 0), "\n")
cat("  - Producer-consumers:", sum(conservative_network_strains$is_validated_producer &
                                   conservative_network_strains$n_usable_pyoverdines > 0), "\n\n")

# ---- Create conservative canonical signatures ----
cat("Creating conservative canonical signatures...\n")

conservative_strains <- conservative_network_strains %>%
  mutate(
    # Conservative production signature: only validated groups
    conservative_production_signature = map_chr(validated_production_groups, function(x) {
      if (length(x) == 0) return("0")
      paste(sort(x), collapse = ",")
    }),
    # Utilization signature remains the same (already uses validated pairs)
    conservative_utilization_signature = map_chr(usable_pyoverdines, function(x) {
      if (length(x) == 0) return("none")
      paste(sort(x), collapse = ",")
    }),
    # Full conservative canonical signature
    conservative_canonical_signature = paste(
      conservative_production_signature,
      conservative_utilization_signature,
      sep = "|"
    )
  )

# ---- Collapse strains by conservative signatures ----
cat("Collapsing strains into functional groups...\n")

functional_groups <- conservative_strains %>%
  group_by(conservative_canonical_signature) %>%
  summarize(
    functional_group_id = cur_group_id(),
    strain_count = n(),
    n_validated_production_groups = first(n_validated_production_groups),
    n_usable_pyoverdines = first(n_usable_pyoverdines),
    is_validated_producer = first(is_validated_producer),
    is_validated_multi_producer = first(is_validated_multi_producer),
    # Store strain names for mapping
    strain_names = list(strainName),
    .groups = 'drop'
  ) %>%
  arrange(functional_group_id)

cat("Collapse results:\n")
cat("  - Total functional groups:", nrow(functional_groups), "\n")
cat("  - Total strains collapsed:", sum(functional_groups$strain_count), "\n")
cat("  - Average collapse ratio:", sprintf("%.2f", mean(functional_groups$strain_count)), "strains per group\n")

# ---- Create strain-to-functional-group mapping ----
cat("\nCreating strain to functional group mapping...\n")

strain_mapping <- conservative_strains %>%
  left_join(
    functional_groups %>% select(conservative_canonical_signature, functional_group_id),
    by = "conservative_canonical_signature"
  ) %>%
  select(strainName, functional_group_id, conservative_canonical_signature)

# ---- Validation checks ----
cat("\nRunning validation checks...\n")

# Check that all strains in mapping have unique functional_group_id per signature
mapping_check <- strain_mapping %>%
  group_by(conservative_canonical_signature) %>%
  summarize(
    unique_groups = n_distinct(functional_group_id),
    .groups = 'drop'
  ) %>%
  filter(unique_groups > 1)

if (nrow(mapping_check) > 0) {
  cat("WARNING: Some signatures map to multiple functional_group_ids!\n")
} else {
  cat("✓ All signatures map to unique functional_group_ids\n")
}

# Check that all functional groups have consistent metadata
metadata_consistency <- functional_groups %>%
  group_by(functional_group_id) %>%
  summarize(
    consistent_prod = n_distinct(n_validated_production_groups) == 1,
    consistent_usable = n_distinct(n_usable_pyoverdines) == 1,
    .groups = 'drop'
  ) %>%
  filter(!consistent_prod | !consistent_usable)

if (nrow(metadata_consistency) > 0) {
  cat("WARNING: Some functional groups have inconsistent metadata!\n")
} else {
  cat("✓ All functional groups have consistent metadata\n")
}

# ---- Save outputs ----
cat("\nSaving outputs...\n")

# Ensure directories exist
dir.create("data/interim", recursive = TRUE, showWarnings = FALSE)
dir.create("docs/phase_02", recursive = TRUE, showWarnings = FALSE)

# Save functional groups and mapping
saveRDS(functional_groups, "data/interim/functional_groups_conservative.rds")
saveRDS(strain_mapping, "data/interim/strain_functional_group_mapping_conservative.rds")

cat("Saved:\n")
cat("  - data/interim/functional_groups_conservative.rds\n")
cat("  - data/interim/strain_functional_group_mapping_conservative.rds\n")

# ---- Generate summary report ----
cat("\nGenerating summary report...\n")

report <- c(
  "Phase 2 - Step 3: Functional Group Collapse Summary (Conservative)",
  paste("Run timestamp:", Sys.time()),
  "",
  "Input Data:",
  paste("  - Total strains analyzed:", nrow(conservative_strains)),
  paste("  - Total strains in dataset:", nrow(strain_repertoires)),
  "",
  "Conservative Network Criteria:",
  paste("  - Strains with validated production OR usable pyoverdines:"),
  paste("  - Validated producers included:", sum(conservative_strains$is_validated_producer)),
  paste("  - Strains with usable pyoverdines included:", sum(conservative_strains$n_usable_pyoverdines > 0)),
  "",
  "Collapse Results:",
  paste("  - Total functional groups created:", nrow(functional_groups)),
  paste("  - Total strains collapsed:", sum(functional_groups$strain_count)),
  paste("  - Average collapse ratio:", sprintf("%.2f", mean(functional_groups$strain_count)), "strains per group"),
  paste("  - Maximum group size:", max(functional_groups$strain_count)),
  paste("  - Minimum group size:", min(functional_groups$strain_count)),
  "",
  "Functional Group Composition:",
  paste("  - Groups with validated production:", sum(functional_groups$is_validated_producer)),
  paste("  - Groups without production (pure consumers):", sum(!functional_groups$is_validated_producer)),
  paste("  - Groups with usable pyoverdines:", sum(functional_groups$n_usable_pyoverdines > 0)),
  paste("  - Groups without usable pyoverdines:", sum(functional_groups$n_usable_pyoverdines == 0)),
  "",
  "Next Steps:",
  "  - Run 04_strategy_classification.R to classify functional groups by strategy",
  "  - Validate group assignments and metadata consistency",
  ""
)

writeLines(report, "docs/phase_02/functional_group_collapse_summary.txt")
cat("Saved summary report to docs/phase_02/functional_group_collapse_summary.txt\n")

# Print summary to console
cat("\n")
cat(paste(report, collapse = "\n"))
cat("\n")

cat("========================================\n")
cat("Functional group collapse complete!\n")
cat("========================================\n")
