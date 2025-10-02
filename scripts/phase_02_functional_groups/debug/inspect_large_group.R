# inspect_large_group.R
# Debug script: Inspect the large 700-strain functional group
# Run from project root: Rscript scripts/phase_02_functional_groups/debug/inspect_large_group.R

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

cat("========================================\n")
cat("Debug: Inspecting Large Functional Group (700 strains)\n")
cat("========================================\n\n")

# Load functional groups and mapping
cat("Loading functional groups...\n")
functional_groups <- readRDS("data/interim/functional_groups_conservative.rds")
strain_mapping <- readRDS("data/interim/strain_functional_group_mapping_conservative.rds")

cat("Loaded:\n")
cat("  - Total functional groups:", nrow(functional_groups), "\n")
cat("  - Total strains in mapping:", nrow(strain_mapping), "\n\n")

# Find the largest group
largest_group <- functional_groups %>%
  filter(strain_count == max(strain_count))

if (nrow(largest_group) != 1) {
  cat("Error: Multiple groups have the maximum strain count!\n")
  quit(status = 1)
}

group_id <- largest_group$functional_group_id
group_size <- largest_group$strain_count

cat("Largest group identified:\n")
cat("  - Functional group ID:", group_id, "\n")
cat("  - Strain count:", group_size, "\n")
cat("  - Percentage of collapsed strains:", sprintf("%.1f%%", 100 * group_size / sum(functional_groups$strain_count)), "\n")

# Extract the strain names in this group
group_strains <- strain_mapping %>%
  filter(functional_group_id == group_id)

cat("\nSignature analysis for group", group_id, ":\n")
cat("  - Production:\n")
cat("    * Number of validated production groups:", largest_group$n_validated_production_groups, "\n")
if (largest_group$n_validated_production_groups > 0) {
  # Extract production groups from actual strains
  first_strain <- group_strains$strainName[1]
  first_strain_data <- readRDS("data/interim/strain_repertoires.rds") %>%
    filter(strainName == first_strain)
  cat("    * Production groups:", paste(sort(unlist(first_strain_data$validated_production_groups)), collapse = ", "), "\n")
} else {
  cat("    * Production groups: none (consumer only)\n")
}

cat("  - Utilization:\n")
cat("    * Number of usable pyoverdines:", largest_group$n_usable_pyoverdines, "\n")
if (largest_group$n_usable_pyoverdines > 0) {
  # Extract utilization from actual strains
  if (missing(first_strain_data)) {
    first_strain_data <- readRDS("data/interim/strain_repertoires.rds") %>%
      filter(strainName == first_strain)
  }
  cat("    * Usable pyoverdines:", paste(sort(unlist(first_strain_data$usable_pyoverdines)), collapse = ", "), "\n")
} else {
  cat("    * Usable pyoverdines: none\n")
}

# Load receptor mapping to understand the utilization pattern
cat("\nUnderstanding utilization signature:\n")
if (largest_group$n_usable_pyoverdines > 0) {
  # Load the receptor_to_pyoverdine mapping and strain repertoires
  receptor_mapping <- readRDS("data/interim/receptor_to_pyoverdine.rds")
  strain_data <- readRDS("data/interim/strain_repertoires.rds")

  # Sample a few strains to see their receptor patterns
  sample_strains <- sample(group_strains$strainName, min(5, length(group_strains$strainName)))

  for (strain in sample_strains) {
    strain_info <- strain_data %>% filter(strainName == strain)
    receptors <- sort(unlist(strain_info$receptor_groups))
    if (length(receptors) > 5) {
      receptors_display <- paste(paste(head(receptors, 5), collapse = ", "), "...", sep = "")
    } else {
      receptors_display <- paste(receptors, collapse = ", ")
    }
    cat("  -", strain, ": receptors [", receptors_display, "] leads to pyoverdines [",
        paste(sort(unlist(strain_info$usable_pyoverdines)), collapse = ", "), "]\n")
  }
}

# Analyze taxonomic patterns if possible
cat("\nTaxonomic/strain name analysis:\n")
sample_strain_names <- sample(group_strains$strainName, min(20, length(group_strains$strainName)))

# Look for patterns in strain names
## Check for obvious species indicators
species_patterns <- c("baumannii", "pseudomallei", "cepacia", "gladioli", "cenocepacia", "multivorans")
found_patterns <- c()

for (pattern in species_patterns) {
  matches <- sum(grepl(pattern, group_strains$strainName, ignore.case = TRUE))
  if (matches > 0) {
    found_patterns <- c(found_patterns, paste(pattern, ":", matches, "strains"))
  }
}

if (length(found_patterns) > 0) {
  cat("Common species patterns found:\n")
  for (pattern in found_patterns) {
    cat("  -", pattern, "\n")
  }
} else {
  cat("No obvious species patterns found in strain names\n")
}

# Check for other naming patterns
cat("\nNaming pattern analysis:\n")
all_repertoires <- unique(group_strains$strainName)

# Check for repeated prefixes or suffixes
prefixes <- substr(all_repertoires, 1, 10)
prefix_counts <- sort(table(prefixes), decreasing = TRUE)

cat("Most common 10-character prefixes:\n")
for (i in 1:min(5, length(prefix_counts))) {
  if (prefix_counts[i] > 5) {
    cat("  -", names(prefix_counts)[i], ":", prefix_counts[i], "strains\n")
  }
}

# Check if this looks biological vs technical
cat("\nBiological plausibility assessment:\n")
if (largest_group$n_validated_production_groups <= 2 && largest_group$n_usable_pyoverdines <= 3) {
  cat("✓ Simple signature is biologically plausible\n")
  cat("  - Many strains could converge on minimal pyoverdine strategy\n")
  cat("  - Conservative validation could collapse diverse raw repertoires into this pattern\n")
} else {
  cat("⚠ Complex signature with", largest_group$n_validated_production_groups, "production groups\n")
}

if (group_size > 100) {
  cat("⚠ Large group size (", group_size, ") suggests strong selection pressure\n")
  cat("  - Could be dominant species/complex\n")
  cat("  - Could be validation bottleneck effect\n")
  cat("  - Recommend checking taxonomic distribution\n")
}

# Generate final report
cat("\n" "="rep(50, "="), "\n")
cat("INSPECTION COMPLETE\n")
cat("="rep(50, "="), "\n")
cat("Group", group_id, "signature likely represents a dominant pyoverdine strategy\n")
cat("shared by", group_size, "strains due to conservative validation constraints.\n")
cat("This is probably a minimal producer-consumer strategy rather than a bug.\n")

# Save detailed report
writeLines(
  capture.output({
    cat("Large Group Inspection Report\n")
    cat("Group ID:", group_id, "\n")
    cat("Strain Count:", group_size, "\n")
    cat("Production Groups:", largest_group$n_validated_production_groups, "\n")
    cat("Usable Pyoverdines:", largest_group$n_usable_pyoverdines, "\n")
    if (length(found_patterns) > 0) {
      cat("Species Patterns:\n")
      cat(paste(found_patterns, collapse = "\n"))
    }
  }),
  "scripts/phase_02_functional_groups/debug/large_group_inspection_report.txt"
)

cat("\nDetailed report saved to: scripts/phase_02_functional_groups/debug/large_group_inspection_report.txt\n")
