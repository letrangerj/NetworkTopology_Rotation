# 01_strain_identification.R
# Phase 2: Functional Group Classification - Step 1
# Extract and validate strain identifiers from syn.rds and rec.rds
# Run from project root: Rscript scripts/phase_02_functional_groups/01_strain_identification.R

library(dplyr)

# Set random seed for reproducibility
set.seed(42)

# Load Phase 1 RDS files
cat("Loading RDS files...\n")
syn <- readRDS("data/interim/rds/syn.rds")
rec <- readRDS("data/interim/rds/rec.rds")

# Extract strain identifiers
cat("Extracting strain identifiers...\n")

# From syn.rds: strainName field
syn_strains <- unique(syn$strainName)
cat("Unique strains in syn.rds:", length(syn_strains), "\n")
cat("Sample syn strains:", paste(head(syn_strains, 5), collapse = ", "), "\n")

# From rec.rds: foldername field (contains strain names)
rec_strains <- unique(rec$foldername)
cat("Unique strains in rec.rds:", length(rec_strains), "\n")
cat("Sample rec strains:", paste(head(rec_strains, 5), collapse = ", "), "\n")

# Validate consistency
cat("\nValidating strain identifier consistency...\n")

# Check overlap
common_strains <- intersect(syn_strains, rec_strains)
only_syn <- setdiff(syn_strains, rec_strains)
only_rec <- setdiff(rec_strains, syn_strains)

cat("Strains in both datasets:", length(common_strains), "\n")
cat("Strains only in syn.rds:", length(only_syn), "\n")
cat("Strains only in rec.rds:", length(only_rec), "\n")

if (length(only_syn) > 0) {
  cat("Sample strains only in syn.rds:", paste(head(only_syn, 3), collapse = ", "), "\n")
}
if (length(only_rec) > 0) {
  cat("Sample strains only in rec.rds:", paste(head(only_rec, 3), collapse = ", "), "\n")
}

# Check if foldername matches strainName pattern
# For rec.rds, foldername should be the strain identifier
# Verify this by checking if foldername values match expected strain format

# Save strain lists for next steps
strain_data <- list(
  syn_strains = syn_strains,
  rec_strains = rec_strains,
  common_strains = common_strains,
  only_syn = only_syn,
  only_rec = only_rec
)

saveRDS(strain_data, "data/interim/strain_identifiers.rds")
cat("\nSaved strain identifier data to data/interim/strain_identifiers.rds\n")

# Summary
cat("\nStrain Identification Summary:\n")
cat("- Total unique strains across both datasets:", length(unique(c(syn_strains, rec_strains))), "\n")
cat("- Strains with synthetase data:", length(syn_strains), "\n")
cat("- Strains with receptor data:", length(rec_strains), "\n")
cat("- Strains with both:", length(common_strains), "\n")

cat("\nStrain identification complete. Proceed to repertoire calculation.\n")
