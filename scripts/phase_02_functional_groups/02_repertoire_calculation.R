# 02_repertoire_calculation.R
# Phase 2: Functional Group Classification - Step 2
# Calculate production and utilization repertoires per strain
# Run from project root: Rscript scripts/phase_02_functional_groups/02_repertoire_calculation.R

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
})

# Set random seed for reproducibility
set.seed(42)

cat("========================================\n")
cat("Phase 2 - Step 2: Repertoire Calculation\n")
cat("========================================\n\n")

# Load Phase 1 RDS files
cat("Loading RDS files...\n")
syn <- readRDS("data/interim/rds/syn.rds")
rec <- readRDS("data/interim/rds/rec.rds")
pairing <- readRDS("data/interim/rds/pairing_result.rds")

cat("Loaded:\n")
cat("  syn: ", length(syn$strainName), "entries\n")
cat("  rec: ", length(rec$foldername), "entries\n")
cat("  pairing dimensions: ", paste(dim(pairing), collapse=" x "), "\n\n")

# ---- Extract validated pairing pairs ----
# pairing is a 26x2 matrix where each cell contains a list with a numeric vector
# Column 1: synthetase groups (pyoverdine groups)
# Column 2: receptor groups

cat("Extracting validated pairs from pairing_result...\n")

# Extract all synthetase groups (column 1) and receptor groups (column 2)
# Each cell pairing[i,j][[1]] contains a numeric vector
synthetase_groups_list <- list()
receptor_groups_list <- list()

for (i in 1:nrow(pairing)) {
  # Defensive extraction: the matrix cells contain lists whose first element is
  # itself a numeric vector (e.g. pairing[i,1][[1]] gives c(53,54,...))
  raw_syn_cell <- pairing[i, 1][[1]]
  raw_rec_cell <- pairing[i, 2][[1]]

  # Unwrap nested lists to atomic vectors (handles list(list(c(...))) and similar)
  syn_vec <- if (is.list(raw_syn_cell)) unlist(raw_syn_cell) else raw_syn_cell
  rec_vec <- if (is.list(raw_rec_cell)) unlist(raw_rec_cell) else raw_rec_cell

  # Defensive coercion to integer; use suppressWarnings to avoid noisy warnings
  syn_int <- if (length(syn_vec) == 0) integer(0) else suppressWarnings(as.integer(syn_vec))
  rec_int <- if (length(rec_vec) == 0) integer(0) else suppressWarnings(as.integer(rec_vec))

  # Diagnostics: report if coercion produced NAs (non-numeric content)
  if (any(is.na(syn_int) & !is.na(syn_vec))) {
    cat(sprintf("Warning: non-numeric synthetase values in pairing row %d: %s\n", i, paste(head(syn_vec, 10), collapse = ",")))
  }
  if (any(is.na(rec_int) & !is.na(rec_vec))) {
    cat(sprintf("Warning: non-numeric receptor values in pairing row %d: %s\n", i, paste(head(rec_vec, 10), collapse = ",")))
  }

  # Store list-columns (possibly empty integer vectors)
  synthetase_groups_list[[i]] <- syn_int
  receptor_groups_list[[i]] <- rec_int
}

# Create lookup table: each receptor group maps to its paired pyoverdine group
# Based on the pairing structure, each row i represents a lock-key pair
valid_pairs <- data.frame(
  pair_id = 1:nrow(pairing),
  pyoverdine_group = I(synthetase_groups_list),  # Use I() to store as list-column
  receptor_group = I(receptor_groups_list)
)

cat("Found", nrow(valid_pairs), "validated pyoverdine-receptor pairs\n")

# Expand to create receptor -> pyoverdine lookup
# For each pair, all receptors in that group can use all pyoverdines in that group
receptor_to_pyoverdine <- valid_pairs %>%
  mutate(
    mappings = map2(pyoverdine_group, receptor_group, ~tidyr::expand_grid(pyoverdine_group = .x, receptor_group = .y))
  ) %>%
  select(pair_id, mappings) %>%
  tidyr::unnest(cols = c(mappings))

cat("Expanded to", nrow(receptor_to_pyoverdine), "receptor-pyoverdine mappings\n")

# Save lookup
saveRDS(receptor_to_pyoverdine, "data/interim/receptor_to_pyoverdine.rds")
cat("Saved receptor->pyoverdine lookup\n\n")

# ---- Build production repertoire per strain ----
cat("Calculating production repertoires per strain...\n")

# For each strain, find all non-zero synthetase groups
production_by_strain <- data.frame(
  strainName = syn$strainName,
  group = as.integer(syn$group),
  stringsAsFactors = FALSE
) %>%
  filter(group != 0 & !is.na(group)) %>%
  group_by(strainName) %>%
  summarize(
    production_groups = list(unique(group)),
    n_production_groups = length(unique(group)),
    .groups = 'drop'
  )

cat("Strains with production (non-zero groups):", nrow(production_by_strain), "\n")

# Check for strains with multiple production groups (should be rare/none per Gu et al.)
multi_prod <- production_by_strain %>% filter(n_production_groups > 1)
if (nrow(multi_prod) > 0) {
  cat("WARNING: Found", nrow(multi_prod), "strains with multiple production groups!\n")
  cat("Sample strains:", paste(head(multi_prod$strainName, 5), collapse = ", "), "\n")
}

# ---- Build receptor repertoire per strain ----
cat("\nCalculating receptor repertoires per strain...\n")

receptors_by_strain <- data.frame(
  strainName = rec$foldername,
  receptor_group = as.integer(rec$group),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(receptor_group)) %>%
  group_by(strainName) %>%
  summarize(
    receptor_groups = list(unique(receptor_group)),
    n_receptor_groups = length(unique(receptor_group)),
    .groups = 'drop'
  )

cat("Strains with receptors:", nrow(receptors_by_strain), "\n\n")

# ---- Combine and translate to usable pyoverdines ----
cat("Combining production and receptor data...\n")

# Get all unique strains (union)
all_strains <- unique(c(syn$strainName, rec$foldername))
cat("Total unique strains:", length(all_strains), "\n")

# Build strain repertoires
strain_repertoires <- data.frame(
  strainName = all_strains,
  stringsAsFactors = FALSE
) %>%
  left_join(production_by_strain, by = "strainName") %>%
  left_join(receptors_by_strain, by = "strainName") %>%
  mutate(
    production_groups = map(production_groups, ~ if(is.null(.x)) integer(0) else .x),
    receptor_groups = map(receptor_groups, ~ if(is.null(.x)) integer(0) else .x),
    n_production_groups = replace_na(n_production_groups, 0),
    n_receptor_groups = replace_na(n_receptor_groups, 0)
  )

cat("Translating receptor groups to usable pyoverdines...\n")

# For each strain, translate receptor groups to usable pyoverdines via lookup
# Only receptors in validated pairs count
strain_repertoires <- strain_repertoires %>%
  mutate(
    usable_pyoverdines = map(receptor_groups, function(rec_groups) {
      if (length(rec_groups) == 0) return(integer(0))
      # Find which pyoverdines these receptors can use
      matched <- receptor_to_pyoverdine %>%
        filter(receptor_group %in% rec_groups) %>%
        pull(pyoverdine_group) %>%
        unique()
      return(sort(matched))
    }),
    n_usable_pyoverdines = map_int(usable_pyoverdines, length)
  )

# ---- Validation checks ----
cat("\nRunning validation checks...\n")

# 1. Strains with multiple production groups
cat("  - Strains with >1 production group:", sum(strain_repertoires$n_production_groups > 1), "\n")

# 2. Receptors excluded (not in validated pairs)
all_receptor_groups <- unique(unlist(strain_repertoires$receptor_groups))
validated_receptor_groups <- unique(receptor_to_pyoverdine$receptor_group)
excluded_receptors <- setdiff(all_receptor_groups, validated_receptor_groups)
cat("  - Total unique receptor groups:", length(all_receptor_groups), "\n")
cat("  - Receptor groups in validated pairs:", length(validated_receptor_groups), "\n")
cat("  - Receptor groups excluded (not validated):", length(excluded_receptors), "\n")

# 3. Strains with receptors but no usable pyoverdines
strains_receptors_no_usable <- strain_repertoires %>%
  filter(n_receptor_groups > 0 & n_usable_pyoverdines == 0)
cat("  - Strains with receptors but NO usable pyoverdines:", nrow(strains_receptors_no_usable), "\n")

# ---- Create canonical signatures ----
cat("\nCreating canonical repertoire signatures...\n")

strain_repertoires <- strain_repertoires %>%
  mutate(
    production_signature = map_chr(production_groups, function(x) {
      if (length(x) == 0) return("0")
      paste(sort(x), collapse = ",")
    }),
    utilization_signature = map_chr(usable_pyoverdines, function(x) {
      if (length(x) == 0) return("none")
      paste(sort(x), collapse = ",")
    }),
    canonical_signature = paste(production_signature, utilization_signature, sep = "|")
  )

# ---- Save outputs ----
cat("\nSaving outputs...\n")
saveRDS(strain_repertoires, "data/interim/strain_repertoires.rds")
cat("Saved strain_repertoires.rds\n")

# ---- Generate summary report ----
cat("\nGenerating summary report...\n")

report <- c(
  "Phase 2 - Repertoire Calculation Summary",
  paste("Run timestamp:", Sys.time()),
  "",
  "Input Data:",
  paste("  - Synthetase entries:", length(syn$strainName)),
  paste("  - Receptor entries:", length(rec$foldername)),
  paste("  - Validated pairing pairs:", nrow(pairing)),
  paste("  - Expanded receptor-pyoverdine mappings:", nrow(receptor_to_pyoverdine)),
  "",
  "Strain Counts:",
  paste("  - Total unique strains:", nrow(strain_repertoires)),
  paste("  - Strains with production (non-zero):", sum(strain_repertoires$n_production_groups > 0)),
  paste("  - Strains with receptors:", sum(strain_repertoires$n_receptor_groups > 0)),
  paste("  - Strains with usable pyoverdines:", sum(strain_repertoires$n_usable_pyoverdines > 0)),
  "",
  "Quality Checks:",
  paste("  - Strains with >1 production group:", sum(strain_repertoires$n_production_groups > 1)),
  paste("  - Receptor groups excluded (not validated):", length(excluded_receptors)),
  paste("  - Strains with receptors but no usable pyoverdines:", nrow(strains_receptors_no_usable)),
  ""
)

if (nrow(multi_prod) > 0) {
  report <- c(report,
    "WARNING - Strains with multiple production groups:",
    paste("  ", head(multi_prod$strainName, 10), collapse = "\n"),
    ""
  )
}

if (length(excluded_receptors) > 0) {
  report <- c(report,
    "Excluded receptor groups (sample):",
    paste("  ", paste(head(excluded_receptors, 20), collapse = ", ")),
    ""
  )
}

if (nrow(strains_receptors_no_usable) > 0) {
  report <- c(report,
    "Strains with receptors but no usable pyoverdines (sample):",
    paste("  ", paste(head(strains_receptors_no_usable$strainName, 10), collapse = ", ")),
    ""
  )
}

report <- c(report,
  "Next Steps:",
  "  - Run 03_functional_group_collapse.R to group strains by canonical signatures",
  "  - Classify functional groups by strategy (Nonproducer/Single/Multi-receptor producer)",
  ""
)

# Create docs/phase_02 if needed
dir.create("docs/phase_02", recursive = TRUE, showWarnings = FALSE)

writeLines(report, "docs/phase_02/repertoire_calculation_summary.txt")
cat("Saved summary report to docs/phase_02/repertoire_calculation_summary.txt\n")

# Print summary to console
cat("\n")
cat(paste(report, collapse = "\n"))
cat("\n")

cat("========================================\n")
cat("Repertoire calculation complete!\n")
cat("========================================\n")
