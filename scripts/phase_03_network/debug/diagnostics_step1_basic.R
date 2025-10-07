#!/usr/bin/env Rscript
#
# diagnostics_step1_basic.R
#
# Basic diagnostics for Phase 3 Step 1 agent sets:
# - Functional group size distribution (top groups, cumulative coverage)
# - Pyoverdine / production-group counts (raw vs validated)
# - Active validated pyoverdine components and pair-row activity
# - Save diagnostic tables into scripts/phase_03_network/debug/
#
# Run from project root:
#   Rscript scripts/phase_03_network/debug/diagnostics_step1_basic.R
#
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(readr)
})

# Set deterministic behavior for sampling (if any)
set.seed(42)

# Paths (relative to project root)
debug_dir <- "scripts/phase_03_network/debug"
dir.create(debug_dir, recursive = TRUE, showWarnings = FALSE)

# Input files (expected)
fg_path <- "data/interim/nodes_functional_groups_conservative.rds"
strain_nodes_path <- "data/interim/nodes_strains_conservative.rds"
strain_rep_path <- "data/interim/strain_repertoires.rds"
receptor_to_pyov_path <- "data/interim/receptor_to_pyoverdine.rds"
# pairing_result may be in one of two locations depending on earlier steps
pairing_paths <- c(
  "data/interim/rds/pairing_result.rds",
  "data/interim/rds/pairing_result_v7.rds",
  "data/interim/rds/pairing_result.rds",
  "data/interim/rds/pairing_result.rds" # duplicate fallback harmless
)

# Helper to attempt reading a file or abort with a friendly message
safe_readRDS <- function(p) {
  if (!file.exists(p)) stop(sprintf("Required file not found: %s", p))
  readRDS(p)
}

# Load required inputs with graceful errors
cat("Loading input RDS files...\n")
if (!file.exists(fg_path) || !file.exists(strain_nodes_path) || !file.exists(strain_rep_path) || !file.exists(receptor_to_pyov_path)) {
  stop("One or more required input files are missing in data/interim/. Please run Phase 1/2 scripts first.")
}

fg <- safe_readRDS(fg_path)
strain_nodes <- safe_readRDS(strain_nodes_path)
strain_rep <- safe_readRDS(strain_rep_path)
receptor_to_pyov <- safe_readRDS(receptor_to_pyov_path)

# Attempt to load pairing_result if present
pairing <- NULL
for (pp in pairing_paths) {
  if (file.exists(pp)) {
    pairing <- readRDS(pp)
    pairing_source <- pp
    break
  }
}
if (is.null(pairing)) {
  cat("Warning: pairing_result RDS not found in expected locations. Pair-level mapping diagnostics will be omitted.\n")
} else {
  cat("Loaded pairing_result from:", pairing_source, "\n")
}

# 1) FG size distribution and top groups
cat("Computing FG size distribution...\n")
fg_summary <- fg %>%
  arrange(desc(strain_count)) %>%
  mutate(rank = row_number()) %>%
  mutate(cumulative_strains = cumsum(strain_count),
         cumulative_fraction = cumulative_strains / sum(strain_count))

# Save only the succinct top-10 FG summary to limit output size
top10 <- fg_summary %>% slice(1:10)
write_csv(top10, file.path(debug_dir, "fg_top10.csv"))
cat("Top 10 FG written to:", file.path(debug_dir, "fg_top10.csv"), "\n")

# Compute fraction of strains in the largest FG (kept for console output only)
largest_fraction <- fg_summary$cumulative_fraction[1]
cat(sprintf("Fraction of strains in largest FG: %.4f\n", largest_fraction))

# 2) Raw production group counts vs validated components
cat("Computing production-group / pyoverdine counts...\n")
# production_groups may be named differently in strain_repertoires; try common names
if ("production_groups" %in% names(strain_rep)) {
  raw_prod_list <- unique(unlist(strain_rep$production_groups))
} else if ("validated_production_groups" %in% names(strain_rep)) {
  raw_prod_list <- unique(unlist(strain_rep$validated_production_groups, recursive = TRUE))
} else {
  raw_prod_list <- integer(0)
  cat("Warning: no production_groups or validated_production_groups field in strain_repertoires.\n")
}

n_raw_prod <- length(unique(raw_prod_list))
cat("Unique raw production group IDs found (across strains):", n_raw_prod, "\n")

# Validated components from receptor_to_pyov (pyoverdine_group)
validated_components <- sort(unique(as.integer(receptor_to_pyov$pyoverdine_group)))
n_validated_components <- length(validated_components)
cat("Unique validated pyoverdine component IDs (from receptor_to_pyov):", n_validated_components, "\n")

# If we have pairing, compute pair-level counts (number of rows)
if (!is.null(pairing)) {
  # pairing likely is a matrix/data.frame with two columns (synthetase list, receptor list)
  n_pair_rows <- nrow(pairing)
  cat("Pairing result contains rows (lock-key pairs):", n_pair_rows, "\n")
} else {
  n_pair_rows <- NA
}

# Determine how many validated components are active in FG-level sets (production OR usable)
fg_components_prod <- unique(unlist(fg$validated_production_set, recursive = TRUE))
fg_components_util <- unique(unlist(fg$usable_pyoverdine_set, recursive = TRUE))
components_active <- sort(unique(c(as.integer(fg_components_prod), as.integer(fg_components_util))))
n_components_active <- length(components_active)
cat("Validated components active in FG sets (prod or util):", n_components_active, "\n")

# Validated-but-inactive components
validated_inactive <- setdiff(validated_components, components_active)
n_validated_inactive <- length(validated_inactive)
cat("Validated components not present in FG-level prod/util sets:", n_validated_inactive, "\n")

# Save concise pyov summary (do not generate per-component inactive file to limit outputs)
pyov_summary <- tibble(
  metric = c("n_raw_production_groups", "n_validated_components", "n_pair_rows", "n_components_active", "n_validated_inactive"),
  value = c(n_raw_prod, n_validated_components, ifelse(is.na(n_pair_rows), -1, n_pair_rows), n_components_active, n_validated_inactive)
)
write_csv(pyov_summary, file.path(debug_dir, "pyov_summary.csv"))

# 3) Map component -> pair_id (if pairing available) and compute active pair rows
if (!is.null(pairing)) {
  cat("Building component -> pair_id map from pairing_result...\n")
  # pairing[i,1] expected to be synthetase component list; defensive unwrapping
  pair_map <- map_df(seq_len(nrow(pairing)), function(i) {
    syn_cell <- pairing[i, 1][[1]]
    syn_vec <- if (is.list(syn_cell)) unlist(syn_cell) else syn_cell
    syn_vec <- syn_vec[!is.na(syn_vec)]
    if (length(syn_vec) == 0) return(tibble(pair_id = integer(0), component_pyov = integer(0)))
    tibble(pair_id = i, component_pyov = as.integer(syn_vec))
  })
  write_csv(pair_map, file.path(debug_dir, "pair_map_component_to_pairid.csv"))

  # Which pair_ids have components that are active in FG sets?
  active_pair_ids <- sort(unique(pair_map$pair_id[pair_map$component_pyov %in% components_active]))
  inactive_pair_ids <- setdiff(unique(pair_map$pair_id), active_pair_ids)

  cat("Active pair rows (count):", length(active_pair_ids), "\n")
  cat("Inactive pair rows (count):", length(inactive_pair_ids), "\n")

  # Save active/inactive pair ids
  write_csv(tibble(active_pair_id = active_pair_ids), file.path(debug_dir, "active_pair_ids.csv"))
  write_csv(tibble(inactive_pair_id = inactive_pair_ids), file.path(debug_dir, "inactive_pair_ids.csv"))
} else {
  cat("Skipping pair-level mapping because pairing_result was not found.\n")
}

# 4) Inspect FG_067 in detail (if present)
if ("FG_067" %in% fg$agent_id) {
  cat("Inspecting FG_067 (large group) in detail...\n")
  fg67 <- fg %>% filter(agent_id == "FG_067")
  members <- unlist(fg67$member_strains)
  # Save member list
  write_csv(tibble(strainName = members), file.path(debug_dir, "FG_067_member_strains.csv"))
  # Pull their repertoires (validated_production_groups, usable_pyoverdines)
  if ("strainName" %in% names(strain_rep)) {
    fg67_reps <- strain_rep %>% filter(strainName %in% members) %>%
      select(strainName, production_groups = production_groups, validated_production_groups = validated_production_groups, usable_pyoverdines = usable_pyoverdines)
    # For readability, collapse list-columns to strings for CSV
    fg67_csv <- fg67_reps %>%
      mutate(
        production_groups = map_chr(production_groups, ~ if (length(.x) == 0) "" else paste(sort(as.integer(.x)), collapse = ";")),
        validated_production_groups = map_chr(validated_production_groups, ~ if (length(.x) == 0) "" else paste(sort(as.integer(.x)), collapse = ";")),
        usable_pyoverdines = map_chr(usable_pyoverdines, ~ if (length(.x) == 0) "" else paste(sort(as.integer(.x)), collapse = ";"))
      )
    write_csv(fg67_csv, file.path(debug_dir, "FG_067_members_repertoires_sample.csv"))
  } else {
    cat("strain_repertoires does not contain field 'strainName'—cannot extract FG_067 member repertoires.\n")
  }

  # Compute frequency of production components and usable pyovs within FG_067
  fg67_prod_freq <- strain_rep %>% filter(strainName %in% members) %>%
    pull(validated_production_groups) %>% flatten_int() %>% table() %>% as.data.frame(stringsAsFactors = FALSE)
  names(fg67_prod_freq) <- c("component", "count")
  write_csv(as_tibble(fg67_prod_freq), file.path(debug_dir, "FG_067_production_component_freq.csv"))

  fg67_util_freq <- strain_rep %>% filter(strainName %in% members) %>%
    pull(usable_pyoverdines) %>% flatten_int() %>% table() %>% as.data.frame(stringsAsFactors = FALSE)
  names(fg67_util_freq) <- c("component", "count")
  write_csv(as_tibble(fg67_util_freq), file.path(debug_dir, "FG_067_util_component_freq.csv"))

} else {
  cat("FG_067 not found in FG node table; skipping FG_067 inspection.\n")
}

# 5) Summary report (plain text) for quick review
report_lines <- list()
report_lines <- c(report_lines, sprintf("Diagnostics run: %s", Sys.time()))
report_lines <- c(report_lines, sprintf("Total functional groups (FG rows): %d", nrow(fg)))
report_lines <- c(report_lines, sprintf("Total conservative strains (sum FG strain_count): %d", sum(fg$strain_count)))
report_lines <- c(report_lines, sprintf("Largest FG fraction (top 1): %.4f", largest_fraction))
report_lines <- c(report_lines, sprintf("Unique raw production groups (approx): %d", n_raw_prod))
report_lines <- c(report_lines, sprintf("Unique validated pyov components (from receptor_to_pyov): %d", n_validated_components))
report_lines <- c(report_lines, sprintf("Validated pyov components active in FG sets: %d", n_components_active))
report_lines <- c(report_lines, sprintf("Validated components not present in FG sets: %d", n_validated_inactive))
if (!is.null(pairing)) {
  report_lines <- c(report_lines, sprintf("Pair-level rows (lock-key pairs) in pairing_result: %d", n_pair_rows))
  report_lines <- c(report_lines, sprintf("Active pair rows (components mapped to FG sets): %d", length(active_pair_ids)))
  report_lines <- c(report_lines, sprintf("Inactive pair rows: %d", length(inactive_pair_ids)))
}
write_lines(report_lines, file.path(debug_dir, "diagnostics_step1_basic_summary.txt"))

cat("Diagnostics complete. Outputs written to:", debug_dir, "\n")
cat("Key output files:\n")
cat(" - fg_top10.csv\n - pyov_summary.csv\n - diagnostics_step1_basic_summary.txt\n")
if ("FG_067" %in% fg$agent_id) cat(" - FG_067_member_strains.csv\n - FG_067_members_repertoires_sample.csv\n")

# End of script
