#!/usr/bin/env Rscript
#
# diagnostics_step1_detailed.R
#
# Detailed diagnostics for Phase 3 Step 1:
# - Build component -> pair_id mapping (from pairing_result)
# - Produce pair-level summaries (components per pair, active components, FG/strain coverage)
# - Deep inspection of FG_067: full member lists, repertoires, receptor frequencies, pair_id mapping
# - Save CSVs and a short plain-text summary for quick review
#
# Usage (project root):
#   Rscript scripts/phase_03_network/debug/diagnostics_step1_detailed.R
#
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(readr)
  library(stringr)
})

# Paths and directories
debug_dir <- "scripts/phase_03_network/debug"
dir.create(debug_dir, recursive = TRUE, showWarnings = FALSE)

# Input paths (expected)
fg_nodes_path <- "data/interim/nodes_functional_groups_conservative.rds"
strain_nodes_path <- "data/interim/nodes_strains_conservative.rds"
strain_rep_path <- "data/interim/strain_repertoires.rds"
receptor_to_pyov_path <- "data/interim/receptor_to_pyoverdine.rds"
pairing_candidates <- c(
  "data/interim/rds/pairing_result.rds",
  "data/interim/rds/pairing_result_v7.rds",
  "data/interim/rds/pairing_result_v6.rds",
  "data/interim/rds/pairing_result.mat.rds",
  "data/interim/pairing_result.rds",
  "data/interim/pairing_result_v7.rds"
)

# Helper to read or stop with clear message
safe_readRDS <- function(p) {
  if (!file.exists(p)) stop(sprintf("Required file not found: %s", p))
  readRDS(p)
}

# Verify inputs exist
required_inputs <- c(fg_nodes_path, strain_nodes_path, strain_rep_path, receptor_to_pyov_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop("Missing required input RDS files:\n", paste(missing_inputs, collapse = "\n"),
       "\nPlease ensure Phase 1/2 outputs exist in data/interim/.")
}

# Load baseline inputs
fg <- safe_readRDS(fg_nodes_path)
strain_nodes <- safe_readRDS(strain_nodes_path)
strain_rep <- safe_readRDS(strain_rep_path)
receptor_to_pyov <- safe_readRDS(receptor_to_pyov_path)

# Attempt to load pairing_result (mandatory for pair-level mapping)
pairing_path <- NULL
for (p in pairing_candidates) {
  if (file.exists(p)) {
    pairing_path <- p
    break
  }
}
if (is.null(pairing_path)) {
  stop("pairing_result RDS not found in expected locations. Cannot perform pair-level diagnostics. Searched:\n",
       paste(pairing_candidates, collapse = "\n"))
}
pairing <- safe_readRDS(pairing_path)

cat("Loaded inputs:\n")
cat(" - FG nodes:", fg_nodes_path, "\n")
cat(" - Strain nodes:", strain_nodes_path, "\n")
cat(" - Strain repertoires:", strain_rep_path, "\n")
cat(" - receptor_to_pyov:", receptor_to_pyov_path, "\n")
cat(" - pairing_result:", pairing_path, "\n\n")

# Defensive: understanding pairing structure
# pairing usually is a matrix/data.frame with nrow = 26, column 1 = synthetase group lists, column 2 = receptor group lists
if (!is.matrix(pairing) && !is.data.frame(pairing)) {
  stop("Unexpected structure for pairing_result: expected matrix or data.frame where each cell is a list.")
}

# Build pair_map: component_pyov -> pair_id
cat("Building pair_map (component pyov -> pair_id)...\n")
pair_map <- map_df(seq_len(nrow(pairing)), function(i) {
  # defensive extraction: pairing[i,1][[1]] may be a vector or list
  syn_cell <- NULL
  # try several access patterns
  try({
    syn_cell <- pairing[i, 1][[1]]
  }, silent = TRUE)
  if (is.null(syn_cell)) {
    # fallback: try column name first if present
    if ("pyoverdine" %in% colnames(pairing)) {
      syn_cell <- pairing[i, "pyoverdine"][[1]]
    } else {
      syn_cell <- pairing[i, 1]
    }
  }
  syn_vec <- if (is.list(syn_cell)) unlist(syn_cell) else syn_cell
  syn_vec <- syn_vec[!is.na(syn_vec)]
  if (length(syn_vec) == 0) {
    return(tibble(pair_id = integer(0), component_pyov = integer(0)))
  }
  tibble(pair_id = i, component_pyov = as.integer(syn_vec))
})

# Save pair_map
# Determine which components are active in FG-level sets (production or utilization)
fg_components_prod <- unique(unlist(fg$validated_production_set, recursive = TRUE))
fg_components_util <- unique(unlist(fg$usable_pyoverdine_set, recursive = TRUE))
components_active <- sort(unique(c(as.integer(fg_components_prod), as.integer(fg_components_util))))
components_active <- components_active[!is.na(components_active)]

# Annotate pair_map with active flag (pair_map exists in memory)
pair_map <- pair_map %>%
  mutate(active_in_FG = component_pyov %in% components_active)

# Pair-level activity summary: how many components active per pair
pair_activity <- pair_map %>%
  group_by(pair_id) %>%
  summarize(
    n_components = n(),
    n_active_components = sum(active_in_FG),
    prop_active = ifelse(n_components > 0, n_active_components / n_components, 0),
    active_components_sample = paste(head(sort(component_pyov[active_in_FG]), 10), collapse = ","),
    .groups = "drop"
  ) %>%
  arrange(desc(prop_active), pair_id)

# Save minimal essential pair-level outputs
pair_activity_file <- file.path(debug_dir, "pair_activity_summary.csv")
write_csv(pair_activity, pair_activity_file)
cat("Saved pair activity summary to:", pair_activity_file, "\n")

# Map components -> pair_id lookup vector for quick mapping (kept in memory)
component_to_pairid <- pair_map %>% distinct(component_pyov, pair_id) %>% arrange(component_pyov)

# Compute per-pair FG-level and strain-level coverage (needed for active pair detection)
cat("Computing per-pair FG and strain coverage (aggregated to pair_id)...\n")
fg_prod_expanded <- fg %>%
  select(agent_id, functional_group_id, strain_count, validated_production_set) %>%
  mutate(validated_production_set = map(validated_production_set, ~ if (length(.x) == 0) integer(0) else as.integer(.x))) %>%
  unnest_longer(validated_production_set) %>%
  rename(component_pyov = validated_production_set) %>%
  filter(!is.na(component_pyov)) %>%
  left_join(component_to_pairid, by = "component_pyov") %>%
  filter(!is.na(pair_id))

fg_util_expanded <- fg %>%
  select(agent_id, functional_group_id, strain_count, usable_pyoverdine_set) %>%
  mutate(usable_pyoverdine_set = map(usable_pyoverdine_set, ~ if (length(.x) == 0) integer(0) else as.integer(.x))) %>%
  unnest_longer(usable_pyoverdine_set) %>%
  rename(component_pyov = usable_pyoverdine_set) %>%
  filter(!is.na(component_pyov)) %>%
  left_join(component_to_pairid, by = "component_pyov") %>%
  filter(!is.na(pair_id))

pair_fg_prod_counts <- fg_prod_expanded %>%
  group_by(pair_id) %>%
  summarize(n_fg_producing = n_distinct(agent_id),
            total_strains_producing = sum(strain_count), .groups = "drop")

pair_fg_util_counts <- fg_util_expanded %>%
  group_by(pair_id) %>%
  summarize(n_fg_utilizing = n_distinct(agent_id),
            total_strains_utilizing = sum(strain_count), .groups = "drop")

# Compose pair_coverage without writing the verbose component-level CSVs
pair_coverage <- pair_activity %>%
  left_join(pair_fg_prod_counts, by = "pair_id") %>%
  left_join(pair_fg_util_counts, by = "pair_id") %>%
  mutate(
    n_fg_producing = replace_na(n_fg_producing, 0),
    total_strains_producing = replace_na(total_strains_producing, 0),
    n_fg_utilizing = replace_na(n_fg_utilizing, 0),
    total_strains_utilizing = replace_na(total_strains_utilizing, 0)
  ) %>%
  arrange(desc(total_strains_utilizing + total_strains_producing))

pair_coverage_file <- file.path(debug_dir, "pair_level_coverage_summary.csv")
write_csv(pair_coverage, pair_coverage_file)
cat("Saved pair-level coverage summary to:", pair_coverage_file, "\n")

# Save active/inactive pair id lists (essential)
active_pair_ids <- pair_coverage %>% filter((n_fg_producing + n_fg_utilizing) > 0) %>% pull(pair_id)
inactive_pair_ids <- setdiff(pair_activity$pair_id, active_pair_ids)
write_csv(tibble(active_pair_id = active_pair_ids), file.path(debug_dir, "active_pair_ids.csv"))
write_csv(tibble(inactive_pair_id = inactive_pair_ids), file.path(debug_dir, "inactive_pair_ids.csv"))
cat("Saved active_pair_ids.csv and inactive_pair_ids.csv\n")

# Compute how many of the 26 pair rows are active (have any FG production or utilization)
active_pair_ids <- pair_coverage %>% filter((n_fg_producing + n_fg_utilizing) > 0) %>% pull(pair_id)
inactive_pair_ids <- setdiff(pair_level_summary$pair_id, active_pair_ids)
cat(sprintf("Active pair rows: %d / %d\n", length(active_pair_ids), nrow(pair_level_summary)))
write_csv(tibble(active_pair_id = active_pair_ids), file.path(debug_dir, "active_pair_ids.csv"))
write_csv(tibble(inactive_pair_id = inactive_pair_ids), file.path(debug_dir, "inactive_pair_ids.csv"))

# Deep inspection of FG_067 (if present)
if ("FG_067" %in% fg$agent_id) {
  cat("Deep inspection: FG_067 (trimmed outputs)\n")
  fg67 <- fg %>% filter(agent_id == "FG_067")
  fg67_members <- unlist(fg67$member_strains)
  # Save member list (essential)
  write_csv(tibble(strainName = fg67_members), file.path(debug_dir, "FG_067_member_strains_full.csv"))

  # Extract repertoires for members
  if (!"strainName" %in% names(strain_rep)) {
    alt_names <- names(strain_rep)
    stop("strain_repertoires does not contain 'strainName' column. Available columns: ", paste(alt_names, collapse = ", "))
  }

  fg67_reps <- strain_rep %>%
    filter(strainName %in% fg67_members) %>%
    select(strainName, production_groups = production_groups, validated_production_groups = validated_production_groups, receptor_groups = receptor_groups, usable_pyoverdines = usable_pyoverdines)

  # Compute receptor group frequency (essential)
  fg67_receptors_long <- fg67_reps %>%
    mutate(receptor_groups = map(receptor_groups, ~ if (length(.x) == 0) integer(0) else as.integer(.x))) %>%
    unnest_longer(receptor_groups) %>%
    filter(!is.na(receptor_groups)) %>%
    rename(receptor_group = receptor_groups) %>%
    group_by(receptor_group) %>%
    summarize(count = n(), .groups = "drop") %>%
    arrange(desc(count))

  write_csv(fg67_receptors_long, file.path(debug_dir, "FG_067_receptor_group_freq.csv"))

  # Production component frequency (validated components) and map to pair_id, then summarize by pair_id (essential)
  fg67_prod_comp_long <- fg67_reps %>%
    mutate(validated_production_groups = map(validated_production_groups, ~ if (length(.x) == 0) integer(0) else as.integer(.x))) %>%
    unnest_longer(validated_production_groups) %>%
    filter(!is.na(validated_production_groups)) %>%
    rename(component_pyov = validated_production_groups) %>%
    group_by(component_pyov) %>%
    summarize(count = n(), .groups = "drop") %>%
    arrange(desc(count))

  fg67_prod_pairs <- fg67_prod_comp_long %>%
    left_join(component_to_pairid, by = c("component_pyov" = "component_pyov")) %>%
    group_by(pair_id) %>%
    summarize(n_comp = n(), total_count = sum(count), components = paste(component_pyov, collapse = ","), .groups = "drop") %>%
    arrange(desc(total_count))

  write_csv(fg67_prod_pairs, file.path(debug_dir, "FG_067_production_pairid_freq.csv"))

  # Utilization component frequency and map to pair_id (essential)
  fg67_util_comp_long <- fg67_reps %>%
    mutate(usable_pyoverdines = map(usable_pyoverdines, ~ if (length(.x) == 0) integer(0) else as.integer(.x))) %>%
    unnest_longer(usable_pyoverdines) %>%
    filter(!is.na(usable_pyoverdines)) %>%
    rename(component_pyov = usable_pyoverdines) %>%
    group_by(component_pyov) %>%
    summarize(count = n(), .groups = "drop") %>%
    arrange(desc(count))

  fg67_util_pairs <- fg67_util_comp_long %>%
    left_join(component_to_pairid, by = c("component_pyov" = "component_pyov")) %>%
    group_by(pair_id) %>%
    summarize(n_comp = n(), total_count = sum(count), components = paste(component_pyov, collapse = ","), .groups = "drop") %>%
    arrange(desc(total_count))

  write_csv(fg67_util_pairs, file.path(debug_dir, "FG_067_util_pairid_freq.csv"))

  cat("FG_067 diagnostics (trimmed) saved to debug directory.\n")
} else {
  cat("FG_067 not present in FG node table; skipping FG_067 deep inspection.\n")
}

# Final plain-text summary
summary_lines <- list()
summary_lines <- c(summary_lines, sprintf("Detailed diagnostics run: %s", Sys.time()))
summary_lines <- c(summary_lines, sprintf("FG node rows: %d", nrow(fg)))
summary_lines <- c(summary_lines, sprintf("Sum strains (from FG rows): %d", sum(fg$strain_count)))
summary_lines <- c(summary_lines, sprintf("Number of pair rows (lock-key pairs): %d", length(unique(pair_map$pair_id))))
summary_lines <- c(summary_lines, sprintf("Unique components in pair_map: %d", nrow(component_to_pairid)))
summary_lines <- c(summary_lines, sprintf("Validated components active in FG sets (prod or util): %d", length(components_active)))
summary_lines <- c(summary_lines, sprintf("Pair rows active (any FG production or utilization): %d", length(active_pair_ids)))
summary_lines <- c(summary_lines, sprintf("Pair rows inactive: %d", length(inactive_pair_ids)))
# Top pairs by total strains util+prod (trimmed to top 5)
top_pairs <- pair_coverage %>% mutate(total = total_strains_producing + total_strains_utilizing) %>% arrange(desc(total)) %>% slice_head(n = 5) %>%
  transmute(pair_id, n_components = n_components, n_active_components = n_active_components, total_strains_producing, total_strains_utilizing, total)
summary_lines <- c(summary_lines, "", "Top pair rows by total strains (prod+util) [top 5]:")
summary_lines <- c(summary_lines, capture.output(print(top_pairs)))

summary_file <- file.path(debug_dir, "diagnostics_step1_detailed_summary.txt")
write_lines(summary_lines, summary_file)
cat("Detailed diagnostics summary written to:", summary_file, "\n")

cat("\nAll detailed diagnostics complete. Essential files written to:", debug_dir, "\n")
cat("Essential outputs:\n")
cat(" - pair_activity_summary.csv\n - pair_level_coverage_summary.csv\n - active_pair_ids.csv\n")
if ("FG_067" %in% fg$agent_id) {
  cat(" - FG_067_member_strains_full.csv\n - FG_067_receptor_group_freq.csv\n - FG_067_production_pairid_freq.csv\n - FG_067_util_pairid_freq.csv\n")
}
# End of script
