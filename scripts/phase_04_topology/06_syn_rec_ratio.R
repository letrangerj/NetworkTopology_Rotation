# Phase 4 Step 6: Simplified Syn/Rec Ratio Analysis
# Focus on interpretable biological metrics without composite scoring
#
# Inputs:
# - data/interim/nodes_functional_groups_conservative.csv
# - data/interim/nodes_pyoverdines_conservative.csv
# - data/interim/edges_functional_groups_conservative.csv
# - results/phase_04/modularity/module_assignments_utilization_sk.csv
# - results/phase_04/nodf/module_nodf_analysis_utilization.csv
#
# Outputs:
# - results/phase_04/syn_rec_ratio/pyov_metrics_simplified.csv
# - results/phase_04/syn_rec_ratio/high_exploitation_candidates.csv
# - results/phase_04/syn_rec_ratio/high_ratio_candidates.csv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

# Paths
nodes_fg_path <- "data/interim/nodes_functional_groups_conservative.csv"
nodes_pyo_path <- "data/interim/nodes_pyoverdines_conservative.csv"
edges_fg_path <- "data/interim/edges_functional_groups_conservative.csv"
util_module_assign_path <- "results/phase_04/modularity/module_assignments_utilization_sk.csv"
util_module_nodf_path <- "results/phase_04/nodf/module_nodf_analysis_utilization.csv"

out_dir <- "results/phase_04/syn_rec_ratio"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading data...")
nodes_fg <- readr::read_csv(nodes_fg_path, show_col_types = FALSE)
nodes_pyo <- readr::read_csv(nodes_pyo_path, show_col_types = FALSE)
edges <- readr::read_csv(edges_fg_path, show_col_types = FALSE)
util_modules <- readr::read_csv(util_module_assign_path, show_col_types = FALSE)
util_mod_nodf <- readr::read_csv(util_module_nodf_path, show_col_types = FALSE)

# Get pyoverdine utilization module assignments
pyo_util_modules <- util_modules %>%
  filter(node_type == "PYOV", layer == "utilization") %>%
  select(node_id, module_id)

# Get module statistics
util_mod_stats <- util_mod_nodf %>%
  select(module_id, n_rows, n_cols, nodf_total, fill_percentage) %>%
  rename(
    module_fg_count = n_rows,
    module_pyo_count = n_cols,
    module_nestedness = nodf_total,
    module_fill_pct = fill_percentage
  )

# Separate production vs utilization edges
prod_edges <- edges %>%
  filter(edge_type == "production", source_type == "functional_group", target_type == "pyoverdine") %>%
  select(fg_id = source, pyo_id = target)

util_edges <- edges %>%
  filter(edge_type == "utilization", source_type == "pyoverdine", target_type == "functional_group") %>%
  select(pyo_id = source, fg_id = target)

# FG lookup for strategy classification
fg_lookup <- nodes_fg %>%
  select(agent_id, strain_count, strategy, has_production, has_utilization)

# Build utilization breakdown by strategy
message("Computing utilization breakdown by functional group strategy...")

# Create producer lookup for each pyoverdine
fg_produces_pyo <- prod_edges %>%
  distinct(fg_id, pyo_id) %>%
  mutate(produces_this_pyo = TRUE)

# Classify utilizers by their production strategy
util_breakdown <- util_edges %>%
  left_join(fg_lookup, by = c("fg_id" = "agent_id")) %>%
  left_join(fg_produces_pyo, by = c("fg_id", "pyo_id")) %>%
  mutate(
    produces_this_pyo = coalesce(produces_this_pyo, FALSE),
    utilizer_type = case_when(
      produces_this_pyo ~ "cooperator",
      has_production ~ "crossfeeder",
      !has_production ~ "nonproducer",
      TRUE ~ "unknown"
    )
  )

# Aggregate metrics per pyoverdine
message("Calculating core metrics per pyoverdine...")

pyo_metrics <- util_breakdown %>%
  group_by(pyo_id) %>%
  summarise(
    # Basic counts (FG level)
    utilizers_fg_total = n_distinct(fg_id),
    cooperators_fg = sum(utilizer_type == "cooperator"),
    crossfeeders_fg = sum(utilizer_type == "crossfeeder"),
    nonproducers_fg = sum(utilizer_type == "nonproducer"),

    # Strain-weighted counts
    utilizers_strains_total = sum(strain_count),
    cooperators_strains = sum(strain_count[utilizer_type == "cooperator"]),
    crossfeeders_strains = sum(strain_count[utilizer_type == "crossfeeder"]),
    nonproducers_strains = sum(strain_count[utilizer_type == "nonproducer"]),
    .groups = "drop"
  ) %>%
  # Add fractions
  mutate(
    # FG-level fractions
    cooperators_fg_frac = cooperators_fg / utilizers_fg_total,
    crossfeeders_fg_frac = crossfeeders_fg / utilizers_fg_total,
    nonproducers_fg_frac = nonproducers_fg / utilizers_fg_total,

    # Strain-level fractions
    cooperators_strain_frac = cooperators_strains / utilizers_strains_total,
    crossfeeders_strain_frac = crossfeeders_strains / utilizers_strains_total,
    nonproducers_strain_frac = nonproducers_strains / utilizers_strains_total
  )

# Get producer counts per pyoverdine
prod_counts <- prod_edges %>%
  left_join(fg_lookup, by = c("fg_id" = "agent_id")) %>%
  group_by(pyo_id) %>%
  summarise(
    producers_fg_total = n_distinct(fg_id),
    producers_strains_total = sum(strain_count),
    .groups = "drop"
  )

# Get basic pyoverdine info
pyo_info <- nodes_pyo %>%
  select(
    pyo_id = node_id,
    pyo_label = label,
    original_pyov_id,
    fg_utilization_count,
    strain_utilization_count,
    fg_production_count,
    strain_production_count
  )

# Combine all metrics
final_metrics <- pyo_info %>%
  left_join(prod_counts, by = "pyo_id") %>%
  left_join(pyo_metrics, by = "pyo_id") %>%
  left_join(pyo_util_modules, by = c("pyo_id" = "node_id")) %>%
  left_join(util_mod_stats, by = "module_id") %>%
  # Fill missing producer counts with 0
  mutate(
    producers_fg_total = coalesce(producers_fg_total, 0L),
    producers_strains_total = coalesce(producers_strains_total, 0L),
    utilizers_fg_total = coalesce(utilizers_fg_total, 0L),
    utilizers_strains_total = coalesce(utilizers_strains_total, 0L)
  ) %>%
  # Calculate syn/rec ratios
  mutate(
    # Core ratios (utilizers/producers)
    fg_util_prod_ratio = ifelse(producers_fg_total > 0,
      utilizers_fg_total / producers_fg_total,
      Inf
    ),
    strain_util_prod_ratio = ifelse(producers_strains_total > 0,
      utilizers_strains_total / producers_strains_total,
      Inf
    ),

    # Alternative: receptors/synthetases (from nodes_pyo data)
    fg_rec_syn_ratio = ifelse(fg_production_count > 0,
      fg_utilization_count / fg_production_count,
      Inf
    ),
    strain_rec_syn_ratio = ifelse(strain_production_count > 0,
      strain_utilization_count / strain_production_count,
      Inf
    ),

    # Module context flags
    is_large_module = coalesce(module_fg_count >= 30, FALSE),
    is_nested_module = coalesce(module_nestedness >= 50, FALSE),
    is_major_module = module_id %in% c(0, 3, 8) # From Phase 4 results
  ) %>%
  # Clean up column order
  select(
    pyo_id, pyo_label, original_pyov_id,

    # Producer counts
    producers_fg_total, producers_strains_total,
    fg_production_count, strain_production_count,

    # Utilizer counts
    utilizers_fg_total, utilizers_strains_total,
    fg_utilization_count, strain_utilization_count,

    # Core ratios
    fg_util_prod_ratio, strain_util_prod_ratio,
    fg_rec_syn_ratio, strain_rec_syn_ratio,

    # Utilizer breakdown (counts)
    cooperators_fg, crossfeeders_fg, nonproducers_fg,
    cooperators_strains, crossfeeders_strains, nonproducers_strains,

    # Utilizer breakdown (fractions)
    cooperators_fg_frac, crossfeeders_fg_frac, nonproducers_fg_frac,
    cooperators_strain_frac, crossfeeders_strain_frac, nonproducers_strain_frac,

    # Module context
    module_id, module_fg_count, module_pyo_count, module_nestedness, module_fill_pct,
    is_large_module, is_nested_module, is_major_module
  ) %>%
  arrange(desc(strain_util_prod_ratio), desc(nonproducers_strain_frac))

# Save main results
out_main <- file.path(out_dir, "pyov_metrics_simplified.csv")
readr::write_csv(final_metrics, out_main)
message("Wrote: ", out_main)

# Identify high exploitation candidates
high_exploit <- final_metrics %>%
  filter(
    nonproducers_strain_frac >= 0.5, # At least 50% nonproducers
    utilizers_strains_total >= 10 # Minimum usage threshold
  ) %>%
  arrange(desc(nonproducers_strain_frac), desc(strain_util_prod_ratio)) %>%
  head(20)

out_exploit <- file.path(out_dir, "high_exploitation_candidates.csv")
readr::write_csv(high_exploit, out_exploit)
message("Wrote: ", out_exploit)

# Identify high ratio candidates (imbalanced usage vs production)
high_ratio <- final_metrics %>%
  filter(
    strain_util_prod_ratio >= 5, # At least 5x more utilizers than producers
    utilizers_strains_total >= 20 # Minimum usage threshold
  ) %>%
  arrange(desc(strain_util_prod_ratio), desc(utilizers_strains_total)) %>%
  head(20)

out_ratio <- file.path(out_dir, "high_ratio_candidates.csv")
readr::write_csv(high_ratio, out_ratio)
message("Wrote: ", out_ratio)

# Summary statistics
message("\n=== SUMMARY STATISTICS ===")
cat("Total pyoverdines analyzed:", nrow(final_metrics), "\n")
cat("Pyoverdines with >50% nonproducer utilizers:", sum(final_metrics$nonproducers_strain_frac >= 0.5, na.rm = TRUE), "\n")
cat("Pyoverdines with >10x utilizer/producer ratio:", sum(final_metrics$strain_util_prod_ratio >= 10, na.rm = TRUE), "\n")
cat("Pyoverdines in major modules (0,3,8):", sum(final_metrics$is_major_module, na.rm = TRUE), "\n")

# Top candidates by different criteria
message("\nTop 10 by strain utilizer/producer ratio:")
print(final_metrics %>%
  select(pyo_id, strain_util_prod_ratio, nonproducers_strain_frac, utilizers_strains_total, module_id) %>%
  head(10), n = 10)

message("\nTop 10 by nonproducer fraction:")
print(final_metrics %>%
  arrange(desc(nonproducers_strain_frac)) %>%
  select(pyo_id, nonproducers_strain_frac, strain_util_prod_ratio, utilizers_strains_total, module_id) %>%
  head(10), n = 10)

message("\nTop 10 by FG utilizer/producer ratio:")
print(final_metrics %>%
  arrange(desc(fg_util_prod_ratio)) %>%
  select(pyo_id, fg_util_prod_ratio, nonproducers_fg_frac, utilizers_fg_total, module_id) %>%
  head(10), n = 10)
