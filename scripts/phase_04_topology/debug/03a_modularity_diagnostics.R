#!/usr/bin/env Rscript

# Debug script: Step 03a - Modularity Diagnostics
# ------------------------------------------------
# Purpose:
#  - Compute per-node and per-module co-assignment stability from the coassignment matrix
#  - Compute pairwise replicate similarity (Adjusted Rand Index)
#  - Produce summary CSVs and diagnostic plots
#
# Inputs (produced by 03_modularity_analysis.R):
#  - results/phase_04/modularity/modules_fg_raw_replicates_memberships.csv
#  - results/phase_04/modularity/coassignment_matrix.rds
#  - results/phase_04/modularity/modules_fg_consensus.csv
#  - results/phase_04/modularity/modularity_scores.csv
#
# Outputs (diagnostics):
#  - results/phase_04/modularity/diagnostics/node_stability_by_consensus_module.csv
#  - results/phase_04/modularity/diagnostics/module_stability_summary.csv
#  - results/phase_04/modularity/diagnostics/replicate_pairwise_ARI.csv
#  - results/phase_04/modularity/diagnostics/replicate_pairwise_ARI_hist.pdf
#  - figures/network_topology/step03_modularity_figures/diagnostics/node_stability_boxplot_by_module.pdf
#  - figures/network_topology/step03_modularity_figures/diagnostics/module_mean_coassignment_bar.pdf
#  - figures/network_topology/step03_modularity_figures/diagnostics/replicate_pairwise_ARI_hist.pdf
#  - results/phase_04/modularity/diagnostics/diagnostics_summary.txt
#
# Usage:
#  Rscript scripts/phase_04_topology/debug/03a_modularity_diagnostics.R
#
# Notes:
#  - This script is intended for interactive debugging and exploratory assessment of
#    modularity stability. It is intentionally conservative about stopping on missing files
#    and reports diagnostic outputs into a dedicated diagnostics folder.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(mclust) # adjustedRandIndex
})

# -----------------------------
# Paths
# -----------------------------
replicate_memberships_csv <- "results/phase_04/modularity/modules_fg_raw_replicates_memberships.csv"
coassign_rds             <- "results/phase_04/modularity/coassignment_matrix.rds"
consensus_csv            <- "results/phase_04/modularity/modules_fg_consensus.csv"
mod_scores_csv           <- "results/phase_04/modularity/modularity_scores.csv"

out_diag_dir  <- "results/phase_04/modularity/diagnostics"
fig_diag_dir  <- "figures/network_topology/step03_modularity_figures/diagnostics"

# Ensure output directories exist
if (!dir.exists(out_diag_dir)) dir.create(out_diag_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(fig_diag_dir)) dir.create(fig_diag_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Sanity checks for inputs
# -----------------------------
required_files <- c(replicate_memberships_csv, coassign_rds, consensus_csv, mod_scores_csv)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files for diagnostics:\n - ", paste(missing_files, collapse = "\n - "),
       "\nRun Step 03 modularity analysis and ensure outputs exist.")
}

# -----------------------------
# Load inputs
# -----------------------------
message("Loading replicate memberships...")
memberships <- read_csv(replicate_memberships_csv, show_col_types = FALSE)
if (!"node_id" %in% colnames(memberships)) {
  stop("Expected first column 'node_id' in replicate memberships CSV.")
}
node_ids <- memberships$node_id
rep_df <- memberships %>% select(-node_id)

message("Loading co-assignment matrix...")
coassign <- readRDS(coassign_rds)
if (!is.matrix(coassign)) {
  stop("Co-assignment object is not a matrix. Check: ", coassign_rds)
}
# Ensure coassign has dimnames
if (is.null(rownames(coassign)) || is.null(colnames(coassign))) {
  stop("coassignment matrix missing row/col names.")
}

message("Loading consensus assignments...")
consensus <- read_csv(consensus_csv, show_col_types = FALSE)
if (!all(c("node_id", "consensus_module") %in% colnames(consensus))) {
  stop("Consensus CSV must contain 'node_id' and 'consensus_module' columns.")
}
consensus_map <- consensus$consensus_module
names(consensus_map) <- consensus$node_id

message("Loading modularity scores...")
mod_scores <- read_csv(mod_scores_csv, show_col_types = FALSE)
if (!"modularity" %in% colnames(mod_scores)) {
  warning("modularity scores CSV does not contain 'modularity' column; proceeding but some summaries will be NA.")
}

# -----------------------------
# Validate node sets align
# -----------------------------
# Nodes present in consensus and coassign should match membership nodes
nodes_in_coassign <- rownames(coassign)
if (!all(node_ids %in% nodes_in_coassign)) {
  missing_nodes <- setdiff(node_ids, nodes_in_coassign)
  warning("Some nodes from membership file are not present in co-assignment matrix. Missing count: ", length(missing_nodes),
          ". These nodes will be skipped in co-assignment diagnostics.")
}
# Work with intersection for safety
common_nodes <- intersect(node_ids, nodes_in_coassign)
if (length(common_nodes) == 0) stop("No overlapping nodes between membership list and co-assignment matrix.")

# -----------------------------
# Node-level stability (within consensus module)
# -----------------------------
message("Computing node-level stability within consensus modules...")
modules_list <- split(consensus$node_id, consensus$consensus_module)
node_stability_records <- list()

for (mod in names(modules_list)) {
  members <- modules_list[[mod]]
  members_present <- intersect(members, common_nodes)
  if (length(members_present) == 0) {
    next
  }
  if (length(members_present) == 1) {
    # Single-member module: cannot compute within-module co-assignment
    node_stability_records[[mod]] <- tibble(
      node_id = members_present,
      consensus_module = as.integer(mod),
      mean_within = NA_real_,
      median_within = NA_real_,
      sd_within = NA_real_
    )
  } else {
    submat <- coassign[members_present, members_present, drop = FALSE]
    for (mnode in members_present) {
      others <- setdiff(members_present, mnode)
      vals <- as.numeric(submat[mnode, others])
      node_stability_records[[mod]] <- bind_rows(
        node_stability_records[[mod]] %||% tibble(),
        tibble(
          node_id = mnode,
          consensus_module = as.integer(mod),
          mean_within = mean(vals, na.rm = TRUE),
          median_within = median(vals, na.rm = TRUE),
          sd_within = sd(vals, na.rm = TRUE)
        )
      )
    }
  }
}

# bind rows
node_stability_df <- bind_rows(node_stability_records)
# If some consensus nodes were not present in coassign, add them with NA
missing_consensus_nodes <- setdiff(consensus$node_id, node_stability_df$node_id)
if (length(missing_consensus_nodes) > 0) {
  node_stability_df <- bind_rows(
    node_stability_df,
    tibble(node_id = missing_consensus_nodes,
           consensus_module = as.integer(consensus_map[missing_consensus_nodes]),
           mean_within = NA_real_,
           median_within = NA_real_,
           sd_within = NA_real_)
  )
}

# Save node-level stability
write_csv(node_stability_df, file.path(out_diag_dir, "node_stability_by_consensus_module.csv"))

# Plot: boxplot of node mean_within by module
p_node_box <- ggplot(node_stability_df, aes(x = factor(consensus_module), y = mean_within)) +
  geom_boxplot(fill = "lightblue", na.rm = TRUE) +
  labs(x = "Consensus module", y = "Node mean within-module co-assignment",
       title = "Node stability (mean within-module co-assignment) by consensus module") +
  theme_minimal()
ggsave(file.path(fig_diag_dir, "node_stability_boxplot_by_module.pdf"), p_node_box, width = 8, height = 4)

# -----------------------------
# Module-level stability
# -----------------------------
message("Computing module-level stability metrics...")
module_stability <- node_stability_df %>%
  group_by(consensus_module) %>%
  summarise(module_mean_within = mean(mean_within, na.rm = TRUE),
            module_median_within = median(mean_within, na.rm = TRUE),
            n_members = n()) %>%
  arrange(desc(module_mean_within))

write_csv(module_stability, file.path(out_diag_dir, "module_stability_summary.csv"))

# Plot: barplot of module mean within
p_module_bar <- ggplot(module_stability, aes(x = reorder(as.factor(consensus_module), -module_mean_within), y = module_mean_within)) +
  geom_col(fill = "orange") +
  labs(x = "Consensus module", y = "Mean within-module co-assignment",
       title = "Module-level mean co-assignment (stability)") +
  theme_minimal()
ggsave(file.path(fig_diag_dir, "module_mean_coassignment_bar.pdf"), p_module_bar, width = 8, height = 4)

# -----------------------------
# Replicate similarity (pairwise ARI)
# -----------------------------
message("Computing pairwise ARI between replicates...")
rep_matrix <- as.matrix(rep_df)
n_reps <- ncol(rep_matrix)
if (n_reps < 2) {
  warning("Less than 2 replicates found; cannot compute pairwise ARI.")
  ari_vals <- numeric(0)
  ari_df <- tibble(rep_pair = character(0), ARI = numeric(0))
} else {
  comb_idx <- combn(n_reps, 2)
  ari_vals <- apply(comb_idx, 2, function(idx) {
    a <- rep_matrix[, idx[1]]
    b <- rep_matrix[, idx[2]]
    mclust::adjustedRandIndex(a, b)
  })
  ari_df <- tibble(rep_pair = apply(comb_idx, 2, function(x) paste0("rep", x[1], "_rep", x[2])), ARI = ari_vals)
  write_csv(ari_df, file.path(out_diag_dir, "replicate_pairwise_ARI.csv"))

  # Plot ARI histogram
  p_ari <- ggplot(ari_df, aes(x = ARI)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    labs(title = "Pairwise ARI between Louvain replicates", x = "Adjusted Rand Index", y = "Count") +
    theme_minimal()
  ggsave(file.path(fig_diag_dir, "replicate_pairwise_ARI_hist.pdf"), p_ari, width = 6, height = 4)
}

# -----------------------------
# Diagnostics summary
# -----------------------------
message("Writing diagnostics summary...")
summary_lines <- c()
summary_lines <- c(summary_lines, paste0("Timestamp: ", Sys.time()))
summary_lines <- c(summary_lines, sprintf("Replicates: %d", n_reps))
if ("modularity" %in% colnames(mod_scores)) {
  summary_lines <- c(summary_lines, sprintf("Mean modularity: %.4f", mean(mod_scores$modularity, na.rm = TRUE)))
  summary_lines <- c(summary_lines, sprintf("SD modularity: %.4f", sd(mod_scores$modularity, na.rm = TRUE)))
} else {
  summary_lines <- c(summary_lines, "Mean modularity: NA (modularity column missing)")
  summary_lines <- c(summary_lines, "SD modularity: NA (modularity column missing)")
}
summary_lines <- c(summary_lines, sprintf("Mean pairwise ARI: %s", ifelse(length(ari_vals) > 0, format(mean(ari_vals, na.rm = TRUE), digits = 4), "NA")))
summary_lines <- c(summary_lines, sprintf("Median pairwise ARI: %s", ifelse(length(ari_vals) > 0, format(median(ari_vals, na.rm = TRUE), digits = 4), "NA")))
summary_lines <- c(summary_lines, "")
summary_lines <- c(summary_lines, "Top module stability summary (descending module_mean_within):")
summary_lines <- c(summary_lines, capture.output(print(head(module_stability, 10))))
summary_lines <- c(summary_lines, "")
summary_lines <- c(summary_lines, "Notes:")
summary_lines <- c(summary_lines, " - Node stability (mean_within) close to 1 indicates core nodes; values < 0.5 indicate weak membership.")
summary_lines <- c(summary_lines, " - Low mean pairwise ARI (< 0.5) suggests high variability across replicates; consider increasing replicate count or switching to Leiden.")
writeLines(summary_lines, file.path(out_diag_dir, "diagnostics_summary.txt"))

message("Diagnostics complete. Outputs saved to: ", out_diag_dir)
message("Figures saved to: ", fig_diag_dir)

# End of script
