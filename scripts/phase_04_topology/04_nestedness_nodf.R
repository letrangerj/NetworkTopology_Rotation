#!/usr/bin/env Rscript

# Phase 04 - Step 04: Nestedness (NODF) Analysis - Minimal Output
# ---------------------------------------------------------------
# Purpose:
#   Calculate nestedness metrics for FG-level and STR-level bipartite networks
#   using primary bipartite::nested() with vegan::nestednodf() cross-validation.
#   Minimal output: text summary and figures only (no RDS/CSV files).
#
# Outputs:
#   - figures/network_topology/step04_nodf_matrix_fg.pdf
#   - figures/network_topology/step04_nodf_matrix_str.pdf
#   - figures/network_topology/step04_nodf_comparison.pdf
#   - docs/phase_04/logs/step04_nestedness_summary.txt
#
# Usage:
#   Rscript scripts/phase_04_topology/04_nestedness_nodf.R

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(bipartite)
  library(vegan)
  library(pheatmap)
})

# -----------------------------
# Helper functions
# -----------------------------
safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

# Deterministic incidence matrix builder
build_incidence_matrix <- function(adj_production, adj_utilization) {
  adj_util_reoriented <- t(adj_utilization)
  incidence_matrix <- adj_production | adj_util_reoriented
  incidence_matrix <- as.matrix(incidence_matrix > 0) * 1

  row_sums <- rowSums(incidence_matrix)
  col_sums <- colSums(incidence_matrix)
  row_order <- order(-row_sums, seq_len(nrow(incidence_matrix)))
  col_order <- order(-col_sums, seq_len(ncol(incidence_matrix)))

  return(incidence_matrix[row_order, col_order, drop = FALSE])
}

# Matrix heatmap visual
plot_nestedness_matrix <- function(incidence_matrix, title, nodf_value = NULL) {
  df <- as.data.frame(as.table(incidence_matrix))
  colnames(df) <- c("Agent", "Pyov", "Interaction")

  p <- ggplot(df, aes(x = Pyov, y = Agent, fill = as.factor(Interaction))) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_manual(values = c("white", "black"), name = "Interaction") +
    labs(title = title, x = "Pyoverdine Communities", y = "Agents") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      axis.text.y = element_text(size = 6),
      legend.position = "bottom"
    )

  if (!is.null(nodf_value)) {
    p <- p + labs(title = paste(title, sprintf("(NODF = %.2f)", nodf_value)))
  }

  return(p)
}

# Primary bipartite NODF
compute_nodf_bipartite <- function(incidence_matrix) {
  nodf_total <- NA_real_
  nodf_rows <- NA
  nodf_cols <- NA

  res_primary <- tryCatch({
    suppressWarnings(bipartite::nested(incidence_matrix, method = "NODF"))
  }, error = function(e) NULL)

  if (!is.null(res_primary)) {
    if (is.numeric(res_primary) && length(res_primary) == 1) {
      nodf_total <- as.numeric(res_primary)
    } else if (is.list(res_primary)) {
      possible_fields <- c("statistic", "NODF", "value", "nestedness", "NODF_total")
      for (nm in possible_fields) {
        if (!is.null(res_primary[[nm]]) && is.numeric(res_primary[[nm]]) && length(res_primary[[nm]]) >= 1) {
          nodf_total <- as.numeric(res_primary[[nm]][1])
          break
        }
      }
      if (!is.null(res_primary$NODF_rows)) nodf_rows <- res_primary$NODF_rows
      if (!is.null(res_primary$NODF_cols)) nodf_cols <- res_primary$NODF_cols
    }
  }

  return(list(
    nodf_total = nodf_total,
    nodf_rows = nodf_rows,
    nodf_cols = nodf_cols,
    n_rows = nrow(incidence_matrix),
    n_cols = ncol(incidence_matrix),
    total_interactions = sum(incidence_matrix),
    matrix_fill = sum(incidence_matrix) / (nrow(incidence_matrix) * ncol(incidence_matrix)),
    method = "bipartite::nested (NODF)"
  ))
}

# Vegan cross-validation
compute_nodf_vegan <- function(incidence_matrix) {
  res <- tryCatch({
    suppressWarnings(vegan::nestednodf(incidence_matrix))
  }, error = function(e) NA_real_)

  nodf_val <- NA_real_
  if (is.numeric(res) && length(res) >= 1 && !is.na(res[1])) {
    nodf_val <- as.numeric(res[1])
  }

  list(
    nodf_total = nodf_val,
    method = "vegan::nestednodf",
    convergence = ifelse(!is.na(nodf_val), "success", "failed")
  )
}

# -----------------------------
# Directory setup
# -----------------------------
safe_dir_create("figures/network_topology")
safe_dir_create("docs/phase_04/logs")

cat("=== Phase 04 - Step 04: Nestedness Analysis (Minimal Output) ===\n")
cat(paste("Timestamp:", timestamp(), "\n\n"))

# -----------------------------
# Load and build matrices
# -----------------------------
cat("Loading adjacency matrices...\n")
adj_prod_fg <- readRDS("data/interim/adj_production_FG_agentsxpyov_conservative.rds")
adj_util_fg <- readRDS("data/interim/adj_utilization_FG_pyovxagents_conservative.rds")
adj_prod_str <- readRDS("data/interim/adj_production_STR_agentsxpyov_conservative.rds")
adj_util_str <- readRDS("data/interim/adj_utilization_STR_pyovxagents_conservative.rds")

fg_nodes <- readRDS("data/interim/nodes_functional_groups_conservative.rds")
str_nodes <- readRDS("data/interim/nodes_strains_conservative.rds")

cat("Building incidence matrices...\n")
incidence_fg <- build_incidence_matrix(adj_prod_fg, adj_util_fg)
incidence_str <- build_incidence_matrix(adj_prod_str, adj_util_str)

cat(sprintf("FG: %d × %d (%.1f%% fill)\n", nrow(incidence_fg), ncol(incidence_fg),
            100 * sum(incidence_fg) / (nrow(incidence_fg) * ncol(incidence_fg))))
cat(sprintf("STR: %d × %d (%.1f%% fill)\n", nrow(incidence_str), ncol(incidence_str),
            100 * sum(incidence_str) / (nrow(incidence_str) * ncol(incidence_str))))

# -----------------------------
# Compute NODF values
# -----------------------------
cat("Computing NODF values...\n")
nodf_fg_bipartite <- compute_nodf_bipartite(incidence_fg)
nodf_fg_vegan <- compute_nodf_vegan(incidence_fg)
nodf_str_bipartite <- compute_nodf_bipartite(incidence_str)
nodf_str_vegan <- compute_nodf_vegan(incidence_str)

# -----------------------------
# Generate visualizations
# -----------------------------
cat("Generating visualizations...\n")

# FG matrix
p_fg_matrix <- plot_nestedness_matrix(
  incidence_fg, "Functional Group Network", nodf_fg_bipartite$nodf_total
)
ggsave("figures/network_topology/step04_nodf_matrix_fg.pdf", p_fg_matrix, width = 10, height = 8)

# STR matrix (subset if large)
if (nrow(incidence_str) <= 200) {
  p_str_matrix <- plot_nestedness_matrix(
    incidence_str, "Strain Network", nodf_str_bipartite$nodf_total
  )
  ggsave("figures/network_topology/step04_nodf_matrix_str.pdf", p_str_matrix, width = 10, height = 12)
} else {
  top_strains <- head(order(rowSums(incidence_str), decreasing = TRUE), 200)
  p_str_matrix <- plot_nestedness_matrix(
    incidence_str[top_strains, ], "Strain Network (Top 200)", nodf_str_bipartite$nodf_total
  )
  ggsave("figures/network_topology/step04_nodf_matrix_str.pdf", p_str_matrix, width = 10, height = 12)
}

# Comparison plot
comparison_data <- data.frame(
  Network = c("Functional Groups", "Strains"),
  NODF = c(nodf_fg_bipartite$nodf_total, nodf_str_bipartite$nodf_total),
  Rows = c(nodf_fg_bipartite$n_rows, nodf_str_bipartite$n_rows),
  Columns = c(nodf_fg_bipartite$n_cols, nodf_str_bipartite$n_cols)
)

p_comparison <- ggplot(comparison_data, aes(x = Network, y = NODF, fill = Network)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.2f", NODF)), vjust = -0.5, size = 4) +
  labs(title = "Nestedness (NODF) Comparison", y = "NODF Score", x = "Network Level") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("figures/network_topology/step04_nodf_comparison.pdf", p_comparison, width = 8, height = 6)

# -----------------------------
# Generate summary report
# -----------------------------
summary_lines <- c(
  "Phase 04 - Step 04: Nestedness (NODF) Analysis Summary",
  paste("Timestamp:", timestamp()),
  "",
  "Method:",
  "  Primary: bipartite::nested() with deterministic tie-breaking",
  "  Cross-validation: vegan::nestednodf()",
  "",
  "Functional Group (FG) Level:",
  sprintf("  Network: %d FGs × %d Pyovs", nodf_fg_bipartite$n_rows, nodf_fg_bipartite$n_cols),
  sprintf("  Interactions: %d (%.1f%% fill)", nodf_fg_bipartite$total_interactions,
          100 * nodf_fg_bipartite$matrix_fill),
  sprintf("  NODF (bipartite): %.2f", nodf_fg_bipartite$nodf_total),
  sprintf("  NODF (vegan cross-check): %.2f", nodf_fg_vegan$nodf_total),
  "",
  "Strain (STR) Level:",
  sprintf("  Network: %d Strains × %d Pyovs", nodf_str_bipartite$n_rows, nodf_str_bipartite$n_cols),
  sprintf("  Interactions: %d (%.1f%% fill)", nodf_str_bipartite$total_interactions,
          100 * nodf_str_bipartite$matrix_fill),
  sprintf("  NODF (bipartite): %.2f", nodf_str_bipartite$nodf_total),
  sprintf("  NODF (vegan cross-check): %.2f", nodf_str_vegan$nodf_total),
  "",
  "Interpretation:",
  "  • NODF < 10: Low nestedness (modular structure)",
  "  • NODF 10-30: Moderate nestedness",
  "  • NODF > 30: High nestedness (hierarchical sharing)",
  "",
  "Outputs:",
  "  • figures/network_topology/step04_nodf_matrix_fg.pdf",
  "  • figures/network_topology/step04_nodf_matrix_str.pdf",
  "  • figures/network_topology/step04_nodf_comparison.pdf",
  "  • docs/phase_04/logs/step04_nestedness_summary.txt",
  "",
  "Note: No RDS/CSV files generated (compact output mode)",
  "",
  "Ready for subsequent analyses."
)

writeLines(summary_lines, "docs/phase_04/logs/step04_nestedness_summary.txt")

# Save session info
writeLines(capture.output(sessionInfo()), "docs/phase_04/logs/step04_session_info.txt")

# Display summary
cat(paste(summary_lines, collapse = "\n"), "\n")

cat("\nStep 04 complete. Nestedness analysis finished.\n")
cat("Summary: docs/phase_04/logs/step04_nestedness_summary.txt\n")
cat(sprintf("Figures: %s\n", paste(list.files("figures/network_topology", pattern = "^step04_.*\\.pdf$", full.names = TRUE), collapse = ", ")))
