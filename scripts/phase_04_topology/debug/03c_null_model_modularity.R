#!/usr/bin/env Rscript

# Phase 04 - Step 03c: Null-model modularity test (FG-level, debug/logs-only)
# -------------------------------------------------------------------------
# Purpose:
#   Test whether observed modularity is significantly higher than expected from
#   degree-preserving null models using an independent-generation strategy
#   (same approach as the NODF null analysis). This variant writes only a single
#   human-readable summary text file into the debug logs directory.
#
# Method:
#   - Matrix to randomize: FG × PYO incidence matrix (binary)
#   - Null generation: vegan::oecosimu(method = "curveball") with independent seeds
#   - Statistic: igraph Louvain modularity on combined bipartite graph
#   - Output: single text summary in scripts/phase_04_topology/debug/logs/
#
# Usage:
#   Rscript scripts/phase_04_topology/debug/03c_null_model_modularity.R

suppressPackageStartupMessages({
  library(Matrix)
  library(igraph)
  library(vegan)
  library(dplyr)    # for case_when
})

# -----------------------------
# Configuration (minimal)
# -----------------------------
N_NULLS <- 200
MASTER_SEED <- 2025
PROGRESS_UPDATE_FREQ <- 10

# Ensure debug logs directory exists
safe_dir_create <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
debug_logs_dir <- "scripts/phase_04_topology/debug/logs"
safe_dir_create(debug_logs_dir)

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

cat("=== Phase 04 - Step 03c: Null-model modularity test (FG-level) ===\n")
cat(paste("Timestamp:", timestamp(), "\n"))
cat(sprintf("Config: N_NULLS=%d, method=oecosimu(curveball), independent nulls\n", N_NULLS))

# -----------------------------
# Load data and build incidence matrix
# -----------------------------
adj_prod_fg_path <- "data/interim/adj_production_FG_agentsxpyov_conservative.rds"
adj_util_fg_path <- "data/interim/adj_utilization_FG_pyovxagents_conservative.rds"

stopifnot(file.exists(adj_prod_fg_path), file.exists(adj_util_fg_path))

cat("Loading adjacency matrices...\n")
adj_prod_fg <- readRDS(adj_prod_fg_path)     # FG x PYO
adj_util_fg <- readRDS(adj_util_fg_path)     # PYO x FG

# Build incidence matrix (production OR utilization)
build_incidence <- function(adj_production, adj_utilization) {
  inc <- (adj_production | t(adj_utilization))
  (as.matrix(inc > 0)) * 1
}

inc_fg <- build_incidence(adj_prod_fg, adj_util_fg)
n_edges <- sum(inc_fg)
fill_pct <- 100 * n_edges / (nrow(inc_fg) * ncol(inc_fg))

cat(sprintf("Incidence matrix: %d x %d, %d edges (%.1f%% fill)\n",
            nrow(inc_fg), ncol(inc_fg), n_edges, fill_pct))

# -----------------------------
# Helper: Compute Louvain modularity from incidence
# -----------------------------
modularity_from_incidence <- function(inc) {
  nr <- nrow(inc); nc <- ncol(inc)
  if (nr == 0 || nc == 0) return(NA_real_)

  which_ones <- which(inc > 0, arr.ind = TRUE)
  if (nrow(which_ones) == 0) return(NA_real_)

  fg_ids <- rownames(inc); pyo_ids <- colnames(inc)
  if (is.null(fg_ids)) fg_ids <- paste0("FG_", seq_len(nr))
  if (is.null(pyo_ids)) pyo_ids <- paste0("PYO_", seq_len(nc))

  edge_from <- paste0("FG:", fg_ids[which_ones[,1]])
  edge_to   <- paste0("PYO:", pyo_ids[which_ones[,2]])
  g <- graph_from_data_frame(data.frame(from = edge_from, to = edge_to), directed = FALSE)

  cl <- cluster_louvain(g)
  Q <- modularity(g, membership = membership(cl))
  as.numeric(Q)
}

# -----------------------------
# Observed modularity
# -----------------------------
cat("\nComputing observed modularity...\n")
obs_Q <- modularity_from_incidence(inc_fg)
cat(sprintf("Observed Louvain modularity: %.4f\n", obs_Q))

# -----------------------------
# Null model generation and testing
# -----------------------------
cat("\nGenerating independent nulls (curveball method)...\n")

modularity_stat_fun <- function(x) {
  x_bin <- (as.matrix(x) > 0) * 1
  modularity_from_incidence(x_bin)
}

null_seeds <- MASTER_SEED + (1:N_NULLS) * 1000
null_Q <- rep(NA_real_, N_NULLS)

progress_start <- Sys.time()

for (i in seq_len(N_NULLS)) {
  set.seed(null_seeds[i])

  res <- tryCatch({
    vegan::oecosimu(inc_fg,
                    nestfun = modularity_stat_fun,
                    method = "curveball",
                    nsimul = 2,
                    alternative = "two.sided")
  }, error = function(e) {
    NULL
  })

  if (!is.null(res) && !is.null(res$oecosimu$simulated)) {
    null_Q[i] <- as.numeric(res$oecosimu$simulated[1])
  }

  if (i %% PROGRESS_UPDATE_FREQ == 0 || i == N_NULLS) {
    elapsed <- as.numeric(difftime(Sys.time(), progress_start, units = "secs"))
    rate <- if (elapsed > 0) i / elapsed else NA
    remaining_rate <- if (!is.na(rate) && rate > 0) (N_NULLS - i) / rate else NA
    cat(sprintf("\rProgress: %d/%d (%.1f%%) | Rate: %.2f nulls/sec | ETA: %.0f sec",
                i, N_NULLS, 100 * i / N_NULLS, ifelse(is.na(rate), 0, rate), ifelse(is.na(remaining_rate), 0, remaining_rate)))
    flush.console()
  }
}
cat("\n") # newline after progress

# -----------------------------
# Analyze null results (no CSV or RDS output; only text log)
# -----------------------------
null_clean <- null_Q[!is.na(null_Q)]
n_valid <- length(null_clean)

if (n_valid == 0) {
  warning("No valid nulls were generated - aborting statistical comparison.")
}

null_mean <- if (n_valid > 0) mean(null_clean) else NA_real_
null_sd   <- if (n_valid > 1) sd(null_clean) else NA_real_
null_median <- if (n_valid > 0) median(null_clean) else NA_real_

# empirical p-value (two-tailed) and SES
p_value <- if (n_valid > 0) (sum(abs(null_clean - null_mean) >= abs(obs_Q - null_mean)) + 1) / (n_valid + 1) else NA_real_
ses <- if (!is.na(null_sd) && null_sd > 0) (obs_Q - null_mean) / null_sd else NA_real_
z_score <- ses

significance <- if (!is.na(p_value) && p_value < 0.05) "SIGNIFICANT" else "Not significant"
effect_interpretation <- dplyr::case_when(
  !is.na(ses) & ses > 2 ~ "Strong positive deviation (more modular than expected)",
  !is.na(ses) & ses < -2 ~ "Strong negative deviation (less modular than expected)",
  TRUE ~ "Within expected range given degree constraints"
)

cat("\n--- RESULTS ---\n")
cat(sprintf("Observed modularity:     %.4f\n", obs_Q))
cat(sprintf("Null mean ± SD:         %.4f ± %.4f\n", null_mean, null_sd))
cat(sprintf("Null median:            %.4f\n", null_median))
cat(sprintf("P-value (two-tailed):   %.4f\n", p_value))
cat(sprintf("Standardized Effect Size (SES): %.2f\n", ses))
cat(sprintf("Significance:           %s\n", significance))
cat(sprintf("Valid nulls:            %d/%d (%.1f%%)\n", n_valid, N_NULLS, 100 * n_valid / N_NULLS))
cat(sprintf("Effect interpretation:  %s\n", effect_interpretation))

# -----------------------------
# Write single text summary to debug logs directory
# -----------------------------
summary_lines <- c(
  "Phase 04 - Step 03c: Null-model modularity test (FG-level)",
  paste("Timestamp:", timestamp()),
  sprintf("Method: vegan::oecosimu(curveball) + igraph Louvain"),
  sprintf("Parameters: N_NULLS=%d, independent seeds, nsimul= workaround", N_NULLS),
  "",
  "Network Summary:",
  sprintf("  Incidence matrix: %d x %d (%d edges, %.1f%% fill)",
          nrow(inc_fg), ncol(inc_fg), n_edges, fill_pct),
  "",
  "Modularity Results:",
  sprintf("  Observed modularity:     %.4f", obs_Q),
  sprintf("  Null mean ± SD:         %.4f ± %.4f", null_mean, null_sd),
  sprintf("  Null range:            %s", if (n_valid>0) paste0(sprintf("%.4f - %.4f", min(null_clean), max(null_clean))) else "NA"),
  sprintf("  P-value (two-tailed):   %.4f", p_value),
  sprintf("  SES (Standardized Effect Size): %.2f", ses),
  sprintf("  Significance:           %s", significance),
  sprintf("  Effect interpretation:  %s", effect_interpretation),
  "",
  "Quality Metrics:",
  sprintf("  Valid nulls:            %d/%d (%.1f%%)", n_valid, N_NULLS, 100 * n_valid / N_NULLS),
  sprintf("  Independence method:    Unique seeds per null"),
  "",
  "Interpretation Guide:",
  "  • SES > 2: Network is more modular than degree constraints predict",
  "  • |SES| < 2: Modularity within expected range of null models",
  "  • p < 0.05: Statistically significant deviation from null expectation",
  "",
  "Next steps:",
  if (!is.na(p_value) && p_value < 0.05 && !is.na(ses) && ses > 2) {
    "  • Observed modularity is significant → Invest in robust community detection"
  } else {
    "  • Observed modularity not significant → Focus on other topological metrics"
  },
  "",
  "Note: This script intentionally writes a single text summary file only (no CSV or RDS).",
  "The full null distribution is not saved by design; re-run with modifications if you need the raw vector."
)

summary_path <- file.path(debug_logs_dir, "step03c_null_modularity_summary.txt")
writeLines(summary_lines, con = summary_path)

cat(sprintf("\nSummary written to: %s\n", summary_path))
cat("\n=== Null-model modularity test complete ===\n")
