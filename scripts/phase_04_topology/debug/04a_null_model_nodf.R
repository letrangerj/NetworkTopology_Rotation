#!/usr/bin/env Rscript

# Phase 04 - Step 04a: Fast Null Model NODF Analysis (FG Only)
# -------------------------------------------------------------
# Purpose:
#   Perform null model analysis of NODF nestedness using vegan::oecosimu() for
#   fast null generation while computing NODF with bipartite::nested().
#
# Modified:
#   - Removed RDS and PDF file saving (compact output)
#   - Reduced terminal verbosity (progress bar + final summary only)
#   - Deterministic tie-breaking for reproducibility
#   - suppressWarnings to preserve bipartite::nested results
#
# Outputs:
#   - debug/logs/04_a_null_summary.txt (only file saved)
#
# Usage:
#   Rscript scripts/phase_04_topology/debug/04a_null_model_nodf.R

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(vegan)
  library(bipartite)
  if (!require(progress, quietly = TRUE)) {
    cat("Installing progress package for progress bar...\n")
    install.packages("progress", repos = "https://cran.r-project.org/")
    library(progress)
  }
})

# -----------------------------
# Configuration (optimized for speed)
# -----------------------------
N_NULLS <- 100          # Reduced for fast debug execution
NULL_METHOD <- "curveball"  # Fast degree-preserving method
MASTER_SEED <- 2025
PROGRESS_UPDATE_FREQ <- 5

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

# Create debug output directory (logs only)
safe_dir_create("scripts/phase_04_topology/debug/logs")

cat("=== Phase 04 - Step 04a: Fast Null Model NODF Analysis (FG Only) ===\n")
cat(paste("Timestamp:", timestamp(), "\n"))
cat(paste("Configuration: N_nulls =", N_NULLS, ", method = vegan::oecosimu +", NULL_METHOD, "\n"))
cat("Compact mode: No RDS/PDF outputs, progress bar + summary only\n\n")

# -----------------------------
# Load observed data and results
# -----------------------------
cat("Loading observed NODF results and adjacency matrices...\n")

# Check that Step 04 has been run
if (!file.exists("results/phase_04/nestedness/nestedness_fg.rds")) {
  stop("Step 04 nestedness analysis not found. Run 04_nestedness_nodf.R first.")
}

# Load observed NODF results
observed_fg <- readRDS("results/phase_04/nestedness/nestedness_fg.rds")

# Load adjacency matrices
adj_prod_fg <- readRDS("data/interim/adj_production_FG_agentsxpyov_conservative.rds")
adj_util_fg <- readRDS("data/interim/adj_utilization_FG_pyovxagents_conservative.rds")

# Extract observed NODF value (use primary bipartite method for consistency)
observed_nodf_fg <- observed_fg$primary$nodf_total

if (is.na(observed_nodf_fg)) {
  stop("Observed NODF value is NA. Check Step 04 output.")
}

cat(sprintf("Observed NODF - FG: %.3f\n", observed_nodf_fg))

# -----------------------------
# Helper functions
# -----------------------------

# Build incidence matrix (deterministic tie-breaking)
build_incidence_matrix <- function(adj_production, adj_utilization) {
  adj_util_reoriented <- t(adj_utilization)
  incidence_matrix <- adj_production | adj_util_reoriented
  incidence_matrix <- as.matrix(incidence_matrix > 0) * 1

  # Order by decreasing degree for NODF calculation
  row_sums <- rowSums(incidence_matrix)
  col_sums <- colSums(incidence_matrix)
  # Deterministic tie-breaking: primary key = decreasing degree, secondary = original index
  row_order <- order(-row_sums, seq_len(nrow(incidence_matrix)))
  col_order <- order(-col_sums, seq_len(ncol(incidence_matrix)))

  return(incidence_matrix[row_order, col_order, drop = FALSE])
}

# Compute NODF using bipartite method (suppressWarnings to preserve results)
compute_nodf_bipartite <- function(incidence_matrix) {
  nodf_total <- NA_real_
  nodf_rows <- NA
  nodf_cols <- NA

  # Try primary method: bipartite::nested with NODF request
  res_primary <- tryCatch({
    suppressWarnings(bipartite::nested(incidence_matrix, method = "NODF"))
  }, error = function(e) {
    NULL
  })

  if (!is.null(res_primary)) {
    # If the result is a single numeric, treat as total NODF
    if (is.numeric(res_primary) && length(res_primary) == 1) {
      nodf_total <- as.numeric(res_primary)
    } else if (is.numeric(res_primary) && length(res_primary) > 1) {
      # Sometimes returns a vector, take first element
      nodf_total <- as.numeric(res_primary[1])
    } else if (is.list(res_primary)) {
      # Common possible field names - check defensively
      possible_fields <- c("statistic", "NODF", "value", "nestedness", "NODF_total")
      for (nm in possible_fields) {
        if (!is.null(res_primary[[nm]]) && is.numeric(res_primary[[nm]]) && length(res_primary[[nm]]) >= 1) {
          nodf_total <- as.numeric(res_primary[[nm]][1])
          break
        }
      }
      # Try to extract row/column components if present
      if (!is.null(res_primary$NODF_rows)) nodf_rows <- res_primary$NODF_rows
      if (!is.null(res_primary$NODF_cols)) nodf_cols <- res_primary$NODF_cols
      if (!is.null(res_primary$rows)) nodf_rows <- res_primary$rows
      if (!is.null(res_primary$cols)) nodf_cols <- res_primary$cols
    }
  }

  # Fallback: vegan::nestednodf if primary failed
  if (is.na(nodf_total)) {
    res_vegan <- tryCatch({
      vegan::nestednodf(incidence_matrix)
    }, error = function(e) {
      NA_real_
    })
    if (is.numeric(res_vegan) && length(res_vegan) >= 1 && !is.na(res_vegan[1])) {
      nodf_total <- as.numeric(res_vegan[1])
    }
  }

  return(nodf_total)
}

# Wrapper function for oecosimu that computes NODF using bipartite
nodf_stat_function <- function(x) {
  # Minimal wrapper - no verbose output
  result <- tryCatch({
    compute_nodf_bipartite(as.matrix(x))
  }, error = function(e) {
    NA_real_
  }, warning = function(w) {
    compute_nodf_bipartite(as.matrix(x))
  })

  return(result)
}

# -----------------------------
# Build incidence matrix and test NODF
# -----------------------------
cat("Building FG incidence matrix...\n")

# Build FG incidence matrix
incidence_fg <- build_incidence_matrix(adj_prod_fg, adj_util_fg)
n_edges_fg <- sum(incidence_fg)

cat(sprintf("FG network: %d x %d matrix, %d edges (%.1f%% fill)\n",
            nrow(incidence_fg), ncol(incidence_fg), n_edges_fg,
            100 * n_edges_fg / (nrow(incidence_fg) * ncol(incidence_fg))))

# Test NODF computation on observed matrix
cat("Testing NODF computation...\n")
test_nodf <- compute_nodf_bipartite(incidence_fg)
cat(sprintf("Test NODF result: %.3f\n", test_nodf))

if (is.na(test_nodf)) {
  stop("NODF computation failed on observed matrix. Check bipartite installation.")
}

# Verify consistency with loaded observed value
if (abs(test_nodf - observed_nodf_fg) > 0.01) {
  cat(sprintf("WARNING: Test NODF (%.3f) differs from loaded observed (%.3f)\n",
              test_nodf, observed_nodf_fg))
  cat("Using test NODF value for consistency.\n")
  observed_nodf_fg <- test_nodf
}

# -----------------------------
# Fast null model analysis using vegan::oecosimu
# -----------------------------
cat("Running null model analysis (progress bar only)...\n")

# Set up progress tracking
pb <- progress_bar$new(
  format = "Nulls [:bar] :percent (:current/:total) ETA: :eta",
  total = N_NULLS,
  clear = FALSE,
  width = 70
)

# Generate unique seeds for each null
cat("Generating independent nulls to avoid serial correlation...\n")
start_time <- Sys.time()
null_seeds <- MASTER_SEED + (1:N_NULLS) * 1000

# Initialize storage
null_nodf_fg <- numeric(N_NULLS)
valid_nulls <- 0

# Generate each null independently to avoid serial correlation
for (i in 1:N_NULLS) {
  set.seed(null_seeds[i])

  # Generate single null with nsimul=2 (nsimul=1 has bug) and take first value
  single_result <- tryCatch({
    vegan::oecosimu(incidence_fg,
                    nestfun = nodf_stat_function,
                    method = NULL_METHOD,
                    nsimul = 2,  # Use 2 to avoid nsimul=1 bug
                    alternative = "two.sided")
  }, error = function(e) {
    NULL
  })

  if (!is.null(single_result) && !is.null(single_result$oecosimu$simulated)) {
    null_nodf_fg[i] <- as.numeric(single_result$oecosimu$simulated[1])  # Take first value
    valid_nulls <- valid_nulls + 1
  } else {
    null_nodf_fg[i] <- NA_real_
  }

  # Update progress
  pb$tick()
}

end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
pb$terminate()

# Extract results
observed_nodf_oecosimu <- observed_nodf_fg  # Use the consistent observed value

# Remove NA values for analysis
null_nodf_fg_clean <- null_nodf_fg[!is.na(null_nodf_fg)]

if (length(null_nodf_fg_clean) < 5) {
  stop("Too few valid null NODF values for analysis.")
}

# -----------------------------
# Calculate statistics
# -----------------------------
null_mean_fg <- mean(null_nodf_fg_clean)
null_sd_fg <- sd(null_nodf_fg_clean)
null_median_fg <- median(null_nodf_fg_clean)

# P-value (two-tailed) using original observed value
p_value_fg <- (sum(abs(null_nodf_fg_clean - null_mean_fg) >= abs(observed_nodf_fg - null_mean_fg)) + 1) / (length(null_nodf_fg_clean) + 1)

# Z-score and SES using original observed value
z_score_fg <- (observed_nodf_fg - null_mean_fg) / null_sd_fg
ses_fg <- z_score_fg

cat(sprintf("\nFG Results: Observed=%.3f, Null_mean=%.3f±%.3f, p=%.4f, SES=%.3f\n",
            observed_nodf_fg, null_mean_fg, null_sd_fg, p_value_fg, ses_fg))

# -----------------------------
# Generate summary report (text only)
# -----------------------------
summary_lines <- c(
  "Phase 04 - Step 04a: Independent Null Model NODF Analysis Summary (FG Only)",
  paste("Timestamp:", timestamp()),
  sprintf("Method: Independent vegan::oecosimu(%s) + bipartite::nested", NULL_METHOD),
  sprintf("Runtime: %.1f seconds (%.1f nulls/sec)", runtime, N_NULLS / runtime),
  sprintf("Independence: Individual generation with unique seeds"),
  sprintf("Algorithm: nsimul=2 per null, take first value (avoids nsimul=1 bug)"),
  "",
  "Functional Group (FG) Level:",
  sprintf("  Observed NODF: %.3f", observed_nodf_fg),
  sprintf("  Valid nulls: %d out of %d (%.1f%% success rate)",
          length(null_nodf_fg_clean), N_NULLS, 100 * length(null_nodf_fg_clean) / N_NULLS),
  sprintf("  Null mean ± SD: %.3f ± %.3f", null_mean_fg, null_sd_fg),
  sprintf("  Null median: %.3f", null_median_fg),
  sprintf("  Null range: %.3f - %.3f", min(null_nodf_fg_clean), max(null_nodf_fg_clean)),
  "",
  "Statistical Tests:",
  sprintf("  P-value (two-tailed): %.4f", p_value_fg),
  sprintf("  Z-score: %.3f", z_score_fg),
  sprintf("  Standardized Effect Size (SES): %.3f", ses_fg),
  sprintf("  Significance: %s", ifelse(p_value_fg < 0.05, "SIGNIFICANT", "Not significant")),
  "",
  "Interpretation:",
  "  SES > 2: Strong positive deviation (more nested than expected)",
  "  SES < -2: Strong negative deviation (less nested than expected)",
  "  |SES| < 2: Within expected range given degree constraints",
  "",
  "Independence Validation:",
  if (length(null_nodf_fg_clean) >= 20) {
    first_half <- null_nodf_fg_clean[1:(length(null_nodf_fg_clean)%/%2)]
    second_half <- null_nodf_fg_clean[(length(null_nodf_fg_clean)%/%2 + 1):length(null_nodf_fg_clean)]
    mean_diff <- abs(mean(first_half) - mean(second_half))
    if (mean_diff < 0.1) {
      sprintf("  ✅ Independence check: PASSED (stable means, diff=%.3f)", mean_diff)
    } else {
      sprintf("  ⚠️  Independence check: WARNING (unstable means, diff=%.3f)", mean_diff)
    }
  } else {
    "  ⚠️  Independence check: Not enough nulls for validation"
  },
  "",
  "Note: Compact mode - only text summary saved (no RDS/PDF outputs)",
  "Re-run with modifications if you need raw values or visualizations."
)

writeLines(summary_lines, "scripts/phase_04_topology/debug/logs/04_a_null_summary.txt")

# Display summary
cat("\n", paste(summary_lines, collapse = "\n"), "\n")

cat("\n=== Independent FG-Only Null Model Analysis Complete ===\n")
cat(sprintf("Runtime: %.1f seconds, Success rate: %.1f%%\n", runtime, 100 * length(null_nodf_fg_clean) / N_NULLS))
cat("Summary saved to: scripts/phase_04_topology/debug/logs/04_a_null_summary.txt\n")
cat("Compact mode: No RDS/PDF files generated by design.\n")
