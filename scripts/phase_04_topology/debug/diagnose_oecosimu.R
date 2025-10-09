#!/usr/bin/env Rscript

# Diagnose oecosimu Return Structure
# ----------------------------------
# Purpose: Simple diagnostic to understand what vegan::oecosimu returns
# for different nsimul values, especially nsimul=1

suppressPackageStartupMessages({
  library(vegan)
  library(bipartite)
})

cat("=== Diagnosing vegan::oecosimu Return Structure ===\n")

# Load data
adj_prod_fg <- readRDS("data/interim/adj_production_FG_agentsxpyov_conservative.rds")
adj_util_fg <- readRDS("data/interim/adj_utilization_FG_pyovxagents_conservative.rds")

# Build matrix
adj_util_reoriented <- t(as.matrix(adj_util_fg))
incidence_fg <- as.matrix(adj_prod_fg) | adj_util_reoriented
incidence_fg <- as.matrix(incidence_fg > 0) * 1
row_sums <- rowSums(incidence_fg)
col_sums <- colSums(incidence_fg)
row_order <- order(row_sums, decreasing = TRUE)
col_order <- order(col_sums, decreasing = TRUE)
incidence_fg <- incidence_fg[row_order, col_order, drop = FALSE]

cat(sprintf("Matrix: %d x %d, %d edges\n", nrow(incidence_fg), ncol(incidence_fg), sum(incidence_fg)))

# Simple NODF function
simple_nodf <- function(x) {
  result <- tryCatch({
    res <- bipartite::nested(as.matrix(x), method = "NODF")
    if (is.numeric(res) && length(res) >= 1) {
      return(as.numeric(res[1]))
    }
    return(NA_real_)
  }, error = function(e) {
    return(NA_real_)
  })
}

# Test with nsimul = 1
cat("\n=== Testing nsimul = 1 ===\n")
set.seed(2025)
result1 <- tryCatch({
  vegan::oecosimu(incidence_fg,
                  nestfun = simple_nodf,
                  method = "curveball",
                  nsimul = 1)
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  return(NULL)
})

if (!is.null(result1)) {
  cat("Result1 class:", class(result1), "\n")
  cat("Result1 names:", paste(names(result1), collapse = ", "), "\n")
  cat("Result1 structure:\n")
  str(result1, max.level = 2)
} else {
  cat("Result1 is NULL\n")
}

# Test with nsimul = 3
cat("\n=== Testing nsimul = 3 ===\n")
set.seed(2025)
result3 <- tryCatch({
  vegan::oecosimu(incidence_fg,
                  nestfun = simple_nodf,
                  method = "curveball",
                  nsimul = 3)
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  return(NULL)
})

if (!is.null(result3)) {
  cat("Result3 class:", class(result3), "\n")
  cat("Result3 names:", paste(names(result3), collapse = ", "), "\n")
  cat("Result3 structure:\n")
  str(result3, max.level = 2)
} else {
  cat("Result3 is NULL\n")
}

# Compare extraction approaches
cat("\n=== Extraction Comparison ===\n")

if (!is.null(result1)) {
  cat("nsimul=1 extraction attempts:\n")

  # Method 1: Direct statistic
  if (!is.null(result1$statistic)) {
    cat("  statistic:", result1$statistic, "\n")
  }

  # Method 2: oecosimu component
  if (!is.null(result1$oecosimu)) {
    cat("  oecosimu present, class:", class(result1$oecosimu), "\n")
    if (is.list(result1$oecosimu)) {
      cat("  oecosimu names:", paste(names(result1$oecosimu), collapse = ", "), "\n")
      if (!is.null(result1$oecosimu$simulated)) {
        cat("  oecosimu$simulated:", result1$oecosimu$simulated, "\n")
      }
    }
  }
}

if (!is.null(result3)) {
  cat("\nnsimul=3 extraction attempts:\n")

  # Method 1: Direct statistic
  if (!is.null(result3$statistic)) {
    cat("  statistic:", result3$statistic, "\n")
  }

  # Method 2: oecosimu component
  if (!is.null(result3$oecosimu)) {
    cat("  oecosimu present, class:", class(result3$oecosimu), "\n")
    if (is.list(result3$oecosimu)) {
      cat("  oecosimu names:", paste(names(result3$oecosimu), collapse = ", "), "\n")
      if (!is.null(result3$oecosimu$simulated)) {
        cat("  oecosimu$simulated length:", length(result3$oecosimu$simulated), "\n")
        cat("  oecosimu$simulated values:", paste(result3$oecosimu$simulated, collapse = ", "), "\n")
      }
    }
  }
}

# Test direct curveball without oecosimu wrapper
cat("\n=== Testing Direct Curveball ===\n")
set.seed(2025)
direct_result <- tryCatch({
  # Try using vegan's curveball directly if available
  null_matrix <- vegan::curveball(incidence_fg)
  nodf_direct <- simple_nodf(null_matrix)
  cat("Direct curveball NODF:", nodf_direct, "\n")
  nodf_direct
}, error = function(e) {
  cat("Direct curveball failed:", e$message, "\n")
  NA_real_
})

cat("\n=== Summary ===\n")
cat("This diagnostic should help identify:\n")
cat("1. What oecosimu returns for nsimul=1 vs nsimul>1\n")
cat("2. How to extract simulated values correctly\n")
cat("3. Whether direct curveball is available as alternative\n")

cat("\nDiagnostic complete.\n")
