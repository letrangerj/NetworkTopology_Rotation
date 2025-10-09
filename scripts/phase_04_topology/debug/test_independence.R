#!/usr/bin/env Rscript

# Test Independence of Null Models
# ---------------------------------
# Purpose: Test whether vegan::oecosimu generates truly independent null models
# or if there's correlation between consecutive nulls that could bias results.

suppressPackageStartupMessages({
  library(Matrix)
  library(vegan)
  library(bipartite)
})

cat("=== Testing Null Model Independence ===\n")

# Load the same data
adj_prod_fg <- readRDS("data/interim/adj_production_FG_agentsxpyov_conservative.rds")
adj_util_fg <- readRDS("data/interim/adj_utilization_FG_pyovxagents_conservative.rds")

# Build incidence matrix (same as main script)
build_incidence_matrix <- function(adj_production, adj_utilization) {
  adj_util_reoriented <- t(adj_utilization)
  incidence_matrix <- adj_production | adj_util_reoriented
  incidence_matrix <- as.matrix(incidence_matrix > 0) * 1
  row_sums <- rowSums(incidence_matrix)
  col_sums <- colSums(incidence_matrix)
  row_order <- order(row_sums, decreasing = TRUE, rownames(incidence_matrix))
  col_order <- order(col_sums, decreasing = TRUE, colnames(incidence_matrix))
  return(incidence_matrix[row_order, col_order, drop = FALSE])
}

# NODF computation function
compute_nodf_bipartite <- function(incidence_matrix) {
  nodf_total <- NA_real_
  res_primary <- tryCatch({
    bipartite::nested(incidence_matrix, method = "NODF")
  }, error = function(e) NULL)

  if (!is.null(res_primary)) {
    if (is.numeric(res_primary) && length(res_primary) >= 1) {
      nodf_total <- as.numeric(res_primary[1])
    } else if (is.list(res_primary)) {
      possible_fields <- c("statistic", "NODF", "value", "nestedness", "NODF_total")
      for (nm in possible_fields) {
        if (!is.null(res_primary[[nm]]) && is.numeric(res_primary[[nm]]) && length(res_primary[[nm]]) >= 1) {
          nodf_total <- as.numeric(res_primary[[nm]][1])
          break
        }
      }
    }
  }

  if (is.na(nodf_total)) {
    res_vegan <- tryCatch({
      vegan::nestednodf(incidence_matrix)
    }, error = function(e) NA_real_)
    if (is.numeric(res_vegan) && length(res_vegan) >= 1 && !is.na(res_vegan[1])) {
      nodf_total <- as.numeric(res_vegan[1])
    }
  }

  return(nodf_total)
}

nodf_stat_function <- function(x) {
  compute_nodf_bipartite(as.matrix(x))
}

# Build incidence matrix
incidence_fg <- build_incidence_matrix(adj_prod_fg, adj_util_fg)
cat(sprintf("Matrix: %d x %d, %d edges\n", nrow(incidence_fg), ncol(incidence_fg), sum(incidence_fg)))

# Test 1: Multiple runs with same N should give similar means
cat("\n=== Test 1: Reproducibility with Same N ===\n")
test_runs <- 5
test_n <- 20
means_same_n <- numeric(test_runs)

for (i in 1:test_runs) {
  set.seed(2025)  # Same seed for each run
  result <- vegan::oecosimu(incidence_fg,
                           nestfun = nodf_stat_function,
                           method = "curveball",
                           nsimul = test_n)
  nulls <- as.numeric(result$oecosimu$simulated)
  means_same_n[i] <- mean(nulls, na.rm = TRUE)
  cat(sprintf("Run %d: mean = %.4f, range = %.4f - %.4f\n",
              i, means_same_n[i], min(nulls, na.rm=TRUE), max(nulls, na.rm=TRUE)))
}

cat(sprintf("Mean of means: %.4f, SD of means: %.6f\n", mean(means_same_n), sd(means_same_n)))
if (sd(means_same_n) > 0.01) {
  cat("WARNING: High variability between runs with same seed - reproducibility issue!\n")
} else {
  cat("GOOD: Low variability between runs - reproducible\n")
}

# Test 2: Different N values to see if mean increases
cat("\n=== Test 2: Mean vs N (Independence Test) ===\n")
n_values <- c(10, 20, 30, 50, 100)
means_vs_n <- numeric(length(n_values))

for (i in seq_along(n_values)) {
  n <- n_values[i]
  set.seed(2025)  # Same starting seed
  result <- vegan::oecosimu(incidence_fg,
                           nestfun = nodf_stat_function,
                           method = "curveball",
                           nsimul = n)
  nulls <- as.numeric(result$oecosimu$simulated)
  means_vs_n[i] <- mean(nulls, na.rm = TRUE)
  cat(sprintf("N=%3d: mean = %.4f, SD = %.4f, min = %.4f, max = %.4f\n",
              n, means_vs_n[i], sd(nulls, na.rm=TRUE), min(nulls, na.rm=TRUE), max(nulls, na.rm=TRUE)))
}

# Check for trend
correlation <- cor(n_values, means_vs_n)
cat(sprintf("\nCorrelation between N and mean: %.4f\n", correlation))

if (abs(correlation) > 0.5) {
  cat("WARNING: Strong correlation between N and mean - NON-INDEPENDENT nulls!\n")
} else {
  cat("GOOD: Weak correlation - nulls appear independent\n")
}

# Test 3: Sequential vs batch generation
cat("\n=== Test 3: Sequential vs Batch Generation ===\n")

# Batch generation (current method)
set.seed(2025)
batch_result <- vegan::oecosimu(incidence_fg,
                               nestfun = nodf_stat_function,
                               method = "curveball",
                               nsimul = 20)
batch_nulls <- as.numeric(batch_result$oecosimu$simulated)
batch_mean <- mean(batch_nulls, na.rm = TRUE)

# Sequential generation (one at a time)
set.seed(2025)
sequential_nulls <- numeric(20)
for (i in 1:20) {
  result <- vegan::oecosimu(incidence_fg,
                           nestfun = nodf_stat_function,
                           method = "curveball",
                           nsimul = 1)
  sequential_nulls[i] <- as.numeric(result$oecosimu$simulated)
}
sequential_mean <- mean(sequential_nulls, na.rm = TRUE)

cat(sprintf("Batch generation mean: %.4f\n", batch_mean))
cat(sprintf("Sequential generation mean: %.4f\n", sequential_mean))
cat(sprintf("Difference: %.6f\n", abs(batch_mean - sequential_mean)))

if (abs(batch_mean - sequential_mean) > 0.01) {
  cat("WARNING: Batch vs sequential gives different results - potential dependence!\n")
} else {
  cat("GOOD: Batch and sequential generation agree\n")
}

# Test 4: Autocorrelation test
cat("\n=== Test 4: Autocorrelation Analysis ===\n")
set.seed(2025)
autocorr_result <- vegan::oecosimu(incidence_fg,
                                  nestfun = nodf_stat_function,
                                  method = "curveball",
                                  nsimul = 100)
autocorr_nulls <- as.numeric(autocorr_result$oecosimu$simulated)

# Test for autocorrelation
if (length(autocorr_nulls) > 10) {
  lag1_corr <- cor(autocorr_nulls[-length(autocorr_nulls)], autocorr_nulls[-1], use = "complete.obs")
  cat(sprintf("Lag-1 autocorrelation: %.4f\n", lag1_corr))

  if (abs(lag1_corr) > 0.1) {
    cat("WARNING: Significant autocorrelation detected - nulls are not independent!\n")
  } else {
    cat("GOOD: Low autocorrelation - nulls appear independent\n")
  }
} else {
  cat("Not enough valid nulls for autocorrelation test\n")
}

# Summary and recommendations
cat("\n=== SUMMARY ===\n")
issues <- 0

if (sd(means_same_n) > 0.01) {
  cat("❌ ISSUE: Poor reproducibility\n")
  issues <- issues + 1
}

if (abs(correlation) > 0.5) {
  cat("❌ ISSUE: Mean increases with N (non-independence)\n")
  issues <- issues + 1
}

if (abs(batch_mean - sequential_mean) > 0.01) {
  cat("❌ ISSUE: Batch vs sequential difference\n")
  issues <- issues + 1
}

if (length(autocorr_nulls) > 10 && abs(lag1_corr) > 0.1) {
  cat("❌ ISSUE: Significant autocorrelation\n")
  issues <- issues + 1
}

if (issues == 0) {
  cat("✅ ALL TESTS PASSED: Null models appear independent and valid\n")
} else {
  cat(sprintf("⚠️  %d ISSUES DETECTED: Null model independence questionable\n", issues))
  cat("\nRECOMMENDATIONS:\n")
  cat("1. Use individual oecosimu calls with separate seeds\n")
  cat("2. Consider alternative null model methods\n")
  cat("3. Test with different random number generators\n")
  cat("4. Validate results against known benchmarks\n")
}

cat("\nTest complete. Check results above for independence validation.\n")
