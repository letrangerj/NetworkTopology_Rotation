#!/usr/bin/env Rscript
#
# Debug: Inspect functional group strategy distributions and sample rows
# Writes human-readable summary and CSV samples to docs/phase_03/
#
# Usage:
#   Rscript scripts/phase_03_network/debug/inspect_strategies.R
#
# Outputs:
#   - docs/phase_03/debug_strategy_inspect.txt
#   - docs/phase_03/producer_only_samples.csv
#   - docs/phase_03/multi_receptor_samples.csv
#   - docs/phase_03/strategy_counts.csv
#
# The script is defensive: if input file missing or columns unexpected it reports and exits.

# ---- Configuration ----
fg_csv <- "data/interim/nodes_functional_groups_conservative.csv"
out_dir <- "docs/phase_03"
out_txt <- file.path(out_dir, "debug_strategy_inspect.txt")
out_counts_csv <- file.path(out_dir, "strategy_counts.csv")
out_prodonly_csv <- file.path(out_dir, "producer_only_samples.csv")
out_multirec_csv <- file.path(out_dir, "multi_receptor_samples.csv")

# ---- Helpers ----
mkdirp <- function(d) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Required input file not found: ", path)
  }
  # Read with base R for portability
  df <- tryCatch({
    read.csv(path, stringsAsFactors = FALSE, check.names = TRUE)
  }, error = function(e) {
    stop("Failed to read CSV: ", path, "\n  ", conditionMessage(e))
  })
  df
}

# ---- Run ----
mkdirp(out_dir)

msg_lines <- character()

msg_lines <- c(msg_lines, paste0("Strategy inspection run: ", Sys.time()))
msg_lines <- c(msg_lines, paste0("Input FG CSV: ", fg_csv))

# Load file
fg <- tryCatch({
  safe_read_csv(fg_csv)
}, error = function(e) {
  msg_lines <<- c(msg_lines, paste0("ERROR: ", conditionMessage(e)))
  cat(paste(msg_lines, collapse = "\n"), "\n")
  quit(status = 1)
})

# Check for expected columns
expected_cols <- c("agent_id", "functional_group_id", "strain_count", "strategy")
missing_cols <- setdiff(expected_cols, colnames(fg))
if (length(missing_cols) > 0) {
  msg_lines <- c(msg_lines, paste0("ERROR: Missing expected columns in FG table: ", paste(missing_cols, collapse = ", ")))
  msg_lines <- c(msg_lines, "Available columns: ", paste(colnames(fg), collapse = ", "))
  writeLines(msg_lines, out_txt)
  stop(paste(msg_lines, collapse = "\n"))
}

# Normalize strategy column (trim whitespace)
fg$strategy <- trimws(as.character(fg$strategy))

# Count strategies
strategy_counts <- as.data.frame(table(fg$strategy), stringsAsFactors = FALSE)
colnames(strategy_counts) <- c("strategy", "count")
strategy_counts <- strategy_counts[order(-strategy_counts$count), ]

msg_lines <- c(msg_lines, "")
msg_lines <- c(msg_lines, "Strategy counts (descending):")
for (i in seq_len(nrow(strategy_counts))) {
  msg_lines <- c(msg_lines, sprintf("  %3d  %s", strategy_counts$count[i], strategy_counts$strategy[i]))
}
# Save counts CSV
write.csv(strategy_counts, out_counts_csv, row.names = FALSE)

# Define categories of interest
cat_prodonly <- "Producer-only (no utilization)"
cat_multi_receptor <- "Multi-receptor producer"
cat_single <- "Single-receptor producer"
cat_multi_producer <- "Multi-producer"
cat_nonproducer <- "Nonproducer"

# Helper to sample rows and write CSV
sample_and_write <- function(df, strategy_name, max_n = 10, out_csv) {
  sel <- df[df$strategy == strategy_name, , drop = FALSE]
  n <- nrow(sel)
  if (n == 0) {
    return(list(n = 0, sample = NULL))
  }
  samp_n <- min(n, max_n)
  set.seed(42)
  samp_idx <- if (n <= samp_n) seq_len(n) else sample(seq_len(n), samp_n)
  sample_rows <- sel[samp_idx, , drop = FALSE]
  # Write CSV for review
  write.csv(sample_rows, out_csv, row.names = FALSE)
  list(n = n, sample = sample_rows)
}

# Sample Producer-only
prodonly_res <- sample_and_write(fg, cat_prodonly, max_n = 10, out_csv = out_prodonly_csv)
msg_lines <- c(msg_lines, "")
msg_lines <- c(msg_lines, sprintf("Producer-only groups: %d (sample written to %s)", prodonly_res$n, out_prodonly_csv))
if (!is.null(prodonly_res$sample)) {
  msg_lines <- c(msg_lines, "Sample rows (Producer-only):")
  # print key columns for sample
  pr_cols <- intersect(c("agent_id", "functional_group_id", "strain_count", "strategy", "n_production_groups", "n_usable_pyoverdines"), colnames(prodonly_res$sample))
  for (j in seq_len(nrow(prodonly_res$sample))) {
    row <- prodonly_res$sample[j, , drop = FALSE]
    msg_lines <- c(msg_lines, paste0("  - ", paste(sprintf("%s=%s", pr_cols, as.character(row[1, pr_cols])), collapse = " ; ")))
  }
} else {
  msg_lines <- c(msg_lines, "  (none found)")
}

# Sample Multi-receptor producers
multirec_res <- sample_and_write(fg, cat_multi_receptor, max_n = 10, out_csv = out_multirec_csv)
msg_lines <- c(msg_lines, "")
msg_lines <- c(msg_lines, sprintf("Multi-receptor producer groups: %d (sample written to %s)", multirec_res$n, out_multirec_csv))
if (!is.null(multirec_res$sample)) {
  msg_lines <- c(msg_lines, "Sample rows (Multi-receptor producer):")
  pr_cols <- intersect(c("agent_id", "functional_group_id", "strain_count", "strategy", "n_production_groups", "n_usable_pyoverdines"), colnames(multirec_res$sample))
  for (j in seq_len(nrow(multirec_res$sample))) {
    row <- multirec_res$sample[j, , drop = FALSE]
    msg_lines <- c(msg_lines, paste0("  - ", paste(sprintf("%s=%s", pr_cols, as.character(row[1, pr_cols])), collapse = " ; ")))
  }
} else {
  msg_lines <- c(msg_lines, "  (none found)")
}

# Also report counts for Multi-producer and Single-receptor producer (completeness)
count_single <- sum(fg$strategy == cat_single, na.rm = TRUE)
count_multi_prod <- sum(fg$strategy == cat_multi_producer, na.rm = TRUE)
count_nonproducer <- sum(fg$strategy == cat_nonproducer, na.rm = TRUE)
msg_lines <- c(msg_lines, "")
msg_lines <- c(msg_lines, sprintf("Single-receptor producer groups: %d", count_single))
msg_lines <- c(msg_lines, sprintf("Multi-producer groups: %d", count_multi_prod))
msg_lines <- c(msg_lines, sprintf("Nonproducer groups: %d", count_nonproducer))

# Write summary to text file
writeLines(msg_lines, out_txt)

# Also print to console
cat(paste(msg_lines, collapse = "\n"), "\n")

# Append session info (for provenance)
session_out <- capture.output(sessionInfo())
cat("\nSession info:\n")
cat(paste(session_out, collapse = "\n"), "\n", file = out_txt, append = TRUE)
