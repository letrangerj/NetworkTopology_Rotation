#!/usr/bin/env Rscript
#
# Phase 04 - Step 01: Preflight & Sanity Checks
#
# Purpose:
#   Perform environment, input presence, structural, and internal consistency checks
#   prior to running topology analyses (modularity, nestedness, degree metrics).
#
# Responsibilities:
#   1. Load parameter manifest (JSON preferred).
#   2. Verify existence of required Phase 03 adjacency & node artifacts.
#   3. Validate binary incidence properties (0/1, no NAs).
#   4. Cross-check dimensional concordance (FG / STR node counts vs matrices).
#   5. Confirm pyoverdine column/row name alignment with node dictionary.
#   6. Report inactive pyoverdines (no prod & no util).
#   7. Produce a structured status object + human-readable summary.
#
# Outputs:
#   - docs/phase_04/logs/phase_04_step01_preflight_summary.txt
#   - results/phase_04/preflight_status.rds
#   - (append / create) docs/phase_04/session_info.txt
#
# Exit codes:
#   - 0: All critical checks passed (only non-blocking warnings may exist).
#   - 1: One or more critical failures (downstream steps should not proceed).
#
# Usage:
#   Rscript scripts/phase_04_topology/01_preflight_and_sanity.R
#
# Dependencies: base R, jsonlite, Matrix, tools, (optional) fs
#
# Notes:
#   - Downstream scripts (02_degree_metrics, 03_modularity, etc.) should read
#     'results/phase_04/preflight_status.rds' and halt if status$critical_failures > 0.
#   - This script does NOT modify any adjacency or node data.
#

suppressPackageStartupMessages({
  library(jsonlite)
  library(Matrix)
})

message_banner <- function(...) {
  cat("\n=== ", paste(..., collapse = " "), " ===\n", sep = "")
}

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

# ------------------------------------------------------------------------------
# 0. Directory scaffolding (idempotent)
# ------------------------------------------------------------------------------
safe_dir_create("docs/phase_04")
safe_dir_create("docs/phase_04/logs")
safe_dir_create("results/phase_04")

log_path    <- "docs/phase_04/logs/phase_04_step01_preflight_summary.txt"
status_rds  <- "results/phase_04/preflight_status.rds"
session_log <- "docs/phase_04/session_info.txt"

# ------------------------------------------------------------------------------
# 1. Load parameter manifest
# ------------------------------------------------------------------------------
manifest_txt  <- "docs/phase_04/parameters_manifest.txt"
manifest_json <- "docs/phase_04/parameters_manifest.json"

manifest <- list()
manifest_source <- NA_character_

if (file.exists(manifest_json)) {
  manifest <- tryCatch(fromJSON(manifest_json), error = function(e) NULL)
  if (!is.null(manifest)) {
    manifest_source <- "json"
  }
}

if (is.null(manifest) || length(manifest) == 0) {
  if (file.exists(manifest_txt)) {
    # Fallback: parse minimal key: value lines from text (simple heuristic)
    raw_lines <- readLines(manifest_txt, warn = FALSE)
    kv <- grep("^[A-Za-z0-9_.-]+:\\s*", raw_lines, value = TRUE)
    if (length(kv) > 0) {
      for (ln in kv) {
        parts <- strsplit(ln, ":", fixed = TRUE)[[1]]
        key <- trimws(parts[1])
        val <- trimws(paste(parts[-1], collapse=":"))
        if (nzchar(key) && nzchar(val) && !(key %in% names(manifest))) {
          manifest[[key]] <- val
        }
      }
    }
    manifest_source <- "txt_partial"
  } else {
    manifest_source <- "absent"
  }
}

if (identical(manifest_source, "absent")) {
  warning("Parameter manifest not found; proceeding with default assumptions.")
}

# Helper accessor with default
mf_get <- function(key, default = NULL) {
  if (!is.null(manifest[[key]])) manifest[[key]] else default
}

# ------------------------------------------------------------------------------
# 2. Define expected input files (allow override via manifest)
# ------------------------------------------------------------------------------
expected_files <- list(
  adj_prod_fg  = "data/interim/adj_production_FG_agentsxpyov_conservative.rds",
  adj_util_fg  = "data/interim/adj_utilization_FG_pyovxagents_conservative.rds",
  adj_prod_str = "data/interim/adj_production_STR_agentsxpyov_conservative.rds",
  adj_util_str = "data/interim/adj_utilization_STR_pyovxagents_conservative.rds",
  nodes_fg     = "data/interim/nodes_functional_groups_conservative.rds",
  nodes_str    = "data/interim/nodes_strains_conservative.rds",
  nodes_pyov   = "data/interim/nodes_pyoverdines_conservative.rds",
  degree_seq   = "data/interim/degree_sequences_conservative.rds"
)

# If manifest defines input_files array (JSON form), just record for provenance
declared_inputs <- mf_get("input_files", character(0))

# ------------------------------------------------------------------------------
# 3. Initialize check containers
# ------------------------------------------------------------------------------
checks <- list()
add_check <- function(name, passed, severity = c("critical","warning","info"), details = "") {
  severity <- match.arg(severity)
  checks[[length(checks) + 1]] <<- list(
    name = name,
    passed = isTRUE(passed),
    severity = severity,
    details = details
  )
}

critical_fail <- function(name, condition, details_ok = "", details_fail = "") {
  if (isTRUE(condition)) {
    add_check(name, TRUE, "critical", details_ok)
    TRUE
  } else {
    add_check(name, FALSE, "critical", details_fail)
    FALSE
  }
}

warn_check <- function(name, condition, details_ok = "", details_fail = "") {
  if (isTRUE(condition)) {
    add_check(name, TRUE, "warning", details_ok)
    TRUE
  } else {
    add_check(name, FALSE, "warning", details_fail)
    FALSE
  }
}

info_check <- function(name, condition, details_ok = "", details_fail = "") {
  if (isTRUE(condition)) {
    add_check(name, TRUE, "info", details_ok)
    TRUE
  } else {
    add_check(name, FALSE, "info", details_fail)
    FALSE
  }
}

# ------------------------------------------------------------------------------
# 4. File existence checks
# ------------------------------------------------------------------------------
message_banner("Preflight: File Existence")
all_present <- TRUE
for (nm in names(expected_files)) {
  f <- expected_files[[nm]]
  exists_flag <- file.exists(f)
  add_check(
    name = paste0("exists_", nm),
    passed = exists_flag,
    severity = if (nm %in% c("adj_prod_fg","adj_util_fg","nodes_fg","nodes_pyov")) "critical" else "warning",
    details = if (exists_flag) f else paste("Missing:", f)
  )
  if (!exists_flag && nm %in% c("adj_prod_fg","adj_util_fg","nodes_fg","nodes_pyov")) {
    all_present <- FALSE
  }
  cat(sprintf("  %-22s : %s\n", f, if (exists_flag) "FOUND" else "MISSING"))
}

if (!all_present) {
  cat("\nOne or more CRITICAL files are missing. Downstream topology steps should NOT proceed.\n")
}

# Early abort possibility (do not abort yet—collect all diagnostics)
# ------------------------------------------------------------------------------
# 5. Load core objects (only if present)
# ------------------------------------------------------------------------------
load_or_null <- function(path) {
  if (file.exists(path)) {
    tryCatch(readRDS(path), error = function(e) {
      warning("Failed to read RDS: ", path, " :: ", e$message)
      NULL
    })
  } else NULL
}

adj_prod_fg  <- load_or_null(expected_files$adj_prod_fg)
adj_util_fg  <- load_or_null(expected_files$adj_util_fg)
adj_prod_str <- load_or_null(expected_files$adj_prod_str)
adj_util_str <- load_or_null(expected_files$adj_util_str)
nodes_fg     <- load_or_null(expected_files$nodes_fg)
nodes_str    <- load_or_null(expected_files$nodes_str)
nodes_pyov   <- load_or_null(expected_files$nodes_pyov)
degree_seq   <- load_or_null(expected_files$degree_seq)

# ------------------------------------------------------------------------------
# 6. Structural matrix checks (FG-level)
# ------------------------------------------------------------------------------
message_banner("Structural Checks: FG-level")

if (!is.null(adj_prod_fg) && !is.null(nodes_fg) && !is.null(nodes_pyov)) {
  # Expected dims: rows = FG count; cols = pyov count
  fg_row_ok <- nrow(adj_prod_fg) == nrow(nodes_fg)
  pyov_col_ok <- ncol(adj_prod_fg) == nrow(nodes_pyov)
  critical_fail("fg_prod_rows_match_nodes_fg", fg_row_ok,
                details_ok = sprintf("Rows in prod matrix: %d == FG nodes: %d", nrow(adj_prod_fg), nrow(nodes_fg)),
                details_fail = sprintf("Mismatch: prod rows=%d vs FG nodes=%d", nrow(adj_prod_fg), nrow(nodes_fg)))
  critical_fail("fg_prod_cols_match_nodes_pyov", pyov_col_ok,
                details_ok = sprintf("Cols in prod matrix: %d == Pyov nodes: %d", ncol(adj_prod_fg), nrow(nodes_pyov)),
                details_fail = sprintf("Mismatch: prod cols=%d vs Pyov nodes=%d", ncol(adj_prod_fg), nrow(nodes_pyov)))

  # Utilization matrix dims: rows = pyov count, cols = FG count
  if (!is.null(adj_util_fg)) {
    util_pyov_rows_ok <- nrow(adj_util_fg) == nrow(nodes_pyov)
    util_fg_cols_ok   <- ncol(adj_util_fg) == nrow(nodes_fg)
    critical_fail("fg_util_rows_match_nodes_pyov", util_pyov_rows_ok,
                  details_ok = "Util rows == pyov nodes",
                  details_fail = sprintf("Mismatch: util rows=%d vs pyov nodes=%d", nrow(adj_util_fg), nrow(nodes_pyov)))
    critical_fail("fg_util_cols_match_nodes_fg", util_fg_cols_ok,
                  details_ok = "Util cols == FG nodes",
                  details_fail = sprintf("Mismatch: util cols=%d vs FG nodes=%d", ncol(adj_util_fg), nrow(nodes_fg)))
  }

  # Name alignment
  if (!is.null(rownames(adj_prod_fg))) {
    name_alignment_fg <- identical(sort(rownames(adj_prod_fg)), sort(nodes_fg$agent_id))
    warn_check("fg_prod_rowname_alignment", name_alignment_fg,
               details_ok = "FG production rownames align with FG agent IDs (unordered match)",
               details_fail = "FG production rownames do not match FG agent IDs")
  } else {
    add_check("fg_prod_has_rownames", FALSE, "warning", "No rownames on FG production matrix")
  }
  if (!is.null(colnames(adj_prod_fg))) {
    pyov_name_alignment <- identical(sort(colnames(adj_prod_fg)), sort(nodes_pyov$node_id))
    warn_check("fg_prod_colname_alignment", pyov_name_alignment,
               details_ok = "Pyov colnames align with node_id",
               details_fail = "Pyov colnames mismatch node_id set")
  } else {
    add_check("fg_prod_has_colnames", FALSE, "warning", "No colnames on FG production matrix")
  }
}

# ------------------------------------------------------------------------------
# 7. Structural matrix checks (STR-level)
# ------------------------------------------------------------------------------
message_banner("Structural Checks: STR-level")

if (!is.null(adj_prod_str) && !is.null(nodes_str) && !is.null(nodes_pyov)) {
  str_prod_rows_ok <- nrow(adj_prod_str) == nrow(nodes_str)
  str_prod_cols_ok <- ncol(adj_prod_str) == nrow(nodes_pyov)
  critical_fail("str_prod_rows_match_nodes_str", str_prod_rows_ok,
                details_ok = "Strain production rows == strain nodes",
                details_fail = sprintf("Mismatch: rows=%d vs nodes=%d", nrow(adj_prod_str), nrow(nodes_str)))
  critical_fail("str_prod_cols_match_nodes_pyov", str_prod_cols_ok,
                details_ok = "Strain production cols == pyov nodes",
                details_fail = sprintf("Mismatch: cols=%d vs pyov=%d", ncol(adj_prod_str), nrow(nodes_pyov)))

  if (!is.null(adj_util_str)) {
    str_util_rows_ok <- nrow(adj_util_str) == nrow(nodes_pyov)
    str_util_cols_ok <- ncol(adj_util_str) == nrow(nodes_str)
    critical_fail("str_util_rows_match_nodes_pyov", str_util_rows_ok,
                  details_ok = "Strain utilization rows == pyov nodes",
                  details_fail = sprintf("Mismatch: rows=%d vs pyov=%d", nrow(adj_util_str), nrow(nodes_pyov)))
    critical_fail("str_util_cols_match_nodes_str", str_util_cols_ok,
                  details_ok = "Strain utilization cols == strain nodes",
                  details_fail = sprintf("Mismatch: cols=%d vs nodes=%d", ncol(adj_util_str), nrow(nodes_str)))
  }

  if (!is.null(rownames(adj_prod_str))) {
    name_alignment_str <- identical(sort(rownames(adj_prod_str)), sort(nodes_str$agent_id))
    warn_check("str_prod_rowname_alignment", name_alignment_str,
               details_ok = "Strain production rownames align with strain agent IDs",
               details_fail = "Strain production rownames mismatch agent IDs")
  }
  if (!is.null(colnames(adj_prod_str))) {
    pyov_name_alignment_str <- identical(sort(colnames(adj_prod_str)), sort(nodes_pyov$node_id))
    warn_check("str_prod_colname_alignment", pyov_name_alignment_str,
               details_ok = "Strain production colnames align with pyov node IDs",
               details_fail = "Strain production colnames mismatch pyov node IDs")
  }
}

# ------------------------------------------------------------------------------
# 8. Binary / NA integrity checks
# ------------------------------------------------------------------------------
check_binary <- function(mat, label) {
  if (is.null(mat)) return()
  vals <- unique(mat@x)
  non01 <- vals[!vals %in% c(0,1)]
  has_na <- anyNA(mat@x)
  add_check(
    name = paste0("binary_", label),
    passed = length(non01) == 0,
    severity = "critical",
    details = if (length(non01) == 0) "All entries in {0,1}" else paste("Non-binary values:", paste(non01, collapse=",")))
  add_check(
    name = paste0("no_NA_", label),
    passed = !has_na,
    severity = "critical",
    details = if (!has_na) "No NA entries" else "NA values detected"
  )
}

message_banner("Binary Integrity")
check_binary(adj_prod_fg,  "adj_prod_fg")
check_binary(adj_util_fg,  "adj_util_fg")
check_binary(adj_prod_str, "adj_prod_str")
check_binary(adj_util_str, "adj_util_str")

# ------------------------------------------------------------------------------
# 9. Pyoverdine inactivity report (FG-level primary perspective)
# ------------------------------------------------------------------------------
inactive_pyov_list <- character(0)
if (!is.null(adj_prod_fg) && !is.null(adj_util_fg)) {
  prod_counts <- Matrix::colSums(adj_prod_fg)
  util_counts <- Matrix::rowSums(adj_util_fg)
  inactive_idx <- which(prod_counts == 0 & util_counts == 0)
  if (length(inactive_idx) > 0) {
    inactive_pyov_list <- colnames(adj_prod_fg)[inactive_idx]
    add_check("inactive_pyoverdines_fg_level", FALSE, "warning",
              details = paste("Inactive pyov node_ids:", paste(inactive_pyov_list, collapse = ",")))
  } else {
    add_check("inactive_pyoverdines_fg_level", TRUE, "info", "No inactive pyoverdines at FG-level")
  }
}

# ------------------------------------------------------------------------------
# 10. Degree sequence reproducibility (optional soft check)
# ------------------------------------------------------------------------------
if (!is.null(degree_seq) && !is.null(adj_prod_fg) && !is.null(adj_util_fg)) {
  # Compare stored vs recomputed
  recomputed_fg_prod_out <- as.integer(Matrix::rowSums(adj_prod_fg))
  stored_fg_prod_out <- degree_seq$fg_network$fg_production_out
  if (!is.null(stored_fg_prod_out)) {
    same_length <- length(recomputed_fg_prod_out) == length(stored_fg_prod_out)
    same_values <- same_length && identical(recomputed_fg_prod_out, stored_fg_prod_out)
    warn_check("degree_sequence_fg_production_out_match", same_values,
               details_ok = "Stored FG production out-degree matches recomputed",
               details_fail = "Stored FG production out-degree mismatch (or length mismatch)")
  }
}

# ------------------------------------------------------------------------------
# 11. Summaries and counts
# ------------------------------------------------------------------------------
summary_items <- list(
  manifest_source        = manifest_source,
  n_fg_nodes             = if (!is.null(nodes_fg)) nrow(nodes_fg) else NA_integer_,
  n_str_nodes            = if (!is.null(nodes_str)) nrow(nodes_str) else NA_integer_,
  n_pyov_nodes           = if (!is.null(nodes_pyov)) nrow(nodes_pyov) else NA_integer_,
  fg_prod_edges          = if (!is.null(adj_prod_fg)) sum(adj_prod_fg) else NA_integer_,
  fg_util_edges          = if (!is.null(adj_util_fg)) sum(adj_util_fg) else NA_integer_,
  str_prod_edges         = if (!is.null(adj_prod_str)) sum(adj_prod_str) else NA_integer_,
  str_util_edges         = if (!is.null(adj_util_str)) sum(adj_util_str) else NA_integer_,
  inactive_pyov_count    = length(inactive_pyov_list),
  inactive_pyov_node_ids = inactive_pyov_list
)

# ------------------------------------------------------------------------------
# 12. Compile status object
# ------------------------------------------------------------------------------
checks_df <- do.call(rbind, lapply(checks, function(x) {
  data.frame(
    name = x$name,
    passed = x$passed,
    severity = x$severity,
    details = x$details,
    stringsAsFactors = FALSE
  )
}))

critical_failures <- sum(checks_df$severity == "critical" & !checks_df$passed)
warning_failures  <- sum(checks_df$severity == "warning" & !checks_df$passed)

status <- list(
  timestamp = timestamp(),
  summary = summary_items,
  checks = checks_df,
  critical_failures = critical_failures,
  warning_failures = warning_failures,
  manifest = manifest
)

saveRDS(status, status_rds)

# ------------------------------------------------------------------------------
# 13. Human-readable log
# ------------------------------------------------------------------------------
log_lines <- c(
  "Phase 04 - Step 01: Preflight & Sanity Check Summary",
  paste("Timestamp:", status$timestamp),
  "",
  sprintf("Manifest source: %s", manifest_source),
  "",
  "Node / Matrix Overview:",
  sprintf("  FG nodes         : %s", summary_items$n_fg_nodes),
  sprintf("  STR nodes        : %s", summary_items$n_str_nodes),
  sprintf("  Pyov nodes       : %s", summary_items$n_pyov_nodes),
  sprintf("  FG prod edges    : %s", summary_items$fg_prod_edges),
  sprintf("  FG util edges    : %s", summary_items$fg_util_edges),
  sprintf("  STR prod edges   : %s", summary_items$str_prod_edges),
  sprintf("  STR util edges   : %s", summary_items$str_util_edges),
  "",
  sprintf("Inactive pyoverdines (FG-level): %d", summary_items$inactive_pyov_count),
  if (length(inactive_pyov_list) > 0) paste0("  node_ids: ", paste(inactive_pyov_list, collapse = ", ")) else "",
  "",
  sprintf("Total checks: %d", nrow(checks_df)),
  sprintf("Critical failures: %d", critical_failures),
  sprintf("Warning failures : %d", warning_failures),
  "",
  "Detailed Checks:",
  paste(
    apply(checks_df, 1, function(r)
      sprintf("  [%s] %-40s : %s (%s)",
              toupper(r["severity"]),
              r["name"],
              if (identical(r["passed"], "TRUE") || r["passed"] == "TRUE") "PASS" else "FAIL",
              r["details"])),
    collapse = "\n"
  ),
  "",
  if (critical_failures == 0)
    "STATUS: PASS (no critical failures). Safe to proceed to Step 02 (degree metrics)."
  else
    "STATUS: FAIL (critical failures present). DO NOT proceed until resolved."
)

writeLines(log_lines, log_path)
cat(paste(log_lines, collapse = "\n"), "\n")

# ------------------------------------------------------------------------------
# 14. Session info (append if exists)
# ------------------------------------------------------------------------------
sess_info <- capture.output(sessionInfo())
if (file.exists(session_log)) {
  write("\n--- Append (Preflight) ---\n", file = session_log, append = TRUE)
  writeLines(sess_info, session_log, sep = "\n", useBytes = TRUE)
} else {
  writeLines(sess_info, session_log)
}

# ------------------------------------------------------------------------------
# 15. Exit code logic
# ------------------------------------------------------------------------------
if (critical_failures > 0) {
  cat("\nExiting with code 1 due to critical failures.\n")
  quit(save = "no", status = 1, runLast = FALSE)
} else {
  cat("\nPreflight completed successfully (no critical failures).\n")
  quit(save = "no", status = 0, runLast = FALSE)
}
