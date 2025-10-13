#!/usr/bin/env Rscript

# 03f_validate_modules_matrix.R
# Purpose:
#   Validate matrix-based module extraction for bipartite::computeModules results
#   (extract module membership by using the module × node numeric matrix and max.col)
#
# Behavior:
#   - Loads production and utilization FG-level adjacency matrices (Phase 3 outputs)
#   - Builds binary FG × PYOV incidence matrices for each layer
#   - Runs computeModules() once per layer (seeded)
#   - Inspects the returned object's 'modules' slot:
#       - verifies it's a numeric matrix
#       - checks row/col names and dimensions
#       - extracts node->module assignment using max.col on the modules matrix
#   - Aligns extracted assignments to the filtered incidence row/col names
#   - Writes a compact log and (if successful) a CSV with assignments
#
# Outputs:
#   - scripts/phase_04_topology/debug/logs/03f_validate_modules_matrix.txt
#   - scripts/phase_04_topology/debug/module_assignments_<layer>_validated.csv  (when extraction succeeds)
#   - scripts/phase_04_topology/debug/modules_matrix_<layer>.rds
#
# Notes:
#   - This is a pure diagnostic / validation script and is safe to run repeatedly.
#   - It does not modify canonical project data.
#   - Ensure phase-03 adjacency RDS files exist.

suppressPackageStartupMessages({
  library(bipartite)
  library(Matrix)
  library(dplyr)
  library(tidyr)
})

# -----------------------------
# Utilities
# -----------------------------
safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

# Simple logger
log_path <- "scripts/phase_04_topology/debug/logs/03f_validate_modules_matrix.txt"
safe_dir_create(dirname(log_path))
log_con <- file(log_path, open = "wt")
wlog <- function(fmt, ...) {
  line <- sprintf(fmt, ...)
  cat(line, "\n", file = log_con, append = TRUE)
  cat(line, "\n")
}

wlog("=== 03f Validate Modules Matrix ===")
wlog("Timestamp: %s", timestamp())
wlog("Working directory: %s", getwd())
wlog("")

# -----------------------------
# Input checks
# -----------------------------
paths <- list(
  adj_prod_fg = "data/interim/adj_production_FG_agentsxpyov_conservative.rds",
  adj_util_fg = "data/interim/adj_utilization_FG_pyovxagents_conservative.rds"
)

missing_inputs <- names(paths)[!file.exists(unlist(paths))]
if (length(missing_inputs) > 0) {
  wlog("ERROR: Missing required input files:")
  for (nm in missing_inputs) wlog("  - %s : %s", nm, paths[[nm]])
  wlog("Aborting validation.")
  close(log_con)
  quit(status = 1)
}

# Load adjacency matrices
adj_prod_fg <- readRDS(paths$adj_prod_fg)
adj_util_fg <- readRDS(paths$adj_util_fg)
wlog(
  "Loaded adjacency matrices: prod %d x %d, util %d x %d",
  nrow(adj_prod_fg), ncol(adj_prod_fg), nrow(adj_util_fg), ncol(adj_util_fg)
)

# -----------------------------
# Helpers: incidence build and prune
# -----------------------------
build_incidence <- function(adj_prod, adj_util, layer) {
  if (layer == "production") {
    inc <- as.matrix(adj_prod > 0) * 1
  } else if (layer == "utilization") {
    inc <- as.matrix(t(adj_util) > 0) * 1
  } else {
    stop("Unknown layer")
  }
  inc
}

prune_zero_degree <- function(inc) {
  row_deg <- rowSums(inc)
  col_deg <- colSums(inc)
  zero_rows <- which(row_deg == 0)
  zero_cols <- which(col_deg == 0)
  zero_row_names <- if (length(zero_rows) > 0) rownames(inc)[zero_rows] else character(0)
  zero_col_names <- if (length(zero_cols) > 0) colnames(inc)[zero_cols] else character(0)
  if (length(zero_rows) > 0) inc <- inc[-zero_rows, , drop = FALSE]
  if (length(zero_cols) > 0) inc <- inc[, -zero_cols, drop = FALSE]
  list(incidence = inc, removed = list(rows = zero_row_names, cols = zero_col_names))
}

compute_mod_safe <- function(inc, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  res <- tryCatch(
    {
      bipartite::computeModules(inc, method = "Beckett", deep = FALSE, deleteOriginalFiles = TRUE)
    },
    error = function(e) {
      list(error = TRUE, message = e$message)
    }
  )
  res
}

# Extract assignments from module matrix (matrix-based extractor)
extract_from_modules_matrix <- function(module_result, incidence) {
  # tolerant retrieval of modules container
  mods_mat <- tryCatch(
    {
      if (isS4(module_result) && "modules" %in% slotNames(module_result)) {
        slot(module_result, "modules")
      } else if (!is.null(module_result$modules)) {
        module_result$modules
      } else {
        NULL
      }
    },
    error = function(e) NULL
  )

  if (is.null(mods_mat)) {
    return(list(success = FALSE, reason = "no_modules_slot"))
  }

  if (!is.matrix(mods_mat) || !is.numeric(mods_mat)) {
    return(list(success = FALSE, reason = "modules_not_numeric_matrix"))
  }

  # Q (likelihood) attempt
  Q <- tryCatch(
    {
      if (isS4(module_result) && "likelihood" %in% slotNames(module_result)) {
        as.numeric(slot(module_result, "likelihood"))
      } else if (!is.null(module_result$likelihood)) {
        as.numeric(module_result$likelihood)
      } else {
        NA_real_
      }
    },
    error = function(e) NA_real_
  )

  # node names in modules matrix
  mod_colnames <- colnames(mods_mat)
  mod_rownames <- rownames(mods_mat)
  n_mod_cols <- ncol(mods_mat)

  # our incidence nodes
  fg_names <- rownames(incidence)
  pyov_names <- colnames(incidence)
  n_fg <- length(fg_names)
  n_pyov <- length(pyov_names)
  n_inc_nodes <- n_fg + n_pyov

  wlog("modules matrix dims: %d x %d (rows x cols)", nrow(mods_mat), n_mod_cols)
  if (!is.null(mod_colnames)) {
    wlog("modules colnames available: %d", length(mod_colnames))
    # report how many map to our FG/PYOV names
    fg_map_count <- sum(mod_colnames %in% fg_names)
    py_map_count <- sum(mod_colnames %in% pyov_names)
    wlog("  mapping: %d match FG, %d match PYOV", fg_map_count, py_map_count)
  } else {
    wlog("modules matrix has no column names")
  }

  # Strategy A: name-based mapping if full coverage
  if (!is.null(mod_colnames)) {
    fg_mask <- mod_colnames %in% fg_names
    py_mask <- mod_colnames %in% pyov_names
    if (sum(fg_mask) == n_fg && sum(py_mask) == n_pyov) {
      assign_vec <- max.col(mods_mat, ties.method = "first")
      df_all <- data.frame(node_name = mod_colnames, module = assign_vec, stringsAsFactors = FALSE)
      df_fg <- df_all[df_all$node_name %in% fg_names, , drop = FALSE][match(fg_names, df_all$node_name[df_all$node_name %in% fg_names]), , drop = FALSE]
      df_py <- df_all[df_all$node_name %in% pyov_names, , drop = FALSE][match(pyov_names, df_all$node_name[df_all$node_name %in% pyov_names]), , drop = FALSE]
      # check alignment
      if (nrow(df_fg) == n_fg && nrow(df_py) == n_pyov && !any(is.na(df_fg$module)) && !any(is.na(df_py$module))) {
        modules_df <- data.frame(
          node_id = c(df_fg$node_name, df_py$node_name),
          node_type = c(rep("FG", nrow(df_fg)), rep("PYOV", nrow(df_py))),
          module_id = c(as.integer(df_fg$module), as.integer(df_py$module)),
          stringsAsFactors = FALSE
        )
        return(list(success = TRUE, modules_df = modules_df, mods_mat = mods_mat, Q = Q))
      } else {
        # fall through to index-based fallback
        wlog("Name-based mapping incomplete, will attempt index-based fallback")
      }
    } else {
      wlog("Name-based mapping does not cover all nodes (fg_map=%d/%d, py_map=%d/%d)", sum(mod_colnames %in% fg_names), n_fg, sum(mod_colnames %in% pyov_names), n_pyov)
    }
  }

  # Strategy B: index-based fallback when column counts match
  if (n_mod_cols == n_inc_nodes) {
    wlog("Index-based fallback: modules columns (%d) == incidence nodes (%d). Assuming FG then PYOV ordering.", n_mod_cols, n_inc_nodes)
    assign_vec <- max.col(mods_mat, ties.method = "first")
    fg_assign <- assign_vec[1:n_fg]
    py_assign <- assign_vec[(n_fg + 1):(n_fg + n_pyov)]
    modules_df <- data.frame(
      node_id = c(fg_names, pyov_names),
      node_type = c(rep("FG", n_fg), rep("PYOV", n_pyov)),
      module_id = c(as.integer(fg_assign), as.integer(py_assign)),
      stringsAsFactors = FALSE
    )
    return(list(success = TRUE, modules_df = modules_df, mods_mat = mods_mat, Q = Q))
  }

  # If we get here, cannot align automatically
  wlog("Unable to align modules matrix automatically: cols_in_mods=%d, nodes_in_inc=%d", n_mod_cols, n_inc_nodes)
  return(list(success = FALSE, reason = "dim_mismatch_or_unaligned", mods_cols = n_mod_cols, inc_nodes = n_inc_nodes, mods_colnames = mod_colnames, mods_rownames = mod_rownames))
}

# -----------------------------
# Run validation for both layers
# -----------------------------
layers <- c("production", "utilization")
master_seed <- 2025

for (layer in layers) {
  wlog("\n--- Layer: %s ---", layer)
  inc_raw <- build_incidence(adj_prod_fg, adj_util_fg, layer)
  wlog("Raw incidence dims: %d x %d; edges=%d", nrow(inc_raw), ncol(inc_raw), sum(inc_raw))
  pr <- prune_zero_degree(inc_raw)
  inc <- pr$incidence
  wlog("Pruned incidence dims: %d x %d; edges=%d; removed rows=%d cols=%d", nrow(inc), ncol(inc), sum(inc), length(pr$removed$rows), length(pr$removed$cols))
  if (sum(inc) == 0 || nrow(inc) == 0 || ncol(inc) == 0) {
    wlog("Skipping layer due to empty incidence after pruning.")
    next
  }

  seed <- master_seed + ifelse(layer == "production", 100, 200)
  wlog("Calling computeModules() with seed=%d", seed)
  mod_res <- compute_mod_safe(inc, seed = seed)
  if (is.list(mod_res) && !is.null(mod_res$error) && mod_res$error) {
    wlog("computeModules failed: %s", mod_res$message)
    next
  }
  wlog("computeModules returned object class: %s", paste(class(mod_res), collapse = ", "))

  # Try extraction
  ext <- extract_from_modules_matrix(mod_res, inc)
  if (isTRUE(ext$success)) {
    wlog("SUCCESS: extracted %d module assignments (unique modules = %d)", nrow(ext$modules_df), length(unique(ext$modules_df$module_id)))
    # save validated assignments
    out_csv <- sprintf("scripts/phase_04_topology/debug/module_assignments_%s_validated.csv", layer)
    safe_dir_create(dirname(out_csv))
    write.csv(ext$modules_df, out_csv, row.names = FALSE)
    wlog("Saved validated assignments: %s", out_csv)
    # save raw modules matrix for inspection
    out_rds <- sprintf("scripts/phase_04_topology/debug/modules_matrix_%s.rds", layer)
    saveRDS(ext$mods_mat, out_rds)
    wlog("Saved raw modules matrix (RDS): %s", out_rds)
    # summarize top modules
    comp <- ext$modules_df %>%
      group_by(module_id, node_type) %>%
      tally(name = "count") %>%
      pivot_wider(names_from = node_type, values_from = count, values_fill = 0) %>%
      arrange(desc(FG + PYOV))
    wlog("Top modules (module_id | FG | PYOV):")
    topk <- head(comp, 10)
    if (nrow(topk) == 0) {
      wlog("  <no modules found in summary>")
    } else {
      for (i in seq_len(nrow(topk))) {
        row <- topk[i, ]
        fg_val <- ifelse(!is.null(row$FG), as.integer(row$FG), 0L)
        py_val <- ifelse(!is.null(row$PYOV), as.integer(row$PYOV), 0L)
        wlog("  %s | %d | %d", as.character(row$module_id), fg_val, py_val)
      }
    }
  } else {
    wlog("FAIL: Could not extract modules for layer '%s'. Reason: %s", layer, ifelse(!is.null(ext$reason), ext$reason, "unknown"))
    if (!is.null(ext$mods_cols)) wlog("  modules cols: %d, incidence nodes: %d", ext$mods_cols, ext$inc_nodes)
    if (!is.null(ext$mods_colnames)) wlog("  modules colnames (sample): %s", paste(head(ext$mods_colnames, 20), collapse = ", "))
    if (!is.null(ext$mods_rownames)) wlog("  modules rownames (sample): %s", paste(head(ext$mods_rownames, 20), collapse = ", "))
    wlog("Recommendation: inspect the saved raw module object via bipartite::listModuleInformation() or adapt extractor to the module object layout.")
  }
}

wlog("")
wlog("Session info:")
wlog("%s", paste(capture.output(sessionInfo()), collapse = "\n"))
close(log_con)
cat("03f validation script finished. See log:", log_path, "\n")
