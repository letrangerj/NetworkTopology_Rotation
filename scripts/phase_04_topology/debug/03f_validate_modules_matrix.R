#!/usr/bin/env Rscript

# 03f_validate_modules_matrix.R
# Purpose:
#   Validate matrix-based module extraction for bipartite::computeModules results
#   using orderA/orderB to align nodes and handling the "+2" modules-matrix columns
#   pattern seen in the current bipartite version.
#
# Behavior:
#   - Loads FG-level production/utilization adjacency RDS files (Phase 3 outputs)
#   - Builds FG × PYOV binary incidence for each layer and prunes zero-degree nodes
#   - Runs computeModules() once per layer (seeded)
#   - Extracts node->module assignments by:
#       * Accessing module_result@modules (numeric matrix: modules × nodes)
#       * Using orderA (FG) and orderB (PYOV) to map column blocks
#       * Handling either ncol(mods) == nA+nB OR nA+nB+2 (drop last 2 columns)
#       * Taking column-wise argmax (apply(..., 2, which.max)) to assign each node to a module
#       * Realigning to original incidence order and saving assignments
#   - Writes a compact log and CSV/RDS debug artifacts
#
# Outputs:
#   - scripts/phase_04_topology/debug/logs/03f_validate_modules_matrix.txt
#   - scripts/phase_04_topology/debug/module_assignments_<layer>_validated.csv
#   - scripts/phase_04_topology/debug/modules_matrix_<layer>.rds
#   - scripts/phase_04_topology/debug/order_slots_<layer>.rds

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

log_path <- "scripts/phase_04_topology/debug/logs/03f_validate_modules_matrix.txt"
safe_dir_create(dirname(log_path))
log_con <- file(log_path, open = "wt")

wlog <- function(fmt, ...) {
  line <- sprintf(fmt, ...)
  cat(line, "\n", file = log_con, append = TRUE)
  cat(line, "\n")
}

wlog("=== 03f Validate Modules Matrix (orderA/orderB; col-wise argmax) ===")
wlog("Timestamp: %s", timestamp())

# -----------------------------
# Inputs
# -----------------------------
paths <- list(
  adj_prod_fg = "data/interim/adj_production_FG_agentsxpyov_conservative.rds",
  adj_util_fg = "data/interim/adj_utilization_FG_pyovxagents_conservative.rds"
)
missing_inputs <- names(paths)[!file.exists(unlist(paths))]
if (length(missing_inputs) > 0) {
  wlog("ERROR: Missing required inputs:")
  for (nm in missing_inputs) wlog("  - %s : %s", nm, paths[[nm]])
  close(log_con)
  quit(status = 1)
}

adj_prod_fg <- readRDS(paths$adj_prod_fg)
adj_util_fg <- readRDS(paths$adj_util_fg)
wlog("Loaded adjacencies: prod %d×%d, util %d×%d", nrow(adj_prod_fg), ncol(adj_prod_fg), nrow(adj_util_fg), ncol(adj_util_fg))

# -----------------------------
# Incidence helpers
# -----------------------------
build_incidence <- function(adj_prod, adj_util, layer) {
  if (layer == "production") {
    return(as.matrix(adj_prod > 0) * 1)
  }
  if (layer == "utilization") {
    return(as.matrix(t(adj_util) > 0) * 1)
  }
  stop("unknown layer")
}

prune_zero_degree <- function(inc) {
  rd <- rowSums(inc)
  cd <- colSums(inc)
  zr <- which(rd == 0)
  zc <- which(cd == 0)
  zr_names <- if (length(zr) > 0) rownames(inc)[zr] else character(0)
  zc_names <- if (length(zc) > 0) colnames(inc)[zc] else character(0)
  if (length(zr) > 0) inc <- inc[-zr, , drop = FALSE]
  if (length(zc) > 0) inc <- inc[, -zc, drop = FALSE]
  list(incidence = inc, removed = list(rows = zr_names, cols = zc_names))
}

compute_mod_safe <- function(inc, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  tryCatch(
    {
      bipartite::computeModules(inc, method = "Beckett", deep = FALSE, deleteOriginalFiles = TRUE)
    },
    error = function(e) list(error = TRUE, message = e$message)
  )
}

# -----------------------------
# Core extractor using modules matrix + orderA/orderB (col-wise argmax)
# -----------------------------
extract_with_order <- function(mod_res, incidence, layer_name, log = TRUE) {
  # Retrieve slots
  mods_mat <- tryCatch(
    {
      if (isS4(mod_res) && "modules" %in% slotNames(mod_res)) {
        slot(mod_res, "modules")
      } else {
        mod_res$modules
      }
    },
    error = function(e) NULL
  )

  orderA <- tryCatch(
    {
      if (isS4(mod_res) && "orderA" %in% slotNames(mod_res)) {
        slot(mod_res, "orderA")
      } else {
        mod_res$orderA
      }
    },
    error = function(e) NULL
  )

  orderB <- tryCatch(
    {
      if (isS4(mod_res) && "orderB" %in% slotNames(mod_res)) {
        slot(mod_res, "orderB")
      } else {
        mod_res$orderB
      }
    },
    error = function(e) NULL
  )

  Q <- tryCatch(
    {
      if (isS4(mod_res) && "likelihood" %in% slotNames(mod_res)) {
        as.numeric(slot(mod_res, "likelihood"))
      } else if (!is.null(mod_res$likelihood)) {
        as.numeric(mod_res$likelihood)
      } else {
        NA_real_
      }
    },
    error = function(e) NA_real_
  )

  if (is.null(mods_mat) || !is.matrix(mods_mat) || !is.numeric(mods_mat)) {
    return(list(success = FALSE, reason = "modules_missing_or_not_numeric_matrix", Q = Q))
  }

  inc_row_names <- rownames(incidence)
  inc_col_names <- colnames(incidence)
  nA <- if (!is.null(orderA)) length(orderA) else nrow(incidence)
  nB <- if (!is.null(orderB)) length(orderB) else ncol(incidence)
  nCols <- ncol(mods_mat)

  if (log) {
    wlog("[%s] modules matrix dims: %d x %d (rows x cols)", layer_name, nrow(mods_mat), ncol(mods_mat))
    wlog("[%s] orderA length: %d (inc rows: %d), orderB length: %d (inc cols: %d)", layer_name, nA, nrow(incidence), nB, ncol(incidence))
  }

  # Determine usable columns: exact match or "+2" trailing meta columns
  usable_cols <- NULL
  if (nCols == (nA + nB)) {
    usable_cols <- 1:(nA + nB)
    if (log) wlog("[%s] Using all %d columns (exact match).", layer_name, length(usable_cols))
  } else if (nCols == (nA + nB + 2)) {
    usable_cols <- 1:(nA + nB) # ignore last 2
    if (log) wlog("[%s] Detected +2 columns pattern; using first %d cols, dropping last 2.", layer_name, length(usable_cols))
  } else if (nCols > (nA + nB)) {
    usable_cols <- 1:(nA + nB)
    if (log) wlog("[%s] Warning: nCols=%d > nA+nB=%d. Using first %d columns.", layer_name, nCols, (nA + nB), length(usable_cols))
  } else {
    if (log) wlog("[%s] Mismatch: modules cols (%d) < required (nA+nB=%d).", layer_name, nCols, (nA + nB))
    return(list(success = FALSE, reason = "modules_cols_mismatch", Q = Q))
  }


  # Column-wise argmax across modules rows -> a vector of length = number of nodes (columns)

  # Exclude baseline row (row 1) for all layers, then add +1 to recover original row index

  assign_vec <- apply(mods_mat[-1, usable_cols, drop = FALSE], 2, which.max) + 1L

  if (length(assign_vec) != (nA + nB)) {
    if (log) wlog("[%s] Assignment length %d != nA+nB %d", layer_name, length(assign_vec), (nA + nB))
    return(list(success = FALSE, reason = "assign_length_mismatch", Q = Q))
  }

  # Map ordered names; if orderA/orderB missing, assume identity
  fg_names_ordered <- if (!is.null(orderA)) inc_row_names[orderA] else inc_row_names
  py_names_ordered <- if (!is.null(orderB)) inc_col_names[orderB] else inc_col_names

  # Sanity: lengths must match
  if (length(fg_names_ordered) != nA || length(py_names_ordered) != nB) {
    if (log) wlog("[%s] orderA/orderB size mismatch with incidence names", layer_name)
    return(list(success = FALSE, reason = "order_size_mismatch", Q = Q))
  }

  fg_assign_ordered <- assign_vec[1:nA]
  py_assign_ordered <- assign_vec[(nA + 1):(nA + nB)]

  # Realign to original incidence order using match
  fg_idx <- match(inc_row_names, fg_names_ordered)
  py_idx <- match(inc_col_names, py_names_ordered)
  if (any(is.na(fg_idx)) || any(is.na(py_idx))) {
    if (log) {
      wlog("[%s] Alignment failed: NA indices present (fg NAs=%d, pyov NAs=%d)", layer_name, sum(is.na(fg_idx)), sum(is.na(py_idx)))
    }
    return(list(success = FALSE, reason = "alignment_failed", Q = Q))
  }

  final_fg_assign <- fg_assign_ordered[fg_idx]
  final_py_assign <- py_assign_ordered[py_idx]

  modules_df <- data.frame(
    node_id = c(inc_row_names, inc_col_names),
    node_type = c(rep("FG", length(inc_row_names)), rep("PYOV", length(inc_col_names))),
    module_id = c(as.integer(final_fg_assign), as.integer(final_py_assign)),
    stringsAsFactors = FALSE
  )

  n_modules <- length(unique(modules_df$module_id))
  if (log) wlog("[%s] Extraction OK: %d nodes, %d unique modules; Q=%s", layer_name, nrow(modules_df), n_modules, ifelse(is.na(Q), "NA", sprintf("%.6f", Q)))

  list(
    success = TRUE, Q = Q, modules_df = modules_df,
    orderA = orderA, orderB = orderB, mods_mat = mods_mat
  )
}

# -----------------------------
# Driver
# -----------------------------
layers <- c("production", "utilization")
master_seed <- 2025

for (layer in layers) {
  wlog("")
  wlog("--- Layer: %s ---", layer)
  inc_raw <- build_incidence(adj_prod_fg, adj_util_fg, layer)
  wlog("Raw incidence: %d×%d; edges=%d", nrow(inc_raw), ncol(inc_raw), sum(inc_raw))
  pr <- prune_zero_degree(inc_raw)
  inc <- pr$incidence
  wlog("Pruned incidence: %d×%d; edges=%d; removed rows=%d, cols=%d", nrow(inc), ncol(inc), sum(inc), length(pr$removed$rows), length(pr$removed$cols))
  if (sum(inc) == 0 || nrow(inc) == 0 || ncol(inc) == 0) {
    wlog("Skipping (empty incidence)")
    next
  }

  seed <- master_seed + ifelse(layer == "production", 100, 200)
  wlog("computeModules(seed=%d)...", seed)
  mod_res <- compute_mod_safe(inc, seed = seed)
  if (is.list(mod_res) && !is.null(mod_res$error) && mod_res$error) {
    wlog("computeModules failed: %s", mod_res$message)
    next
  }
  wlog("computeModules class: %s", paste(class(mod_res), collapse = ","))

  ext <- extract_with_order(mod_res, inc, layer_name = layer, log = TRUE)
  if (isTRUE(ext$success)) {
    # Save validated assignments
    out_csv <- sprintf("scripts/phase_04_topology/debug/module_assignments_%s_validated.csv", layer)
    write.csv(ext$modules_df, out_csv, row.names = FALSE)
    wlog("Saved validated assignments: %s", out_csv)

    # Save raw modules matrix and order slots
    out_mods <- sprintf("scripts/phase_04_topology/debug/modules_matrix_%s.rds", layer)
    saveRDS(ext$mods_mat, out_mods)
    wlog("Saved modules matrix (RDS): %s", out_mods)

    out_orders <- sprintf("scripts/phase_04_topology/debug/order_slots_%s.rds", layer)
    saveRDS(list(orderA = ext$orderA, orderB = ext$orderB), out_orders)
    wlog("Saved orderA/orderB slots (RDS): %s", out_orders)

    # Composition summary (top 10 modules by total nodes)
    comp <- ext$modules_df %>%
      group_by(module_id, node_type) %>%
      tally(name = "count") %>%
      pivot_wider(names_from = node_type, values_from = count, values_fill = 0) %>%
      arrange(desc(FG + PYOV))
    wlog("Top modules (module_id | FG | PYOV):")
    if (nrow(comp) == 0) {
      wlog("  <no modules detected>")
    } else {
      topk <- head(comp, 10)
      for (i in seq_len(nrow(topk))) {
        row <- topk[i, ]
        fg_val <- ifelse(!is.null(row$FG), as.integer(row$FG), 0L)
        py_val <- ifelse(!is.null(row$PYOV), as.integer(row$PYOV), 0L)
        wlog("  %s | %d | %d", as.character(row$module_id), fg_val, py_val)
      }
    }
  } else {
    wlog("FAIL: %s", ifelse(!is.null(ext$reason), ext$reason, "unknown"))
  }
}

# -----------------------------
# Session info
# -----------------------------
wlog("")
wlog("Session info:")
wlog("%s", paste(capture.output(sessionInfo()), collapse = "\n"))
close(log_con)
cat("03f validation completed. See log:", log_path, "\n")
