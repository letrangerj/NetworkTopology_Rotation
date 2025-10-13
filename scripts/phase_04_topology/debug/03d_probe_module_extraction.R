#!/usr/bin/env Rscript

# Debug: Probe module extraction from bipartite::computeModules across versions
# File: scripts/phase_04_topology/debug/03d_probe_module_extraction.R
#
# Purpose:
#   - Run computeModules() on FG-level production/utilization incidence matrices
#   - Systematically inspect the returned object structure
#   - Try multiple extraction strategies for module memberships:
#       1) Direct slot/list "modules" as list of two vectors (rows, cols)
#       2) Direct slot/list "modules" as one concatenated vector (rows + cols)
#       3) Direct slot/list "modules" as flat vector of n_rows + n_cols
#       4) Fallback to bipartite::listModuleInformation()
#   - Validate that the extracted membership lengths match nrow and ncol
#   - Report which strategy worked, counts of modules, and composition
#   - Save a compact probe log and (if successful) a CSV of module assignments
#
# Outputs:
#   - scripts/phase_04_topology/debug/logs/03d_module_extraction_probe.txt
#   - scripts/phase_04_topology/debug/module_assignments_<layer>_probe.csv (if membership extracted)
#
# Notes:
#   - This script does NOT alter any project data.
#   - It is safe to run repeatedly.

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(bipartite)
})

# -----------------------------
# Utilities
# -----------------------------
safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

capture_str <- function(x, max.level = 2) {
  # Return str() as character vector
  out <- utils::capture.output(str(x, max.level = max.level))
  if (length(out) == 0) out <- "<empty>"
  out
}

# Build layer-specific incidence matrix (FG × PYO)
build_layer_incidence <- function(adj_production, adj_utilization, layer = "production") {
  if (layer == "production") {
    incidence <- as.matrix(adj_production > 0) * 1
  } else if (layer == "utilization") {
    # adj_utilization is expected to be PYO × FG; transpose to FG × PYO
    incidence <- as.matrix(t(adj_utilization) > 0) * 1
  } else {
    stop("Layer must be 'production' or 'utilization'")
  }
  incidence
}

# Remove zero-degree rows/cols and return names removed (pre-prune space)
remove_zero_degree <- function(incidence) {
  row_deg <- rowSums(incidence)
  col_deg <- colSums(incidence)
  zero_r_idx <- which(row_deg == 0)
  zero_c_idx <- which(col_deg == 0)
  zero_r_names <- if (length(zero_r_idx) > 0) rownames(incidence)[zero_r_idx] else character(0)
  zero_c_names <- if (length(zero_c_idx) > 0) colnames(incidence)[zero_c_idx] else character(0)
  if (length(zero_r_idx) > 0) incidence <- incidence[-zero_r_idx, , drop = FALSE]
  if (length(zero_c_idx) > 0) incidence <- incidence[, -zero_c_idx, drop = FALSE]
  list(
    incidence = incidence,
    removed = list(
      zero_rows = zero_r_names,
      zero_cols = zero_c_names
    )
  )
}

# Compute modules with error handling
compute_modularity_safe <- function(incidence, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  res <- tryCatch(
    {
      bipartite::computeModules(
        incidence,
        method = "Beckett",
        deep = FALSE,
        deleteOriginalFiles = TRUE
      )
    },
    error = function(e) {
      message("computeModules failed: ", e$message)
      NULL
    }
  )
  res
}

# Attempt multiple extraction strategies; return standardized membership df or NULL
extract_membership <- function(module_result, incidence) {
  if (is.null(module_result)) {
    return(list(df = NULL, strategy = "none", notes = "module_result=NULL"))
  }

  # Helpers
  n_rows <- nrow(incidence)
  n_cols <- ncol(incidence)
  row_names <- rownames(incidence)
  col_names <- colnames(incidence)

  # Get Q (S4/S3 tolerant)
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

  # Get "modules" container (S4/S3 tolerant)
  mods <- tryCatch(
    {
      if (isS4(module_result) && "modules" %in% slotNames(module_result)) {
        slot(module_result, "modules")
      } else {
        module_result$modules
      }
    },
    error = function(e) NULL
  )

  notes <- character(0)

  # Strategy 1: list-of-two vectors (row, col)
  if (!is.null(mods) && is.list(mods) && length(mods) >= 2) {
    m1 <- mods[[1]]
    m2 <- mods[[2]]
    if (is.atomic(m1) && is.atomic(m2) && length(m1) == n_rows && length(m2) == n_cols) {
      df <- data.frame(
        node_id = c(row_names, col_names),
        node_type = c(rep("FG", n_rows), rep("PYOV", n_cols)),
        module_id = c(as.integer(m1), as.integer(m2)),
        stringsAsFactors = FALSE
      )
      return(list(df = df, strategy = "mods_list_two_vectors", notes = paste(notes, collapse = "; "), Q = Q))
    } else {
      notes <- c(notes, sprintf(
        "mods list-of-two not matching: len(m1)=%s, len(m2)=%s vs rows=%d cols=%d",
        length(m1), length(m2), n_rows, n_cols
      ))
    }
  }

  # Strategy 2: mods[[1]] is single concatenated vector
  if (!is.null(mods) && is.list(mods) && length(mods) >= 1 && is.atomic(mods[[1]])) {
    part <- mods[[1]]
    if (length(part) == (n_rows + n_cols)) {
      row_modules <- as.integer(part[1:n_rows])
      col_modules <- as.integer(part[(n_rows + 1):(n_rows + n_cols)])
      df <- data.frame(
        node_id = c(row_names, col_names),
        node_type = c(rep("FG", n_rows), rep("PYOV", n_cols)),
        module_id = c(row_modules, col_modules),
        stringsAsFactors = FALSE
      )
      return(list(df = df, strategy = "mods_concat_vector", notes = paste(notes, collapse = "; "), Q = Q))
    } else {
      notes <- c(notes, sprintf("mods[[1]] length=%s != %d", length(mods[[1]]), n_rows + n_cols))
    }
  }

  # Strategy 3: mods is a flat vector of length n_rows + n_cols
  if (!is.null(mods) && is.atomic(mods) && length(mods) == (n_rows + n_cols)) {
    row_modules <- as.integer(mods[1:n_rows])
    col_modules <- as.integer(mods[(n_rows + 1):(n_rows + n_cols)])
    df <- data.frame(
      node_id = c(row_names, col_names),
      node_type = c(rep("FG", n_rows), rep("PYOV", n_cols)),
      module_id = c(row_modules, col_modules),
      stringsAsFactors = FALSE
    )
    return(list(df = df, strategy = "mods_flat_concat", notes = paste(notes, collapse = "; "), Q = Q))
  }

  # Strategy 4: fallback to listModuleInformation()
  mi <- tryCatch(bipartite::listModuleInformation(module_result), error = function(e) NULL)
  if (!is.null(mi)) {
    # We don't know the exact structure; we will log and try to infer common shapes.
    # Common pattern: mi$modules is a data.frame with row/col module vectors or per-node mapping.
    if (!is.null(mi$modules) && is.data.frame(mi$modules)) {
      # Try columns that look like row/col modules
      cn <- colnames(mi$modules)
      # Heuristics: look for integer columns with lengths matching rows or cols
      candidates <- lapply(mi$modules, function(col) as.integer(col))
      # If there is a single vector of length n_rows + n_cols, try to align by names if present
      # Since this is version-dependent, we only extract structure summaries and return NULL here.
      notes <- c(notes, "listModuleInformation() available but structure is version-dependent; see logs for details")
      return(list(df = NULL, strategy = "listModuleInformation_only", notes = paste(notes, collapse = "; "), Q = Q, mi = mi))
    } else {
      notes <- c(notes, "listModuleInformation() returned non-data.frame modules or missing expected fields")
      return(list(df = NULL, strategy = "listModuleInformation_unhelpful", notes = paste(notes, collapse = "; "), Q = Q, mi = mi))
    }
  }

  list(df = NULL, strategy = "failed_all_strategies", notes = paste(notes, collapse = "; "), Q = Q)
}

summarize_membership <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(list(n_modules = 0L, composition = NULL))
  }
  n_modules <- length(unique(df$module_id))
  comp <- df %>%
    group_by(module_id, node_type) %>%
    summarise(count = n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = node_type, values_from = count, values_fill = 0L) %>%
    arrange(desc(FG + PYOV))
  list(n_modules = n_modules, composition = comp)
}

# -----------------------------
# Setup and I/O
# -----------------------------
safe_dir_create("scripts/phase_04_topology/debug/logs")

log_file <- "scripts/phase_04_topology/debug/logs/03d_module_extraction_probe.txt"
con <- file(log_file, open = "wt")

writeln <- function(...) {
  cat(paste0(..., "\n"), file = con, append = TRUE)
}

writeln("=== 03d Probe: Module Extraction Debug ===")
writeln(sprintf("Timestamp: %s", timestamp()))
writeln("")

# -----------------------------
# Load inputs
# -----------------------------
paths <- list(
  adj_prod_fg = "data/interim/adj_production_FG_agentsxpyov_conservative.rds",
  adj_util_fg = "data/interim/adj_utilization_FG_pyovxagents_conservative.rds"
)

missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing) > 0) {
  writeln("ERROR: Missing required input files:")
  for (nm in missing) writeln(sprintf("  - %s: %s", nm, paths[[nm]]))
  close(con)
  stop("Required inputs missing; see log for details.")
}

adj_prod_fg <- readRDS(paths$adj_prod_fg)
adj_util_fg <- readRDS(paths$adj_util_fg)

# -----------------------------
# Probe both layers
# -----------------------------
layers <- c("production", "utilization")
master_seed <- 2025

for (layer in layers) {
  writeln(sprintf("\n--- Layer: %s ---", layer))

  # Build and filter incidence
  inc_raw <- build_layer_incidence(adj_prod_fg, adj_util_fg, layer)
  pre_dims <- dim(inc_raw)
  pre_edges <- sum(inc_raw)
  filt <- remove_zero_degree(inc_raw)
  inc <- filt$incidence
  post_dims <- dim(inc)
  post_edges <- sum(inc)
  fill <- ifelse(nrow(inc) * ncol(inc) > 0, post_edges / (nrow(inc) * ncol(inc)), 0)

  writeln(sprintf("Matrix pre-filter: %d × %d, edges=%d", pre_dims[1], pre_dims[2], pre_edges))
  writeln(sprintf("Zero-degree removed: rows=%d, cols=%d", length(filt$removed$zero_rows), length(filt$removed$zero_cols)))
  writeln(sprintf("Matrix post-filter: %d × %d, edges=%d, fill=%.3f", post_dims[1], post_dims[2], post_edges, fill))

  if (sum(inc) == 0 || any(dim(inc) == 0)) {
    writeln("SKIP: No interactions after filtering.")
    next
  }

  # Compute modules (single replicate for probe)
  seed <- master_seed + ifelse(layer == "production", 100, 200)
  writeln(sprintf("Running computeModules (seed=%d)...", seed))
  mod_res <- compute_modularity_safe(inc, seed = seed)

  if (is.null(mod_res)) {
    writeln("computeModules returned NULL. Cannot proceed for this layer.")
    next
  }

  # Inspect structure
  writeln("module_result class/slots:")
  writeln(paste0("  class: ", paste(class(mod_res), collapse = ", ")))
  if (isS4(mod_res)) {
    writeln(paste0("  slotNames: ", paste(slotNames(mod_res), collapse = ", ")))
  } else {
    writeln(paste0("  names: ", paste(names(mod_res), collapse = ", ")))
  }

  # Try to get "modules" container and log its structure
  mods <- tryCatch(
    {
      if (isS4(mod_res) && "modules" %in% slotNames(mod_res)) slot(mod_res, "modules") else mod_res$modules
    },
    error = function(e) NULL
  )

  writeln("Structure of 'modules' container (str):")
  for (ln in capture_str(mods, max.level = 1)) writeln(paste0("  ", ln))

  # Also try listModuleInformation()
  mi <- tryCatch(bipartite::listModuleInformation(mod_res), error = function(e) NULL)
  if (!is.null(mi)) {
    writeln("Structure of listModuleInformation(module_result) (str):")
    for (ln in capture_str(mi, max.level = 1)) writeln(paste0("  ", ln))
  } else {
    writeln("listModuleInformation(): not available or failed.")
  }

  # Extract membership
  ext <- extract_membership(mod_res, inc)
  writeln(sprintf("Extraction strategy: %s", ext$strategy))
  if (!is.null(ext$notes) && nzchar(ext$notes)) writeln(sprintf("Notes: %s", ext$notes))
  if (!is.null(ext$Q) && !is.na(ext$Q)) writeln(sprintf("Q (likelihood): %.6f", ext$Q))

  if (!is.null(ext$df)) {
    sm <- summarize_membership(ext$df)
    writeln(sprintf("n_modules: %d", sm$n_modules))
    writeln("Top modules by size (FG + PYOV):")
    if (!is.null(sm$composition) && nrow(sm$composition) > 0) {
      # Print up to top 10
      topk <- head(sm$composition, 10)
      writeln(paste0(capture.output(print(topk)), collapse = "\n"))
    } else {
      writeln("  <no composition>")
    }

    # Save extracted membership CSV to debug folder
    out_csv <- sprintf("scripts/phase_04_topology/debug/module_assignments_%s_probe.csv", layer)
    suppressWarnings(write.csv(ext$df, out_csv, row.names = FALSE))
    writeln(sprintf("Saved module assignments (probe): %s", out_csv))
  } else {
    writeln("Failed to extract module assignments using all strategies.")
  }
}

writeln("\nSession Info:")
writeln(paste0(capture.output(sessionInfo()), collapse = "\n"))

close(con)

cat("03d probe complete. See log: scripts/phase_04_topology/debug/logs/03d_module_extraction_probe.txt\n")
