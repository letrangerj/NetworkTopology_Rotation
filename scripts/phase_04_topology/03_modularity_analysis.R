#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(bipartite)
  library(stringr)
  library(tibble)
  library(pbapply)
  library(rlang)
  library(parallel)
})

pboptions(type = "txt")

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

build_matrix <- function(edges, row_col, col_col) {
  row_sym <- sym(row_col)
  col_sym <- sym(col_col)

  edges %>%
    distinct(!!row_sym, !!col_sym) %>%
    mutate(value = 1L) %>%
    pivot_wider(names_from = !!col_sym, values_from = value, values_fill = 0L) %>%
    arrange(!!row_sym) %>%
    column_to_rownames(var = row_col) %>%
    as.matrix()
}

module_to_vector <- function(res) {
  if (!inherits(res, "moduleWeb")) {
    stop("Unsupported object returned by computeModules")
  }

  mods_mat <- slot(res, "modules")
  if (is.null(mods_mat) || !is.matrix(mods_mat)) {
    stop("computeModules result lacks modules matrix")
  }

  data_mat <- slot(res, "data")
  fg_names <- rownames(data_mat)
  py_names <- colnames(data_mat)

  orderA <- slot(res, "orderA")
  orderB <- slot(res, "orderB")

  nA <- if (length(orderA) > 0) length(orderA) else length(fg_names)
  nB <- if (length(orderB) > 0) length(orderB) else length(py_names)
  expected_cols <- nA + nB

  if (ncol(mods_mat) < expected_cols) {
    stop("Module matrix has fewer columns than expected nodes")
  }

  membership <- mods_mat[-1, 1:expected_cols, drop = FALSE]
  assignment_idx <- apply(membership, 2, which.max)

  fg_order <- if (length(orderA) > 0) fg_names[orderA] else fg_names
  py_order <- if (length(orderB) > 0) py_names[orderB] else py_names

  fg_assign <- assignment_idx[seq_len(nA)]
  py_assign <- assignment_idx[nA + seq_len(nB)]

  fg_assign <- fg_assign[match(fg_names, fg_order)]
  py_assign <- py_assign[match(py_names, py_order)]

  setNames(c(fg_assign, py_assign), c(fg_names, py_names))
}

parallel_apply <- function(X, cfg, fun, base_seed = NULL) {
  if ("pbmclapply" %in% getNamespaceExports("pbapply") && .Platform$OS.type != "windows") {
    pbapply::pbmclapply(
      X,
      fun,
      mc.cores = cfg$worker_count,
      mc.preschedule = FALSE,
      mc.set.seed = TRUE
    )
  } else {
    cl <- parallel::makeCluster(cfg$worker_count)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    if (!is.null(base_seed)) {
      parallel::clusterSetRNGStream(cl, iseed = base_seed)
    } else if (!is.null(cfg$master_seed)) {
      parallel::clusterSetRNGStream(cl, iseed = cfg$master_seed)
    }
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages({
      library(bipartite)
      library(tibble)
      library(dplyr)
      library(purrr)
      library(stringr)
      library(rlang)
    }))
    pbapply::pblapply(X, fun, cl = cl)
  }
}

run_modules <- function(mat, seeds, cfg, base_seed = NULL) {
  results <- parallel_apply(
    seeds,
    cfg,
    base_seed = base_seed,
    fun = function(seed) {
      tryCatch({
        res <- computeModules(mat, method = cfg$method, steps = cfg$replicate_steps)
        modules_vec <- module_to_vector(res)
        tibble(
          replicate = seed,
          Q = as.numeric(slot(res, "likelihood")),
          modules = list(modules_vec)
        )
      }, error = function(e) {
        warning(sprintf("computeModules failed for seed %s: %s", seed, e$message), call. = FALSE)
        tibble(replicate = seed, Q = NA_real_, modules = list(NULL))
      })
    }
  )
  bind_rows(results)
}

null_distribution <- function(mat, seeds, cfg, base_seed = NULL) {
  unlist(parallel_apply(
    seeds,
    cfg,
    base_seed = base_seed,
    fun = function(seed) {
      tryCatch({
        null_mat <- swap.web(mat, N = cfg$swap_iterations)
        res <- computeModules(null_mat, method = cfg$method, steps = cfg$null_steps)
        as.numeric(slot(res, "likelihood"))
      }, error = function(e) {
        warning(sprintf("Null model failed for seed %s: %s", seed, e$message), call. = FALSE)
        NA_real_
      })
    }
  ))
}

summarise_nulls <- function(q_obs, q_null) {
  q_valid <- q_null[!is.na(q_null)]
  mu <- mean(q_valid)
  sigma <- sd(q_valid)
  z <- if (sigma > 0) (q_obs - mu) / sigma else NA_real_
  p <- mean(q_valid >= q_obs)
  tibble(Q_obs = q_obs, Q_null_mean = mu, Q_null_sd = sigma, z_score = z, p_value = p)
}

save_membership <- function(mod_vec, layer_name, path) {
  df <- tibble(node_id = names(mod_vec), module_id = as.integer(mod_vec)) %>%
    mutate(
      node_type = case_when(
        str_starts(node_id, "FG_") ~ "FG",
        str_starts(node_id, "PYO_") ~ "PYOV",
        TRUE ~ "UNKNOWN"
      ),
      layer = layer_name
    )
  saveRDS(df, path)
  df
}

# configuration
config <- list(
  edges_file = "data/interim/edges_functional_groups_conservative.csv",
  replicate_seeds = 1:2,
  null_seeds = 1:10,
  swap_iterations = 1e4,
  method = "Beckett",
  replicate_steps = 5e3,
  null_steps = 5e3,
  master_seed = 2025,
  worker_count = 2,
  output_dir = "results/phase_04/modularity"
)

safe_dir_create(config$output_dir)
safe_dir_create("docs/phase_04/logs")
cat("=== Phase 04 Step 03: Modularity via pollination-style workflow ===\n")
cat("Timestamp:", timestamp(), "\n\n")

if (!file.exists(config$edges_file)) stop("Missing edges file: ", config$edges_file)

edges <- read_csv(config$edges_file, show_col_types = FALSE)

production_mat <- edges %>%
  filter(edge_type == "production") %>%
  transmute(row = source, col = target)

utilization_mat <- edges %>%
  filter(edge_type == "utilization") %>%
  transmute(row = target, col = source)

if (nrow(production_mat) == 0 || nrow(utilization_mat) == 0) stop("Edges missing for one or both layers")

prod_incidence <- build_matrix(production_mat, "row", "col")
util_incidence <- build_matrix(utilization_mat, "row", "col")

layers <- list(
  production = prod_incidence,
  utilization = util_incidence
)

analyze_layer <- function(layer_name, mat, cfg, layer_index) {
  work_mat <- mat
  if (any(rowSums(work_mat) == 0) || any(colSums(work_mat) == 0)) {
    work_mat <- work_mat[rowSums(work_mat) > 0, colSums(work_mat) > 0, drop = FALSE]
  }

  cat("Layer:", layer_name, "| matrix dims:", nrow(work_mat), "x", ncol(work_mat), "| density:", signif(sum(work_mat) / (nrow(work_mat) * ncol(work_mat)), 3), "\n")

  base_seed_rep <- if (!is.null(cfg$master_seed)) {
    RNGkind("L'Ecuyer-CMRG")
    cfg$master_seed + layer_index
  } else {
    NULL
  }
  if (!is.null(base_seed_rep)) set.seed(base_seed_rep)

  replicate_tbl <- run_modules(work_mat, cfg$replicate_seeds, cfg, base_seed = base_seed_rep) %>%
    filter(!is.na(Q))
  if (nrow(replicate_tbl) == 0) {
    stop("No successful module detection for layer: ", layer_name)
  }

  best_row <- replicate_tbl %>% arrange(desc(Q)) %>% slice(1)
  modules_vec <- best_row$modules[[1]]
  if (is.null(modules_vec) || length(modules_vec) == 0) {
    stop("Empty module assignment for layer: ", layer_name)
  }

  assignments_path <- file.path(cfg$output_dir, paste0("module_assignments_", layer_name, ".rds"))
  membership_df <- save_membership(modules_vec, layer_name, assignments_path)

  replicate_out <- replicate_tbl %>% select(replicate, Q) %>% mutate(layer = layer_name)
  replicate_file <- file.path(cfg$output_dir, paste0("replicate_Q_", layer_name, ".csv"))
  write_csv(replicate_out, replicate_file)

  base_seed_null <- if (!is.null(cfg$master_seed)) cfg$master_seed + 1000 + layer_index else NULL
  if (!is.null(base_seed_null)) set.seed(base_seed_null)

  q_null <- null_distribution(work_mat, cfg$null_seeds, cfg, base_seed = base_seed_null)
  null_summary <- summarise_nulls(best_row$Q, q_null)
  null_file <- file.path(cfg$output_dir, paste0("null_summary_", layer_name, ".csv"))
  write_csv(null_summary, null_file)

  list(
    Q_best = best_row$Q,
    replicate_stats = replicate_tbl %>% summarise(mean_Q = mean(Q), sd_Q = sd(Q), min_Q = min(Q), max_Q = max(Q)),
    membership = membership_df,
    null_summary = null_summary,
    assignments_path = assignments_path,
    replicate_file = replicate_file,
    null_file = null_file
  )
}

layer_results <- list()
layer_names <- names(layers)
for (i in seq_along(layer_names)) {
  nm <- layer_names[i]
  layer_results[[nm]] <- analyze_layer(nm, layers[[nm]], config, i)
}

summary_lines <- c(
  "Phase 04 Step 03 summary",
  paste("Timestamp:", timestamp()),
  paste("Replicate seeds:", paste(range(config$replicate_seeds), collapse = "-")),
  paste("Null replicates:", length(config$null_seeds)),
  paste("Swap iterations per null:", format(config$swap_iterations, scientific = FALSE)),
  ""
)

for (layer_name in names(layer_results)) {
  res <- layer_results[[layer_name]]
  qs <- res$replicate_stats
  ns <- res$null_summary
  summary_lines <- c(
    summary_lines,
    str_to_title(layer_name),
    sprintf("  Best Q: %.4f", res$Q_best),
    sprintf("  Replicate Q mean±sd: %.4f ± %.4f", qs$mean_Q, qs$sd_Q),
    sprintf("  Null mean±sd: %.4f ± %.4f", ns$Q_null_mean, ns$Q_null_sd),
    sprintf("  z-score: %.2f | p(one-tailed): %.3f", ns$z_score, ns$p_value),
    sprintf("  Modules saved: %s", res$assignments_path),
    sprintf("  Replicate Q file: %s", res$replicate_file),
    sprintf("  Null summary file: %s", res$null_file),
    ""
  )
}

summary_file <- "docs/phase_04/logs/step03_modularity_summary.txt"
writeLines(summary_lines, summary_file)
cat(paste(summary_lines, collapse = "\n"), "\n")
