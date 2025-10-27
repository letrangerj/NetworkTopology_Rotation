#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
})

# Inputs
nodes_pyo_csv <- "data/interim/nodes_pyoverdines_conservative.csv"
pairing_rds <- "data/interim/rds/pairing_result.rds"
syn_rds <- "data/interim/rds/syn.rds"
rec_rds <- "data/interim/rds/rec.rds"
rec_to_pyov <- "data/interim/receptor_to_pyoverdine.rds" # optional mapping (rec_group -> original_pyov_id)

# Outputs
out_dir <- "results/phase_04/syn_rec_ratio"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_map_summary <- file.path(out_dir, "pyov_group_mapping.csv")
out_map_syn <- file.path(out_dir, "pyov_to_syn_detail.csv")
out_map_rec <- file.path(out_dir, "pyov_to_rec_detail.csv")
out_log <- file.path(out_dir, "pyov_group_mapping_log.txt")

cat("Linking pyoverdines to original syn/rec RDS sources...\n")


pyo_nodes <- readr::read_csv(nodes_pyo_csv, show_col_types = FALSE) %>%
  transmute(
    pyo_id = node_id,
    pyov_label = label,
    original_pyov_id = original_pyov_id,
    fg_production_count = fg_production_count,
    strain_production_count = strain_production_count,
    fg_utilization_count = fg_utilization_count,
    strain_utilization_count = strain_utilization_count
  )


# Helper to safely readRDS with message
safe_readRDS <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  readRDS(path)
}

pairing <- safe_readRDS(pairing_rds)
syn_obj <- safe_readRDS(syn_rds)
rec_obj <- safe_readRDS(rec_rds)
rec_pyov_map <- safe_readRDS(rec_to_pyov) # may be NULL

# ---- Parse pairing_result ----
# Goal: we want a data.frame with columns:
#   pair_index, syn_groups (list<int>), rec_groups (list<int>)
#
# This accommodates common forms:
# 1) list(final_pairs = data.frame/list-columns with col1=syn_groups, col2=rec_groups)
# 2) final_pairs as a list of length K, each element length 2 (syn, rec)
# 3) direct data.frame/matrix with 2 columns
#
extract_pairs <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.null(x$final_pairs)) x <- x$final_pairs
  # case: list of pairs
  if (is.list(x) && is.null(dim(x))) {
    # each element is a pair: list(syn_groups, rec_groups) or named
    df <- tibble(
      pair_index = seq_along(x),
      syn_groups = map(x, ~ {
        el <- .
        if (is.list(el) && length(el) >= 1) {
          v <- el[[1]]
          as.integer(unlist(v))
        } else if (!is.null(el$syn) || !is.null(el$syn_groups)) {
          as.integer(unlist(el$syn %||% el$syn_groups))
        } else {
          integer()
        }
      }),
      rec_groups = map(x, ~ {
        el <- .
        if (is.list(el) && length(el) >= 2) {
          v <- el[[2]]
          as.integer(unlist(v))
        } else if (!is.null(el$rec) || !is.null(el$rec_groups)) {
          as.integer(unlist(el$rec %||% el$rec_groups))
        } else {
          integer()
        }
      })
    )
    return(df)
  }
  # case: data.frame/matrix with two list-cols
  if (is.data.frame(x) || is.matrix(x)) {
    xdf <- as.data.frame(x, stringsAsFactors = FALSE)
    nm <- names(xdf)
    if (is.null(nm) || length(nm) < 2) {
      names(xdf) <- c("X1", "X2")[seq_len(ncol(xdf))]
      nm <- names(xdf)
    }
    # normalize list columns to integer vectors
    as_int_list <- function(col) {
      if (is.list(col)) {
        map(col, ~ as.integer(unlist(.)))
      } else {
        list(as.integer(unlist(col)))
      }
    }
    df <- tibble(
      pair_index = seq_len(nrow(xdf)),
      syn_groups = as_int_list(xdf[[1]]),
      rec_groups = as_int_list(xdf[[2]])
    )
    return(df)
  }
  NULL
}

`%||%` <- function(a, b) if (is.null(a)) b else a

pairs_df <- extract_pairs(pairing)

if (is.null(pairs_df)) {
  stop("Unable to parse pairing_result.rds into syn/rec group pairs.")
}

# ---- Build mapping of pyoverdine group id (original_pyov_id) ----
# Strategy:
# - If receptor_to_pyoverdine.rds exists, use it to assign each rec_group -> original_pyov_id.
# - Then for each pair (syn_groups, rec_groups), assign pyov_id = modal original_pyov_id across its rec_groups.
# - Finally aggregate per pyov_id: union of syn_groups and rec_groups.
#
# This is robust when 'original_pyov_id' is not sequential and aligns to curated IDs (as seen in nodes CSV).
#
rec_to_pyov_tbl <- NULL
if (!is.null(rec_pyov_map)) {
  # Expect a data.frame with at least: rec_group, original_pyov_id (or pyov_group/pyov_id)
  rec_to_pyov_tbl <- tryCatch(
    {
      as_tibble(rec_pyov_map)
    },
    error = function(e) NULL
  )
  # Normalize column names
  if (!is.null(rec_to_pyov_tbl)) {
    cn <- names(rec_to_pyov_tbl)
    if (!"rec_group" %in% cn) {
      # try to guess: group, recGroup, etc.
      cand <- cn[grepl("rec", tolower(cn)) & grepl("group", tolower(cn))]
      if (length(cand) >= 1) rec_to_pyov_tbl <- rec_to_pyov_tbl %>% rename(rec_group = !!cand[1])
    }
    if (!"original_pyov_id" %in% cn) {
      cand <- cn[grepl("pyov", tolower(cn)) | grepl("pyo", tolower(cn))]
      if (length(cand) >= 1) rec_to_pyov_tbl <- rec_to_pyov_tbl %>% rename(original_pyov_id = !!cand[1])
    }
  }
}


# Assign pyov_id to each pair_index using rec groups (expanded to all PYOs in the pair)

pairs_long <- pairs_df %>%
  mutate(rec_group = rec_groups) %>%
  unnest(rec_group)

if (!is.null(rec_to_pyov_tbl) && "rec_group" %in% names(rec_to_pyov_tbl) && "original_pyov_id" %in% names(rec_to_pyov_tbl)) {
  pairs_long <- pairs_long %>%
    left_join(rec_to_pyov_tbl %>% distinct(rec_group, original_pyov_id), by = "rec_group")
} else {
  cat("NOTE: receptor_to_pyoverdine mapping not found; pair_index assignment may be incomplete.\n", file = out_log, append = TRUE)
}

# Compute per pair and original_pyov_id the support (number of rec_groups mapped)
pair_pyov_support <- pairs_long %>%
  group_by(pair_index, original_pyov_id) %>%
  summarise(
    n_rec_groups_support = n_distinct(rec_group[!is.na(rec_group)]),
    .groups = "drop"
  ) %>%
  filter(!is.na(original_pyov_id))

# Resolve conflicts: some original_pyov_id may appear in multiple pair_index.
# Choose the pair_index with max support; break ties by smallest pair_index.
original_to_pair <- pair_pyov_support %>%
  group_by(original_pyov_id) %>%
  arrange(desc(n_rec_groups_support), pair_index) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(original_pyov_id, pair_index)

# Log conflicts for transparency
conflicts <- pair_pyov_support %>%
  semi_join(
    pair_pyov_support %>% count(original_pyov_id) %>% filter(n > 1),
    by = "original_pyov_id"
  )

if (nrow(conflicts) > 0) {
  cat("WARNING: original_pyov_id mapped to multiple pair_index. Resolved by max support then smallest pair_index.\n", file = out_log, append = TRUE)
  suppressWarnings(
    write.table(conflicts, file = out_log, append = TRUE, sep = ",", row.names = FALSE, col.names = TRUE)
  )
}

# Map syn/rec groups per pair for later attachment
pairs_kv <- pairs_df %>%
  transmute(
    pair_index,
    syn_groups = syn_groups,
    rec_groups = rec_groups
  )

# Assign pair_index to all PYO nodes via original_pyov_id
pyov_map <- pyo_nodes %>%
  left_join(original_to_pair, by = "original_pyov_id")

# Unresolved PYOs (no rec->pyov map available)
unresolved <- pyov_map %>% filter(is.na(pair_index))

if (nrow(unresolved) > 0) {
  cat(sprintf("WARNING: %d pyoverdines had no assigned pair_index via receptor mapping.\n", nrow(unresolved)), file = out_log, append = TRUE)

  cat(paste("Unresolved original_pyov_id:", paste(head(unresolved$original_pyov_id, 20), collapse = ", "), if (nrow(unresolved) > 20) "..."), file = out_log, append = TRUE, sep = "\n")
}



# Build final aggregation per PYO: use the pair's syn/rec groups for every PYO in that pair

pyov_summary <- pyo_nodes %>%
  left_join(pyov_map %>% select(pyo_id, pair_index), by = "pyo_id") %>%
  left_join(pairs_kv, by = "pair_index") %>%
  mutate(
    n_syn_groups = map_int(syn_groups, ~ length(unique(.x %||% integer()))),
    n_rec_groups = map_int(rec_groups, ~ length(unique(.x %||% integer()))),
    syn_groups_csv = map_chr(syn_groups, ~ paste(sort(unique(.x %||% integer())), collapse = ";")),
    rec_groups_csv = map_chr(rec_groups, ~ paste(sort(unique(.x %||% integer())), collapse = ";"))
  ) %>%
  select(pyo_id, pyov_label, original_pyov_id, pair_index, n_syn_groups, n_rec_groups, syn_groups_csv, rec_groups_csv)


# ---- Attach syn/rec original metadata detail tables ----
# syn_obj and rec_obj are expected to be lists/dfs from .mat with at least a 'group' column.
# We will try to coerce and summarise:
#
# syn_detail: pyov x syn_group -> counts, example strains/regions
# rec_detail: pyov x rec_group -> counts, example rec names/tags
#
syn_tbl <- NULL
rec_tbl <- NULL

coerce_to_tibble <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.data.frame(x)) {
    return(as_tibble(x))
  }
  # Try to convert list of equal-length vectors to tibble
  if (is.list(x)) {
    lens <- lengths(x)
    if (length(unique(lens)) == 1) {
      return(as_tibble(x))
    }
  }
  NULL
}

syn_tbl <- coerce_to_tibble(syn_obj)
rec_tbl <- coerce_to_tibble(rec_obj)

# Normalise group columns
if (!is.null(syn_tbl)) {
  if (!"group" %in% names(syn_tbl)) {
    # try alternative names
    cand <- names(syn_tbl)[grepl("group", tolower(names(syn_tbl)))]
    if (length(cand) >= 1) syn_tbl <- syn_tbl %>% rename(group = !!cand[1])
  }
}
if (!is.null(rec_tbl)) {
  if (!"group" %in% names(rec_tbl)) {
    cand <- names(rec_tbl)[grepl("group", tolower(names(rec_tbl)))]
    if (length(cand) >= 1) rec_tbl <- rec_tbl %>% rename(group = !!cand[1])
  }
}

# Pick a few descriptive columns if present
pick_cols <- function(tbl, preferred) {
  cols <- intersect(preferred, names(tbl))
  unique(c("group", cols))
}

syn_cols <- if (!is.null(syn_tbl)) pick_cols(syn_tbl, c("strainName", "regionName", "regionIdentifier", "assemblyDefinition", "location")) else character()
rec_cols <- if (!is.null(rec_tbl)) pick_cols(rec_tbl, c("recname", "fragmentname", "foldername", "tag", "domloc")) else character()

# Expand pyov -> syn_group mapping
pyov_syn_detail <- pyov_summary %>%
  select(pyo_id, original_pyov_id, syn_groups = syn_groups_csv) %>%
  mutate(syn_group = strsplit(syn_groups, ";")) %>%
  unnest(syn_group, keep_empty = TRUE) %>%
  mutate(syn_group = suppressWarnings(as.integer(syn_group))) %>%
  filter(!is.na(syn_group))

if (!is.null(syn_tbl) && length(syn_cols) > 0) {
  pyov_syn_detail <- pyov_syn_detail %>%
    left_join(syn_tbl %>% select(all_of(syn_cols)), by = c("syn_group" = "group")) %>%
    group_by(pyo_id, original_pyov_id, syn_group) %>%
    summarise(
      n_entries = n(),
      ex_strain = first(na.omit(strainName %||% NA_character_)),
      ex_region = first(na.omit(regionName %||% NA_character_)),
      .groups = "drop"
    )
}

# Expand pyov -> rec_group mapping
pyov_rec_detail <- pyov_summary %>%
  select(pyo_id, original_pyov_id, rec_groups = rec_groups_csv) %>%
  mutate(rec_group = strsplit(rec_groups, ";")) %>%
  unnest(rec_group, keep_empty = TRUE) %>%
  mutate(rec_group = suppressWarnings(as.integer(rec_group))) %>%
  filter(!is.na(rec_group))

if (!is.null(rec_tbl) && length(rec_cols) > 0) {
  pyov_rec_detail <- pyov_rec_detail %>%
    left_join(rec_tbl %>% select(all_of(rec_cols)), by = c("rec_group" = "group")) %>%
    group_by(pyo_id, original_pyov_id, rec_group) %>%
    summarise(
      n_entries = n(),
      ex_recname = first(na.omit(recname %||% NA_character_)),
      ex_fragment = first(na.omit(fragmentname %||% NA_character_)),
      .groups = "drop"
    )
}


# Attach counts from nodes (if present) and compute strain-level counts from RDS as fallback
# 1) Bring in node-provided counts (may be NA if not present in nodes file)
pyo_counts_from_nodes <- pyo_nodes %>% select(
  pyo_id, fg_production_count, strain_production_count, fg_utilization_count, strain_utilization_count
)
pyov_summary <- pyov_summary %>%
  left_join(pyo_counts_from_nodes, by = "pyo_id")


# 2) Compute strain producer/utilizer counts from syn/rec RDS when nodes file lacks counts
syn_counts <- NULL
rec_counts <- NULL

# Build producer strain counts from syn_tbl (choose best available strain identifier column)
if (!is.null(syn_tbl)) {
  syn_id_col <- intersect(c("strainName", "assemblyDefinition", "regionIdentifier", "foldername"), names(syn_tbl))
  if (length(syn_id_col) >= 1) {
    syn_ids <- syn_tbl %>%
      select(group, !!all_of(syn_id_col[1])) %>%
      rename(strain_id = !!all_of(syn_id_col[1]))
    syn_counts <- pyov_summary %>%
      select(pyo_id, syn_groups_csv) %>%
      mutate(syn_group = strsplit(syn_groups_csv, ";")) %>%
      unnest(syn_group, keep_empty = TRUE) %>%
      mutate(syn_group = suppressWarnings(as.integer(syn_group))) %>%
      filter(!is.na(syn_group)) %>%
      left_join(syn_ids, by = c("syn_group" = "group")) %>%
      group_by(pyo_id) %>%
      summarise(strain_production_count_rds = n_distinct(strain_id), .groups = "drop")
  }
}

# Build utilizer strain counts from rec_tbl (choose best available strain identifier column)
if (!is.null(rec_tbl)) {
  rec_id_col <- intersect(c("strainName", "foldername", "fragmentname", "tag"), names(rec_tbl))
  if (length(rec_id_col) >= 1) {
    rec_ids <- rec_tbl %>%
      select(group, !!all_of(rec_id_col[1])) %>%
      rename(strain_id = !!all_of(rec_id_col[1]))
    rec_counts <- pyov_summary %>%
      select(pyo_id, rec_groups_csv) %>%
      mutate(rec_group = strsplit(rec_groups_csv, ";")) %>%
      unnest(rec_group, keep_empty = TRUE) %>%
      mutate(rec_group = suppressWarnings(as.integer(rec_group))) %>%
      filter(!is.na(rec_group)) %>%
      left_join(rec_ids, by = c("rec_group" = "group")) %>%
      group_by(pyo_id) %>%
      summarise(strain_utilization_count_rds = n_distinct(strain_id), .groups = "drop")
  }
}


# Ensure join targets are data frames (not NULL) to avoid auto_copy() errors
if (is.null(syn_counts)) {
  syn_counts <- data.frame(
    pyo_id = character(0),
    strain_production_count_rds = integer(0),
    stringsAsFactors = FALSE
  )
}
if (is.null(rec_counts)) {
  rec_counts <- data.frame(
    pyo_id = character(0),
    strain_utilization_count_rds = integer(0),
    stringsAsFactors = FALSE
  )
}

pyov_summary <- pyov_summary %>%
  left_join(syn_counts, by = "pyo_id") %>%
  left_join(rec_counts, by = "pyo_id") %>%
  mutate(
    strain_production_count = dplyr::coalesce(strain_production_count, strain_production_count_rds),
    strain_utilization_count = dplyr::coalesce(strain_utilization_count, strain_utilization_count_rds)
  ) %>%
  select(
    pyo_id, pyov_label, original_pyov_id, pair_index,
    n_syn_groups, n_rec_groups, syn_groups_csv, rec_groups_csv,
    fg_production_count, strain_production_count,
    fg_utilization_count, strain_utilization_count
  )


# Save outputs
readr::write_csv(pyov_summary, out_map_summary)
readr::write_csv(pyov_syn_detail, out_map_syn)
readr::write_csv(pyov_rec_detail, out_map_rec)


cat("Wrote:\n")
cat(" - ", out_map_summary, "\n", sep = "")
cat(" - ", out_map_syn, "\n", sep = "")
cat(" - ", out_map_rec, "\n", sep = "")
