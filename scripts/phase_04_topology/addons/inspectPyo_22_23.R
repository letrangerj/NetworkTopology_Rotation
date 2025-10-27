suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(purrr); library(stringr)
})

syn_rds <- "data/interim/rds/syn.rds"
rec_rds <- "data/interim/rds/rec.rds"

safe_read <- function(p) tryCatch(readRDS(p), error = function(e) NULL)

syn_obj <- safe_read(syn_rds)
rec_obj <- safe_read(rec_rds)

# Try multiple coercion strategies
to_tbl <- function(x) {
  if (is.null(x)) return(NULL)
  # 1) Already a data.frame or tibble
  if (is.data.frame(x)) return(as_tibble(x))
  # 2) List of equal-length vectors (field = vector)
  if (is.list(x)) {
    lens <- lengths(x)
    if (length(lens) > 0 && length(unique(lens)) == 1 && is.atomic(x[[1]])) {
      # Every element is a vector of same length: bind into a tibble
      out <- tryCatch(as_tibble(x), error = function(e) NULL)
      if (!is.null(out)) return(out)
    }
    # 3) List of records (each element is a list with the same names)
    if (all(map_lgl(x, ~ is.list(.x)))) {
      nm_list <- map(x, names)
      # if many are NULL, abort
      if (all(map_lgl(nm_list, ~ !is.null(.x)))) {
        common <- Reduce(intersect, nm_list)
        if (length(common) > 0) {
          # Extract common fields across all rows
          rows <- map(x, ~ .x[common])
          # Some fields may still be nested; simplify atomics when possible
          rows <- map(rows, ~ setNames(.x, common))
          out <- tryCatch(bind_rows(rows), error = function(e) NULL)
          if (!is.null(out)) return(as_tibble(out))
        }
      }
    }
  }
  NULL
}

syn_tbl <- to_tbl(syn_obj)
rec_tbl <- to_tbl(rec_obj)

# Normalize group column name if needed
normalize_group <- function(tbl) {
  if (is.null(tbl)) return(NULL)
  if ("group" %in% names(tbl)) return(tbl)
  gcol <- names(tbl)[grepl("\\bgroup\\b", tolower(names(tbl)))]
  if (length(gcol) >= 1) {
    tbl %>% rename(group = !!gcol[1])
  } else {
    tbl
  }
}

syn_tbl <- normalize_group(syn_tbl)
rec_tbl <- normalize_group(rec_tbl)

syn_groups <- c(22L, 23L)
rec_groups <- 5421:5428

cat("=== SYNTHETASE SIDE (groups 22, 23) ===\n")
if (!is.null(syn_tbl)) {
  cand <- tryCatch(syn_tbl %>% filter(group %in% syn_groups), error = function(e) NULL)
  if (!is.null(cand) && nrow(cand) > 0) {
    # Grab likely informative columns if present
    cols <- intersect(c("clusterblast","regionName","regionIdentifier","assemblyDefinition","location"), names(cand))
    # Make a small human-readable digest
    if ("clusterblast" %in% cols) {
      cb <- cand %>%
        mutate(
          clusterblast_txt = if (is.list(clusterblast))
            map_chr(clusterblast, ~ paste(as.character(.x)[1:min(3L, length(.x))], collapse=" | "))
          else as.character(clusterblast)
        ) %>%
        transmute(group, clusterblast_snippet = substr(clusterblast_txt, 1, 200)) %>%
        distinct(group, clusterblast_snippet)
      cat("Top clusterblast snippets:\n")
      print(cb, n = Inf)
    }
    if ("regionName" %in% cols) {
      cat("Example regionName by group:\n")
      print(cand %>% group_by(group) %>% summarize(ex_region = first(na.omit(regionName))), n = Inf)
    }
    if ("regionIdentifier" %in% cols) {
      cat("Example regionIdentifier by group:\n")
      print(cand %>% group_by(group) %>% summarize(ex_id = first(na.omit(regionIdentifier))), n = Inf)
    }
  } else {
    cat("Could not filter syn_tbl by groups 22/23 (structure may be nested more deeply).\n")
  }
} else {
  cat("syn.rds could not be flattened into a table; consider inspecting str(readRDS('.../syn.rds')).\n")
}

cat("\n=== RECEPTOR SIDE (groups 5421–5428) ===\n")
if (!is.null(rec_tbl)) {
  cand <- tryCatch(rec_tbl %>% filter(group %in% rec_groups), error = function(e) NULL)
  if (!is.null(cand) && nrow(cand) > 0) {
    if ("recname" %in% names(cand)) {
      cat("Common receptor names (recname):\n")
      print(cand %>% count(recname, sort = TRUE) %>% head(12))
    } else {
      cat("recname missing; showing group counts by foldername/fragmentname if available:\n")
      print(cand %>% count(group, foldername, fragmentname, sort = TRUE) %>% head(12))
    }
  } else {
    cat("Could not filter rec_tbl by groups 5421–5428 (structure may be nested more deeply).\n")
  }
} else {
  cat("rec.rds could not be flattened into a table; consider inspecting str(readRDS('.../rec.rds')).\n")
}
