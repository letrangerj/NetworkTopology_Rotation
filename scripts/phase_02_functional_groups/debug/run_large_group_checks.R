#!/usr/bin/env Rscript
# run_large_group_checks.R
# Debug utilities: run three checks for the large functional group
# 1) GCF prefix counts for the largest group
# 2) Signature and representative strain details for the largest group
# 3) Receptor group frequencies that map to the group's usable pyoverdines
#
# Usage (project root):
#   Rscript scripts/phase_02_functional_groups/debug/run_large_group_checks.R
#
# Outputs (in same debug folder):
#   - large_group_prefix_counts.csv
#   - large_group_signature_details.csv
#   - large_group_receptor_frequencies.csv
#   - large_group_checks_summary.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
})

cat("=== Running large-group checks ===\n")

# Paths (adjust if your repo layout differs)
fg_path <- "data/interim/functional_groups_conservative.rds"
mapping_path <- "data/interim/strain_functional_group_mapping_conservative.rds"
strain_reps_path <- "data/interim/strain_repertoires.rds"
recmap_path <- "data/interim/receptor_to_pyoverdine.rds"

# Output debug paths
out_dir <- "scripts/phase_02_functional_groups/debug"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_prefix_csv <- file.path(out_dir, "large_group_prefix_counts.csv")
out_signature_csv <- file.path(out_dir, "large_group_signature_details.csv")
out_rec_freq_csv <- file.path(out_dir, "large_group_receptor_frequencies.csv")
out_text <- file.path(out_dir, "large_group_checks_summary.txt")

# Load inputs with basic checks
if (!file.exists(fg_path)) stop("Missing file: ", fg_path)
if (!file.exists(mapping_path)) stop("Missing file: ", mapping_path)
if (!file.exists(strain_reps_path)) stop("Missing file: ", strain_reps_path)
if (!file.exists(recmap_path)) stop("Missing file: ", recmap_path)

functional_groups <- readRDS(fg_path)
strain_mapping <- readRDS(mapping_path)
strain_reps <- readRDS(strain_reps_path)
receptor_map <- readRDS(recmap_path)

cat("Loaded inputs:\n")
cat(" - functional_groups:", nrow(functional_groups), "groups\n")
cat(" - strain_mapping:", nrow(strain_mapping), "rows\n")
cat(" - strain_repertoires:", nrow(strain_reps), "strains\n")
cat(" - receptor_to_pyoverdine entries:", nrow(receptor_map), "\n\n")

# Identify the largest functional group
max_size <- max(functional_groups$strain_count, na.rm = TRUE)
largest_fg <- functional_groups %>% filter(strain_count == max_size)
if (nrow(largest_fg) > 1) {
  # If multiple groups tie for max, pick the first but warn
  warning("Multiple functional groups tie for max size; selecting the first one.")
  largest_fg <- largest_fg[1, , drop = FALSE]
}
group_id <- largest_fg$functional_group_id
group_size <- largest_fg$strain_count

cat(sprintf("Largest group: ID=%s, size=%d\n\n", group_id, group_size))

# Get all strains in that group
group_strains <- strain_mapping %>%
  filter(functional_group_id == group_id) %>%
  pull(strainName)

# 1) GCF prefix counts (first 10 chars)
prefixes <- substr(group_strains, 1, 10)
prefix_tbl <- sort(table(prefixes), decreasing = TRUE)
prefix_df <- data.frame(prefix = names(prefix_tbl), count = as.integer(prefix_tbl), stringsAsFactors = FALSE)

# Save prefix counts
write.csv(prefix_df, out_prefix_csv, row.names = FALSE)
cat("Wrote prefix counts to", out_prefix_csv, "\n")

# 2) Signature details and representative strains
# We want the validated production group list and usable_pyoverdines list for the group.
# functional_groups stores counts; the repertoires hold lists per strain. We'll pick a few representative strains
rep_sample <- head(group_strains, 10)  # up to 10 representative strains
rep_info <- strain_reps %>%
  filter(strainName %in% rep_sample) %>%
  transmute(
    strainName,
    production_groups = sapply(production_groups, function(x) if (length(x)==0) "" else paste(sort(x), collapse = ",")),
    validated_production_groups = sapply(validated_production_groups, function(x) if (length(x)==0) "" else paste(sort(x), collapse = ",")),
    receptor_groups = sapply(receptor_groups, function(x) if (length(x)==0) "" else paste(sort(x), collapse = ",")),
    usable_pyoverdines  = sapply(usable_pyoverdines, function(x) if (length(x)==0) "" else paste(sort(x), collapse = ","))
  )

# Also produce a single-row summary for the group using the first representative strain when list fields are needed
representative <- rep_info$strainName[1]
rep_row <- rep_info %>% filter(strainName == representative)

# Capture the group's production/usable values based on representative
group_validated_prod <- rep_row$validated_production_groups
group_usable_pyovs <- rep_row$usable_pyoverdines

# Build signature table (representatives)
signature_df <- rep_info
write.csv(signature_df, out_signature_csv, row.names = FALSE)
cat("Wrote signature details for representative strains to", out_signature_csv, "\n")

# 3) Receptor group frequencies that map to the group's usable pyoverdines
# Determine target pyoverdine ids
target_pyovs <- if (is.na(group_usable_pyovs) || group_usable_pyovs == "") {
  character(0)
} else {
  strsplit(group_usable_pyovs, ",", fixed = TRUE)[[1]] %>% trimws() %>% as.integer()
}

cat("Group usable pyoverdines (parsed):", paste(target_pyovs, collapse = ", "), "\n")

if (length(target_pyovs) == 0) {
  cat("No usable pyoverdines found for representative strain; skipping receptor-frequency analysis.\n")
  rec_freq_df <- data.frame()
} else {
  # Find receptor groups that map to these pyovs in receptor_map
  recs_for_pyov <- receptor_map %>%
    filter(pyoverdine_group %in% target_pyovs) %>%
    pull(receptor_group) %>%
    unique() %>%
    as.integer()

  # Count occurrences of these receptor groups across all strains
  # Gather all receptor group occurrences from strain_reps
  all_rec_list <- unlist(strain_reps$receptor_groups)
  # Some receptor lists may be empty; ensure integer
  all_rec_list <- as.integer(all_rec_list)
  all_rec_list <- all_rec_list[!is.na(all_rec_list)]

  # Count frequency for recs_for_pyov
  if (length(recs_for_pyov) == 0) {
    rec_freq_df <- data.frame(receptor_group = integer(0), count = integer(0))
  } else {
    rec_counts <- table(factor(all_rec_list, levels = sort(unique(all_rec_list))))
    rec_counts_subset <- rec_counts[as.character(recs_for_pyov)]
    # Convert to DF and fill NA with 0
    rec_freq_df <- data.frame(
      receptor_group = as.integer(names(rec_counts_subset)),
      count = as.integer(rec_counts_subset),
      stringsAsFactors = FALSE
    )
    rec_freq_df$count[is.na(rec_freq_df$count)] <- 0
    # Sort descending
    rec_freq_df <- rec_freq_df %>% arrange(desc(count))
  }

  write.csv(rec_freq_df, out_rec_freq_csv, row.names = FALSE)
  cat("Wrote receptor frequency table to", out_rec_freq_csv, "\n")
}

# Summary text report
report_lines <- c()
report_lines <- c(report_lines, "Large Group Checks Summary")
report_lines <- c(report_lines, paste0("Timestamp: ", Sys.time()))
report_lines <- c(report_lines, "")
report_lines <- c(report_lines, paste0("Largest functional group ID: ", group_id))
report_lines <- c(report_lines, paste0("Largest group size: ", group_size, " strains"))
report_lines <- c(report_lines, "")
report_lines <- c(report_lines, "Prefix distribution (top 20):")
top_prefix <- head(prefix_df, 20)
for (i in seq_len(nrow(top_prefix))) {
  report_lines <- c(report_lines, paste0("  - ", top_prefix$prefix[i], ": ", top_prefix$count[i]))
}
report_lines <- c(report_lines, "")

report_lines <- c(report_lines, "Representative strain signature details (up to 10 strains):")
for (i in seq_len(nrow(signature_df))) {
  r <- signature_df[i, ]
  report_lines <- c(report_lines, paste0("  - ", r$strainName,
                                         " | validated_production_groups: [", r$validated_production_groups, "]",
                                         " | usable_pyoverdines: [", r$usable_pyoverdines, "]",
                                         " | receptor_groups: [", r$receptor_groups, "]"))
}
report_lines <- c(report_lines, "")

if (length(target_pyovs) == 0) {
  report_lines <- c(report_lines, "No usable pyoverdines detected for the representative strain; receptor frequency analysis skipped.")
} else {
  report_lines <- c(report_lines, paste0("Target usable pyoverdines for group (representative): ", paste(target_pyovs, collapse = ", ")))
  report_lines <- c(report_lines, "Top receptor groups that map to these pyoverdines (by occurrences across all strains):")
  if (nrow(rec_freq_df) == 0) {
    report_lines <- c(report_lines, "  - (no matching receptor groups found)")
  } else {
    top_rec <- head(rec_freq_df, 50)
    for (i in seq_len(nrow(top_rec))) {
      report_lines <- c(report_lines, paste0("  - receptor_group ", top_rec$receptor_group[i], ": ", top_rec$count[i], " occurrences"))
    }
  }
}

# Save text summary
writeLines(report_lines, out_text)
cat("Wrote summary report to", out_text, "\n")

# Also print main summary to console
cat("\n=== Summary ===\n")
cat(paste(report_lines, collapse = "\n"))
cat("\n=== End ===\n")

# Exit normally
invisible(NULL)
