# scripts/phase_01_explore_validate/convert_core_files_to_rds.R
# Simplified version: Convert core files to RDS format
# Only processes: pairing_result_v7.mat, Syn_PWA_inGroup.h5, Rec_MSA_inGroup.h5
# Run from project root: Rscript scripts/phase_01_explore_validate/convert_core_files_to_rds.R

library(rhdf5)     # For HDF5 files
library(R.matlab)  # For MATLAB v7 files

# Helper to convert uint8 matrix from HDF5 to string vector
convert_h5_uint8_to_strings <- function(raw_data) {
    if (is.matrix(raw_data)) {
        char_matrix <- apply(raw_data, c(1,2), function(x) rawToChar(as.raw(x)))
        string_vector <- apply(char_matrix, 1, paste, collapse = "")
        return(trimws(string_vector))  # Remove trailing spaces
    } else {
        return(raw_data)
    }
}

# Set project paths
project_root <- "/home/wulong/RotationWorks/PhylogenicTree"
interim_dir <- file.path(project_root, "data", "interim")
rds_dir <- file.path(interim_dir, "rds")
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# Function to load MATLAB v7 pairing result file
load_pairing_result <- function() {
    v7_path <- file.path(interim_dir, "pairing_result_v7.mat")
    if (!file.exists(v7_path)) {
        stop("MATLAB v7 file not found: ", v7_path)
    }

    cat("Loading pairing_result_v7.mat...\n")

    # Load MATLAB file
    data <- tryCatch({
        readMat(v7_path)
    }, error = function(e) {
        stop("Error reading pairing_result_v7.mat: ", e$message)
    })

    if (!is.null(data$final.pairs)) {
        cat("Successfully loaded final.pairs with", nrow(data$final.pairs), "pairs\n")
    }

    return(data)
}

# Function to load clustering HDF5 files (Syn_PWA_inGroup.h5 and Rec_MSA_inGroup.h5)
load_clustering_file <- function(h5_file) {
    h5_path <- file.path(interim_dir, h5_file)
    if (!file.exists(h5_path)) {
        stop("HDF5 file not found: ", h5_path)
    }

    cat("Loading", h5_file, "...\n")

    data <- list()

    tryCatch({
        # Read distance matrix
        data$distm <- h5read(h5_path, "distm")
        cat("  distm dimensions:", paste(dim(data$distm), collapse = "x"), "\n")

        # Read dataname and convert strings
        raw_dataname <- h5read(h5_path, "dataname")
        data$dataname <- convert_h5_uint8_to_strings(raw_dataname)
        cat("  dataname length:", length(data$dataname), "\n")

    }, error = function(e) {
        stop("Error reading ", h5_file, ": ", e$message)
    })

    return(data)
}

# Load the three core files
cat("=== Loading Core Files ===\n")

# 1. Load pairing result
pairing_data <- load_pairing_result()

# 2. Load syn clustering data
syn_cluster_data <- load_clustering_file("Syn_PWA_inGroup.h5")

# 3. Load rec clustering data
rec_cluster_data <- load_clustering_file("Rec_MSA_inGroup.h5")

# Save as RDS files
cat("\n=== Saving as RDS Files ===\n")

# Save pairing result
if (!is.null(pairing_data)) {
    rds_path <- file.path(rds_dir, "pairing_result_v7.rds")
    saveRDS(pairing_data, rds_path)
    cat("Saved pairing_result_v7.rds to", rds_path, "\n")
}

# Save syn clustering data
if (!is.null(syn_cluster_data)) {
    rds_path <- file.path(rds_dir, "Syn_PWA_inGroup.rds")
    saveRDS(syn_cluster_data, rds_path)
    cat("Saved Syn_PWA_inGroup.rds to", rds_path, "\n")
}

# Save rec clustering data
if (!is.null(rec_cluster_data)) {
    rds_path <- file.path(rds_dir, "Rec_MSA_inGroup.rds")
    saveRDS(rec_cluster_data, rds_path)
    cat("Saved Rec_MSA_inGroup.rds to", rds_path, "\n")
}

# Close HDF5 handles
H5close()

cat("\nConversion complete! RDS files saved to:", rds_dir, "\n")
cat("Files converted:\n")
cat("- pairing_result_v7.rds\n")
cat("- Syn_PWA_inGroup.rds\n")
cat("- Rec_MSA_inGroup.rds\n")
