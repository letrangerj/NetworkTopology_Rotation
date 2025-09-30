# inspect_rds.R
# Script to load and inspect all RDS files in data/interim/rds/
# Prints str() summary and key head() checks for each object
# Run from project root: Rscript scripts/phase_01_explore_validate/inspect_rds.R

root_path <- '/home/wulong/RotationWorks/PhylogenicTree'
rds_dir <- file.path(root_path, 'data', 'interim', 'rds')

# Get list of RDS files
files <- list.files(rds_dir, pattern='\\.rds$', full.names=TRUE)

if (length(files) == 0) {
  cat("No RDS files found in", rds_dir, "\n")
  quit(status=1)
}

for (f in files) {
  cat('=== Loading:', basename(f), '===\n')

  # Load the RDS object
  obj <- tryCatch({
    readRDS(f)
  }, error = function(e) {
    cat("Error loading", basename(f), ":", e$message, "\n\n")
    next
  })

  # Print str() summary
  cat('str() summary:\n')
  print(str(obj, max.level=1))

  cat('\nKey head() checks:\n')

  # Check for common fields and print head()
  if ('strainName' %in% names(obj)) {
    cat('head(strainName): ', paste(head(obj$strainName), collapse=', '), '\n')
  }
  if ('recname' %in% names(obj)) {
    cat('head(recname): ', paste(head(obj$recname), collapse=', '), '\n')
  }
  if ('group' %in% names(obj)) {
    cat('head(group): ', paste(head(obj$group), collapse=', '), '\n')
  }
  if ('dataname' %in% names(obj)) {
    cat('head(dataname): ', paste(head(obj$dataname), collapse=', '), '\n')
  }
  if ('distm' %in% names(obj)) {
    cat('dim(distm): ', paste(dim(obj$distm), collapse=' x '), '\n')
  }
  if ('final.pairs' %in% names(obj)) {
    cat('dim(final.pairs): ', paste(dim(obj$final.pairs), collapse=' x '), '\n')
  }

  cat('\n')
}

cat("Inspection complete.\n")
