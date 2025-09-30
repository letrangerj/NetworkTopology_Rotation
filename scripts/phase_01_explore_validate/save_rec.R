library(rhdf5)

# load rec
rec_path <- '/home/wulong/RotationWorks/PhylogenicTree/data/interim/rec.h5'

# Function to convert uint8 matrix to vector of strings
uint8_matrix_to_vector <- function(uint8_mat) {
  char_mat <- intToUtf8(uint8_mat, multiple = TRUE)
  char_mat <- matrix(char_mat, nrow = nrow(uint8_mat), ncol = ncol(uint8_mat))
  apply(char_mat, 1, function(x) trimws(paste(x, collapse="")))
}

# Read all fields
rec <- list()

rec$count <- h5read(rec_path, "/rec/count")
rec$foldername <- uint8_matrix_to_vector(h5read(rec_path, "/rec/foldername"))
rec$fragmentname <- uint8_matrix_to_vector(h5read(rec_path, "/rec/fragmentname"))
rec$boarders <- h5read(rec_path, "/rec/boarders")
rec$seqs <- uint8_matrix_to_vector(h5read(rec_path, "/rec/seqs"))
rec$tag <- uint8_matrix_to_vector(h5read(rec_path, "/rec/tag"))
rec$recname <- uint8_matrix_to_vector(h5read(rec_path, "/rec/recname"))
rec$domloc <- h5read(rec_path, "/rec/domloc")
rec$domseq <- uint8_matrix_to_vector(h5read(rec_path, "/rec/domseq"))
rec$group <- h5read(rec_path, "/rec/group")

# Check the structure
str(rec)
head(rec$foldername)
head(rec$recname)

# Save as RDS
saveRDS(rec, file = '/home/wulong/RotationWorks/PhylogenicTree/data/interim/rec.rds')
