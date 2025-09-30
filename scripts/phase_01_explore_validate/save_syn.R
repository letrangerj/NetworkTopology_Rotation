library(rhdf5)
library(jsonlite)

syn_path <- '/home/wulong/RotationWorks/PhylogenicTree/data/interim/syn.h5'
syn_output_rds <- '/home/wulong/RotationWorks/PhylogenicTree/data/interim/rds/syn.rds'

# Function to convert uint8 matrix to vector of strings
uint8_matrix_to_vector <- function(uint8_mat) {
    char_mat <- intToUtf8(uint8_mat, multiple = TRUE)
    char_mat <- matrix(char_mat, nrow = nrow(uint8_mat), ncol = ncol(uint8_mat))
    apply(char_mat, 1, function(x) trimws(paste(x, collapse="")))
}

# Read all fields
syn <- list()

syn$clusterblast <- uint8_matrix_to_vector(h5read(syn_path, "/syn/clusterblast"))
syn$strainName <- uint8_matrix_to_vector(h5read(syn_path, "/syn/strainName"))
syn$regionName <- uint8_matrix_to_vector(h5read(syn_path, "/syn/regionName"))
syn$assemblyDefinition <- uint8_matrix_to_vector(h5read(syn_path, "/syn/assemblyDefinition"))
syn$location <- h5read(syn_path, "/syn/location")

# Read JSON fields
wholeCDS_json_uint8 <- h5read(syn_path, "/syn/wholeCDS_json")
wholeCDS_json_str <- rawToChar(as.raw(wholeCDS_json_uint8))
syn$wholeCDS <- fromJSON(wholeCDS_json_str, simplifyVector = FALSE)

biosynCDS_json_uint8 <- h5read(syn_path, "/syn/biosynCDS_json")
biosynCDS_json_str <- rawToChar(as.raw(biosynCDS_json_uint8))
syn$biosynCDS <- fromJSON(biosynCDS_json_str, simplifyVector = FALSE)

syn$group <- h5read(syn_path, "/syn/group")

syn$regionIdentifier <- uint8_matrix_to_vector(h5read(syn_path, "/syn/regionIdentifier"))

# Check the structure
str(syn, max.level = 2)
head(syn$strainName)
head(syn$regionIdentifier)

# Save as RDS
saveRDS(syn, file = syn_output_rds)
