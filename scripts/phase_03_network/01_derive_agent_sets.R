
# Step 1: Derive and Persist Validated Agent Sets

# Purpose: Build validated agent sets (functional groups and strains) with persistent repertoire data

# for network construction in Phase 3



# Load required libraries

library(dplyr)

library(tidyr)

library(purrr)

# Set random seed for reproducibility

set.seed(42)



# Define paths

interim_dir <- "data/interim"

output_dir <- "data/interim"

script_id <- "scripts/phase_03_network/01_derive_agent_sets.R"



# Ensure docs dir for this phase exists
dir.create("docs/phase_03", recursive = TRUE, showWarnings = FALSE)

# Load required datasets

cat("Loading Phase 2 datasets...\n")



# Load functional groups with their classifications

fg_conservative <- readRDS(file.path(interim_dir, "functional_groups_conservative.rds"))

fg_classified <- readRDS(file.path(interim_dir, "functional_groups_classified_conservative.rds"))



# Load strain-to-functional-group mapping

strain_to_fg <- readRDS(file.path(interim_dir, "strain_functional_group_mapping_conservative.rds"))



# Load strain repertoires with validated production/usable pyoverdines

strain_rep <- readRDS(file.path(interim_dir, "strain_repertoires.rds"))



# Load receptor-to-pyoverdine mapping

receptor_to_pyov <- readRDS(file.path(interim_dir, "receptor_to_pyoverdine.rds"))



cat("Data loading complete.\n")

cat("Strain repertoires:", nrow(strain_rep), "strains\n")

cat("Functional groups:", nrow(fg_conservative), "groups\n")



# Helper to normalize list vectors to comparable strings
normalize_set <- function(x, empty_token = "∅") {
  x <- as.integer(x)
  if (length(x) == 0 || all(is.na(x))) return(empty_token)
  paste(sort(unique(x)), collapse = ",")
}

# Step 1.1: Build Functional Group Agent Sets

cat("\nBuilding functional group agent sets...\n")



# For each functional group, collect member strains and validate repertoires

fg_sets <- strain_to_fg %>%

  group_by(functional_group_id) %>%

  summarise(

    member_strains = list(as.character(strainName)),

    strain_count = n(),

    .groups = "drop"

  ) %>%

  mutate(

    # Collect per-strain validated production and usable sets for each FG
    strain_production_sets = map(member_strains, function(strains) {

      strain_rep %>%
        filter(strainName %in% strains) %>%
        select(strainName, validated_production_groups)
    }),

    strain_usable_sets = map(member_strains, function(strains) {

      strain_rep %>%
        filter(strainName %in% strains) %>%
        select(strainName, usable_pyoverdines)
    }),

    # Normalized signatures per FG

    prod_signatures = map(strain_production_sets, ~ vapply(.x$validated_production_groups, normalize_set, character(1))),
    usable_signatures = map(strain_usable_sets, ~ vapply(.x$usable_pyoverdines, normalize_set, character(1))),
    n_unique_prod_signatures = map_int(prod_signatures, ~ length(unique(.x))),
    n_unique_usable_signatures = map_int(usable_signatures, ~ length(unique(.x))),
    production_consistent = n_unique_prod_signatures == 1,
    utilization_consistent = n_unique_usable_signatures == 1,
    # Canonical FG-level sets (use the first member's set; groups should be consistent by construction)
    validated_production_set = map(strain_production_sets, ~ if (nrow(.x) == 0) integer(0) else as.integer(sort(unique(.x$validated_production_groups[[1]])))),

    usable_pyoverdine_set   = map(strain_usable_sets, ~ if (nrow(.x) == 0) integer(0) else as.integer(sort(unique(.x$usable_pyoverdines[[1]]))))
  )

# Validation: Check that all member strains in each group have identical signatures

cat("Validating functional group consistency...\n")



validation_results <- fg_sets %>%
  transmute(
    functional_group_id,
    strain_count,
    production_consistent,
    utilization_consistent,
    n_unique_production = n_unique_prod_signatures,
    n_unique_usable = n_unique_usable_signatures
  )

# Report validation results
n_issues <- validation_results %>% filter(!production_consistent | !utilization_consistent) %>% nrow()

cat("Validation complete. Issues found:", n_issues, "\n")



if (n_issues > 0) {

  cat("WARNING: Inconsistent functional groups detected!\n")

  write.csv(validation_results %>% filter(!production_consistent | !utilization_consistent),

            file.path(output_dir, "functional_group_validation_issues.csv"), row.names = FALSE)

}



# Create final functional group agent table

fg_agents <- fg_sets %>%

  select(functional_group_id, strain_count, member_strains,

         validated_production_set, usable_pyoverdine_set) %>%
  mutate(
    # Assign stable IDs
    agent_id = sprintf("FG_%03d", row_number()),

    # Calculate counts

    n_production_groups = map_int(validated_production_set, length),

    n_usable_pyoverdines = map_int(usable_pyoverdine_set, length),

    # Add strategy from classification

    strategy = fg_classified$strategy[match(functional_group_id, fg_classified$functional_group_id)],

    # Add metadata
    has_production = n_production_groups > 0,
    has_utilization = n_usable_pyoverdines > 0
  )

# Step 1.2: Build Strain Agent Sets

cat("\nBuilding strain agent sets...\n")



# Get conservative strains (validated producers OR consumers)

strain_agents <- strain_rep %>%

  mutate(

    is_conservative = (n_validated_production_groups > 0) | (n_usable_pyoverdines > 0)

  ) %>%

  filter(is_conservative) %>%

  mutate(

    # Assign stable IDs

    agent_id = sprintf("STR_%04d", row_number()),

    # Keep essential fields

    strainName = as.character(strainName),

    # Ensure sets are integer vectors

    validated_production_set = map(validated_production_groups, ~ as.integer(sort(unique(.x)))),

    usable_pyoverdine_set = map(usable_pyoverdines, ~ as.integer(sort(unique(.x)))),

    has_production = n_validated_production_groups > 0,
    has_utilization = n_usable_pyoverdines > 0
  ) %>%

  select(

    agent_id, strainName,
    validated_production_set, usable_pyoverdine_set,

    n_validated_production_groups, n_usable_pyoverdines,

    has_production, has_utilization
  )

cat("Strain agents created:", nrow(strain_agents), "strains\n")



# Step 1.3: Non-circular validation sampling

cat("\nPerforming non-circular validation...\n")



validation_sample <- sample(nrow(fg_agents), min(10, nrow(fg_agents)), replace = FALSE)



for (fg_idx in validation_sample) {

  fg_id <- fg_agents$functional_group_id[fg_idx]

  fg_row <- fg_agents %>% filter(functional_group_id == fg_id)
  member_strains <- fg_row$member_strains[[1]]
  expected_prod <- sort(as.integer(fg_row$validated_production_set[[1]]))

  expected_usable <- sort(as.integer(fg_row$usable_pyoverdine_set[[1]]))



  # Recompute from raw strain data

  recomputed <- strain_rep %>% filter(strainName %in% member_strains)
  prod_matches <- all(map_lgl(recomputed$validated_production_groups, ~ identical(sort(as.integer(.x)), expected_prod)))
  usable_matches <- all(map_lgl(recomputed$usable_pyoverdines, ~ identical(sort(as.integer(.x)), expected_usable)))

  if (!prod_matches || !usable_matches) {
    cat("MISMATCH in FG", fg_id, "- Prod all-match:", prod_matches, "Usable all-match:", usable_matches, "\n")

  }

}



cat("Non-circular validation complete.\n")



# Step 1.4: Extract Pyoverdine Set Summary

cat("\nExtracting pyoverdine set summary...\n")



all_pyoverdines <- sort(unique(unlist(strain_rep$usable_pyoverdines)))

all_production_groups <- sort(unique(unlist(strain_rep$validated_production_groups)))



cat("Total unique pyoverdines in dataset:", length(all_pyoverdines), "\n")

cat("Total production groups:", length(all_production_groups), "\n")



# Step 1.5: Save all outputs

cat("\nSaving outputs...\n")



# Save functional group agents

saveRDS(fg_agents, file.path(output_dir, "nodes_functional_groups_conservative.rds"))

write.csv(
  fg_agents %>%
    select(agent_id, functional_group_id, strain_count, strategy,

           has_production, has_utilization, n_production_groups, n_usable_pyoverdines),

  file.path(output_dir, "nodes_functional_groups_conservative.csv"),
  row.names = FALSE
)



# Save strain agents

saveRDS(strain_agents, file.path(output_dir, "nodes_strains_conservative.rds"))

write.csv(
  strain_agents %>%
    select(agent_id, strainName, n_validated_production_groups, n_usable_pyoverdines,

           has_production, has_utilization),

  file.path(output_dir, "nodes_strains_conservative.csv"),
  row.names = FALSE
)



# Save validation results

saveRDS(validation_results, file.path(output_dir, "functional_group_validation_results.rds"))


# Create and save session info

session_info <- capture.output(sessionInfo())

writeLines(session_info, "docs/phase_03/session_info.txt")



# Create summary report

summary_report <- sprintf("

## Phase 3 Step 1: Agent Sets Summary



### Functional Groups

- Total functional groups: %d

- Average strain count per group: %.2f

- Groups with production: %d (%.1f%%)

- Groups with utilization: %d (%.1f%%)



### Strains

- Total conservative strains: %d

- Strains with production: %d (%.1f%%)

- Strains with utilization: %d (%.1f%%)



### Validation

- Consistency issues found: %d

- Non-circular validation: Performed on %d random groups



### Dataset Statistics

- Unique pyoverdines in dataset: %d

- Unique production groups: %d

- Conservative production groups (validated pairs): 26



Files created:

- nodes_functional_groups_conservative.rds/csv

- nodes_strains_conservative.rds/csv

- functional_group_validation_results.rds
",
nrow(fg_agents), mean(fg_agents$strain_count),

sum(fg_agents$has_production), sum(fg_agents$has_production)/nrow(fg_agents)*100,

sum(fg_agents$has_utilization), sum(fg_agents$has_utilization)/nrow(fg_agents)*100,

nrow(strain_agents),

sum(strain_agents$has_production), sum(strain_agents$has_production)/nrow(strain_agents)*100,

sum(strain_agents$has_utilization), sum(strain_agents$has_utilization)/nrow(strain_agents)*100,

n_issues, length(validation_sample),

length(all_pyoverdines), length(all_production_groups)

)


writeLines(summary_report, "docs/phase_03/agent_sets_summary.txt")



cat("\nPhase 3 Step 1 complete!\n")

cat("Files saved to:", output_dir, "\n")
