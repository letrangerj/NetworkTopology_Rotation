<file_path>
PhylogenicTree/docs/phase_02/phase_02_plan.md
</file_path>

<edit_description>
Create a summary document for Phase 2 plan
</edit_description>

# Phase 2 Plan: Functional Group Classification

## Overview
Phase 2 implements the strain-level collapse algorithm to create functional groups based on identical pyoverdine production and utilization repertoires, following Gu et al. (2025). This phase corrects the initial plan by properly translating receptor groups to usable pyoverdine types via the validated pairing data.

## Key Objectives
- Implement strain-level collapse based on identical production/utilization repertoires
- Classify functional groups into three strategic categories
- Create mapping tables and validate against biological expectations
- Document collapse criteria and statistics

## Implementation Steps

### 1. Data Loading and Strain Identification
- Load `syn.rds`, `rec.rds`, and `pairing_result.rds`
- Extract consistent strain identifiers from `syn$strainName` and `rec$foldername`
- Validate identifier consistency across datasets

### 2. Receptor-to-Pyoverdine Translation
- Convert `pairing_result.rds` to lookup table: receptor_group → pyoverdine_group
- Filter to only validated pairs (26 total)

### 3. Strain Repertoire Calculation
For each strain:
- **Production repertoire**: pyoverdine group from `syn$group` (0 = non-producer)
- **Receptor groups**: list of receptor groups from `rec$group`
- **Utilization repertoire**: translate receptors to usable pyoverdines via pairing lookup

### 4. Functional Group Collapse
- Group strains with identical production + utilization repertoires
- Create unique functional group IDs
- Generate strain-to-functional-group mapping

### 5. Strategy Classification
Classify each functional group:
- **Nonproducer**: production_group = 0, ≥1 usable pyoverdines
- **Single-receptor producer**: production_group ≠ 0, exactly 1 usable pyoverdine
- **Multi-receptor producer**: production_group ≠ 0, >1 usable pyoverdines

### 6. Validation and Documentation
- Cross-validate with pairing data
- Generate summary statistics
- Document edge cases and assumptions

## Script Organization
```
scripts/phase_02_functional_groups/
├── 01_strain_identification.R
├── 02_repertoire_calculation.R
├── 03_functional_group_collapse.R
├── 04_strategy_classification.R
└── 05_validation_summary.R
```

## Expected Outputs
- `data/interim/functional_groups.rds`: Functional group assignments
- `data/interim/strain_functional_group_mapping.rds`: Strain-level mapping
- `docs/phase_02/functional_group_methodology.md`: Detailed methodology
- `docs/phase_02/functional_group_summary.txt`: Statistics and validation

## Validation Considerations
- Ensure no strains have multiple production groups (per Gu et al.)
- Handle receptors not in validated pairs (exclude from utilization)
- Validate strain identifier consistency
- Document any data quality issues

This plan ensures biological accuracy and methodological rigor for the bipartite network construction in Phase 3.