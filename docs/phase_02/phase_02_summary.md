# Phase 2 Summary: Functional Group Classification (Conservative)

## Overview
Phase 2 implemented the strain-level collapse into functional groups using conservative, validated pyoverdine-receptor pairs (26 lock-key groups). We:
- Calculated per-strain production and utilization repertoires with validated-only logic for production (conservative).
- Collapsed strains into functional groups based on identical conservative production and utilization signatures.
- Classified functional groups into biologically meaningful strategies.

All outputs are reproducible, documented, and ready for Phase 3 (network construction). This phase closely follows the Gu et al. (2025) methodology, with explicit documentation of conservative vs. inclusive choices.

---

## Goals Accomplished
- Added validated-only production annotations and counts (conservative labeling).
- Translated receptors → usable pyoverdines via the validated pairing table.
- Collapsed 1,809 strains into 177 conservative functional groups (average ≈ 10.22 strains/group).
- Classified groups by conservative strategies (producer types vs. nonproducers).
- Investigated data features driving large group sizes and multi-producer behavior.

---

## Scripts and What They Do
- `scripts/phase_02_functional_groups/02_repertoire_calculation.R`
  - Loads `syn.rds`, `rec.rds`, and `pairing_result.rds`
  - Builds `receptor_to_pyoverdine.rds` (validated receptor→pyoverdine lookup; 1,149 mappings)
  - Computes per-strain:
    - `production_groups` (raw `syn$group`)
    - `validated_production_groups` and `n_validated_production_groups` (conservative)
    - `usable_pyoverdines` and `n_usable_pyoverdines` (via validated mapping only)
  - Saves: `data/interim/strain_repertoires.rds`
  - Reports: `docs/phase_02/repertoire_calculation_summary.txt`
  - Audit: `docs/phase_02/production_validation_audit.csv`

- `scripts/phase_02_functional_groups/03_functional_group_collapse.R`
  - Filters strains for conservative network: `is_validated_producer` OR `n_usable_pyoverdines > 0`
  - Builds conservative canonical signatures:
    - Production = `validated_production_groups`
    - Utilization = `usable_pyoverdines`
  - Collapses to functional groups by identical signatures
  - Saves:
    - `data/interim/functional_groups_conservative.rds`
    - `data/interim/strain_functional_group_mapping_conservative.rds`
  - Report: `docs/phase_02/functional_group_collapse_summary.txt`

- `scripts/phase_02_functional_groups/04_strategy_classification.R`
  - Classifies conservative functional groups using:
    - Multi-producer: >1 validated production group
    - Single-receptor producer: 1 production, 1 usable
    - Multi-receptor producer: 1 production, >1 usable
    - Nonproducer: 0 production, >0 usable
    - Producer-only (no utilization): >0 production, 0 usable
    - Unconnected: 0 production, 0 usable
  - Saves:
    - `data/interim/functional_groups_classified_conservative.rds`
    - `data/interim/strategy_classification_conservative.csv`
  - Report: `docs/phase_02/strategy_classification_summary.txt`

- Debug utilities (selected):
  - `scripts/phase_02_functional_groups/debug/count_conservative_producers.R`
    - Counts conservative producer categories (validated-only)
  - `scripts/phase_02_functional_groups/debug/inspect_large_group.R`
    - Inspects the largest functional group (signature, strains, receptor→pyoverdine drivers)
  - `scripts/phase_02_functional_groups/debug/run_large_group_checks.R`
    - Summarizes GCF prefix distribution and receptor frequencies mapping to the group’s usable pyoverdines

---

## Key Results

### Repertoire Calculation (Conservative)
- Total unique strains: 2,808
- Raw producers (any non-zero `syn$group`): 1,938
- Conservative producers (≥1 validated production group): 1,584
- Strains with usable pyoverdines (via validated pairs): 1,805
- Debug (conservative producers):
  - Category 1 (multi-producer conservative, ≥2 validated): 347 strains
  - Category 2 (raw multi-producer but ≤1 validated): 241 strains
  - Category 3 (only unvalidated syn groups): 354 strains

Outputs:
- `data/interim/strain_repertoires.rds`
- `docs/phase_02/repertoire_calculation_summary.txt`
- `docs/phase_02/production_validation_audit.csv`

### Functional Group Collapse (Conservative)
- Strains included in conservative collapse: 1,809 (validated producers OR consumers)
- Functional groups created: 177
- Average collapse ratio: 10.22 strains/group
- Group sizes: min = 1, max = 700
- Composition:
  - Groups with validated production: 134
  - Groups without production (pure consumers): 43
  - Groups with usable pyoverdines: 173
  - Groups without usable pyoverdines: 4

Outputs:
- `data/interim/functional_groups_conservative.rds`
- `data/interim/strain_functional_group_mapping_conservative.rds`
- `docs/phase_02/functional_group_collapse_summary.txt`

### Strategy Classification (Conservative)
- Strategy distribution (groups | strains):
  - Multi-producer: 33 groups | 347 strains
  - Multi-receptor producer: 95 groups | 1,221 strains
  - Single-receptor producer: 2 groups | 12 strains
  - Nonproducer: 43 groups | 225 strains
  - Producer-only (no utilization): 4 groups | 4 strains

Outputs:
- `data/interim/functional_groups_classified_conservative.rds`
- `data/interim/strategy_classification_conservative.csv`
- `docs/phase_02/strategy_classification_summary.txt`

---

## Two Key Questions in Phase 2

### 1) Multi-producer Problem
- Observation: Many strains had multiple raw production groups (588 strains), but only a subset were validated.
- Strategy used (conservative default):
  - Annotate `validated_production_groups` and use `n_validated_production_groups` for final labeling.
  - A strain is a conservative producer if `n_validated_production_groups > 0`.
  - A conservative multi-producer if `n_validated_production_groups > 1`.
- Alternative (not applied, candidate for later sensitivity analysis):
  - Enforce a single production group per strain (deterministic selection rule), or
  - Include unvalidated synthetases (all `syn$group`) in production for an inclusive analysis.
  - We deferred these to preserve a high-confidence conservative network now.

### 2) Strange Large Group (700 strains)
- Finding:
  - Largest group: Functional group ID 67 with signature:
    - Production: 1 validated group (22)
    - Utilization: pyoverdines 22, 23
  - Receptor mapping analysis:
    - A single validated receptor group (`5425`) is extremely common (816 occurrences overall) and maps to pyoverdines 22/23.
  - GCF prefix distribution indicates the group spans many projects (e.g., `GCF_001987*`, `GCF_001985*`, `GCF_001986*`, etc.), not a single batch artifact.
- Reason it exists:
  - Conservative validation compresses diverse receptor genotypes into a common functional outcome (usable pyoverdines 22/23).
  - A dominant validated synthetase group (22) and a very common validated receptor group (5425) together create a high-frequency functional signature.
  - The result is a biologically plausible dominant producer–consumer strategy under conservative assumptions, not a code error.

---

## What This Phase Achieved
- Delivered conservative, validated strain repertoires and functional groups suitable for robust network construction.
- Implemented clear, reproducible logic to separate validated vs. unvalidated production for downstream analyses.
- Produced classification of functional strategies, documenting the prevalence of multi-receptor producers and quantifying conservative multi-producers.
- Investigated and explained key data features (multi-producers and the large functional group), laying the groundwork for transparent interpretation of network metrics.

---

## Ready for Phase 3: Network Construction
With conservative functional groups and classifications completed, the next steps are to:
- Build the bipartite, directed network:
  - Nodes: functional groups and validated pyoverdine groups
  - Edges: production (group → pyoverdine) and utilization (pyoverdine → group)
- Validate graph properties and compute initial summary metrics.

All inputs for Phase 3 are in `data/interim/` and documented under `docs/phase_02/`.
