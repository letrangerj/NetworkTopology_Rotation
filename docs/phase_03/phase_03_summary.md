# Phase 3 Summary: Dual Conservative Network Construction & Visualization

## Overview
Phase 3 implemented the construction of conservative, validated-only bipartite siderophore interaction networks and produced exploratory visualizations modeled on the Gu et al. (2025) Fig. 3A style. The implementation produced two network layers:

- Functional-group-level network (nodes = functional groups collapsed by identical repertoires)
- Strain-level network (nodes = conservative, validated strains)

The pipeline persisted canonical list-columns at the agent level, constructed contiguous pyoverdine indices, built oriented adjacency matrices (explicit orientation semantics), exported binary edge lists, computed degree sequences, ran consistency checks, and produced a publication-style bipartite network visualization.

This document summarizes what was done, scripts created, key results, diagnostics, deviations from the original plan, and recommended next steps.

---

## Goals of Phase 3 (as planned)
- Build bipartite directed networks with:
  - Agents (functional groups & strains)
  - Validated pyoverdine groups (lock–key)
  - Production edges (agent → pyov)
  - Utilization edges (pyov → agent)
- Persist contiguous pyoverdine index mapping for reproducible adjacency construction
- Keep edge lists binary (no weights)
- Provide validation artifacts: adjacency ↔ edge consistency, inactive pyov list, self-producer/consumer examples, component summaries
- Produce a visualization closely matching the style and encodings in Gu et al. Fig 3A
- Prepare degree sequences and diagnostics for Phase 4 topology analysis

---

## Scripts (Phase 3 pipeline)
The following scripts were created under `scripts/phase_03_network/`:

- `01_derive_agent_sets.R`  
  Derive and persist canonical agent tables:
  - `data/interim/nodes_functional_groups_conservative.rds/.csv`
  - `data/interim/nodes_strains_conservative.rds/.csv`  
  Performs non-circular validation (sampled checks).

- `02_pyov_index_mapping.R`  
  Create contiguous pyoverdine index mapping (matrix indices) and node dictionary:
  - `data/interim/nodes_pyoverdines_conservative.rds/.csv`  
  Adds `matrix_index`, `node_id` (`PYO_01`, `PYO_02`, ...) and active flags/counts where available.

- `03_build_adjacencies.R`  
  Construct oriented adjacency matrices and binary edge lists (FG & STR):
  - Production adjacency (Agents × Pyov): `data/interim/adj_production_*_agentsxpyov_conservative.rds`
  - Utilization adjacency (Pyov × Agents): `data/interim/adj_utilization_*_pyovxagents_conservative.rds`
  - Edge lists (RDS + CSV): `data/interim/edges_functional_groups_conservative.*`, `data/interim/edges_strains_conservative.*`
  - Degree sequences: `data/interim/degree_sequences_conservative.rds`
  Also produces diagnostic CSVs and component summaries.

- `04_visualize_network_*` (variants: `fixed`, `minimal`, `base`)  
  Produce the bipartite network visualization matching Fig 3A encodings, plus basic network metrics and image files under `figures/network/`.

All scripts include logging and save session info for provenance.

---

## Key results & artifacts produced

### Agent- and pyov-level summaries
- Functional groups: **177** (collapsed from conservative strains)
- Conservative strains included: **1,809**
- Functional groups with validated production: **134**
- Functional groups with usable pyov(s): **173**
- Unique (validated) synthetase/pyov identifiers expanded from validated pairs: **93** (see Decision section)
- Conservative validated pairs (pairing_result rows): **26** (source grouping)

Files:
- `data/interim/nodes_functional_groups_conservative.rds/.csv`
- `data/interim/nodes_strains_conservative.rds/.csv`
- `data/interim/nodes_pyoverdines_conservative.rds/.csv`

### Adjacency / Edge statistics (current implementation)
- FG network (177 agents × 93 pyovs)
  - Production edges (FG → PYO): **175**
  - Utilization edges (PYO → FG): **2,238**
  - Total FG-level unique edges: **2,413**
- Strain network (1,809 strains × 93 pyovs)
  - Production edges (STR → PYO): **1,948**
  - Utilization edges (PYO → STR): **14,198**
  - Total STR-level unique edges: **16,146**
- Pyov activity: **0 inactive pyovs** (all 93 have production and/or utilization)
- Weak component summary: **2** components (giant component + smaller piece)

Files:
- `data/interim/adj_production_FG_agentsxpyov_conservative.rds`
- `data/interim/adj_utilization_FG_pyovxagents_conservative.rds`
- `data/interim/adj_production_STR_agentsxpyov_conservative.rds`
- `data/interim/adj_utilization_STR_pyovxagents_conservative.rds`
- `data/interim/edges_functional_groups_conservative.rds/.csv`
- `data/interim/edges_strains_conservative.rds/.csv`
- `data/interim/degree_sequences_conservative.rds`
- Diagnostics:
  - `docs/phase_03/step_03_adjacencies/adjacency_counts_check.csv`
  - `docs/phase_03/step_03_adjacencies/inactive_pyovs.csv`
  - `docs/phase_03/step_03_adjacencies/self_producer_consumer_examples_fg.csv`
  - `docs/phase_03/step_03_adjacencies/fg_component_summary.csv`
  - `docs/phase_03/step_03_adjacencies/str_component_summary.csv`

### Visualization outputs (Step 4)
A bipartite network plot (Fig 3A style) was produced from the FG-level network:
- `figures/network/fg_pyov_network_panel3A.pdf`
- `figures/network/fg_pyov_network_panel3A.png`
- Visualization metrics: `docs/phase_03/visualization_metrics.csv`
- Visualization summary: `docs/phase_03/visualization_summary.txt`

The plot implements:
- Circles (functional groups) sized by `strain_count`
- Hexagon-like markers (pyov nodes) sized by usage counts (proxied by per-pyov counts)
- Circle color = strategy class (single-receptor / multi-receptor / nonproducer)
- Edge coloring and arrowheads approximated the FG-origin strategy (edges visualized with directional arrows)
- Force-directed layout to reveal hub structure; transparency and thin lines used to reduce clutter

---

## Diagnostics & QA results
- Non-circular validation: sampled checks on 10 FGs — **no inconsistencies found**
- Adjacency ↔ edge consistency checks: sums of adjacency matrices match exported edge counts (saved)
- Inactive pyovs: **none** (per FG-level checks)
- Self-producer/consumer cases: identified and saved (`self_producer_consumer_examples_fg.csv`)
- Edge deduplication: edges are binary and deduplicated before saving
- Provenance: session info and summaries saved in `docs/phase_03/step_01_agent_sets/session_info.txt`, `docs/phase_03/step_02_pyov_index/session_info_pyov_index.txt`, and `docs/phase_03/step_03_adjacencies/session_info_adj_build.txt`

---

## Decision point — 93-node vs 26-node pyov representation
During implementation we encountered an important representation choice:

- The original plan specified pyov node set = 26 validated lock–key groups (one node per pairing row).
- The pipeline as executed expanded the synthetase lists embedded in the 26 pairing rows and created **one pyov node per synthetase group id**, producing **93 pyov nodes**.
  - Rationale for 93-node choice: preserves within-pair synthetase diversity and is more granular.
  - Rationale for 26-node choice: matches a common literature presentation of "one lock–key group per pair" and simplifies comparison to Gu et al.

Status:
- Both representations are defensible biologically. The current artifacts reflect the **93-node** variant.
- A collapsed **26-node** variant (mapping every synthetase id to its `pair_id` → single pyov node) can be produced quickly by re-running the mapping + adjacency steps with a collapse-by-`pair_id` mapping.

Recommendation:
- For reproducibility and sensitivity analysis, produce both variants and compare key topology metrics (modularity, NODF nestedness, degree distributions). If results are similar, present the 26-node network for direct literature comparison and the 93-node as a finer-resolution supplement.

---

## Items still outstanding (Phase 3 plan checklist)
The majority of Phase 3 computations and validations were implemented; a few documentation and comparison artifacts remain to fully satisfy the plan:

- [ ] `docs/phase_03/gu_comparison_required.txt` — mandatory comparison to Gu et al. (2025) reference metrics (placeholders currently)
- [ ] `docs/phase_03/phase_03_runtime_estimates.txt` — runtime and memory estimates for Phase 4 modularity and null-model runs
- [ ] `docs/phase_03/phase_03_network_construction_log.txt` — narrative log of scripts run, parameters, anomalies (partial logs exist in step-specific session_info files)
- [ ] (Optional) Produce a 26-node collapsed pyov network for sensitivity / literature comparison

These are primarily documentation and decision artifacts and can be created without reprocessing the full datasets.

---

## Recommended next steps (prepare for Phase 4)
1. Decide pyov representation:
   - (A) Use the current **93-node** networks for Phase 4 analyses; OR
   - (B) Produce the **26-node** collapsed variant and run Phase 4 on both (recommended).
2. Produce `gu_comparison_required.txt` with reference metrics from Gu et al. (2025) and insert our summary metrics for comparison; flag substantial deviations for sensitivity checks.
3. Create `phase_03_runtime_estimates.txt` indicating:
   - matrix sizes, edge counts, and recommended memory/cores for `computeModules()` and `swap.web()` null-models
   - suggested initial randomization replicate counts (e.g., start with 100, scale to 1000 as feasible)
4. Prepare Phase 4 scripts (modularity, NODF nestedness, null models) using the `degree_sequences_conservative.rds` and adjacency matrices as inputs:
   - `scripts/phase_04_topology/01_compute_modules.R`
   - `scripts/phase_04_topology/02_nestedness_nodf.R`
   - `scripts/phase_04_topology/03_null_models_swap_web.R`
   - Include reproducible seeds and parallelization options.
5. (Planned) Metadata integration for habitat/pathogenicity subsetting in Phase 6: prepare a strain → metadata join key to allow stratified analyses.

---

## How to reproduce the Phase 3 run (quick)
From the project root (`PhylogenicTree`):

1. Derive agent sets (Step 1):
   - `Rscript scripts/phase_03_network/01_derive_agent_sets.R`

2. Create pyov mapping (Step 2):
   - `Rscript scripts/phase_03_network/02_pyov_index_mapping.R`

3. Build adjacencies and edge lists (Step 3):
   - `Rscript scripts/phase_03_network/03_build_adjacencies.R`

4. Visualize network (Step 4):
   - `Rscript scripts/phase_03_network/04_visualize_network_fixed.R`  
   (alternative lighter variants exist: `04_visualize_network_minimal.R`, `04_visualize_network_base.R`)

All scripts save provenance (session info) and diagnostic CSVs under `docs/phase_03/`.

---

## Files & locations (phase outputs)
- Node tables:
  - `data/interim/nodes_functional_groups_conservative.rds`
  - `data/interim/nodes_functional_groups_conservative.csv`
  - `data/interim/nodes_strains_conservative.rds`
  - `data/interim/nodes_strains_conservative.csv`
  - `data/interim/nodes_pyoverdines_conservative.rds`
  - `data/interim/nodes_pyoverdines_conservative.csv`
- Adjacencies:
  - `data/interim/adj_production_FG_agentsxpyov_conservative.rds`
  - `data/interim/adj_utilization_FG_pyovxagents_conservative.rds`
  - `data/interim/adj_production_STR_agentsxpyov_conservative.rds`
  - `data/interim/adj_utilization_STR_pyovxagents_conservative.rds`
- Edge lists:
  - `data/interim/edges_functional_groups_conservative.rds/.csv`
  - `data/interim/edges_strains_conservative.rds/.csv`
- Degree sequences:
  - `data/interim/degree_sequences_conservative.rds`
- Diagnostics & summaries:
  - `docs/phase_03/step_03_adjacencies/adjacency_counts_check.csv`
  - `docs/phase_03/step_03_adjacencies/inactive_pyovs.csv`
  - `docs/phase_03/step_03_adjacencies/self_producer_consumer_examples_fg.csv`
  - `docs/phase_03/step_03_adjacencies/fg_component_summary.csv`
  - `docs/phase_03/step_03_adjacencies/str_component_summary.csv`
  - `docs/phase_03/visualization_metrics.csv`
  - `docs/phase_03/visualization_summary.txt`
- Visual outputs:
  - `figures/network/fg_pyov_network_panel3A.pdf`
  - `figures/network/fg_pyov_network_panel3A.png`

---

## Final remarks
Phase 3 has been implemented end-to-end and the conservative bipartite networks (functional-group and strain-level) and visualizations were generated successfully. The main open item is the representational choice between the expanded synthetase-level pyov set (93 nodes) and the strict lock–key pair collapse (26 nodes). For scientific transparency and reproducibility, I recommend producing both variants (quick re-run of Steps 2–3 with a pair_id collapse) and then proceeding to Phase 4 network topology analyses with both (or selecting one after sensitivity checks). I can prepare the Phase 4 scripts and the 26-node collapse variant as the next actions.

If you want, I will:
- produce the 26-node collapsed variant now and rebuild adjacencies, then regenerate the visualization for direct comparison, and
- generate `gu_comparison_required.txt` and `phase_03_runtime_estimates.txt` to finish Phase 3 documentation.

Which of these would you like me to do next?