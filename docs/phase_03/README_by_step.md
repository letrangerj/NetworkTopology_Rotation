# Phase 03 — README by Step

This document summarizes Phase 3 (Network Construction) artifacts, explains the key pyoverdine mapping decision (93 vs 26), and organizes all produced files by the step that created them. Use this README to locate outputs, understand provenance, and choose the appropriate pyoverdine-node granularity for downstream analysis.

---

## High-level summary

- Phase 3 constructs two conservative bipartite directed networks:
  - Functional-group-level network: functional groups (agents) ↔ pyoverdines
  - Strain-level network: conservative strains (agents) ↔ pyoverdines
- The project uses validated pairing data from Phase 2 as the biological basis for edges.
- A mapping decision controls pyoverdine node granularity:
  - "93-node" mapping (current run): one pyoverdine node per validated synthetase group ID (expanded from pairing rows). This preserves synthetase-level detail.
  - "26-node" collapse (optional): one pyoverdine node per lock–key row (pair_id). This produces a single node per curated lock–key pair and matches the 1:1 lock–key representation used in some literature.

Recommendation: keep both variants for sensitivity analysis. Use the 26-node collapsed mapping for direct comparisons to publications that report one node per validated pair, and use the 93-node mapping to preserve synthetase-level diversity.

---

## Mapping decision: 93 vs 26 (detailed)

- Source structure:
  - The validated pairing table is organized as ~26 rows. Each row represents a validated lock–key relationship but may contain multiple synthetase (pyoverdine) group IDs in the synthetase column.
  - Expanding those rows yields a larger set of synthetase IDs (93 in the current run).
- 93-node variant:
  - Each validated synthetase ID becomes a separate pyoverdine node.
  - Pros: maximal retention of synthetase diversity, finer-grain topology.
  - Cons: not one-to-one with pairing rows; may complicate direct comparisons with literature using a single node per pair.
- 26-node variant (collapse-by-pair_id):
  - All synthetase IDs from one pairing row are collapsed to the same pyoverdine node (the `pair_id` acts as the pyov node id).
  - Pros: one node per validated lock–key pair, simpler comparisons with publications and conceptual lock–key analyses.
  - Cons: within-pair diversity is hidden (synthetase variants represented as a single node).

When reporting results, always document which pyov mapping variant was used and include a sensitivity check between both variants when possible.

---

## Files produced by step (organized)

Below the files are grouped by the Phase 3 step that produced them. Paths are relative to the project root.

### Step 1 — Agent sets (`scripts/phase_03_network/01_derive_agent_sets.R`)
Purpose: derive canonical functional-group and strain agent tables with per-agent validated production and usable pyoverdine sets.

- `data/interim/nodes_functional_groups_conservative.rds`  
  (tabular object; includes list-columns `validated_production_set` and `usable_pyoverdine_set`, `strain_count`, `strategy`)
- `data/interim/nodes_functional_groups_conservative.csv`  
  (flattened CSV summary of FG node table)
- `data/interim/nodes_strains_conservative.rds`  
  (strain node table with `validated_production_set`, `usable_pyoverdine_set`, counts, and `agent_id` like `STR_0001`)
- `data/interim/nodes_strains_conservative.csv`  
  (flattened CSV summary of strain nodes)
- `data/interim/functional_group_validation_results.rds`  
  (detailed validation diagnostics; empty/benign if no issues)
- `docs/phase_03/step_01_agent_sets/agent_sets_summary.txt`  
  (human-readable summary and provenance)
- `docs/phase_03/step_01_agent_sets/session_info.txt`  
  (R session info at runtime)

Notes:
- Step 1 performed non-circular validation by sampling groups and verifying that members share canonical signatures.

### Step 2 — Pyoverdine index mapping (`scripts/phase_03_network/02_pyov_index_mapping.R`)
Purpose: create a contiguous pyoverdine mapping (matrix indices) used to order adjacency matrices deterministically.

- `data/interim/nodes_pyoverdines_conservative.rds`  
  (pyov node dictionary with `original_pyov_id`, `matrix_index`, `node_id` e.g. `PYO_01`, and active flags)
- `data/interim/nodes_pyoverdines_conservative.csv`  
  (flattened CSV of pyov node dictionary)
- `docs/phase_03/step_02_pyov_index/pyov_index_mapping_summary.txt`  
  (mapping summary, provenance, and the 93 vs 26 discussion)
- `docs/phase_03/step_02_pyov_index/session_info_pyov_index.txt`  
  (R session info for this step)

Notes:
- The current mapping produced 93 pyov nodes (one per validated synthetase id found in the `receptor_to_pyoverdine` expansion).
- The mapping includes `active_production` and `active_utilization` flags computed from FG/strain node tables where available.

### Step 3 — Adjacency construction (`scripts/phase_03_network/03_build_adjacencies.R`)
Purpose: build binary adjacency matrices and edge lists for FG and strain networks using the pyov mapping.

- FG-level adjacencies:
  - `data/interim/adj_production_FG_agentsxpyov_conservative.rds`  
    (sparse matrix, rows = FGs, cols = pyovs; agent → pyov production)
  - `data/interim/adj_utilization_FG_pyovxagents_conservative.rds`  
    (sparse matrix, rows = pyovs, cols = FGs; pyov → agent utilization)
  - `data/interim/edges_functional_groups_conservative.rds` / `.csv`  
    (binary edge list for FG-level network; columns: `source`, `target`, `edge_type`, `source_type`, `target_type`)
- Strain-level adjacencies:
  - `data/interim/adj_production_STR_agentsxpyov_conservative.rds`  
    (sparse matrix, rows = strains, cols = pyovs)
  - `data/interim/adj_utilization_STR_pyovxagents_conservative.rds`  
    (sparse matrix, rows = pyovs, cols = strains)
  - `data/interim/edges_strains_conservative.rds` / `.csv`  
    (strain-level binary edge list)
- Degree sequences:
  - `data/interim/degree_sequences_conservative.rds`  
    (structured list with FG and strain degree vectors for production/utilization)
- Diagnostics and summaries:
  - `docs/phase_03/step_03_adjacencies/phase_03_network_summary.txt`  
    (run summary; includes the 93 vs 26 mapping note and counts)
  - `docs/phase_03/step_03_adjacencies/adjacency_counts_check.csv`  
    (adjacency vs edge-list counts)
  - `docs/phase_03/step_03_adjacencies/inactive_pyovs.csv`  
    (pyovs with zero production & zero utilization at FG-level — empty if none)
  - `docs/phase_03/step_03_adjacencies/self_producer_consumer_examples_fg.csv`  
    (FGs that both produce and utilize the same pyov(s))
  - `docs/phase_03/step_03_adjacencies/fg_component_summary.csv`  
    (weak component summary for FG-level undirected graph)
  - `docs/phase_03/step_03_adjacencies/str_component_summary.csv`  
    (weak component summary for strain-level undirected graph)
  - `docs/phase_03/step_03_adjacencies/session_info_adj_build.txt`  
    (R session info for adjacency construction)

Notes:
- Adjacency semantics are explicit:
  - Production adjacency matrix is Agents (rows) × Pyov (cols); entry [i,j]=1 maps to an edge agent_i → pyov_j.
  - Utilization adjacency matrix is Pyov (rows) × Agents (cols); entry [r,c]=1 maps to an edge pyov_r → agent_c.
- The Step 3 outputs in this run used the 93-node pyov mapping (see above).

---

## Next steps and recommended actions

1. If you want the strict 26-node lock–key representation, produce the collapsed mapping:
   - Collapsing approach (recommended): map each original synthetase ID to its `pair_id` (one `pair_id` per pairing row), then assign `matrix_index` 1..26 in that order. Re-run the adjacency builder with the collapsed mapping to produce a 26-column pyov matrix.
   - Alternative (not recommended): pick a representative synthetase per pairing row (arbitrary).
   - I recommend keeping both 93-node and 26-node networks and comparing metrics as a sensitivity analysis.

2. Phase 4 preparation:
   - Use `data/interim/degree_sequences_conservative.rds` to plan null-model runs and resource estimates.
   - Suggested initial null-models: `swap.web()`-style edge swaps preserving degree sequences (start with 100 randomizations for quick runtime testing, then scale to 1,000+ for final inference).
   - For modularity: use bipartite-aware modularity detection with multiple replicates (n≥10) and store seeds for reproducibility.

3. Reporting:
   - When reporting network metrics (modularity, nestedness/NODF, degree distributions), always state which pyov mapping was used (93 vs 26).
   - If the largest FG is disproportionately large, run a sensitivity check excluding it or stratifying null models to preserve its size.

---

## Short usage examples (what scripts to run)
- Derive agent sets (Step 1):
  - `Rscript scripts/phase_03_network/01_derive_agent_sets.R`
- Create pyov mapping (Step 2):
  - `Rscript scripts/phase_03_network/02_pyov_index_mapping.R`
- Build adjacencies and edge lists (Step 3):
  - `Rscript scripts/phase_03_network/03_build_adjacencies.R`

Be sure to inspect:
- `docs/phase_03/step_01_agent_sets/agent_sets_summary.txt`
- `docs/phase_03/step_02_pyov_index/pyov_index_mapping_summary.txt`
- `docs/phase_03/step_03_adjacencies/phase_03_network_summary.txt`

These documents provide provenance and important notes about mapping and activity flags.

---

## Provenance and reproducibility
- Each step saves `session_info` to `docs/phase_03/step_*/*session_info*.txt`.
- All intermediate R objects are stored under `data/interim/` with descriptive filenames.
- Document the pyov mapping choice (93 vs 26) clearly in any downstream reports or publications.

---

If you want, I can:
- Produce the 26-node collapsed mapping and re-run Step 3 to create the 26-node variant (recommended for sensitivity),
- Or generate Phase 4 analysis templates (modularity, nestedness, null-model generation) tailored to the chosen pyov mapping.
