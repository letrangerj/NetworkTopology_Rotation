# Phase 3 Plan: Dual Conservative Network Construction (Updated)

## Overview
Phase 3 constructs TWO conservative bipartite, directed siderophore interaction networks and prepares robust diagnostics and null-model inputs for Phase 4.

Networks (both conservative, validated-only):
- Functional-group-level network: 177 functional groups × validated pyoverdines
- Strain-level network: 1,809 conservative strains × validated pyoverdines

Primary principles:
- Use only validated pyoverdine–receptor lock–key pairs (pairing_result: 26 groups) to define production and utilization edges.
- Persist actual production and utilization sets at the agent level (list-columns) so edge construction is deterministic and auditable.
- Standardize adjacency orientation semantics and file naming to avoid ambiguity.
- Provide mandatory diagnostics including direct comparison with Gu et al. (2025) before any topology analysis.
- Prepare degree sequences and runtime estimates to enable efficient null-model generation and computational feasibility planning.

---

## Inputs (From Phase 2)
- `data/interim/functional_groups_conservative.rds` (177 groups; includes counts, strategies)
- `data/interim/strain_repertoires.rds` (per-strain validated production & usable pyoverdines)
- `data/interim/receptor_to_pyoverdine.rds` (validated receptor → pyov mapping)
- `data/interim/rds/pairing_result.rds` (26 validated pyoverdine clusters)
- `data/interim/functional_groups_classified_conservative.rds` (strategy metadata)

---

## Key Design Decisions (concise)
- Dual networks to assess collapse effects and permit direct comparison to reference literature.
- Production adjacency stored as Agents × Pyoverdines (rows = agents; cols = pyov); entry [i,j] = 1 ⇒ agent_i produces pyov_j (edge: agent_i → pyov_j).
- Utilization adjacency stored as Pyoverdines × Agents (rows = pyov; cols = agents); entry [j,i] = 1 ⇒ pyov_j → agent_i (edge: pyov_j → agent_i). This avoids inverting semantics during edge export.
- Create and persist a contiguous pyoverdine index mapping (original pyov id → matrix index 1..26 → node_id `PYO_01`...).
- Functional-group nodes include `strain_count` as `node_weight`; weighted degree sequences are generated and saved alongside binary degrees.
- Edge lists will be binary (no `weight` column) to avoid premature quantitative assumptions.

---

## Network Definitions
- Pyoverdine node set:
  - `U_all`: all validated pyoverdine group ids (26 total)
  - `U_active`: subset of `U_all` that appear in ≥1 production or utilization edge (reported)
- Agent node sets:
  - `V_fg`: functional groups (FG_001 … FG_177)
  - `V_str`: strains (STR_0001 … STR_1809)
- Edge semantics:
  - Production edges: Agent → Pyoverdine
  - Utilization edges: Pyoverdine → Agent

---

## Detailed Implementation Steps

### Step 0 — Preflight & Conventions
- Define file & naming conventions (examples):
  - Adjacencies:
    - `data/interim/adj_production_FG_agentsxpyov_conservative.rds` (rows FG; cols PYO matrix_index)
    - `data/interim/adj_utilization_FG_pyovxagents_conservative.rds` (rows PYO; cols FG)
    - Equivalent filenames for strain-level: `..._STR_...`
  - Node tables:
    - `data/interim/nodes_functional_groups_conservative.{rds,csv}`
    - `data/interim/nodes_strains_conservative.{rds,csv}`
    - `data/interim/nodes_pyoverdines_conservative.{rds,csv}`
  - Degree sequences:
    - `data/interim/degree_sequences_conservative.rds`
  - Summary & diagnostics:
    - `docs/phase_03/phase_03_network_summary.txt`
    - `docs/phase_03/gu_comparison_required.txt` (mandatory)
- Save session info and package versions in `docs/phase_03/session_info.txt`.

### Step 1 — Derive & Persist Agent Sets (Non-circular validation)
- Functional-group sets:
  - For each FG, collect member strains via `strain_functional_group_mapping_conservative.rds`.
  - Recompute each member strain's validated production set and usable pyoverdine set from `strain_repertoires.rds`.
  - Assert equality across members; if any group fails, log error and save diagnostics (do not silently coerce).
  - Persist list-columns in FG table:
    - `validated_production_set` (integer vector)
    - `usable_pyoverdine_set` (integer vector)
  - Include `strain_count` (= node weight) and `strategy` fields.
- Strain agent table:
  - Filter `strain_repertoires.rds` to conservative strains (validated producer OR usable pyov).
  - Assign stable IDs `STR_0001...` and persist validated sets and counts.

**Non-circular validation (explicit):**
- Randomly sample 5–10 FG entries; for each:
  - Recompute per-strain sets directly from `strain_repertoires.rds`.
  - Compare to stored FG-level sets and log mismatches with sample-heavy diagnostic output (filename: `scripts/phase_03/debug/group_set_verification.txt`).

### Step 2 — Pyoverdine Index Mapping (contiguous indexing)
- Create a pyov lookup table mapping:
  - `original_pyov_id` → `matrix_index` (1..|U_all|) → `node_id` (`PYO_01`, `PYO_02`, ...)
- Use `matrix_index` for adjacency matrix column/row ordering to ensure contiguous matrices and simple indexing.
- Persist `nodes_pyoverdines_conservative.{rds,csv}` with mapping and placeholder fields `active_production`, `active_utilization` (filled post-edge construction).

### Step 3 — Build Adjacency Matrices (explicit orientations)
For both FG and STR networks:

3.1 Production adjacency (Agents × Pyov)
- Rows = Agents (FG or STR)
- Columns = Pyov (matrix_index order from Step 2)
- Entry [i, j] = 1 if pyov_j ∈ agent's validated production set
- Save with explicit object name indicating agent type.

3.2 Utilization adjacency (Pyov × Agents)
- Rows = Pyov (matrix_index order)
- Columns = Agents (same ordering used in production adjacency rows)
- Entry [r, c] = 1 if pyov_r ∈ agent_c's usable pyoverdine set
- Save separately; orientation directly corresponds to edges pyov → agent.

**Rationale**: This removes orientation inversion and makes edge semantics direct from matrices.

### Step 4 — Edge List Construction (binary edges)
- Convert production adjacency (Agents × Pyov) to edges:
  - For each [i,j]==1 → `source = agent_id`, `target = pyov_id`, `edge_type = production`
- Convert utilization adjacency (Pyov × Agents) to edges:
  - For each [r,c]==1 → `source = pyov_id`, `target = agent_id`, `edge_type = utilization`
- Do NOT include a `weight` column; edges are binary.
- Add `source_type` and `target_type` metadata for clarity.
- Deduplicate edges (should be unique by construction).
- Save edges as `.rds` and `.csv`.

### Step 5 — Validation & Consistency Checks (both networks)
- **Bipartite integrity**: ensure agent IDs and pyov IDs have disjoint namespaces.
- **Set integrity**: if Step 1 verification sampled groups, include results here.
- **Edge ↔ adjacency consistency**: sum adjacency rows/cols match edge counts.
- **Duplicate edge detection**: zero duplicates expected.
- **Inactive pyoverdines**: list pyov in `U_all` with zero total degree (prod + util); save `inactive_pyov_list.csv`.
- **Self-loops**: agents where the same pyov appears in both production and usable sets (report count & examples).
- **Connectivity**:
  - Weak components: treat entire directed network as undirected; report number of components, size distribution, proportion in largest component.
  - Directed reachability metrics (preserve directionality):
    - For each pyov: number of agents directly reachable (consumers).
    - For each agent: number of pyov directly reachable (produced).
    - Optionally report two-step reachability (pyov→agent→pyov) as exploratory metric.
- **Spot checks**: randomly sample 10 agents and assert adjacency row/col equals their stored sets.

### Step 6 — Degree Sequence Extraction & Structure
- Compute and persist degree sequences for both networks with the structured, unambiguous layout:

```r
list(
  fg_network = list(
    fg_production_out = integer_vector,    # FG -> PYO (row sums of production adjacency)
    fg_utilization_in = integer_vector,    # FG (consumers) <- PYO (col sums of utilization adjacency)
    pyo_production_in = integer_vector,    # PYO <- FG (col sums of production adjacency)
    pyo_utilization_out = integer_vector   # PYO -> FG (row sums of utilization adjacency)
  ),
  strain_network = list(
    str_production_out = integer_vector,
    str_utilization_in = integer_vector,
    pyo_production_in = integer_vector,
    pyo_utilization_out = integer_vector
  )
)
```

- Also compute weighted versions (functional-group `node_weight = strain_count`):
  - Weighted agent degrees (sum of incident edges weighted by node weight when appropriate — document formula)
- Save as `data/interim/degree_sequences_conservative.rds`.

### Step 7 — Network Summary & Gu et al. Comparison (MANDATORY)
- For each network, produce a summary table:
  - |Agents|, |Pyov_all|, |Pyov_active|
  - Production edges, Utilization edges, Total edges
  - Utilization:Production edge ratio
  - Density (production & utilization)
  - Top-5 hubs (agents & pyov) by degree
  - Self-loop count & examples
  - Weak component counts and giant component proportions
- Mandatory comparison:
  - Create `docs/phase_03/gu_comparison_required.txt` that reports our metrics alongside Gu et al. reference values (if available).
  - If Avail. Gu-values are not in the repo, include placeholders and a short guidance text for how to populate them.
  - If utilization:production ratio deviates substantially from Gu et al. (e.g., our ratio << 6.0), flag the result and recommend sensitivity analyses (relax pairing or alternative collapse).
- Save combined FG + STR comparison table as RDS + CSV.

### Step 8 — Computational Feasibility Estimates & Guidance
- Estimate memory and runtime complexity for core Phase 4 metrics on the constructed networks:
  - Heuristics: `computeModules()` and null-model runs scale with edges and number of nodes; NODF scales ~O(n^2) where n ~ max(rows, cols).
  - Produce `docs/phase_03/phase_03_runtime_estimates.txt` with:
    - Matrix sizes, edge counts, recommended memory (RAM) and cores for typical runs
    - Recommended null-model parallelization parameters and an initial adaptive randomization workflow
- If estimates exceed practical limits, recommend:
  - Subsampling strategies (stratified by FG size or taxon)
  - Limiting null-model runs initially (e.g., 100 → 1000 adaptive)
  - Using cluster / multi-core resources for `swap.web()` randomizations

### Step 9 — Save All Outputs & Provenance
- Persist all outputs with clear naming and versioning:
  - Adjacencies, edge lists, node dictionaries, degree sequences, network summaries, runtime estimates, session info.
- Save a `docs/phase_03/phase_03_network_construction_log.txt` listing scripts run, parameters, and a short narrative of any anomalies encountered.

### Step 10 — Prepare for Phase 4
- Confirm degree sequences available and consistent.
- Ensure Gu et al. comparison completed with interpretation notes.
- Provide recommended Phase 4 parameters (null-model types, number of replicates, seeds, and any subsampling rules).

---

## Expected Outputs (Concise)
- Adjacency matrices (production: agents×pyov; utilization: pyov×agents) for FG and STR networks
- Edge lists (binary) for production & utilization (FG & STR)
- Node dictionaries for FG, STR, PYO (with contiguous pyov `matrix_index`)
- Degree sequences saved in structured, explicit format (including weighted variants)
- Network summary tables and mandatory Gu et al. comparison report
- Runtime/memory estimate report
- Diagnostic artifacts: set-verification logs, inactive-pyov lists, self-loop report

---

## Quality Assurance / Risk Mitigation Checklist
- [ ] Contiguous pyov matrix indexing table created and used consistently
- [ ] FG-level list-columns persisted (`validated_production_set`, `usable_pyoverdine_set`) and validated against raw strain repertoires (non-circular checks)
- [ ] Utilization adjacency stored as pyov × agents (direct mapping to pyov→agent edges)
- [ ] No `weight` column in initial edge lists (binary edges only)
- [ ] Functional-group node table includes `strain_count` (node weight) and weighted-degree sequences saved
- [ ] Mandatory Gu et al. comparison produced; flag if discrepancies exceed thresholds
- [ ] Runtime estimate file produced and reviewed to select Phase 4 strategy (subsampling / parallelization if required)
- [ ] Self-loop counts reported and interpreted biologically
- [ ] Degree sequence object uses the clear nested layout required by Phase 4

---

## Future Extensions (Scaffolded but not executed here)
- Relaxed pairing inclusion (similarity-based) and sensitivity networks
- Weighted edge models (receptor copy numbers, expression data)
- Metadata-driven subnetwork extraction (habitat, pathogenicity)
- Additional null-model variants and larger-scale randomization strategies

---

## Ready for Next Action
This plan is finalized and incorporates:
- adjacency orientation fixes,
- contiguous pyov indexing,
- non-circular validation,
- mandatory Gu et al. comparison,
- runtime feasibility checks,
- degree sequence structure clarity,
- self-loop detection,
- removal of premature edge weights, and
- node-weight infrastructure.

If you approve, I will begin implementing the updated plan starting with Step 1 (derive & persist group-level sets, pyov lookup, and node dictionaries) and the runtime estimator script.
