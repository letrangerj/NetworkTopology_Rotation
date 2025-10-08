# Phase 4 Plan: Bipartite Network Topology Analysis & Optional Inference

## Overview
Phase 4 performs quantitative topological characterization of the conservative siderophore (pyoverdine) interaction networks produced in Phase 3. Two bipartite, directed networks are analyzed:

- Functional-group network (FG-level): functional groups × pyoverdine groups  
- Strain network (STR-level): strains × pyoverdine groups

Primary goals:
- Compute core bipartite metrics (modularity, nestedness, degree structure, entropy, Syn/Rec ratios).
- Produce compact, reproducible outputs and figures that let the analyst decide whether deeper null-model inference or sensitivity tests are required.
- Treat null-model significance testing, robustness quantification, and sensitivity analyses as OPTIONAL follow-ups (performed only when warranted by first-pass metrics).

All outputs must be reproducible, versioned, and accompanied by minimal provenance (parameters and session info).

---

## Objectives
1. Derive bipartite modular structure and report module assignments and summary statistics.  
2. Quantify nestedness (NODF) with at least one reproducible implementation and cross-check when feasible.  
3. Characterize degree distributions for producers and consumers on both network levels.  
4. Compute Shannon entropy of functional-group frequencies (strain-weighted).  
5. Calculate Syn/Rec ratios per pyov node and globally.  
6. Provide clear, consolidated outputs and visualizations to determine whether optional inference (nulls, sensitivity) is needed.

OPTIONAL (run only if first-pass results are surprising or require formal testing):
- Generate null-model distributions and test metric significance.  
- Quantify modularity robustness and produce sensitivity analyses (largest-group removal, alternative pyov aggregation).

---

## Inputs (from Phase 3)
Required:
- `data/interim/adj_production_FG_agentsxpyov_conservative.rds`
- `data/interim/adj_utilization_FG_pyovxagents_conservative.rds`
- `data/interim/adj_production_STR_agentsxpyov_conservative.rds`
- `data/interim/adj_utilization_STR_pyovxagents_conservative.rds`
- `data/interim/nodes_functional_groups_conservative.rds`
- `data/interim/nodes_strains_conservative.rds`
- `data/interim/nodes_pyoverdines_conservative.rds`
- `data/interim/degree_sequences_conservative.rds`

Optional (only if sensitivity tests run):
- Alternative/variant mapping artifacts (same file forms under a variant folder)

---

## Core Metric Definitions
- Modularity (bipartite): community detection on the bipartite incidence matrix; report modularity scores and module memberships.  
- Nestedness (NODF): nested ordered overlap metric; record implementation details.  
- Degree distributions: producer (agent → pyov) out-degree, pyov in-degree; consumer (pyov → agent) out-degree, agent in-degree.  
- Shannon entropy: H = -Σ p_i log p_i for FG strain-count distribution; report normalized entropy.  
- Syn/Rec ratio: per-pyov synthetase count divided by receptor-mediated consumer count; include global summary statistics.  
- (OPTIONAL) Null model metrics: sampled metric vectors from degree-preserving randomized webs.

---

## Design Principles
- Reproducibility: fix master seed and record all derived seeds; save session info.  
- Minimality: produce consolidated result artifacts (one RDS per network containing numeric arrays, replicates, and figure pointers) to reduce file proliferation.  
- Binary-first: compute everything on presence/absence matrices; defer weights to optional analyses.  
- Optional inference: null-model and sensitivity steps are available but not executed by default. They are run only upon request after inspecting first-pass results.  
- Documentation: each primary output accompanied by a small JSON or text metadata file describing input files, parameters, and seed.

---

## Implementation: Step-by-step

### Step 0 — Environment & Preflight
- Create scaffolding:
  - `scripts/phase_04_topology/`
  - `results/phase_04/`
  - `figures/network_topology/`
  - `docs/phase_04/logs/`
- Write parameter manifest: `docs/phase_04/parameters_manifest.txt` (master seed, replicate counts, null settings (optional)).
- Save `docs/phase_04/session_info.txt` (package versions / runtime info).

### Step 1 — Data Loading & Sanity Checks
- Load production and utilization adjacency matrices. Confirm:
  - Dimensions and name consistency with node tables
  - Binary / integer values expected
  - No duplicated row/column names
- Construct canonical incidence matrix used by subsequent metrics (document orientation).

### Step 2 — Degree & Distribution Summaries
- Compute degree vectors for FG and STR networks:
  - FG production out-degree (per-FG producing pyovs)
  - FG utilization in-degree (per-FG consuming pyovs)
  - Pyov production in-degree and utilization out-degree
  - The same for STR-level
- Produce:
  - Summary statistics (min, median, mean, max, SD)
  - Degree histograms and CCDF (saved under `figures/network_topology/`)
- Save consolidated outputs:
  - `results/phase_04/degree_metrics_fg.rds`
  - `results/phase_04/degree_metrics_str.rds`

### Step 3 — Modularity (replicate analysis)
- Build a binary incidence matrix (FG × pyov and STR × pyov as applicable).
- Run modularity detection replicates:
  - Pilot: R = 20 replicates; Final: R = 100 (adjust if variance or runtime requires)
- Save a single consolidated object:
  - `results/phase_04/modularity_fg.rds` (contains replicate partitions, Q values, replicate seeds, and consensus partition)
- Produce figure showing consensus modules and module-size distribution: `figures/network_topology/modularity_consensus_fg.pdf`.
- Repeat for STR network if resources permit.

### Step 4 — Nestedness (NODF)
- Compute NODF using one stable implementation. Optionally cross-check with a second implementation and report both values and methods.
- Save:
  - `results/phase_04/nestedness_fg.rds`
  - `results/phase_04/nestedness_str.rds`
- Save visualization: `figures/network_topology/nodf_fg.pdf`

### Step 5 — Shannon Entropy
- Compute FG entropy from `strain_count` (H and normalized H/H_max).
- Save:
  - `results/phase_04/entropy_fg.rds`
- Save figure: `figures/network_topology/fg_entropy.pdf`

### Step 6 — Syn/Rec Ratios
- For each pyov node compute:
  - synthetase count (producer count, FG- or strain-level as chosen)
  - receptor-based consumer count
  - ratio and handling policy for zero denominators
- Save:
  - `results/phase_04/syn_rec_fg.rds`
- Save distribution figure: `figures/network_topology/syn_rec_distribution.pdf`

### OPTIONAL Step 7 — Null Model Framework (execute only on request)
- Null strategy: degree-preserving edge swaps on the incidence matrix.
- Pilot: N = 100; Final: N = 1000 (run only if needed)
- Compute chosen metrics per null replicate (e.g., NODF and, optionally, modularity).
- Save compact metric vectors to `results/phase_04/nulls_fg.rds` (and `_str.rds` when run).
- Produce null-vs-observed plots as needed.

### OPTIONAL Step 8 — Statistical Testing & Multiple Comparisons
- When nulls are available:
  - Compute empirical p-values using standard empirical formula with +1 correction.
  - Compute z-scores.
  - Apply FDR correction across tested metrics.
- Save a single significance table:
  - `results/phase_04/significance_table.csv`

### OPTIONAL Step 9 — Sensitivity Analyses (execute only on request)
- Largest FG removal: recompute core metrics and record percentage changes.
- Alternative pyov aggregation (if desired): rebuild incidence and recompute core metrics; report deltas.
- Save consolidated sensitivity summary:
  - `results/phase_04/sensitivity_summary.rds`

### Step 10 — Aggregated Reporting
- Consolidate primary results into a small number of artifacts:
  - `results/phase_04_topology_summary.rds` — canonical RDS bundle for programmatic inspection (metrics, replicates, pointers to figures and input versions).
  - `results/phase_04_summary_metrics.csv` — flat tabular summary for reporting and plotting.
- Create `docs/phase_04/phase_04_metrics_summary.txt` (human-readable high-level summary).

---

## File & Naming Conventions (consolidated)
- Primary programmatic artifact per network: `results/phase_04_topology_summary.rds` (contains FG & STR subsets).  
- Flat summary table for reports: `results/phase_04_summary_metrics.csv`.  
- Visual outputs in `figures/network_topology/` (PDF/PNG).  
- Optional heavy artifacts (nulls, sensitivity) stored under `results/phase_04/` and referenced inside the main RDS bundle when present.  
- Metadata and seeds: `docs/phase_04/parameters_manifest.txt`, `docs/phase_04/session_info.txt`, `docs/phase_04/phase_04_run_log.txt`.

---

## Default Parameters (recommended)
- Master seed: `2025` (derive replicate seeds deterministically)
- Modularity replicates (pilot/final): FG 20 / 100; STR 10 / 50
- Null-run N (pilot/final): 100 / 1000 (OPTIONAL)
- Swaps per null: heuristic `10 × edge count` (tunable)
- Convergence threshold: increase N if SD of null metric > 5% of null mean

---

## Quality Assurance Checklist
- [ ] Adjacency dimensions match node tables and Phase 3 degree sequences  
- [ ] Degree vectors reproduce Phase 3 stored sequences  
- [ ] Modularity replicates logged with seeds; consensus partition reported  
- [ ] NODF computed and method recorded; cross-check stored if performed  
- [ ] Entropy and Syn/Rec calculations include provenance and handling rules for zeros  
- [ ] Optional null runs preserve degree sequences (sample validation) — performed only if requested  
- [ ] All outputs accompanied by metadata and session info

---

## Risk Mitigation
- High runtime for STR nulls: run pilot sample; parallelize where available; use sparse matrices.  
- Modularity instability: increase replicates; present consensus with uncertainty.  
- NODF discrepancies: confirm matrix orientation and method args.  
- Dominant FG effects: report metrics both with and without major FG(s) if influence suspected.  
- Storage pressure: store only metric vectors for nulls; avoid storing full randomized webs.

---

## Extension Hooks (future)
- Weighted modularity using strain_count as node weight.  
- Alternative null models (probabilistic configuration models).  
- Module enrichment analysis using strain metadata (Phase 6).  
- Phylogenetic/temporal overlays if lineage/time data available.
