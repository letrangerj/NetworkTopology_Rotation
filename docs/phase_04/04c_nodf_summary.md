# Phase 04 — Step 04c: NODF Nestedness Summary

## Run metadata

| Field | Value |
|---|---|
| Script | `scripts/phase_04_topology/04c_nodf_analysis_only.py` |
| Inputs (edges) | `data/interim/edges_functional_groups_conservative.csv` (conservative FG set; 93 PYO nodes representation) |
| Module assignments | `results/phase_04/modularity/module_assignments_production_sk.csv` / `results/phase_04/modularity/module_assignments_utilization_sk.csv` |
| Master seed | `42` (layer-specific offsets applied by script) |
| Null model type | Configuration-model (stub-matching) — randomized bipartite matrices |
| Null replicates | `100` per layer |
| Output directories | `results/phase_04/nodf` and `figures/network_topology/step04_nestedness` |
| Summary log | `docs/phase_04/logs/04c_nodf_analysis_summary.txt` |

---

## Per-layer summary (observed NODF vs null)

| Layer | FG × PYO | Fill (%) | NODF_total (obs) | NODF_rows | NODF_cols | Null_mean | Null_sd | z-score | p (two-tailed) | Modules analyzed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Production | 134 × 93 | 1.404 | 1.976 | 2.915 | 1.037 | 2.575 | 0.141 | −4.237 | 0.0000 | 15 |
| Utilization | 173 × 93 | 13.910 | 19.471 | 23.117 | 15.826 | 25.716 | 0.404 | −15.464 | 0.0000 | 10 |

Notes:
- Fill = 100 × (number of ones) / (rows × cols).
- NODF values and null summaries are taken from run artifacts in `results/phase_04/nodf`.
- Z-score computed as (NODF_obs − Null_mean) / Null_sd; p-value is empirical two‑tailed based on null replicates.

---

## Modularity — module-level nestedness (summary table)

Below are per-module statistics computed from Phase 03 module assignments. The module tables now show only `NODF_total` (rounded to 2 decimal places) along with module sizes and fill percentages.

### Production layer — module-level NODF

| Module ID | n_rows (FG) | n_cols (PYO) | Module size | Fill (%) | NODF_total |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0  | 17 | 8  | 25 | 18.38 | 35.90 |
| 1  | 11 | 10 | 21 | 13.64 | 21.98 |
| 2  | 14 | 9  | 23 | 14.29 | 23.85 |
| 3  | 5  | 4  | 9  | 25.00 | 0.00 |
| 4  | 8  | 5  | 13 | 20.00 | 0.00 |
| 5  | 13 | 10 | 23 | 19.23 | 35.02 |
| 6  | 6  | 4  | 10 | 25.00 | 0.00 |
| 7  | 7  | 6  | 13 | 19.05 | 16.67 |
| 8  | 19 | 11 | 30 | 12.92 | 11.48 |
| 9  | 5  | 4  | 9  | 25.00 | 0.00 |
| 10 | 3  | 3  | 6  | 33.33 | 0.00 |
| 11 | 7  | 4  | 11 | 25.00 | 0.00 |
| 12 | 8  | 7  | 15 | 14.29 | 0.00 |
| 13 | 6  | 4  | 10 | 25.00 | 0.00 |
| 14 | 5  | 4  | 9  | 25.00 | 0.00 |

Key takeaways (production):
- Several mid-to-large modules (M0, M5, M2) show substantial internal nestedness (NODF_total > 20).
- Many smaller modules have NODF_total = 0; these modules often have equal-degree or complete submatrix structures that yield no hierarchical ordering under NODF.

---

### Utilization layer — module-level NODF

| Module ID | n_rows (FG) | n_cols (PYO) | Module size | Fill (%) | NODF_total |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 31 | 16 | 47 | 66.13 | 85.71 |
| 1 | 4  | 2  | 6  | 100.00 | 0.00 |
| 2 | 10 | 9  | 19 | 33.33 | 0.00 |
| 3 | 55 | 21 | 76 | 67.53 | 84.09 |
| 4 | 5  | 3  | 8  | 100.00 | 0.00 |
| 5 | 6  | 6  | 12 | 100.00 | 0.00 |
| 6 | 4  | 4  | 8  | 100.00 | 0.00 |
| 7 | 5  | 4  | 9  | 100.00 | 0.00 |
| 8 | 48 | 25 | 73 | 37.67 | 54.58 |
| 9 | 3  | 2  | 5  | 100.00 | 0.00 |

Key takeaways (utilization):
- Very large, dense modules (M0, M3) exhibit extremely high internal nestedness (NODF_total > 80), indicating strong hierarchical core structures.
- Several small modules report NODF_total = 0 despite high fill; such patterns indicate equal-degree complete subgraphs lacking hierarchical order by NODF.
- Module 8 is large and moderately nested (NODF_total ≈ 54.58), indicating an additional nested core beyond M0 and M3.

---

## Overall comparison: global NODF vs module-level NODF

The run summary includes module-level summary statistics and a direct comparison of the observed overall NODF to the distribution of module NODF values. These comparison metrics are reported below (values taken from the run log: `docs/phase_04/logs/04c_nodf_analysis_summary.txt`).

### Production — overall vs modules

- Number of modules considered: 15  
- Module NODF mean ± sd: 9.661 ± 13.569  
- Observed overall NODF vs module mean: z = −0.566, p = 0.5712

Interpretation: The observed overall NODF (1.98) is slightly below the mean of module-level NODF values, but this difference is not statistically significant (p ≈ 0.57). Module-level nestedness is heterogeneous: a few modules are strongly nested while many small modules have zero nestedness, resulting in a moderate module mean and large variance.

### Utilization — overall vs modules

- Number of modules considered: 10  
- Module NODF mean ± sd: 22.438 ± 37.061  
- Observed overall NODF vs module mean: z = −0.080, p = 0.9362

Interpretation: The observed overall NODF (19.47) is effectively indistinguishable from the mean of module-level NODF values (p ≈ 0.94). Large, highly nested modules (e.g., M0, M3) substantially elevate the module-level mean and variance; when aggregated, the global NODF aligns with the central tendency of module-level nestedness in this layer.

---

## Methods at-a-glance

| Component | Details |
|---|---|
| Incidence construction | Build binary FG × PYO matrices from deduplicated edges; convert sparse → dense for NODF and plotting. |
| Sorting | Rows and columns deterministically sorted by degree (descending) with index tie-breakers before visualization & NODF. |
| NODF computation | Custom Almeida‑Neto (2008) style: compute pairwise overlaps; skip equal-degree pairs; average row & column components → NODF_total. |
| Null model used | Configuration-model (stub-matching) — preserve endpoint stub counts, shuffle, instantiate randomized bipartite matrix; binarize (>0) for NODF. |
| Within-module analysis | Use Phase 03 module assignments; analyze modules with ≥ 2 rows and ≥ 2 columns; compute module-level NODF and fill%. |
| Reproducibility | Master seed = 42; per-layer offsets; null replicates n = 100. |

---

## Notes on provenance and reproducibility

- Script run: `scripts/phase_04_topology/04c_nodf_analysis_only.py`  
- Inputs: conservative FG set from Phase 03; expanded 93-node PYO representation was used.  
- Seed: master seed = 42 (layer-specific offsets applied by the script).  
- Null replicates: 100 per layer (configuration-model nulls).  
- Full run summary and per-module breakdown are available in `docs/phase_04/logs/04c_nodf_analysis_summary.txt`.

---

End of Phase 04 Step 04c NODF summary.