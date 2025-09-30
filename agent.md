# Siderophore Network Analysis Project - Agent Instructions

## Project Management Guidelines

### Folder Structure and Organization
The project must maintain a clear, standardized directory structure to ensure reproducibility and collaboration:
- data/
  - raw/ (immutable; original inputs). Place all source .mat files here. Current inputs are located under Burkholderiaceae/; recommended location: data/raw/Burkholderiaceae/
  - interim/ (per-phase intermediate outputs)
  - processed/ (final analysis-ready artifacts)
- scripts/
  - phase_01_explore_validate/ (with debug/)
  - phase_02_functional_groups/ (with debug/)
  - phase_03_network/ (with debug/)
  - phase_04_topology/ (with debug/)
  - phase_05_visualization/ (with debug/)
  - phase_06_comparative/ (with debug/)
- results/
  - metrics/
  - null_models/
  - stats/
- figures/
  - network/
  - distributions/
  - comparative/
- reports/
  - methods/
  - data_quality/
  - analysis/
- docs/ (software versions, seeds, run logs)

### Development and Debugging Guidelines

1. **Modular Organization**
   - Each phase gets its own subfolder in `scripts/`
   - Create `debug/` subfolder within task folders when debugging is needed
   - Use descriptive file names: `01_load_data.R`, `02_validate_groups.R`

2. **Basic Quality Control**
   - Document random seeds and key parameters in script headers
   - Keep `data/raw/` unchanged after initial loading
   - Save intermediate results at each phase completion


## Project Overview
This project constructs and analyzes siderophore-mediated iron interaction networks from bacterial genomic data, following the methodological framework established in Gu et al. (2025) Science Advances. The goal is to understand microbial social interactions through pyoverdine (siderophore) production and receptor networks using bipartite network topology analysis.

## Data Structure

### Input Files (.mat format)
1. **syn.mat** - Biosynthetic Gene Cluster (BGC) information. Contains a `syn` struct with the following fields: `clusterblast`, `strainName`, `regionName`, `assemblyDefinition`, `location`, `wholeCDS`, `biosynCDS`, `group`, `regionIdentifier`.
   - Contains potential pyoverdine BGCs (where group ≠ 0).
   - Groups correspond to clustering in `Syn_PWA_inGroup.mat`.
   - Pre-clustered using 0.2 threshold based on sequence distance.

2. **rec.mat** - TonB-Plug receptor information. Contains a `rec` struct with the following fields: `count`, `foldername`, `fragmentname`, `boarders`, `seqs`, `tag`, `recname`, `domloc`, `domseq`, `group`.
   - All potential TonB-Plug receptors from genomes.
   - Groups correspond to clustering in `Rec_MSA_inGroup.mat`.
   - Pre-clustered using 0.2 threshold based on sequence distance.

3. **pairing_result.mat** - Refined synthetase-receptor pairing. Contains a `final_pairs` cell array.
   - Based on FpvA receptor fine-grouping from SA paper.
   - Column 1: Synthetase genes in same group.
   - Column 2: Receptor genes in same group.
   - ~26 pyoverdine-receptor groups total.
   - Other combinations likely not true pyoverdine siderophores.

4. **Syn_PWA_inGroup.mat** - Synthetase gene clustering data. Contains a `cluster_info` struct with the fields `dataname` and `distm`.

5. **Rec_MSA_inGroup.mat** - Receptor gene clustering data. Contains a `cluster_info` struct with the fields `dataname` and `distm`.

## Primary Goals

### 1. Network Construction
Build **bipartite directed network** with:
- **Microbial functional groups**: Strains collapsed by identical pyoverdine production and receptor utilization repertoires (following Gu et al. strain-level collapse rule)
- **Lock-key pyoverdine groups**: The 26 verified pyoverdine-receptor pairs from pairing_result.mat
- **Production edges**: Functional group → pyoverdine group (directional)
- **Utilization edges**: Pyoverdine group → functional group (directional)

### 2. Network Topology Analysis
Calculate bipartite-specific network metrics:
- **Modularity**: Using `bipartite::computeModules()` with multiple replicates for stochastic stability
- **Nestedness**: NODF metric using `bipartite::nested()` (specify version for reproducibility)
- **Degree distributions**: In-degree (production) and out-degree (utilization) for pyoverdine nodes
- **Shannon entropy**: Of functional group frequency distribution H = -Σ(pi × log(pi))
- **Syn/Rec ratios**: Both per lock-key group and global average (synthetase count / receptor count)

### 3. Comparative Analysis
Compare network properties across conditions when metadata available:
- Different habitats (requires strain-to-habitat mapping)
- Pathogenic vs. non-pathogenic organisms (requires lifestyle annotation)
- Functional group strategies (single vs. multi-receptor producers vs. nonproducers)

## Technical Requirements

### Software Stack
**Primary: R**
```r
# Required packages
library(R.matlab)      # Read .mat files
library(bipartite)     # Bipartite network analysis (specify: computeModules vs metaComputeModules)
library(igraph)        # General network metrics
library(ggplot2)       # Visualization
library(dplyr)         # Data manipulation
library(vegan)         # Alternative nestedness calculation (nestednodf)
```

**Secondary: MATLAB**
```matlab
% MATLAB is installed and accessible via the 'matlab' command.
% For data preprocessing and coevolution algorithm validation; matrix ops and distance matrix analysis.
% Preferred batch/background usage examples:
%   matlab -batch "yourFunction"  OR  matlab -nodisplay -nosplash -r "try; yourFunction; catch; exit(1); end; exit"
```

### Critical Implementation Specifications

1. **Functional Group Definition**
   - Strain-level collapse: Group strains with identical production + utilization profiles
   - Use exact repertoire matching (not similarity-based)
   - Document collapse criteria for reproducibility

2. **Distance Matrix Integration**
   - Validate the 26 lock-key groups using Syn_PWA_inGroup.mat and Rec_MSA_inGroup.mat
   - Implement coevolution strength verification (correlation between distance matrices)
   - Do not treat pairing_result.mat as hard constraints without validation

3. **Null Model Specification**
   - Use `bipartite::swap.web()` for edge-swapping null models
   - Preserve bipartite structure and degree sequences
   - Generate 1000+ randomizations for statistical power
   - Document exact null model algorithm for reproducibility

4. **Metadata Integration Plan**
   - Define strain identifier mapping strategy
   - Specify external database joins (NCBI, PATRIC, etc.)
   - Handle missing metadata systematically
   - Document metadata source and version
   - Usage and timing: Metadata are primarily used in Phase 6 for comparative analyses (habitat- and lifestyle-based subsetting and metric comparison). Phases 1–5 do not require metadata and can proceed now

## Phase-Based Implementation

### Phase 1: Data Exploration and Validation
- [ ] Load and inspect all .mat files structure and dimensions
- [ ] Map strain/organism identifiers across all files
- [ ] Verify correspondence between syn.mat groups and Syn_PWA_inGroup.mat
- [ ] Verify correspondence between rec.mat groups and Rec_MSA_inGroup.mat
- [ ] Validate the 26 pyoverdine-receptor pair assignments in pairing_result.mat
- [ ] Create comprehensive data summary statistics
- [ ] Identify and document any missing or inconsistent data
- [ ] Test coevolution algorithm on distance matrices

### Phase 2: Functional Group Classification
- [ ] Implement strain-level collapse algorithm based on identical repertoires
- [ ] Classify organisms into three categories:
  - Single-receptor producers: 1 pyoverdine synthesis + 1 receptor
  - Multi-receptor producers: 1 pyoverdine synthesis + multiple receptors
  - Nonproducers: 0 pyoverdine synthesis + ≥1 receptor
- [ ] Create functional group mapping tables
- [ ] Validate functional groups against biological expectations
- [ ] Document group collapse criteria and statistics

### Phase 3: Network Construction
- [ ] Build bipartite adjacency matrix (functional groups × pyoverdine groups)
- [ ] Create directed edge lists for production and utilization
- [ ] Implement network validation checks (connectivity, isolated components)
- [ ] Verify network structure matches expected bipartite properties
- [ ] Cross-validate with Gu et al. network statistics where possible
- [ ] Generate network summary statistics (nodes, edges, density)

### Phase 4: Network Topology Analysis
- [ ] Calculate basic network metrics with multiple random seeds
- [ ] Implement bipartite modularity analysis with `metaComputeModules()`
- [ ] Calculate nestedness (NODF) using both `bipartite` and `vegan` packages
- [ ] Analyze degree distributions for pyoverdine nodes
- [ ] Calculate Shannon entropy of functional group frequencies
- [ ] Compute Syn/Rec ratios per lock-key group and globally
- [ ] Generate null model distributions for all metrics
- [ ] Perform statistical significance testing

### Phase 5: Visualization and Network Representation
- [ ] Create bipartite network layout visualizations
- [ ] Generate degree distribution plots (log-scale, heavy-tail analysis)
- [ ] Create network property comparison plots
- [ ] Develop interactive network exploration interface
- [ ] Create summary figures matching Gu et al. style
- [ ] Generate publication-quality network diagrams

### Phase 6: Comparative Analysis (Metadata-Dependent)
- [ ] Integrate strain metadata (habitat, pathogenicity, taxonomy)
- [ ] Subset networks by habitat categories
- [ ] Subset networks by pathogenic vs. non-pathogenic lifestyle
- [ ] Calculate network metrics for each subset
- [ ] Perform statistical comparisons between network properties
- [ ] Implement multiple testing corrections
- [ ] Analyze functional group strategy distributions across conditions
- [ ] Generate comparative visualization panels

## Critical Considerations

### Biological Constraints
- Only BGCs with group ≠ 0 in syn.mat represent validated pyoverdine clusters
- Only the 26 pairs in pairing_result.mat are experimentally validated
- Distance matrices must be used to validate, not bypass, the coevolution algorithm
- Network modules must have biological interpretability

### Technical Specifications
- **Bipartite modularity**: Use stochastic algorithms with multiple replicates (n≥10)
- **NODF calculation**: Document which package/version for consistency
- **Null models**: Specify exact algorithm (edge-swapping with degree preservation)
- **Statistical testing**: Apply FDR correction for multiple comparisons
- **Reproducibility**: Set random seeds for all stochastic processes

### Data Quality Requirements
- Cross-validate distance matrix clustering with reported groups
- Verify strain identifier consistency across files
- Handle edge cases (strains with no receptors, orphan synthetases)
- Document all data filtering and quality control steps

## Expected Outputs

### Data Products
- Processed bipartite adjacency matrices
- Functional group classification tables
- Network topology metric calculations
- Statistical significance test results
- Null model distributions

### Visualizations
- Bipartite network diagrams with biological annotations
- Degree distribution plots with power-law fitting
- Network property comparison panels
- Interactive network exploration tools
- Comparative analysis dashboards

### Analysis Reports
- Network topology characterization with statistical validation
- Comparative analysis results across conditions
- Biological interpretation of network modules and hubs
- Methodological validation against Gu et al. reference
- Data quality assessment and limitations documentation

## Success Criteria
- Successfully construct validated bipartite siderophore interaction network
- Reproduce key network topology patterns from Gu et al. (2025)
- Generate statistically significant network modules with biological interpretation
- Deliver robust comparative analysis across available metadata categories
- Provide reproducible analysis pipeline with documented methodology

## Risk Mitigation Strategies
- **Distance matrix validation**: Implement coevolution algorithm before accepting pairing_result.mat
- **Null model reliability**: Test multiple randomization methods for consistency
- **Metadata integration**: Plan for missing data scenarios with alternative analyses
- **Reproducibility**: Document all software versions, random seeds, and parameter choices
- **Biological validation**: Cross-reference results with known siderophore biology literature

---
*This project implements the network topology analysis framework from Gu et al. (2025) "Siderophore synthetase-receptor gene coevolution reveals habitat- and pathogen-specific bacterial iron interaction networks" Science Advances 11:eadq5038.*
