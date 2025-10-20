#!/usr/bin/env python3

import os
from datetime import datetime
from typing import Dict, List, Optional, Tuple
from tqdm import tqdm
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from numpy.random import Generator, default_rng
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr, norm
import warnings
warnings.filterwarnings('ignore')


def _bipartite_config_model(biadj: csr_matrix, random_state: Optional[int] = None) -> csr_matrix:
    """
    Create a bipartite configuration model null network.
    Preserves exact row and column degree sequences while randomizing connections.
    """
    n_rows, n_cols = biadj.shape
    
    # Extract degree sequences
    row_degrees = np.asarray(biadj.sum(axis=1)).ravel().astype(int)
    col_degrees = np.asarray(biadj.sum(axis=0)).ravel().astype(int)
    
    # Quick validation
    if row_degrees.sum() != col_degrees.sum():
        raise ValueError("Row and column degree sums must match")
    
    if row_degrees.sum() == 0:
        # Empty network
        return csr_matrix((n_rows, n_cols), dtype=np.float64)
    
    # Create stub lists (one for each edge endpoint)
    row_stubs = np.repeat(np.arange(n_rows), row_degrees)
    col_stubs = np.repeat(np.arange(n_cols), col_degrees)
    
    # Randomly match row stubs to column stubs
    rng = default_rng(random_state)
    rng.shuffle(col_stubs)
    
    # Create new bipartite adjacency matrix
    data = np.ones(len(row_stubs), dtype=np.float64)
    randomized = csr_matrix((data, (row_stubs, col_stubs)), shape=(n_rows, n_cols))
    
    # Remove potential self-loops (shouldn't exist in bipartite but check)
    randomized.setdiag(0)
    randomized.eliminate_zeros()
    
    return randomized


def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S %Z")


def safe_dir(path: str) -> None:
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)


def build_biadjacency(edges: pd.DataFrame, src_col: str, tgt_col: str) -> Tuple[csr_matrix, List[str], List[str]]:
    """Build bipartite adjacency matrix, consistent with 03a_modularity_analysis_sk.py"""
    edges = edges[[src_col, tgt_col]].drop_duplicates().copy()
    rows = sorted(edges[src_col].unique())
    cols = sorted(edges[tgt_col].unique())
    r_index = {v: i for i, v in enumerate(rows)}
    c_index = {v: i for i, v in enumerate(cols)}
    r_idx = edges[src_col].map(r_index).to_numpy()
    c_idx = edges[tgt_col].map(c_index).to_numpy()
    data = np.ones(len(edges), dtype=np.float64)
    biadj = csr_matrix((data, (r_idx, c_idx)), shape=(len(rows), len(cols)))
    return biadj, rows, cols


def build_incidence_matrix(adj_production: csr_matrix, adj_utilization: csr_matrix) -> np.ndarray:
    """
    Build incidence matrix with deterministic sorting consistent with R version
    """
    # Ensure both matrices have the same shape
    if adj_production.shape != adj_utilization.shape:
        raise ValueError(f"Matrix shapes don't match: production {adj_production.shape} vs utilization {adj_utilization.shape}")
    
    # For utilization network, we need to transpose to match production orientation
    # This creates a bipartite incidence matrix where rows are agents and columns are pyoverdines
    incidence_matrix = adj_production + adj_utilization.transpose()
    incidence_matrix = (incidence_matrix > 0).astype(int)
    
    # Convert to dense matrix for sorting
    incidence_dense = incidence_matrix.toarray()
    
    # Sort by degree with deterministic tie-breaking
    row_sums = np.sum(incidence_dense, axis=1)
    col_sums = np.sum(incidence_dense, axis=0)
    row_order = np.lexsort((np.arange(len(row_sums)), -row_sums))
    col_order = np.lexsort((np.arange(len(col_sums)), -col_sums))
    
    return incidence_dense[row_order, :][:, col_order]


def compute_nodf(incidence_matrix: np.ndarray) -> Dict[str, float]:
    """
    Compute NODF nestedness metric
    Based on Almeida-Neto et al. 2008 algorithm
    """
    n_rows, n_cols = incidence_matrix.shape
    
    if n_rows < 2 or n_cols < 2:
        return {'nodf_total': 0.0, 'nodf_rows': 0.0, 'nodf_cols': 0.0}
    
    # Calculate row degrees
    row_degrees = np.sum(incidence_matrix, axis=1)
    col_degrees = np.sum(incidence_matrix, axis=0)
    
    # Calculate row NODF
    nodf_rows = 0.0
    row_pairs = 0
    for i in range(n_rows):
        for j in range(i + 1, n_rows):
            di, dj = row_degrees[i], row_degrees[j]
            # skip equal-degree pairs and zero-denominator cases to avoid 0/0
            if di == dj:
                continue
            denom = min(di, dj)
            if denom <= 0:
                continue
            overlap = int(np.sum(incidence_matrix[i] & incidence_matrix[j]))
            nodf_rows += (overlap / denom) * 100
            row_pairs += 1
    
    nodf_rows = nodf_rows / row_pairs if row_pairs > 0 else 0.0
    
    # Calculate column NODF
    nodf_cols = 0.0
    col_pairs = 0
    for i in range(n_cols):
        for j in range(i + 1, n_cols):
            di, dj = col_degrees[i], col_degrees[j]
            # skip equal-degree pairs and zero-denominator cases
            if di == dj:
                continue
            denom = min(di, dj)
            if denom <= 0:
                continue
            overlap = int(np.sum(incidence_matrix[:, i] & incidence_matrix[:, j]))
            nodf_cols += (overlap / denom) * 100
            col_pairs += 1
    
    nodf_cols = nodf_cols / col_pairs if col_pairs > 0 else 0.0
    
    # Total NODF
    nodf_total = (nodf_rows + nodf_cols) / 2
    
    return {
        'nodf_total': nodf_total,
        'nodf_rows': nodf_rows,
        'nodf_cols': nodf_cols,
        'n_rows': n_rows,
        'n_cols': n_cols,
        'fill_percentage': 100 * np.sum(incidence_matrix) / (n_rows * n_cols)
    }



def _compute_null_summary(nodf_obs: float, nodf_null: np.ndarray) -> Dict[str, float]:
    # filter out NaN/Inf from null values
    valid = nodf_null[np.isfinite(nodf_null)] if nodf_null.size else np.array([], dtype=float)
    if valid.size == 0:
        return {
            'NODF_obs': nodf_obs,
            'NODF_null_mean': float('nan'),
            'NODF_null_sd': float('nan'),
            'z_score': float('nan'),
            'p_value': float('nan'),
            'n_valid_nulls': 0
        }
    mean_nodf = float(np.mean(valid))
    sd_nodf = float(np.std(valid, ddof=1)) if valid.size > 1 else float('nan')
    if np.isfinite(sd_nodf) and sd_nodf > 0.0:
        z_score = float((nodf_obs - mean_nodf) / sd_nodf)
    else:
        z_score = float('nan')
    
    # Two-tailed test for nestedness analysis
    if nodf_obs >= mean_nodf:
        p_value = float(np.mean(valid >= nodf_obs))
    else:
        p_value = float(np.mean(valid <= nodf_obs))
    p_value = min(p_value * 2, 1.0)  # Two-tailed adjustment, maximum is 1.0
    
    return {
        'NODF_obs': nodf_obs,
        'NODF_null_mean': mean_nodf,
        'NODF_null_sd': sd_nodf,
        'z_score': z_score,
        'p_value': p_value,
        'n_valid_nulls': int(valid.size)
    }


def analyze_within_module_nestedness(incidence_matrix: np.ndarray,
                                   module_assignments: pd.DataFrame,
                                   row_nodes: List[str],
                                   col_nodes: List[str]) -> Dict:
    """Analyze nestedness within each module"""
    
    # Separate row and column module assignments
    fg_modules = module_assignments[module_assignments['node_type'] == 'FG']
    pyo_modules = module_assignments[module_assignments['node_type'] == 'PYOV']
    
    # Create node to module mapping
    fg_module_map = dict(zip(fg_modules['node_id'], fg_modules['module_id']))
    pyo_module_map = dict(zip(pyo_modules['node_id'], pyo_modules['module_id']))
    
    # Get all modules
    all_modules = set(fg_module_map.values()) | set(pyo_module_map.values())
    
    module_results = []
    
    for module_id in all_modules:
        # Get rows and columns for this module
        module_rows = [i for i, node in enumerate(row_nodes) if node in fg_module_map and fg_module_map[node] == module_id]
        module_cols = [j for j, node in enumerate(col_nodes) if node in pyo_module_map and pyo_module_map[node] == module_id]
        
        if len(module_rows) < 2 or len(module_cols) < 2:
            continue  # Too small to compute NODF
        
        # Extract submatrix
        submatrix = incidence_matrix[np.ix_(module_rows, module_cols)]
        
        # Compute NODF
        nodf_result = compute_nodf(submatrix)
        
        module_results.append({
            'module_id': module_id,
            'n_rows': len(module_rows),
            'n_cols': len(module_cols),
            'nodf_total': nodf_result['nodf_total'],
            'nodf_rows': nodf_result['nodf_rows'],
            'nodf_cols': nodf_result['nodf_cols'],
            'fill_percentage': nodf_result['fill_percentage']
        })
    
    return {
        'module_results': module_results,
        'n_modules': len(module_results)
    }


def plot_nestedness_matrix(incidence_matrix: np.ndarray,
                          title: str,
                          nodf_value: float = None,
                          output_path: str = None) -> None:
    """Generate nestedness matrix heatmap"""
    
    plt.figure(figsize=(12, 8))
    
    # Create heatmap
    sns.heatmap(incidence_matrix,
                cmap='binary',
                cbar=False,
                square=True,
                xticklabels=False,
                yticklabels=False)
    
    if nodf_value is not None:
        title = f"{title}\n(NODF = {nodf_value:.2f})"
    
    plt.title(title, fontsize=14, pad=20)
    plt.xlabel('Pyoverdine', fontsize=12)
    plt.ylabel('Agents', fontsize=12)
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        png_path = os.path.splitext(output_path)[0] + '.png'
        plt.savefig(png_path, dpi=300, bbox_inches='tight')
    
    plt.close()


def plot_module_nestedness_matrix(incidence_matrix: np.ndarray,
                                 module_assignments: pd.DataFrame,
                                 row_nodes: List[str],
                                 col_nodes: List[str],
                                 title: str,
                                 nodf_value: float = None,
                                 output_path: str = None) -> None:
    """Generate module-sorted nestedness matrix heatmap with module annotations.
    Annotation bars are drawn directly on the same axis as the heatmap to remove
    the large gap caused by a separate annotation subplot.
    """
    from matplotlib.patches import Rectangle  # local import to avoid top-level changes

    # Extract functional group module assignments
    fg_modules = module_assignments[module_assignments['node_type'] == 'FG']
    fg_module_map = dict(zip(fg_modules['node_id'], fg_modules['module_id']))

    # Get module assignments for each row node in the correct order
    row_modules = []
    for i, node in enumerate(row_nodes):
        if node in fg_module_map:
            row_modules.append(fg_module_map[node])
        else:
            row_modules.append(-1)  # Unknown/unassigned module

    # Group rows by module (preserving original internal degree order within each module)
    module_order = sorted(list(set(row_modules)), key=lambda x: (x == -1, x))
    row_indices_by_module = []

    for module_id in module_order:
        module_indices = [i for i, m in enumerate(row_modules) if m == module_id]

        # Sort within module by degree to maintain nestedness visualization
        module_degrees = [np.sum(incidence_matrix[i, :]) for i in module_indices]
        sorted_indices = [module_indices[i] for i in np.lexsort((np.arange(len(module_degrees)), -np.array(module_degrees)))]
        row_indices_by_module.extend(sorted_indices)

    # Reorder matrix and module labels
    module_sorted_matrix = incidence_matrix[row_indices_by_module, :]
    sorted_row_modules = [row_modules[i] for i in row_indices_by_module]

    # Create module annotation colors
    unique_modules = sorted(list(set(sorted_row_modules)))
    n_modules = len(unique_modules)
    colors = plt.cm.Set3(np.linspace(0, 1, max(n_modules, 1)))
    module_colors = {module_id: colors[i] for i, module_id in enumerate(unique_modules)}

    # Single axis figure (annotations drawn on same axis as heatmap)
    fig, ax_heatmap = plt.subplots(figsize=(14, 9))

    # Create heatmap
    sns.heatmap(module_sorted_matrix,
                cmap='binary',
                cbar=False,
                square=True,
                xticklabels=False,
                yticklabels=False,
                ax=ax_heatmap)

    if nodf_value is not None:
        title = f"{title}\n(NODF = {nodf_value:.2f})"

    ax_heatmap.set_title(title, fontsize=14, pad=20)
    ax_heatmap.set_xlabel('Pyoverdine', fontsize=12)
    ax_heatmap.set_ylabel('', fontsize=12)  # Keep y-label empty; annotations are visual

    # Determine geometry
    n_rows = module_sorted_matrix.shape[0]
    n_cols = module_sorted_matrix.shape[1]

    # Bar width in data coordinates (small fraction of heatmap width)
    bar_width = max(0.5, n_cols * 0.03)

    # Obtain current y-limits to place annotation bars correctly regardless of orientation
    y0, y1 = ax_heatmap.get_ylim()
    y_min, y_max = min(y0, y1), max(y0, y1)
    # y coordinates in heatmap go from 0..n_rows; start at y_min
    current_pos = 0

    for module_id in unique_modules:
        module_size = sorted_row_modules.count(module_id)
        color = module_colors[module_id]

        # Rectangle positioned to the left of the heatmap (negative x coordinates)
        rect = Rectangle((-bar_width, y_min + current_pos), bar_width, module_size,
                         transform=ax_heatmap.transData, facecolor=color, alpha=0.8,
                         edgecolor='white', linewidth=0.5)
        ax_heatmap.add_patch(rect)

        # Add module label if space allows
        if module_size > 5:
            ax_heatmap.text(-bar_width * 0.05, y_min + current_pos + module_size / 2,
                            f'M{module_id}' if module_id != -1 else 'Unassigned',
                            ha='right', va='center', fontsize=8,
                            bbox=dict(boxstyle="round,pad=0.15", facecolor="white", alpha=0.8),
                            transform=ax_heatmap.transData)

        current_pos += module_size

    # Expand x-limits to make left annotation visible
    ax_heatmap.set_xlim(-bar_width * 1.5, n_cols)
    # Ensure y-limits cover all rows (heatmap may have inverted limits)
    ax_heatmap.set_ylim(y_min, y_max)
    
    # Add left y-axis title for functional groups
    ax_heatmap.text(-bar_width * 2.75, (y_min + y_max) / 2, 'Functional group',
                    rotation=90, ha='center', va='center', fontsize=12,
                    transform=ax_heatmap.transData)
    
    # Remove ticks on the left annotation area
    ax_heatmap.set_xticks([])
    ax_heatmap.set_yticks([])

    # Add legend for modules if not too many
    if n_modules <= 10:
        legend_elements = [Rectangle((0, 0), 1, 1, facecolor=module_colors[mid],
                                     alpha=0.8, label=(f'M{mid}' if mid != -1 else 'Unassigned'))
                           for mid in unique_modules]
        # Place legend just inside the left margin of the axis
        ax_heatmap.legend(handles=legend_elements, loc='center left',
                          bbox_to_anchor=(-0.02, 0.5), fontsize=8)

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        png_path = os.path.splitext(output_path)[0] + '.png'
        plt.savefig(png_path, dpi=300, bbox_inches='tight')
    
    plt.close()


def plot_modularity_nodf_relationship(modularity_data: pd.DataFrame,
                                     nodf_results: Dict,
                                     output_path: str = None) -> None:
    """Generate modularity-nestedness relationship scatter plots"""
    
    if 'module_results' not in nodf_results or not nodf_results['module_results']:
        return
    
    module_data = pd.DataFrame(nodf_results['module_results'])
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Modularity-Nestedness Relationship Analysis', fontsize=16)
    
    # Module size vs NODF
    axes[0, 0].scatter(module_data['n_rows'] + module_data['n_cols'],
                      module_data['nodf_total'], alpha=0.7)
    axes[0, 0].set_xlabel('Module Size (nodes)')
    axes[0, 0].set_ylabel('NODF')
    axes[0, 0].set_title('Module Size vs NODF')
    
    # Calculate correlation
    if len(module_data) > 2:
        corr, p_val = pearsonr(module_data['n_rows'] + module_data['n_cols'],
                              module_data['nodf_total'])
        axes[0, 0].text(0.05, 0.95, f'r = {corr:.3f}\np = {p_val:.3f}',
                       transform=axes[0, 0].transAxes,
                       bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # Fill percentage vs NODF
    axes[0, 1].scatter(module_data['fill_percentage'],
                      module_data['nodf_total'], alpha=0.7)
    axes[0, 1].set_xlabel('Fill Percentage (%)')
    axes[0, 1].set_ylabel('NODF')
    axes[0, 1].set_title('Fill Percentage vs NODF')
    
    if len(module_data) > 2:
        corr, p_val = pearsonr(module_data['fill_percentage'],
                              module_data['nodf_total'])
        axes[0, 1].text(0.05, 0.95, f'r = {corr:.3f}\np = {p_val:.3f}',
                       transform=axes[0, 1].transAxes,
                       bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # Rows vs columns
    axes[1, 0].scatter(module_data['n_rows'],
                      module_data['n_cols'],
                      c=module_data['nodf_total'],
                      cmap='viridis', alpha=0.7)
    axes[1, 0].set_xlabel('Number of Rows')
    axes[1, 0].set_ylabel('Number of Columns')
    axes[1, 0].set_title('Module Structure (colored by NODF)')
    plt.colorbar(axes[1, 0].collections[0], ax=axes[1, 0], label='NODF')
    
    # NODF distribution
    axes[1, 1].hist(module_data['nodf_total'], bins=10, alpha=0.7, edgecolor='black')
    axes[1, 1].set_xlabel('NODF')
    axes[1, 1].set_ylabel('Frequency')
    axes[1, 1].set_title('NODF Distribution Across Modules')
    axes[1, 1].axvline(np.mean(module_data['nodf_total']),
                       color='red', linestyle='--',
                       label=f'Mean: {np.mean(module_data["nodf_total"]):.2f}')
    axes[1, 1].legend()
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        png_path = os.path.splitext(output_path)[0] + '.png'
        plt.savefig(png_path, dpi=300, bbox_inches='tight')
    
    plt.close()


def analyze_layer_nodf(layer_name: str,
                      edges: pd.DataFrame,
                      src_col: str,
                      tgt_col: str,
                      outdir: str,
                      figures_dir: str,
                      modularity_assignments: pd.DataFrame,
                      random_state: Optional[int] = None,
                      n_null_replicates: int = 100) -> Dict:
    """Analyze NODF nestedness for a single layer"""
    
    print(f"[{ts()}] Layer: {layer_name}")
    
    # Build bipartite adjacency matrix
    biadj, rows, cols = build_biadjacency(edges, src_col, tgt_col)
    
    # For NODF analysis, we work directly with the bipartite adjacency matrix
    # Convert to incidence matrix (binary matrix)
    incidence_matrix = (biadj > 0).astype(int).toarray()
    
    # Sort by degree with deterministic tie-breaking (consistent with R version)
    row_sums = np.sum(incidence_matrix, axis=1)
    col_sums = np.sum(incidence_matrix, axis=0)
    row_order = np.lexsort((np.arange(len(row_sums)), -row_sums))
    col_order = np.lexsort((np.arange(len(col_sums)), -col_sums))
    
    incidence_matrix = incidence_matrix[row_order, :][:, col_order]
    
    # Update row and column node lists to match the sorted order
    rows = [rows[i] for i in row_order]
    cols = [cols[j] for j in col_order]
    
    print(f"Network: {len(rows)} x {len(cols)} matrix, {np.sum(incidence_matrix)} edges "
          f"({np.sum(incidence_matrix) / (len(rows) * len(cols)) * 100:.1f}% fill)")
    
    # Compute observed NODF
    nodf_observed = compute_nodf(incidence_matrix)
    print(f"[{ts()}] {layer_name} | Observed NODF: {nodf_observed['nodf_total']:.3f}")
    
    # Null model analysis
    null_nodf_values = []
    valid_nulls = 0
    
    if n_null_replicates > 0:
        print(f"[{ts()}] {layer_name} | Running {n_null_replicates} null models...")
        
        master_rng = default_rng(random_state)
        
        for null_idx in tqdm(range(n_null_replicates), desc=f'{layer_name} nulls'):
            config_seed = None if random_state is None else int(master_rng.integers(0, 2**31))
            
            try:
                # Generate null network
                null_biadj = _bipartite_config_model(biadj, random_state=config_seed)
                
                # Build null incidence matrix (same process as observed)
                null_incidence = (null_biadj > 0).astype(int).toarray()
                
                # Apply same sorting as observed matrix
                null_row_sums = np.sum(null_incidence, axis=1)
                null_col_sums = np.sum(null_incidence, axis=0)
                null_row_order = np.lexsort((np.arange(len(null_row_sums)), -null_row_sums))
                null_col_order = np.lexsort((np.arange(len(null_col_sums)), -null_col_sums))
                null_incidence = null_incidence[null_row_order, :][:, null_col_order]
                
                # Compute null NODF
                null_nodf = compute_nodf(null_incidence)
                null_nodf_values.append(null_nodf['nodf_total'])
                valid_nulls += 1
                
            except Exception as exc:
                print(f"[{ts()}] {layer_name} | Null replicate {null_idx + 1} failed: {exc}")
    
    # Compute null model summary
    null_summary = _compute_null_summary(nodf_observed['nodf_total'],
                                      np.array(null_nodf_values))
    
    print(f"[{ts()}] {layer_name} | Null NODF: {null_summary['NODF_null_mean']:.3f} ± {null_summary['NODF_null_sd']:.3f}")
    print(f"[{ts()}] {layer_name} | Z-score: {null_summary['z_score']:.3f}, p-value: {null_summary['p_value']:.4f}")
    
    # Within-module nestedness analysis
    module_analysis = analyze_within_module_nestedness(incidence_matrix,
                                                     modularity_assignments,
                                                     rows, cols)
    
    print(f"[{ts()}] {layer_name} | Analyzed {module_analysis['n_modules']} modules for nestedness")
    
    # Compute module-level summary statistics (mean, sd) and compare to overall NODF
    module_summary = {
        'mean_module_nodf': float('nan'),
        'sd_module_nodf': float('nan'),
        'z_vs_module_mean': float('nan'),
        'p_vs_module_mean': float('nan'),
        'n_modules_for_stats': 0
    }
    if module_analysis and module_analysis.get('module_results'):
        module_vals = np.array([m['nodf_total'] for m in module_analysis['module_results']], dtype=float)
        if module_vals.size > 0:
            mean_module = float(np.mean(module_vals))
            sd_module = float(np.std(module_vals, ddof=1)) if module_vals.size > 1 else float('nan')
            module_summary['mean_module_nodf'] = mean_module
            module_summary['sd_module_nodf'] = sd_module
            module_summary['n_modules_for_stats'] = int(module_vals.size)
            # z-score comparing observed overall NODF to module mean (using module SD)
            if np.isfinite(sd_module) and sd_module > 0.0:
                z_mod = float((nodf_observed['nodf_total'] - mean_module) / sd_module)
                p_mod = float(min(2.0 * (1.0 - norm.cdf(abs(z_mod))), 1.0))
                module_summary['z_vs_module_mean'] = z_mod
                module_summary['p_vs_module_mean'] = p_mod

    # Save results
    results = {
        'layer': layer_name,
        'nodf_observed': nodf_observed,
        'null_summary': null_summary,
        'module_analysis': module_analysis,
        'module_summary': module_summary,
        'null_nodf_values': null_nodf_values
    }
    
    # Save data files
    safe_dir(outdir)
    
    # Observed NODF results
    obs_df = pd.DataFrame([{
        'layer': layer_name,
        'nodf_total': nodf_observed['nodf_total'],
        'nodf_rows': nodf_observed['nodf_rows'],
        'nodf_cols': nodf_observed['nodf_cols'],
        'n_rows': nodf_observed['n_rows'],
        'n_cols': nodf_observed['n_cols'],
        'fill_percentage': nodf_observed['fill_percentage']
    }])
    obs_path = os.path.join(outdir, f"nodf_results_{layer_name}.csv")
    obs_df.to_csv(obs_path, index=False)
    
    # Null model results
    if null_nodf_values:
        null_df = pd.DataFrame({
            'replicate': range(1, len(null_nodf_values) + 1),
            'nodf': null_nodf_values
        })
        null_path = os.path.join(outdir, f"null_replicates_{layer_name}.csv")
        null_df.to_csv(null_path, index=False)
        
        # Null model summary
        summary_df = pd.DataFrame([null_summary])
        summary_path = os.path.join(outdir, f"null_summary_{layer_name}.csv")
        summary_df.to_csv(summary_path, index=False)
    
    # Module analysis results
    if module_analysis['module_results']:
        module_df = pd.DataFrame(module_analysis['module_results'])
        module_path = os.path.join(outdir, f"module_nodf_analysis_{layer_name}.csv")
        module_df.to_csv(module_path, index=False)
    
    # Generate visualizations
    safe_dir(figures_dir)
    
    # Nestedness matrix heatmap
    matrix_path = os.path.join(figures_dir, f"04c_nodf_matrix_{layer_name}.pdf")
    plot_nestedness_matrix(incidence_matrix,
                          f"{layer_name.title()} Network Nestedness",
                          nodf_observed['nodf_total'],
                          matrix_path)
    
    # Module-sorted nestedness matrix plots
    module_matrix_path = os.path.join(figures_dir, f"04c_module_nestedness_matrix_{layer_name}.pdf")
    plot_module_nestedness_matrix(incidence_matrix,
                                modularity_assignments,
                                rows, cols,
                                f"{layer_name.title()} Network Nestedness (Module Sorted)",
                                nodf_observed['nodf_total'],
                                module_matrix_path)
    
    # Modularity-nestedness relationship plot
    if module_analysis['module_results']:
        relationship_path = os.path.join(figures_dir, f"04c_modularity_nodf_relationship_{layer_name}.pdf")
        plot_modularity_nodf_relationship(modularity_assignments,
                                        module_analysis,
                                        relationship_path)
    
    results.update({
        'obs_path': obs_path,
        'null_path': os.path.join(outdir, f"null_replicates_{layer_name}.csv") if null_nodf_values else None,
        'summary_path': os.path.join(outdir, f"null_summary_{layer_name}.csv"),
        'module_path': os.path.join(outdir, f"module_nodf_analysis_{layer_name}.csv") if module_analysis['module_results'] else None,
        'matrix_path': matrix_path,
        'module_matrix_path': module_matrix_path,
        'relationship_path': os.path.join(figures_dir, f"04c_modularity_nodf_relationship_{layer_name}.pdf") if module_analysis['module_results'] else None
    })
    
    return results


def main():
    cfg = {
        'edges_file': 'data/interim/edges_functional_groups_conservative.csv',
        'production_modularity_file': 'results/phase_04/modularity/module_assignments_production_sk.csv',
        'utilization_modularity_file': 'results/phase_04/modularity/module_assignments_utilization_sk.csv',
        'output_dir': 'results/phase_04/nodf',
        'figures_dir': 'figures/network_topology/step04_nestedness',
        'log_dir': 'docs/phase_04/logs',
        'random_state': 42,
        'n_null_replicates': 100,
    }
    
    safe_dir(cfg['output_dir'])
    safe_dir(cfg['figures_dir'])
    safe_dir(cfg['log_dir'])
    
    print("=== Phase 04 Step 04c: NODF Nestedness Analysis (Modularity Integration) ===")
    print(f"Timestamp: {ts()}")
    print(f"Null replicates per layer: {cfg['n_null_replicates']}")
    print()
    
    # Check required files
    required_files = [
        cfg['edges_file'],
        cfg['production_modularity_file'],
        cfg['utilization_modularity_file']
    ]
    
    for file_path in required_files:
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"Missing required file: {file_path}")
    
    # Load data
    print("Loading data...")
    edges = pd.read_csv(cfg['edges_file'])
    prod_modularity = pd.read_csv(cfg['production_modularity_file'])
    util_modularity = pd.read_csv(cfg['utilization_modularity_file'])
    
    print(f"Loaded {len(edges)} edges")
    print(f"Loaded {len(prod_modularity)} production module assignments")
    print(f"Loaded {len(util_modularity)} utilization module assignments")
    
    # Prepare edge data
    prod_edges = edges.query("edge_type == 'production'")[["source", "target"]]
    util_edges = edges.query("edge_type == 'utilization'")[["source", "target"]]
    
    if prod_edges.empty or util_edges.empty:
        raise RuntimeError("Edges missing for one or both layers")
    
    # Define layers
    layers = {
        'production': (prod_edges, 'source', 'target', prod_modularity),
        'utilization': (util_edges.rename(columns={'source': 'col', 'target': 'row'}).rename(columns={'row': 'source', 'col': 'target'}),
                       'source', 'target', util_modularity)
    }
    
    # Analyze each layer
    results = {}
    for layer_idx, (layer_name, (edge_data, src_col, tgt_col, modularity_data)) in enumerate(layers.items()):
        layer_seed = cfg['random_state'] + layer_idx * 10000000 if cfg['random_state'] else None
        results[layer_name] = analyze_layer_nodf(
            layer_name,
            edge_data,
            src_col,
            tgt_col,
            cfg['output_dir'],
            cfg['figures_dir'],
            modularity_data,
            random_state=layer_seed,
            n_null_replicates=cfg['n_null_replicates']
        )
    
    # Generate summary report
    summary_lines = [
        'Phase 04 Step 04c: NODF Nestedness Analysis Summary (Modularity Integration)',
        f"Timestamp: {ts()}",
        ''
    ]
    
    for lname, res in results.items():
        obs = res['nodf_observed']
        null = res['null_summary']
        mod = res['module_analysis']
        
        # Paths (PDF + PNG)
        matrix_pdf = res.get('matrix_path')
        matrix_png = os.path.splitext(matrix_pdf)[0] + '.png' if matrix_pdf else None
        module_pdf = res.get('module_matrix_path')
        module_png = os.path.splitext(module_pdf)[0] + '.png' if module_pdf else None
        relationship_pdf = res.get('relationship_path')
        relationship_png = os.path.splitext(relationship_pdf)[0] + '.png' if relationship_pdf else None
        
        summary_lines.extend([
            f"{lname.title()} Layer:",
            f"  Network: {obs['n_rows']} × {obs['n_cols']} ({obs['fill_percentage']:.1f}% fill)",
            f"  Observed NODF: {obs['nodf_total']:.3f} (rows: {obs['nodf_rows']:.3f}, cols: {obs['nodf_cols']:.3f})",
            f"  Null NODF: {null['NODF_null_mean']:.3f} ± {null['NODF_null_sd']:.3f}",
            f"  Z-score: {null['z_score']:.3f} | p-value: {null['p_value']:.4f}",
            f"  Modules analyzed: {mod['n_modules']}",
            f"  Results saved (observed metrics): {res['obs_path']}",
            f"  Nestedness matrix: {matrix_pdf}" + (f" (png: {matrix_png})" if matrix_png else ""),
            f"  Module-sorted matrix: {module_pdf}" + (f" (png: {module_png})" if module_png else ""),
            ''
        ])
        
        # Relationship plot (if generated)
        if relationship_pdf:
            summary_lines.append(f"  Modularity-NODF relationship plot: {relationship_pdf}" + (f" (png: {relationship_png})" if relationship_png else ""))
            summary_lines.append('')
        
        # Detailed per-module summary so conclusions can be drawn without images
        if mod and mod.get('module_results'):
            summary_lines.append("  Module breakdown (ModuleID: rows, cols, NODF, Fill%):")
            for m in mod['module_results']:
                summary_lines.append(f"    M{m['module_id']}: {m['n_rows']} rows, {m['n_cols']} cols, NODF={m['nodf_total']:.2f}, Fill%={m['fill_percentage']:.1f}")
            # Provide path to module analysis CSV if available
            if res.get('module_path'):
                summary_lines.append(f"  Module analysis CSV: {res['module_path']}")
            summary_lines.append('')
            # Add aggregate module-level statistics for easy interpretation without figures
            module_summary = res.get('module_summary', {})
            if module_summary:
                mean_m = module_summary.get('mean_module_nodf', float('nan'))
                sd_m = module_summary.get('sd_module_nodf', float('nan'))
                z_m = module_summary.get('z_vs_module_mean', float('nan'))
                p_m = module_summary.get('p_vs_module_mean', float('nan'))
                n_mods = module_summary.get('n_modules_for_stats', 0)
                summary_lines.append("  Module-level summary:")
                summary_lines.append(f"    Number of modules considered: {n_mods}")
                summary_lines.append(f"    Module NODF mean ± sd: {mean_m:.3f} ± {sd_m:.3f}")
                summary_lines.append(f"    Observed overall NODF vs module mean: z = {z_m:.3f}, p = {p_m:.4f}")
                summary_lines.append('')
        
        # Add plot generation details (parameters) to the summary to allow interpretation without images
        summary_lines.append("  Plot generation settings:")
        summary_lines.append("    - Heatmap colormap: 'binary'")
        summary_lines.append("    - Heatmap square: True; ticks suppressed")
        summary_lines.append("    - Module annotations: drawn on same axis as heatmap (colored rectangles)")
        summary_lines.append("    - Module color map: Matplotlib Set3")
        summary_lines.append("    - Module annotation bar width: ~3% of heatmap width (adaptive)")
        summary_lines.append('')
    
    # Save summary
    summary_file = os.path.join(cfg['log_dir'], '04c_nodf_analysis_summary.txt')
    with open(summary_file, 'w') as f:
        f.write('\n'.join(summary_lines))
    
    print('\n'.join(summary_lines))
    print(f"\nSummary saved to: {summary_file}")
    print("=== NODF Analysis Complete ===")


if __name__ == '__main__':
    main()
