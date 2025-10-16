#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple


def load_module_data(production_path: str, utilization_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load module assignment data for both layers."""
    prod_df = pd.read_csv(production_path)
    util_df = pd.read_csv(utilization_path)
    return prod_df, util_df


def create_bipartite_matrix(df: pd.DataFrame) -> Tuple[np.ndarray, List[str], List[str]]:
    """Create bipartite matrix from module assignments (FG vs PYOV)."""
    fg_nodes = df[df['node_type'] == 'FG'].sort_values('node_id')
    pyo_nodes = df[df['node_type'] == 'PYOV'].sort_values('node_id')
    
    # Create matrix: FG nodes (rows) x PYO nodes (columns)
    # Value 1 if they share the same module, 0 otherwise
    fg_modules = fg_nodes['module_id'].values
    pyo_modules = pyo_nodes['module_id'].values
    
    matrix = np.zeros((len(fg_nodes), len(pyo_nodes)))
    for i, fg_module in enumerate(fg_modules):
        for j, pyo_module in enumerate(pyo_modules):
            matrix[i, j] = 1 if fg_module == pyo_module else 0
    
    return matrix, fg_nodes['node_id'].tolist(), pyo_nodes['node_id'].tolist()


def create_module_cooccurrence_matrix(df: pd.DataFrame) -> Tuple[np.ndarray, List[int]]:
    """Create module co-occurrence matrix showing weights between modules."""
    # Get all unique modules
    modules = sorted(df['module_id'].unique())
    
    # Count co-occurrences between modules
    n_modules = len(modules)
    matrix = np.zeros((n_modules, n_modules))
    
    # Count nodes in each module
    fg_nodes = df[df['node_type'] == 'FG']
    pyo_nodes = df[df['node_type'] == 'PYOV']
    
    for i, module_i in enumerate(modules):
        for j, module_j in enumerate(modules):
            if i == j:
                # Diagonal: count nodes in this module
                n_nodes = len(df[df['module_id'] == module_i])
                matrix[i, j] = n_nodes
            else:
                # Off-diagonal: count if modules coexist (for bipartite structure)
                # This could be enhanced with actual network connections
                matrix[i, j] = len(fg_nodes[fg_nodes['module_id'] == module_i]) * len(pyo_nodes[pyo_nodes['module_id'] == module_j]) / 100
    
    return matrix, modules


def plot_bipartite_heatmap(matrix: np.ndarray, fg_labels: List[str], pyo_labels: List[str], 
                          title: str, output_path: str, max_display: int = 50):
    """Plot bipartite heatmap with reasonable label limits."""
    if len(fg_labels) > max_display or len(pyo_labels) > max_display:
        print(f"Skipping {title} - too many nodes ({len(fg_labels)} x {len(pyo_labels)})")
        return
    
    plt.figure(figsize=(12, 8))
    sns.heatmap(matrix, 
                xticklabels=pyo_labels, 
                yticklabels=fg_labels,
                cmap='Reds', 
                cbar_kws={'label': 'Same Module (1) / Different Module (0)'},
                square=True)
    plt.title(title)
    plt.xlabel('PYO Nodes')
    plt.ylabel('FG Nodes')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_module_heatmap(matrix: np.ndarray, module_labels: List[int], 
                       title: str, output_path: str):
    """Plot module-level heatmap."""
    plt.figure(figsize=(10, 8))
    sns.heatmap(matrix, 
                xticklabels=[f'Module {m}' for m in module_labels], 
                yticklabels=[f'Module {m}' for m in module_labels],
                cmap='Reds', 
                cbar_kws={'label': 'Connection Weight'},
                square=True)
    plt.title(title)
    plt.xlabel('Modules')
    plt.ylabel('Modules')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_module_distribution(prod_df: pd.DataFrame, util_df: pd.DataFrame, output_dir: str):
    """Plot module size distributions."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Production layer plots
    prod_fg_counts = prod_df[prod_df['node_type'] == 'FG']['module_id'].value_counts().sort_index()
    prod_pyo_counts = prod_df[prod_df['node_type'] == 'PYOV']['module_id'].value_counts().sort_index()
    
    axes[0, 0].bar(range(len(prod_fg_counts)), prod_fg_counts.values)
    axes[0, 0].set_title('Production Layer - FG Nodes per Module')
    axes[0, 0].set_xlabel('Module ID')
    axes[0, 0].set_ylabel('Number of FG Nodes')
    axes[0, 0].set_xticks(range(len(prod_fg_counts)))
    axes[0, 0].set_xticklabels([f'M{m}' for m in prod_fg_counts.index], rotation=45)
    
    axes[0, 1].bar(range(len(prod_pyo_counts)), prod_pyo_counts.values, color='orange')
    axes[0, 1].set_title('Production Layer - PYO Nodes per Module')
    axes[0, 1].set_xlabel('Module ID')
    axes[0, 1].set_ylabel('Number of PYO Nodes')
    axes[0, 1].set_xticks(range(len(prod_pyo_counts)))
    axes[0, 1].set_xticklabels([f'M{m}' for m in prod_pyo_counts.index], rotation=45)
    
    # Utilization layer plots
    util_fg_counts = util_df[util_df['node_type'] == 'FG']['module_id'].value_counts().sort_index()
    util_pyo_counts = util_df[util_df['node_type'] == 'PYOV']['module_id'].value_counts().sort_index()
    
    axes[1, 0].bar(range(len(util_fg_counts)), util_fg_counts.values)
    axes[1, 0].set_title('Utilization Layer - FG Nodes per Module')
    axes[1, 0].set_xlabel('Module ID')
    axes[1, 0].set_ylabel('Number of FG Nodes')
    axes[1, 0].set_xticks(range(len(util_fg_counts)))
    axes[1, 0].set_xticklabels([f'M{m}' for m in util_fg_counts.index], rotation=45)
    
    axes[1, 1].bar(range(len(util_pyo_counts)), util_pyo_counts.values, color='orange')
    axes[1, 1].set_title('Utilization Layer - PYO Nodes per Module')
    axes[1, 1].set_xlabel('Module ID')
    axes[1, 1].set_ylabel('Number of PYO Nodes')
    axes[1, 1].set_xticks(range(len(util_pyo_counts)))
    axes[1, 1].set_xticklabels([f'M{m}' for m in util_pyo_counts.index], rotation=45)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'module_distribution.png'), dpi=300, bbox_inches='tight')
    plt.close()


def plot_module_comparison(prod_df: pd.DataFrame, util_df: pd.DataFrame, output_dir: str):
    """Plot comparison of module assignments between layers."""
    # Get common nodes
    prod_nodes = prod_df[['node_id', 'module_id', 'node_type']].rename(columns={'module_id': 'prod_module'})
    util_nodes = util_df[['node_id', 'module_id', 'node_type']].rename(columns={'module_id': 'util_module'})
    
    combined = prod_nodes.merge(util_nodes, on=['node_id', 'node_type'], how='inner')
    
    if len(combined) == 0:
        print("No common nodes found between production and utilization layers")
        return
    
    # Create cross-tabulation
    cross_tab = pd.crosstab(combined['prod_module'], combined['util_module'])
    
    # Plot heatmap of module assignments
    plt.figure(figsize=(10, 8))
    sns.heatmap(cross_tab, 
                cmap='Blues', 
                cbar_kws={'label': 'Number of Nodes'},
                annot=True, fmt='d', 
                square=True)
    plt.title('Module Assignment Comparison: Production vs Utilization')
    plt.xlabel('Utilization Module')
    plt.ylabel('Production Module')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'module_comparison_heatmap.png'), dpi=300, bbox_inches='tight')
    plt.close()


def main():
    # Configuration
    production_file = 'results/phase_04/modularity/module_assignments_production_sk.csv'
    utilization_file = 'results/phase_04/modularity/module_assignments_utilization_sk.csv'
    output_dir = 'figures/network_topology/step03_modularity'
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    print("Loading module assignment data...")
    prod_df, util_df = load_module_data(production_file, utilization_file)
    
    print(f"Production layer: {len(prod_df)} nodes ({prod_df['node_type'].value_counts().to_dict()})")
    print(f"Utilization layer: {len(util_df)} nodes ({util_df['node_type'].value_counts().to_dict()})")
    
    # Plot module distributions
    print("Creating module distribution plots...")
    plot_module_distribution(prod_df, util_df, output_dir)
    
    # Create bipartite matrices and heatmaps
    print("Creating bipartite heatmaps...")
    prod_matrix, prod_fg_labels, prod_pyo_labels = create_bipartite_matrix(prod_df)
    util_matrix, util_fg_labels, util_pyo_labels = create_bipartite_matrix(util_df)
    
    plot_bipartite_heatmap(prod_matrix, prod_fg_labels, prod_pyo_labels,
                          'Production Layer - FG-PYO Module Co-membership',
                          os.path.join(output_dir, 'production_bipartite_heatmap.png'))
    
    plot_bipartite_heatmap(util_matrix, util_fg_labels, util_pyo_labels,
                          'Utilization Layer - FG-PYO Module Co-membership',
                          os.path.join(output_dir, 'utilization_bipartite_heatmap.png'))
    
    # Create module-level heatmaps
    print("Creating module-level heatmaps...")
    prod_module_matrix, prod_modules = create_module_cooccurrence_matrix(prod_df)
    util_module_matrix, util_modules = create_module_cooccurrence_matrix(util_df)
    
    plot_module_heatmap(prod_module_matrix, prod_modules,
                       'Production Layer - Module Connectivity',
                       os.path.join(output_dir, 'production_module_heatmap.png'))
    
    plot_module_heatmap(util_module_matrix, util_modules,
                       'Utilization Layer - Module Connectivity',
                       os.path.join(output_dir, 'utilization_module_heatmap.png'))
    
    # Plot comparison between layers
    print("Creating module comparison plot...")
    plot_module_comparison(prod_df, util_df, output_dir)
    
    print(f"Visualizations saved to {output_dir}")
    print("Generated files:")
    print("- module_distribution.png")
    print("- production_bipartite_heatmap.png")
    print("- utilization_bipartite_heatmap.png") 
    print("- production_module_heatmap.png")
    print("- utilization_module_heatmap.png")
    print("- module_comparison_heatmap.png")


if __name__ == '__main__':
    main()
