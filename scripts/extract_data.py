#!/usr/bin/env python3
"""
Extract data from H5AD files for web visualization.
Outputs JSON files containing cell coordinates, metadata, and gene expression data.
"""

import json
import os
import sys
from pathlib import Path

import numpy as np
import scanpy as sc

# Configuration
H5AD_FILES = [
    "/home/jiaxi/AG_C00841F3_cellbin.h5ad",
]
OUTPUT_DIR = Path(__file__).parent.parent / "data"

# Metadata columns to extract
METADATA_COLS = ['domain', 'domain.fine', 'subclass', 'cluster']

# Number of top variable genes to include for gene expression
N_TOP_GENES = 100


def extract_dataset(h5ad_path: str, output_dir: Path) -> dict:
    """Extract data from a single H5AD file."""
    print(f"Processing: {h5ad_path}")
    
    adata = sc.read_h5ad(h5ad_path)
    dataset_name = Path(h5ad_path).stem
    
    print(f"  Shape: {adata.shape[0]} cells x {adata.shape[1]} genes")
    
    # Extract cell coordinates
    coords = {
        'x': adata.obs['x'].tolist(),
        'y': adata.obs['y'].tolist()
    }
    
    # Extract metadata columns
    metadata = {}
    available_meta = []
    for col in METADATA_COLS:
        if col in adata.obs.columns:
            # Convert to string to handle categorical data
            values = adata.obs[col].astype(str).tolist()
            metadata[col] = values
            available_meta.append(col)
            
            # Get unique categories and counts
            unique, counts = np.unique(values, return_counts=True)
            print(f"  {col}: {len(unique)} categories")
    
    # Get gene list (all genes)
    gene_names = adata.var['name'].tolist() if 'name' in adata.var.columns else adata.var_names.tolist()
    
    # Select top variable genes for expression data
    # Use a memory-efficient approach for sparse matrices
    print("  Calculating gene variances (memory-efficient)...")
    
    from scipy import sparse
    if sparse.issparse(adata.X):
        # Calculate variance for sparse matrix without converting to dense
        # Var(X) = E[X^2] - E[X]^2
        n_cells = adata.X.shape[0]
        
        # E[X] - mean of each column
        mean_gene = np.array(adata.X.mean(axis=0)).flatten()
        
        # E[X^2] - mean of squared values
        mean_sq = np.array(adata.X.power(2).mean(axis=0)).flatten()
        
        # Variance = E[X^2] - E[X]^2
        gene_vars = mean_sq - np.square(mean_gene)
    else:
        gene_vars = np.var(adata.X, axis=0)
    
    if adata.shape[1] > N_TOP_GENES:
        top_gene_indices = np.argsort(gene_vars)[-N_TOP_GENES:][::-1]
        top_genes = [gene_names[i] for i in top_gene_indices]
    else:
        top_genes = gene_names
        top_gene_indices = list(range(len(gene_names)))
    
    print(f"  Selected {len(top_genes)} top variable genes")
    
    # Extract expression for top genes (one column at a time to save memory)
    print("  Extracting gene expression data...")
    gene_expression = {}
    for i, gene_idx in enumerate(top_gene_indices):
        gene_name = top_genes[i]
        col = adata.X[:, gene_idx]
        if sparse.issparse(col):
            expr = np.array(col.toarray()).flatten()
        else:
            expr = np.array(col).flatten()
        gene_expression[gene_name] = [round(float(v), 3) for v in expr]
        
        if (i + 1) % 20 == 0:
            print(f"    Processed {i + 1}/{len(top_genes)} genes")
    
    # Build output data structure
    output = {
        'name': dataset_name,
        'n_cells': adata.shape[0],
        'n_genes': adata.shape[1],
        'coordinates': coords,
        'metadata': metadata,
        'metadata_columns': available_meta,
        'gene_names': gene_names,
        'top_genes': top_genes,
    }
    
    # Save main metadata file (without expression - smaller)
    metadata_file = output_dir / f"{dataset_name}_metadata.json"
    with open(metadata_file, 'w') as f:
        json.dump(output, f)
    print(f"  Saved: {metadata_file}")
    
    # Save expression data separately (larger file)
    expression_file = output_dir / f"{dataset_name}_expression.json"
    with open(expression_file, 'w') as f:
        json.dump({'genes': gene_expression}, f)
    print(f"  Saved: {expression_file}")
    
    return {
        'name': dataset_name,
        'n_cells': adata.shape[0],
        'n_genes': adata.shape[1],
        'metadata_columns': available_meta,
        'n_top_genes': len(top_genes)
    }


def main():
    """Main entry point."""
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}")
    
    # Process each H5AD file
    datasets = []
    for h5ad_path in H5AD_FILES:
        if not os.path.exists(h5ad_path):
            print(f"Warning: File not found: {h5ad_path}")
            continue
        
        try:
            info = extract_dataset(h5ad_path, OUTPUT_DIR)
            datasets.append(info)
        except Exception as e:
            print(f"Error processing {h5ad_path}: {e}")
            raise
    
    # Save dataset index
    index_file = OUTPUT_DIR / "datasets.json"
    with open(index_file, 'w') as f:
        json.dump({'datasets': datasets}, f, indent=2)
    print(f"\nSaved dataset index: {index_file}")
    
    print("\n✅ Data extraction complete!")
    for ds in datasets:
        print(f"  - {ds['name']}: {ds['n_cells']} cells, {ds['n_top_genes']} genes with expression data")


if __name__ == "__main__":
    main()
