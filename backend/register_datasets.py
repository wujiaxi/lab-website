#!/usr/bin/env python3
"""
Batch register datasets from a directory of H5AD files
Usage: python register_datasets.py /path/to/h5ad/files "Article Name"
"""

import sys
import json
import os
from pathlib import Path
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Configuration
BASE_DIR = Path(__file__).parent
REGISTRY_FILE = BASE_DIR / 'datasets.json'
THUMBNAIL_DIR = BASE_DIR / 'static' / 'thumbnails'

THUMBNAIL_DIR.mkdir(parents=True, exist_ok=True)


def load_registry():
    if REGISTRY_FILE.exists():
        with open(REGISTRY_FILE) as f:
            return json.load(f)
    return {'articles': [], 'datasets': {}}


def save_registry(registry):
    with open(REGISTRY_FILE, 'w') as f:
        json.dump(registry, f, indent=2, ensure_ascii=False)


def generate_thumbnail(dataset_id, adata):
    """Generate thumbnail from spatial data"""
    thumbnail_path = THUMBNAIL_DIR / f"{dataset_id}.png"

    fig, ax = plt.subplots(figsize=(4, 4), dpi=100)

    x = adata.obs['x'].values
    y = adata.obs['y'].values

    # Color by first available metadata
    color_col = None
    for col in ['cluster', 'domain', 'subclass', 'cell_type']:
        if col in adata.obs.columns:
            color_col = col
            break

    if color_col:
        categories = adata.obs[color_col].astype('category')
        colors = plt.cm.tab20(np.linspace(0, 1, min(20, len(categories.cat.categories))))
        color_map = dict(zip(categories.cat.categories[:20], colors))
        c = [color_map.get(cat, [0.5, 0.5, 0.5, 1]) for cat in categories]
        ax.scatter(x, y, c=c, s=0.5, alpha=0.8)
    else:
        ax.scatter(x, y, s=0.5, alpha=0.8, c='steelblue')

    ax.set_aspect('equal')
    ax.axis('off')
    plt.tight_layout(pad=0)

    fig.savefig(thumbnail_path, bbox_inches='tight', pad_inches=0.1,
                facecolor='white', dpi=100)
    plt.close(fig)

    print(f"  ✅ Thumbnail saved: {thumbnail_path}")


def register_dataset(h5ad_path, article_name, dataset_id=None):
    """Register a single dataset"""
    h5ad_path = Path(h5ad_path)

    if not h5ad_path.exists():
        print(f"  ❌ File not found: {h5ad_path}")
        return False

    if dataset_id is None:
        dataset_id = h5ad_path.stem

    print(f"\n📁 Processing: {dataset_id}")
    print(f"   Path: {h5ad_path}")

    # Load dataset
    print("   Loading H5AD...")
    adata = sc.read_h5ad(h5ad_path)
    print(f"   Shape: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Get metadata columns
    metadata_cols = [col for col in adata.obs.columns
                     if col in ['cluster', 'domain', 'domain.fine', 'subclass', 'cell_type']]
    print(f"   Metadata columns: {metadata_cols}")

    # Load registry
    registry = load_registry()

    # Add article if new
    if article_name not in registry['articles']:
        registry['articles'].append(article_name)
        print(f"   Added new article: {article_name}")

    # Add dataset
    registry['datasets'][dataset_id] = {
        'path': str(h5ad_path.absolute()),
        'name': dataset_id.replace('_', ' ').replace('-', ' '),
        'article': article_name,
        'n_cells': adata.shape[0],
        'n_genes': adata.shape[1],
        'metadata_columns': metadata_cols,
        'description': f"Spatial transcriptomics data from {article_name}"
    }

    # Save registry
    save_registry(registry)
    print(f"   ✅ Registered in datasets.json")

    # Generate thumbnail
    print("   Generating thumbnail...")
    generate_thumbnail(dataset_id, adata)

    return True


def register_directory(directory, article_name):
    """Register all H5AD files in a directory"""
    directory = Path(directory)

    if not directory.exists():
        print(f"❌ Directory not found: {directory}")
        return

    h5ad_files = list(directory.glob('*.h5ad'))
    print(f"Found {len(h5ad_files)} H5AD files in {directory}")

    success = 0
    for h5ad_path in sorted(h5ad_files):
        if register_dataset(h5ad_path, article_name):
            success += 1

    print(f"\n✅ Successfully registered {success}/{len(h5ad_files)} datasets")


def main():
    if len(sys.argv) < 3:
        print("Usage:")
        print("  Register single file:  python register_datasets.py /path/to/file.h5ad 'Article Name'")
        print("  Register directory:    python register_datasets.py /path/to/directory/ 'Article Name'")
        print("")
        print("Example:")
        print("  python register_datasets.py /data/cerebral_cortex/ 'Nature Communications 2025 - Cerebral Cortex'")
        sys.exit(1)

    path = sys.argv[1]
    article_name = sys.argv[2]

    if os.path.isdir(path):
        register_directory(path, article_name)
    else:
        register_dataset(path, article_name)


if __name__ == '__main__':
    main()
