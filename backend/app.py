#!/usr/bin/env python3
"""
Flask Backend for Spatial Data Explorer
Serves H5AD data on-demand with caching and thumbnail generation
Features:
- Time-based cache with automatic garbage collection
- Memory management for large datasets
- On-demand data loading
"""

import json
import os
import time
import threading
from pathlib import Path
from datetime import datetime

from flask import Flask, jsonify, request, send_file, render_template
from flask_cors import CORS
import numpy as np
import scanpy as sc
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

# Configuration
app = Flask(__name__, static_folder='static', template_folder='templates')
CORS(app)  # Enable CORS for frontend

# Paths
BASE_DIR = Path(__file__).parent
DATA_DIR = Path(os.environ.get('DATA_DIR', '/home/jiaxi/spatial_data'))
REGISTRY_FILE = BASE_DIR / 'datasets.json'
THUMBNAIL_DIR = BASE_DIR / 'static' / 'thumbnails'

# Cache Configuration
CACHE_EXPIRY_HOURS = 6  # Clear data after 6 hours of no use
CACHE_CHECK_INTERVAL = 300  # Check every 5 minutes

# Ensure directories exist
THUMBNAIL_DIR.mkdir(parents=True, exist_ok=True)


# ============== Time-based Cache with GC ==============

class TimedCache:
    """
    Cache with automatic garbage collection based on last access time.
    Data is cleared after CACHE_EXPIRY_HOURS of no access.
    """

    def __init__(self, expiry_hours=6):
        self.cache = {}  # {dataset_id: {'data': adata, 'last_access': timestamp}}
        self.expiry_seconds = expiry_hours * 3600
        self.lock = threading.Lock()

        # Start background GC thread
        self.gc_thread = threading.Thread(target=self._gc_loop, daemon=True)
        self.gc_thread.start()
        print(f"🧹 Cache GC started (expiry: {expiry_hours} hours)")

    def get(self, dataset_id: str):
        """Get data from cache, update last access time"""
        with self.lock:
            if dataset_id in self.cache:
                self.cache[dataset_id]['last_access'] = time.time()
                print(f"📦 Cache HIT: {dataset_id}")
                return self.cache[dataset_id]['data']
        return None

    def set(self, dataset_id: str, data):
        """Store data in cache"""
        with self.lock:
            self.cache[dataset_id] = {
                'data': data,
                'last_access': time.time()
            }
            print(f"💾 Cache SET: {dataset_id} (total cached: {len(self.cache)})")

    def remove(self, dataset_id: str):
        """Manually remove data from cache"""
        with self.lock:
            if dataset_id in self.cache:
                del self.cache[dataset_id]
                print(f"🗑️ Cache REMOVE: {dataset_id}")

    def clear(self):
        """Clear all cached data"""
        with self.lock:
            count = len(self.cache)
            self.cache.clear()
            print(f"🧹 Cache CLEARED: {count} datasets removed")

    def _gc_loop(self):
        """Background garbage collection loop"""
        while True:
            time.sleep(CACHE_CHECK_INTERVAL)
            self._run_gc()

    def _run_gc(self):
        """Run garbage collection - remove expired entries"""
        now = time.time()
        expired = []

        with self.lock:
            for dataset_id, entry in self.cache.items():
                age_seconds = now - entry['last_access']
                if age_seconds > self.expiry_seconds:
                    expired.append(dataset_id)

            for dataset_id in expired:
                age_hours = (now - self.cache[dataset_id]['last_access']) / 3600
                del self.cache[dataset_id]
                print(f"🧹 GC: Removed {dataset_id} (unused for {age_hours:.1f} hours)")

        if expired:
            print(f"🧹 GC complete: {len(expired)} datasets cleared, {len(self.cache)} remaining")

    def get_stats(self):
        """Get cache statistics"""
        with self.lock:
            stats = {
                'cached_datasets': len(self.cache),
                'datasets': {}
            }
            now = time.time()
            for dataset_id, entry in self.cache.items():
                age_hours = (now - entry['last_access']) / 3600
                stats['datasets'][dataset_id] = {
                    'last_access': datetime.fromtimestamp(entry['last_access']).isoformat(),
                    'age_hours': round(age_hours, 2),
                    'expires_in_hours': round(self.expiry_seconds/3600 - age_hours, 2)
                }
            return stats


# Initialize cache
adata_cache = TimedCache(expiry_hours=CACHE_EXPIRY_HOURS)


def load_adata(dataset_id: str):
    """Load AnnData object with caching"""
    # Check cache first
    cached = adata_cache.get(dataset_id)
    if cached is not None:
        return cached

    # Load from disk
    registry = load_registry()
    if dataset_id not in registry['datasets']:
        return None

    h5ad_path = registry['datasets'][dataset_id]['path']
    print(f"📂 Loading dataset from disk: {dataset_id}")
    print(f"   Path: {h5ad_path}")

    adata = sc.read_h5ad(h5ad_path)
    print(f"   Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Store in cache
    adata_cache.set(dataset_id, adata)

    return adata


def load_registry():
    """Load dataset registry"""
    if REGISTRY_FILE.exists():
        with open(REGISTRY_FILE) as f:
            return json.load(f)
    return {'articles': [], 'datasets': {}}


def save_registry(registry):
    """Save dataset registry"""
    with open(REGISTRY_FILE, 'w') as f:
        json.dump(registry, f, indent=2)


def generate_thumbnail(dataset_id: str, adata, force=False):
    """Generate thumbnail image from spatial data"""
    thumbnail_path = THUMBNAIL_DIR / f"{dataset_id}.png"

    if thumbnail_path.exists() and not force:
        return str(thumbnail_path)

    # Create spatial plot
    fig, ax = plt.subplots(figsize=(4, 4), dpi=100)

    x = adata.obs['x'].values
    y = adata.obs['y'].values

    # Color by first available metadata column
    color_col = None
    for col in ['cluster', 'domain', 'subclass', 'cell_type']:
        if col in adata.obs.columns:
            color_col = col
            break

    if color_col:
        categories = adata.obs[color_col].astype('category')
        colors = plt.cm.tab20(np.linspace(0, 1, len(categories.cat.categories)))
        color_map = dict(zip(categories.cat.categories, colors))
        c = [color_map[cat] for cat in categories]
        ax.scatter(x, y, c=c, s=0.5, alpha=0.8)
    else:
        ax.scatter(x, y, s=0.5, alpha=0.8, c='steelblue')

    ax.set_aspect('equal')
    ax.axis('off')
    plt.tight_layout(pad=0)

    fig.savefig(thumbnail_path, bbox_inches='tight', pad_inches=0.1,
                facecolor='white', dpi=100)
    plt.close(fig)

    return str(thumbnail_path)


# ============== API Routes ==============

@app.route('/')
def index():
    """Serve main page"""
    return render_template('index.html')


@app.route('/api/articles')
def get_articles():
    """Get list of all articles/categories"""
    registry = load_registry()
    return jsonify(registry.get('articles', []))


@app.route('/api/datasets')
def get_datasets():
    """Get all datasets, optionally filtered by article"""
    registry = load_registry()
    article = request.args.get('article')

    datasets = []
    for ds_id, ds_info in registry.get('datasets', {}).items():
        if article and ds_info.get('article') != article:
            continue

        datasets.append({
            'id': ds_id,
            'name': ds_info.get('name', ds_id),
            'article': ds_info.get('article', 'Unknown'),
            'n_cells': ds_info.get('n_cells', 0),
            'n_genes': ds_info.get('n_genes', 0),
            'thumbnail': f"/api/datasets/{ds_id}/thumbnail",
            'description': ds_info.get('description', '')
        })

    return jsonify(datasets)


@app.route('/api/datasets/<dataset_id>')
def get_dataset(dataset_id):
    """Get specific dataset metadata"""
    registry = load_registry()

    if dataset_id not in registry.get('datasets', {}):
        return jsonify({'error': 'Dataset not found'}), 404

    ds_info = registry['datasets'][dataset_id]
    return jsonify({
        'id': dataset_id,
        'name': ds_info.get('name', dataset_id),
        'article': ds_info.get('article', 'Unknown'),
        'n_cells': ds_info.get('n_cells', 0),
        'n_genes': ds_info.get('n_genes', 0),
        'metadata_columns': ds_info.get('metadata_columns', []),
        'description': ds_info.get('description', '')
    })


@app.route('/api/datasets/<dataset_id>/thumbnail')
def get_thumbnail(dataset_id):
    """Get or generate thumbnail for dataset"""
    registry = load_registry()

    if dataset_id not in registry.get('datasets', {}):
        return jsonify({'error': 'Dataset not found'}), 404

    thumbnail_path = THUMBNAIL_DIR / f"{dataset_id}.png"

    # Generate if not exists
    if not thumbnail_path.exists():
        adata = load_adata(dataset_id)
        if adata is None:
            return jsonify({'error': 'Could not load dataset'}), 500
        generate_thumbnail(dataset_id, adata)

    return send_file(thumbnail_path, mimetype='image/png')


@app.route('/api/datasets/<dataset_id>/spatial')
def get_spatial_data(dataset_id):
    """Get spatial coordinates and metadata for visualization"""
    adata = load_adata(dataset_id)
    if adata is None:
        return jsonify({'error': 'Dataset not found'}), 404

    # Get requested metadata column
    color_by = request.args.get('color_by', 'cluster')

    response = {
        'n_cells': adata.shape[0],
        'n_genes': adata.shape[1],
        'coordinates': {
            'x': adata.obs['x'].tolist(),
            'y': adata.obs['y'].tolist()
        },
        'metadata_columns': [col for col in ['cluster', 'domain', 'domain.fine', 'subclass']
                            if col in adata.obs.columns],
    }

    # Add requested metadata
    if color_by in adata.obs.columns:
        response['metadata'] = {
            color_by: adata.obs[color_by].astype(str).tolist()
        }

    return jsonify(response)


@app.route('/api/datasets/<dataset_id>/genes')
def search_genes(dataset_id):
    """Search genes in dataset"""
    adata = load_adata(dataset_id)
    if adata is None:
        return jsonify({'error': 'Dataset not found'}), 404

    query = request.args.get('q', '').upper()
    limit = int(request.args.get('limit', 50))

    # Get gene names
    if 'name' in adata.var.columns:
        gene_names = adata.var['name'].tolist()
    else:
        gene_names = adata.var_names.tolist()

    # Filter by query
    if query:
        matches = [g for g in gene_names if query in g.upper()][:limit]
    else:
        # Return top variable genes
        matches = gene_names[:limit]

    return jsonify({'genes': matches})


@app.route('/api/datasets/<dataset_id>/expression/<gene>')
def get_expression(dataset_id, gene):
    """Get expression values for a specific gene"""
    adata = load_adata(dataset_id)
    if adata is None:
        return jsonify({'error': 'Dataset not found'}), 404

    # Find gene index
    if 'name' in adata.var.columns:
        gene_names = adata.var['name'].tolist()
    else:
        gene_names = adata.var_names.tolist()

    if gene not in gene_names:
        return jsonify({'error': f'Gene {gene} not found'}), 404

    gene_idx = gene_names.index(gene)

    # Get expression values
    from scipy import sparse
    col = adata.X[:, gene_idx]
    if sparse.issparse(col):
        expr = np.array(col.toarray()).flatten()
    else:
        expr = np.array(col).flatten()

    return jsonify({
        'gene': gene,
        'expression': [round(float(v), 3) for v in expr]
    })


@app.route('/api/cache')
def get_cache_stats():
    """Get cache statistics"""
    stats = adata_cache.get_stats()
    stats['config'] = {
        'expiry_hours': CACHE_EXPIRY_HOURS,
        'check_interval_seconds': CACHE_CHECK_INTERVAL
    }
    return jsonify(stats)


@app.route('/api/cache/clear', methods=['POST'])
def clear_cache():
    """Clear all cached data"""
    adata_cache.clear()
    return jsonify({'success': True, 'message': 'Cache cleared'})


@app.route('/api/cache/<dataset_id>', methods=['DELETE'])
def remove_from_cache(dataset_id):
    """Remove specific dataset from cache"""
    adata_cache.remove(dataset_id)
    return jsonify({'success': True, 'message': f'Removed {dataset_id} from cache'})


@app.route('/api/register', methods=['POST'])
def register_dataset():
    """Register a new dataset (for admin use)"""
    data = request.json

    required = ['id', 'path', 'article']
    if not all(k in data for k in required):
        return jsonify({'error': f'Missing required fields: {required}'}), 400

    h5ad_path = Path(data['path'])
    if not h5ad_path.exists():
        return jsonify({'error': f'File not found: {h5ad_path}'}), 400

    # Load to get metadata
    adata = sc.read_h5ad(h5ad_path)

    registry = load_registry()

    # Add article if new
    if data['article'] not in registry['articles']:
        registry['articles'].append(data['article'])

    # Add dataset
    registry['datasets'][data['id']] = {
        'path': str(h5ad_path),
        'name': data.get('name', data['id']),
        'article': data['article'],
        'n_cells': adata.shape[0],
        'n_genes': adata.shape[1],
        'metadata_columns': [col for col in adata.obs.columns
                            if col in ['cluster', 'domain', 'domain.fine', 'subclass']],
        'description': data.get('description', '')
    }

    save_registry(registry)

    # Generate thumbnail
    generate_thumbnail(data['id'], adata)

    return jsonify({'success': True, 'id': data['id']})


if __name__ == '__main__':
    print("Starting Spatial Data Explorer Backend...")
    print(f"Data directory: {DATA_DIR}")
    print(f"Registry file: {REGISTRY_FILE}")
    app.run(host='0.0.0.0', port=5000, debug=True)
