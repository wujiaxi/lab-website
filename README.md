# 🔬 Lab Website - Spatial Transcriptomics Explorer

A modern lab website featuring an interactive spatial transcriptomics data explorer built with vanilla JavaScript and Plotly.

**Live Demo**: https://wujiaxi.github.io/lab-website/

## ✨ Features

- **Interactive Data Explorer**: Visualize spatial transcriptomics datasets with 24,993+ cells
- **Gene Expression Viewer**: Search and visualize expression of 100+ top variable genes
- **Categorical Metadata**: Explore cells by domain, subclass, and 126 clusters
- **Responsive Design**: Modern, mobile-friendly interface
- **Static Website**: No backend required, fast and easy to deploy

## 🚀 Quick Start

### Local Development

```bash
# Start a local web server
python3 -m http.server 8000

# Open in browser
http://localhost:8000/index.html
```

### View the Explorer

```bash
http://localhost:8000/explorer.html
```

## 📁 Project Structure

```
lab-website/
├── index.html              # Homepage with team and publications
├── explorer.html           # Spatial data explorer
├── explorer.css            # Explorer styles
├── explorer.js             # Explorer functionality
├── data/                   # Processed datasets
│   ├── datasets.json       # Dataset index
│   ├── *_metadata.json     # Cell coordinates and metadata
│   └── *_expression.json   # Gene expression data
├── scripts/
│   └── extract_data.py     # Script to process H5AD files
└── page-templates/         # Original template files
```

## 👥 Team

- **Dr. Songren Wei** - Principal Investigator
- **Dr. Meng Luo** - Doctor
- **Wenxuan Li** - Researcher
- **Chenyang Li** - Researcher
- **Jiaxi Wu** - Researcher

## 📚 Featured Publication

**Charting the spatial transcriptome of the human cerebral cortex at single-cell resolution**
Songren Wei, et al.
*Nature Communications* (2025)
[Read Paper](https://www.nature.com/articles/s41467-025-62793-9)

## 🔧 Adding New Datasets

### 1. Prepare your H5AD file

Ensure your H5AD file has:
- `obs['x']` and `obs['y']` - spatial coordinates
- `obs` columns for metadata (e.g., 'cluster', 'domain', 'subclass')
- `var['name']` or `var_names` - gene names

### 2. Update the extraction script

Edit `scripts/extract_data.py` and add your H5AD file path:

```python
H5AD_FILES = [
    "/path/to/your/dataset.h5ad",
]
```

### 3. Run the extraction script

```bash
# Activate virtual environment (if needed)
source h5ad-env/bin/activate

# Run extraction
python3 scripts/extract_data.py
```

This will create:
- `data/your_dataset_metadata.json` (coordinates + metadata)
- `data/your_dataset_expression.json` (top 100 variable genes)
- `data/datasets.json` (updated index)

### 4. Commit and push

```bash
git add data/
git commit -m "Add new dataset: your_dataset"
git push origin main
```

GitHub Pages will automatically rebuild your site.

## 🌐 Deployment

### GitHub Pages (Current Setup)

1. Push code to GitHub
2. Go to Settings → Pages
3. Select **main** branch → Save
4. Site will be live at: `https://wujiaxi.github.io/lab-website/`

### Alternative Deployment Options

- **Netlify**: Drag and drop deployment
- **Vercel**: Git integration
- **AWS S3**: Static website hosting

## 🛠️ Technical Details

### Technologies Used

- **Plotly.js**: Interactive plotting library
- **Vanilla JavaScript**: No framework dependencies
- **Python + Scanpy**: Data processing pipeline
- **GitHub Pages**: Free static hosting

### Browser Support

- Chrome/Edge (recommended)
- Firefox
- Safari

### Dataset Format

The explorer expects JSON files in this format:

```json
{
  "name": "dataset_name",
  "n_cells": 24993,
  "n_genes": 31208,
  "coordinates": {
    "x": [0.1, 0.2, ...],
    "y": [0.3, 0.4, ...]
  },
  "metadata": {
    "cluster": ["C1", "C2", ...],
    "domain": ["D1", "D2", ...]
  },
  "top_genes": ["GENE1", "GENE2", ...]
}
```

## 📝 License

This project is open source and available for research purposes.

## 🤝 Contributing

To update the website:

1. Make changes locally
2. Test with local server
3. Commit and push to GitHub
4. GitHub Pages will auto-deploy

## 📧 Contact

For questions or collaborations, contact the lab through the website.
