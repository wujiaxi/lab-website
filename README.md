# 🔬 Spatial Data Explorer

A Flask-based web application for interactive visualization of spatial transcriptomics data (H5AD files).

## ✨ Features

- **Dataset Browser**: Gallery view with thumbnails, filterable by article/category
- **Interactive Explorer**: Visualize spatial coordinates with Plotly.js
- **Gene Expression**: Search and visualize any gene's expression
- **Smart Caching**: Automatic memory management with 6-hour expiry
- **On-demand Loading**: Data loaded from H5AD files only when requested
- **100+ Datasets Support**: Designed for large-scale data exploration

## 📁 Project Structure

```
lab-website/
├── backend/
│   ├── app.py                 # Flask application
│   ├── datasets.json          # Dataset registry
│   ├── register_datasets.py   # Batch registration script
│   ├── requirements.txt       # Python dependencies
│   ├── run.sh                 # Startup script
│   ├── static/
│   │   └── thumbnails/        # Generated thumbnails
│   └── templates/
│       └── index.html         # Frontend UI
├── h5ad-env/                  # Python virtual environment
└── README.md
```

## 🚀 Quick Start

### 1. Install Dependencies

```bash
cd backend
pip install -r requirements.txt
```

### 2. Register Datasets

```bash
# Register all H5AD files in a directory
python register_datasets.py /path/to/h5ad/files/ "Article Name"

# Example
python register_datasets.py /data/cerebral_cortex/ "Nature Communications 2025 - Cerebral Cortex"
```

### 3. Start Server

```bash
# Development mode
bash run.sh

# Production mode (with gunicorn)
bash run.sh prod
```

### 4. Access

Open `http://localhost:5000/` in your browser.

## 👥 Team

- **Dr. Songren Wei** - Principal Investigator
- **Dr. Meng Luo** - Doctor
- **Wenxuan Li** - Researcher
- **Chenyang Li** - Researcher
- **Jiaxi Wu** - Researcher

## 📚 Featured Publication

**Charting the spatial transcriptome of the human cerebral cortex at single-cell resolution**
*Nature Communications* (2025)
[Read Paper](https://www.nature.com/articles/s41467-025-62793-9)

## 🔧 API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/articles` | GET | List all article categories |
| `/api/datasets` | GET | List all datasets |
| `/api/datasets?article=xxx` | GET | Filter datasets by article |
| `/api/datasets/<id>` | GET | Get dataset metadata |
| `/api/datasets/<id>/thumbnail` | GET | Get/generate thumbnail |
| `/api/datasets/<id>/spatial` | GET | Get spatial coordinates |
| `/api/datasets/<id>/genes?q=xxx` | GET | Search genes |
| `/api/datasets/<id>/expression/<gene>` | GET | Get gene expression |
| `/api/cache` | GET | View cache statistics |
| `/api/cache/clear` | POST | Clear all cached data |
| `/api/cache/<id>` | DELETE | Remove dataset from cache |

## 💾 Memory Management

The application uses a time-based cache with automatic garbage collection:

- Data is cached in memory after first access
- Cache expires after **6 hours** of no use
- Background thread checks every 5 minutes
- Manual cache control via API endpoints

```bash
# Check cache status
curl http://localhost:5000/api/cache

# Clear all cache
curl -X POST http://localhost:5000/api/cache/clear

# Remove specific dataset
curl -X DELETE http://localhost:5000/api/cache/dataset_id
```

## 🌐 Deployment

### Production with Nginx

1. Install nginx: `sudo apt install nginx`
2. Configure nginx as reverse proxy to Flask
3. Run Flask with gunicorn: `bash run.sh prod`

Example nginx config:
```nginx
server {
    listen 80;
    server_name your-domain.com;

    location / {
        proxy_pass http://127.0.0.1:5000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
    }
}
```

## 🛠️ Technical Details

- **Backend**: Python Flask
- **Frontend**: Vanilla JavaScript + Plotly.js
- **Data Format**: H5AD (AnnData)
- **Caching**: Custom time-based LRU cache
- **Thumbnail Generation**: Matplotlib

## 📧 Contact

For questions or collaborations, contact the lab.
