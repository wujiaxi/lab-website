# Two-Server Architecture Deployment Guide

This guide explains how to deploy the lab website using a two-server architecture for better performance and scalability.

## Architecture Overview

```
┌──────────────┐
│    User      │
│   Browser    │
└──────┬───────┘
       │ HTTPS
       ▼
┌──────────────────────────────────────┐
│  Server 1: Static/Transfer Server    │
│  IP: 175.41.197.121                  │
│  Domain: pinolilab.org               │
│  ──────────────────────────────────  │
│  ┌─────────────────────────────┐    │
│  │         Nginx               │    │
│  │  • Serves static HTML/CSS   │    │
│  │  • Reverse proxy /api/*     │    │
│  │  • SSL termination          │    │
│  └─────────────────────────────┘    │
└──────────────┬───────────────────────┘
               │ /api/* and /explorer
               │ (HTTP proxy)
               ▼
┌──────────────────────────────────────┐
│  Server 2: Calculation Server        │
│  IP: <YOUR_CALC_SERVER_IP>           │
│  ──────────────────────────────────  │
│  ┌─────────────────────────────┐    │
│  │   Flask + Gunicorn          │    │
│  │  • Process H5AD files       │    │
│  │  • Generate thumbnails      │    │
│  │  • Cache management         │    │
│  └─────────────────────────────┘    │
│  ┌─────────────────────────────┐    │
│  │   H5AD Data Repository      │    │
│  │  /home/jiaxi/spatial_data   │    │
│  └─────────────────────────────┘    │
└──────────────────────────────────────┘
```

## Server Roles

### Server 1: Static/Transfer Server
- **Purpose**: Front-facing web server
- **Technology**: Nginx
- **Responsibilities**:
  - Serve static HTML pages (/, /publications, /research, /contact)
  - Serve static assets (images, CSS, JavaScript)
  - SSL/TLS termination
  - Reverse proxy to calculation server for dynamic content
  
### Server 2: Calculation Server
- **Purpose**: Backend computation and data processing
- **Technology**: Python 3 + Flask + Gunicorn
- **Responsibilities**:
  - Load and process H5AD files
  - Generate spatial data visualizations
  - Handle gene expression queries
  - Cache management for performance

## Deployment Instructions

### Prerequisites

- Two servers (can be same machine for testing)
- Ubuntu/Debian Linux (or similar)
- Root/sudo access
- Domain name pointing to static server (pinolilab.org)

### Step 1: Deploy Static Server

On your static server (175.41.197.121):

```bash
# Clone repository
git clone https://github.com/wujiaxi/lab-website.git
cd lab-website

# Edit the script to set your calculation server IP
nano scripts/deploy_static_server.sh
# Change: CALC_SERVER="localhost:5000" to your actual IP
# Example: CALC_SERVER="192.168.1.100:5000"

# Run deployment script
sudo bash scripts/deploy_static_server.sh
```

The script will:
- ✓ Install Nginx
- ✓ Clone the repository
- ✓ Extract static HTML from Flask templates
- ✓ Configure Nginx for static serving + reverse proxy
- ✓ Optionally set up SSL with Let's Encrypt

### Step 2: Deploy Calculation Server

On your calculation server:

```bash
# Clone repository
git clone https://github.com/wujiaxi/lab-website.git
cd lab-website

# Edit the script to set your data directory
nano scripts/deploy_calculation_server.sh
# Verify: DATA_DIR="/home/jiaxi/spatial_data"

# Run deployment script
sudo bash scripts/deploy_calculation_server.sh
```

The script will:
- ✓ Install Python 3 and scientific stack
- ✓ Set up virtual environment
- ✓ Install Flask dependencies
- ✓ Create systemd service
- ✓ Configure firewall (optional)
- ✓ Start the Flask backend

### Step 3: Register H5AD Datasets

On the calculation server:

```bash
cd /var/www/lab-website/backend
source venv/bin/activate

# Register your H5AD files
python register_datasets.py

# Or manually via API
curl -X POST http://localhost:5000/api/register \
  -H "Content-Type: application/json" \
  -d '{
    "id": "dataset1",
    "path": "/home/jiaxi/spatial_data/sample.h5ad",
    "article": "Spatial Transcriptomics Study",
    "name": "Sample Dataset",
    "description": "Description of the dataset"
  }'
```

### Step 4: Verify Deployment

Test the static server:
```bash
# Test static pages
curl http://pinolilab.org/
curl http://pinolilab.org/publications
curl http://pinolilab.org/research
```

Test the calculation server:
```bash
# From calculation server
curl http://localhost:5000/health
curl http://localhost:5000/api/datasets
```

Test integration:
```bash
# From any machine
curl https://pinolilab.org/api/datasets
curl https://pinolilab.org/api/cache
```

Visit in browser:
- https://pinolilab.org (homepage)
- https://pinolilab.org/explorer (data explorer)
- https://pinolilab.org/publications
- https://pinolilab.org/research

## File Structure

```
lab-website/
├── backend/
│   ├── app.py                    # Original Flask app (single server)
│   ├── app_calc_server.py        # Modified for calculation server
│   ├── requirements.txt          # Python dependencies
│   ├── register_datasets.py      # Dataset registration tool
│   ├── datasets.json             # Dataset registry
│   ├── static/                   # Static assets
│   │   ├── images/
│   │   └── thumbnails/
│   └── templates/                # Jinja2 templates
│       ├── index.html
│       ├── explorer.html
│       ├── publications.html
│       ├── research.html
│       └── contact.html
├── scripts/
│   ├── deploy_static_server.sh   # Static server deployment
│   ├── deploy_calculation_server.sh  # Calc server deployment
│   └── extract_static_pages.py   # Template → HTML converter
└── README_TWO_SERVER.md          # This file
```

## Service Management

### Static Server (Nginx)

```bash
# Check Nginx status
sudo systemctl status nginx

# Reload configuration
sudo nginx -t  # Test config first
sudo systemctl reload nginx

# View logs
sudo tail -f /var/log/nginx/access.log
sudo tail -f /var/log/nginx/error.log
```

### Calculation Server (Flask)

```bash
# Check service status
sudo systemctl status lab-website-calc

# Restart service
sudo systemctl restart lab-website-calc

# View logs
sudo journalctl -u lab-website-calc -f

# Check cache status
curl http://localhost:5000/api/cache
```

## Updating the Application

### Update Static Content

On the static server:
```bash
cd /var/www/lab-website
git pull origin main

# Re-extract static pages
python3 scripts/extract_static_pages.py \
  backend/templates \
  /var/www/lab-website/static-pages

sudo systemctl reload nginx
```

### Update Calculation Backend

On the calculation server:
```bash
cd /var/www/lab-website
git pull origin main

# Update Python dependencies if needed
cd backend
source venv/bin/activate
pip install -r requirements.txt

# Restart service
sudo systemctl restart lab-website-calc
```

## Configuration Files

### Nginx Configuration
Location: `/etc/nginx/sites-available/pinolilab.org`

Key sections:
- Static page serving: `location = /`
- Static assets: `location /static`
- API proxy: `location /api`
- Explorer proxy: `location /explorer`

### Systemd Service
Location: `/etc/systemd/system/lab-website-calc.service`

Key settings:
- Workers: 4 (adjust based on CPU cores)
- Timeout: 300s (for large dataset processing)
- Memory limit: 8G (adjust based on available RAM)

## Troubleshooting

### Static pages not loading
```bash
# Check Nginx logs
sudo tail -f /var/log/nginx/error.log

# Verify static files exist
ls -la /var/www/lab-website/static-pages/

# Re-extract templates
cd /var/www/lab-website
python3 scripts/extract_static_pages.py \
  backend/templates \
  /var/www/lab-website/static-pages
```

### API requests failing
```bash
# Check if calculation server is running
sudo systemctl status lab-website-calc

# Test connection from static server
curl http://<CALC_SERVER_IP>:5000/health

# Check firewall
sudo ufw status

# View Flask logs
sudo journalctl -u lab-website-calc -n 100
```

### Explorer page errors
```bash
# Check browser console for errors
# Verify API endpoints are accessible
curl https://pinolilab.org/api/datasets
curl https://pinolilab.org/api/articles

# Check CORS settings in app_calc_server.py
```

### Large dataset loading slowly
```bash
# Check cache status
curl http://localhost:5000/api/cache

# Increase Gunicorn timeout in systemd service
sudo nano /etc/systemd/system/lab-website-calc.service
# Change: --timeout 300 to --timeout 600

sudo systemctl daemon-reload
sudo systemctl restart lab-website-calc
```

## Performance Tuning

### Nginx Optimization
```nginx
# Add to nginx.conf
worker_processes auto;
worker_connections 1024;

# Enable gzip compression
gzip on;
gzip_types text/plain text/css application/json application/javascript;
```

### Flask/Gunicorn Optimization
```bash
# Adjust workers based on CPU cores
# Rule of thumb: (2 × CPU cores) + 1
# Edit systemd service:
--workers 8  # For 4-core CPU
```

### Cache Configuration
```python
# In app.py, adjust cache settings:
CACHE_EXPIRY_HOURS = 12  # Keep data longer
```

## Security Considerations

1. **Firewall**: Calculation server should only accept connections from static server
2. **SSL**: Use Let's Encrypt for free SSL certificates
3. **API Authentication**: Consider adding API keys for `/api/register` endpoint
4. **CORS**: Configure CORS to only allow your domain
5. **Updates**: Regularly update system packages and Python dependencies

## Monitoring

Set up monitoring for:
- Nginx access/error logs
- Flask application logs
- Disk space (H5AD files can be large)
- Memory usage (dataset caching)
- CPU usage during data processing

## Backup Strategy

Regular backups of:
- H5AD data files: `/home/jiaxi/spatial_data/`
- Dataset registry: `/var/www/lab-website/backend/datasets.json`
- Thumbnails: `/var/www/lab-website/backend/static/thumbnails/`
- Nginx configuration: `/etc/nginx/sites-available/`

## Cost Optimization

- Use the static server for low-cost hosting
- Use calculation server with higher specs only when needed
- Consider auto-scaling for calculation server
- Use CDN for static assets

## Migration from Single Server

If you're currently running the monolithic Flask app:

1. Back up current deployment
2. Deploy static server first (will serve as reverse proxy to old Flask)
3. Deploy calculation server
4. Update Nginx to proxy to calculation server
5. Test thoroughly
6. Shut down old Flask instance

---

## Support

For issues or questions:
- GitHub: https://github.com/wujiaxi/lab-website
- Check logs first (Nginx + Flask)
- Review this guide's troubleshooting section
