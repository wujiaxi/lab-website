# Quick Reference: Two-Server Deployment

## 🚀 Quick Start

### 1. Deploy Calculation Server First (192.168.25.2)
```bash
ssh root@192.168.25.2

git clone https://github.com/wujiaxi/lab-website.git
cd lab-website
sudo bash scripts/deploy_calculation_server.sh
```

### 2. Deploy Transfer/Static Server (43.207.99.62)
```bash
ssh root@43.207.99.62

git clone https://github.com/wujiaxi/lab-website.git
cd lab-website
sudo bash scripts/deploy_static_server.sh
```

### 3. Verify Deployment
```bash
# Test calculation server
curl http://192.168.25.2:5000/health

# Test static server
curl https://pinolilab.org/
curl https://pinolilab.org/api/datasets
```

---

## 📋 Server Configuration

**Pre-configured server IPs:**

| Server | IP Address | Role |
|--------|------------|------|
| Calculation | 192.168.25.2:5000 | Flask + H5AD processing |
| Transfer/Static | 43.207.99.62 | Nginx + SSL |

**Data directory:** `/home/jiaxi/spatial_data`

---

## 🔍 Quick Tests

### Test Static Server
```bash
curl http://pinolilab.org/
curl http://pinolilab.org/publications
curl http://pinolilab.org/static/images/logo.png
```

### Test Calculation Server (192.168.25.2)
```bash
curl http://192.168.25.2:5000/health
curl http://192.168.25.2:5000/api/datasets
curl http://192.168.25.2:5000/api/cache
```

### Test Integration (API via Static Server)
```bash
curl https://pinolilab.org/api/datasets
curl https://pinolilab.org/api/articles
```

---

## 🔧 Common Commands

### Nginx (Static Server)
```bash
# Test configuration
sudo nginx -t

# Reload
sudo systemctl reload nginx

# Logs
sudo tail -f /var/log/nginx/error.log
```

### Flask (Calculation Server)
```bash
# Status
sudo systemctl status lab-website-calc

# Restart
sudo systemctl restart lab-website-calc

# Logs
sudo journalctl -u lab-website-calc -f

# Cache stats
curl http://192.168.25.2:5000/api/cache
```

---

## 📁 File Locations

### Static Server
- Nginx config: `/etc/nginx/sites-available/pinolilab.org`
- Static HTML: `/var/www/lab-website/static-pages/`
- Static assets: `/var/www/lab-website/backend/static/`
- Repository: `/var/www/lab-website/`

### Calculation Server
- Systemd service: `/etc/systemd/system/lab-website-calc.service`
- Flask app: `/var/www/lab-website/backend/`
- Virtual env: `/var/www/lab-website/backend/venv/`
- H5AD data: `/home/jiaxi/spatial_data/`
- Dataset registry: `/var/www/lab-website/backend/datasets.json`

---

## 🔄 Update Procedures

### Update Static Content
```bash
cd /var/www/lab-website
git pull origin main
python3 scripts/extract_static_pages.py backend/templates /var/www/lab-website/static-pages
sudo systemctl reload nginx
```

### Update Calculation Backend
```bash
cd /var/www/lab-website
git pull origin main
sudo systemctl restart lab-website-calc
```

---

## 🐛 Troubleshooting

### Static pages show 404
```bash
ls -la /var/www/lab-website/static-pages/
sudo tail -f /var/log/nginx/error.log
```

### API requests fail
```bash
# Check Flask is running
sudo systemctl status lab-website-calc

# Test from static server
curl http://192.168.25.2:5000/health

# Check firewall
sudo ufw status
```

### Explorer not loading data
```bash
# Browser console will show errors
# Check CORS in app_calc_server.py
# Verify API is accessible
curl https://pinolilab.org/api/datasets
```

---

## 📊 Architecture Diagram

```
User → Nginx (Static Server) → Flask (Calc Server) → H5AD Files
       ├─ / (static HTML)
       ├─ /publications (static HTML)
       ├─ /static/* (images, CSS)
       └─ /api/* (proxied) ────────────→ Process & Return Data
       └─ /explorer (proxied) ──────────→ Render Template
```

---

## 🎯 Key Differences from Single Server

| Aspect | Single Server | Two Server |
|--------|---------------|------------|
| **Flask routes** | All pages + API | Only /explorer + /api/* |
| **Static pages** | Flask templates | Pre-rendered HTML on Nginx |
| **SSL** | Flask or Nginx | Nginx only (static server) |
| **Port 5000** | Public | Private (calc server only) |
| **CORS** | Not needed | Required |
| **Deployment** | One script | Two separate scripts |

---

## ⚡ Performance Notes

- Static pages served directly by Nginx (faster)
- No Flask overhead for simple pages
- Calculation server can be scaled independently
- Long-running H5AD processing doesn't block static content
- Gunicorn workers optimized for CPU-bound tasks

---

## 🔒 Security Notes

- Static server (Nginx) is public-facing with SSL
- Calculation server can be on private network
- Firewall blocks direct access to Flask (port 5000)
- All traffic routes through Nginx reverse proxy

---

For full documentation, see [README_TWO_SERVER.md](./README_TWO_SERVER.md)
