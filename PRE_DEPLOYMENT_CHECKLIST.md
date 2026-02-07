# Pre-Deployment Checklist

Complete this checklist before running deployment scripts.

## 📋 Pre-Deployment Tasks

### Configuration (Already Done ✅)
- [x] `CALC_SERVER` set to `192.168.25.2:5000`
- [x] `STATIC_SERVER_IP` set to `43.207.99.62`
- [x] `DATA_DIR` set to `/home/jiaxi/spatial_data`
- [ ] Verify `DATA_DIR` path exists on calculation server
- [ ] Provide `SSL_EMAIL` for Let's Encrypt certificate (optional)

### Server Access
- [ ] SSH access to transfer server: `ssh root@43.207.99.62`
- [ ] SSH access to calculation server: `ssh root@192.168.25.2`
- [ ] Root/sudo privileges on both servers
- [ ] Servers can reach GitHub (for git clone)

### DNS Configuration
- [ ] Domain `pinolilab.org` points to static server IP
- [ ] Domain `www.pinolilab.org` points to static server IP
- [ ] DNS propagation complete (check with `dig pinolilab.org`)

### Network Configuration
- [ ] Transfer server (43.207.99.62) can reach calculation server (192.168.25.2) on port 5000
- [ ] Test with: `curl http://192.168.25.2:5000` (after Flask is running)
- [ ] Firewall rules allow traffic between servers

### Data Preparation
- [ ] H5AD files copied to `DATA_DIR` on calculation server
- [ ] Data directory permissions set correctly (`chmod 755`)
- [ ] Enough disk space for datasets and cache
- [ ] Backup of existing data (if any)

### GitHub Repository
- [ ] All new files committed to repository
- [ ] Scripts pushed to GitHub
- [ ] Repository accessible from both servers

### Backup Existing Setup (if applicable)
- [ ] Backup current Flask application
- [ ] Backup Nginx configuration
- [ ] Backup systemd services
- [ ] Backup datasets and database
- [ ] Document current deployment for rollback

## 🚀 Deployment Steps

### Step 1: Deploy Calculation Server First (192.168.25.2)

```bash
# SSH to calculation server
ssh root@192.168.25.2

# Clone repository
git clone https://github.com/wujiaxi/lab-website.git
cd lab-website

# Run deployment
sudo bash scripts/deploy_calculation_server.sh
```

**Expected Time:** 10-15 minutes (includes Python package installation)

**Verification:**
- [ ] Flask service running: `sudo systemctl status lab-website-calc`
- [ ] Health check passes: `curl http://192.168.25.2:5000/health`
- [ ] API responds: `curl http://192.168.25.2:5000/api/datasets`
- [ ] No errors in logs: `sudo journalctl -u lab-website-calc -n 50`

---

### Step 2: Deploy Transfer/Static Server (43.207.99.62)

```bash
# SSH to transfer server
ssh root@43.207.99.62

# Clone repository
git clone https://github.com/wujiaxi/lab-website.git
cd lab-website

# Run deployment
sudo bash scripts/deploy_static_server.sh
```

**Expected Time:** 5-10 minutes

**Verification:**
- [ ] Nginx is running: `sudo systemctl status nginx`
- [ ] No configuration errors: `sudo nginx -t`
- [ ] Static pages accessible: `curl http://localhost/`
- [ ] SSL installed (if enabled): `curl https://pinolilab.org/`

---

### Step 2: Deploy Calculation Server

```bash
# SSH to calculation server
ssh user@<CALC_SERVER_IP>

# Clone repository
git clone https://github.com/wujiaxi/lab-website.git
cd lab-website

# Configure deployment
nano DEPLOYMENT_CONFIG.sh
# Update DATA_DIR, STATIC_SERVER_IP, etc.

# Source configuration
source DEPLOYMENT_CONFIG.sh

# Run deployment
sudo -E bash scripts/deploy_calculation_server.sh
```

**Expected Time:** 10-15 minutes (includes Python package installation)

**Verification:**
- [ ] Flask service running: `sudo systemctl status lab-website-calc`
- [ ] Health check passes: `curl http://localhost:5000/health`
- [ ] API responds: `curl http://localhost:5000/api/datasets`
- [ ] No errors in logs: `sudo journalctl -u lab-website-calc -n 50`

### Step 3: Register Datasets (on 192.168.25.2)

### Step 3: Register Datasets

```bash
# On calculation server
cd /var/www/lab-website/backend
source venv/bin/activate

# Option 1: Use registration script
python register_datasets.py

# Option 2: Manual registration via API
curl -X POST http://localhost:5000/api/register \
  -H "Content-Type: application/json" \
  -d '{
    "id": "sample_dataset",
    "path": "/home/jiaxi/spatial_data/sample.h5ad",
    "article": "Research Article Title",
    "name": "Sample Spatial Dataset",
    "description": "Description of the dataset"
  }'

# Verify registration
curl http://192.168.25.2:5000/api/datasets | jq
```

**Verification:**
- [ ] Datasets appear in API response
- [ ] Thumbnails generated in `backend/static/thumbnails/`
- [ ] No errors in Flask logs

---

### Step 4: End-to-End Testing

### Step 4: Update Static Server Nginx Config

If you deployed calculation server on a different machine, update the static server:

```bash
# On static server
sudo nano /etc/nginx/sites-available/pinolilab.org

# Update upstream block:
# upstream calc_server {
#     server <ACTUAL_CALC_SERVER_IP>:5000;
# }

# Test and reload
sudo nginx -t
sudo systemctl reload nginx
```

**Verification:**
- [ ] Nginx config valid: `sudo nginx -t`
- [ ] No errors after reload

---

### Step 5: End-to-End Testing

**Static Pages:**
- [ ] Homepage loads: https://pinolilab.org/
- [ ] Publications page: https://pinolilab.org/publications
- [ ] Research page: https://pinolilab.org/research
- [ ] Contact page: https://pinolilab.org/contact
- [ ] Images load correctly
- [ ] CSS styles applied
- [ ] No console errors

**API Integration:**
- [ ] API accessible through static server: `curl https://pinolilab.org/api/datasets`
- [ ] Articles endpoint: `curl https://pinolilab.org/api/articles`
- [ ] Cache endpoint: `curl https://pinolilab.org/api/cache`
- [ ] No CORS errors in browser console

**Explorer Page:**
- [ ] Explorer loads: https://pinolilab.org/explorer
- [ ] Datasets list displays
- [ ] Can select dataset
- [ ] Spatial plot renders
- [ ] Gene search works
- [ ] Gene expression visualization works
- [ ] No JavaScript errors

**Performance:**
- [ ] Static pages load in < 100ms
- [ ] API responses reasonable (< 5s for large datasets)
- [ ] No timeout errors
- [ ] Memory usage acceptable on calc server

---

## 🔧 Post-Deployment Tasks

### Monitoring Setup
- [ ] Set up log rotation for Nginx logs
- [ ] Set up log rotation for Flask logs
- [ ] Monitor disk space on calc server (for cache)
- [ ] Monitor memory usage on calc server
- [ ] Set up uptime monitoring for both servers

### Security Hardening
- [ ] Verify firewall rules: `sudo ufw status`
- [ ] Ensure Flask port 5000 not publicly accessible
- [ ] Review Nginx security headers
- [ ] Set up fail2ban (optional)
- [ ] Review file permissions

### Documentation
- [ ] Document actual server IPs used
- [ ] Document any configuration changes made
- [ ] Update internal wiki/docs with new architecture
- [ ] Create runbook for common operations

### Backup Strategy
- [ ] Set up automated backups of datasets
- [ ] Set up automated backups of `datasets.json`
- [ ] Set up configuration backups
- [ ] Test restore procedure

---

## 🐛 Troubleshooting

### Static pages show 404 Not Found
```bash
# Check if static files were extracted
ls -la /var/www/lab-website/static-pages/

# Re-run extraction
cd /var/www/lab-website
python3 scripts/extract_static_pages.py \
  backend/templates \
  /var/www/lab-website/static-pages

# Check Nginx logs
sudo tail -f /var/log/nginx/error.log
```

### API requests timeout or fail
```bash
# Check Flask service
sudo systemctl status lab-website-calc
sudo journalctl -u lab-website-calc -n 50

# Test Flask directly
curl http://localhost:5000/health

# Check connectivity from transfer server
# (run on 43.207.99.62)
curl http://192.168.25.2:5000/health

# Check firewall
sudo ufw status
```

### Explorer page won't load data
```bash
# Check browser console for errors
# Common issues:
# 1. CORS errors → Check CORS settings in app_calc_server.py
# 2. 404 errors → Check API proxying in Nginx config
# 3. Timeout → Check Gunicorn timeout settings

# Verify API accessibility
curl https://pinolilab.org/api/datasets

# Check Flask logs
sudo journalctl -u lab-website-calc -f
```

### Out of memory on calculation server
```bash
# Check memory usage
free -h
htop

# Check cache
curl http://192.168.25.2:5000/api/cache

# Clear cache if needed
curl -X POST http://192.168.25.2:5000/api/cache/clear

# Adjust memory limit in systemd service
sudo nano /etc/systemd/system/lab-website-calc.service
# Change MemoryLimit=8G to higher value

sudo systemctl daemon-reload
sudo systemctl restart lab-website-calc
```

---

## 📊 Success Criteria

Deployment is successful when ALL of these are true:

- ✅ Static pages load in browser (/, /publications, /research, /contact)
- ✅ All images and styles load correctly
- ✅ Explorer page loads and displays datasets
- ✅ Can search and visualize genes
- ✅ API endpoints return valid JSON
- ✅ No errors in browser console
- ✅ No errors in Nginx logs
- ✅ No errors in Flask logs
- ✅ SSL certificate valid (green lock in browser)
- ✅ Page load times acceptable
- ✅ Services auto-start on reboot

---

## 🔄 Rollback Plan

If deployment fails and you need to rollback:

### Static Server Rollback
```bash
# Restore old Nginx config
sudo cp /etc/nginx/sites-available/pinolilab.org.backup.<timestamp> \
       /etc/nginx/sites-available/pinolilab.org

# Test and reload
sudo nginx -t
sudo systemctl reload nginx

# Start old Flask if it was stopped
sudo systemctl start lab-website  # Old service name
```

### Calculation Server Rollback
```bash
# Stop new service
sudo systemctl stop lab-website-calc
sudo systemctl disable lab-website-calc

# Remove service file
sudo rm /etc/systemd/system/lab-website-calc.service
sudo systemctl daemon-reload
```

---

## 📞 Support

If you encounter issues:

1. Check the troubleshooting section above
2. Review logs on both servers
3. Consult README_TWO_SERVER.md for detailed information
4. Check GitHub issues
5. Review deployment scripts for configuration errors

---

**Last Updated:** 2026-02-06  
**Version:** 1.0  
**Author:** Deployment automation script
