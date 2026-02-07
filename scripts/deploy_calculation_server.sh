#!/bin/bash

#############################################
# Calculation Server Deployment Script
# Server Role: Flask backend for H5AD data processing
# Technology: Python + Flask + Gunicorn
# Handles: /api/* and /explorer endpoints
#############################################

set -e

echo "=========================================="
echo "Deploying Calculation Server (Flask Backend)"
echo "=========================================="
echo ""

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

# Check if running as root
if [ "$EUID" -ne 0 ]; then
    echo -e "${RED}Please run with sudo: sudo bash deploy_calculation_server.sh${NC}"
    exit 1
fi

# Configuration
APP_DIR="/var/www/lab-website"
REPO_URL="https://github.com/wujiaxi/lab-website.git"
PYTHON_VERSION="python3"
DATA_DIR="/home/jiaxi/spatial_data"
SERVICE_NAME="lab-website-calc"

# Allow connections from static/transfer server
ALLOWED_IP="175.41.197.121"

echo -e "${YELLOW}Configuration:${NC}"
echo "App Directory: $APP_DIR"
echo "Data Directory: $DATA_DIR"
echo "Python Version: $PYTHON_VERSION"
echo "Allowed IP: $ALLOWED_IP (for firewall)"
echo ""
read -p "Press Enter to continue or Ctrl+C to abort..."

echo -e "${YELLOW}Step 1: Installing system dependencies...${NC}"
apt update
apt install -y python3 python3-pip python3-venv git ufw
echo -e "${GREEN}✓ System dependencies installed${NC}"
echo ""

echo -e "${YELLOW}Step 2: Cloning/updating repository...${NC}"
if [ -d "$APP_DIR" ]; then
    echo "Directory exists, pulling latest changes..."
    cd "$APP_DIR"
    git pull origin main
else
    echo "Cloning repository..."
    git clone "$REPO_URL" "$APP_DIR"
    cd "$APP_DIR"
fi
echo -e "${GREEN}✓ Repository ready${NC}"
echo ""

echo -e "${YELLOW}Step 3: Setting up Python virtual environment...${NC}"
cd "$APP_DIR/backend"

# Remove old venv if exists
if [ -d "venv" ]; then
    echo "Removing old virtual environment..."
    rm -rf venv
fi

$PYTHON_VERSION -m venv venv
source venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install dependencies
pip install -r requirements.txt
pip install gunicorn

deactivate
echo -e "${GREEN}✓ Python environment ready${NC}"
echo ""

echo -e "${YELLOW}Step 4: Creating systemd service...${NC}"
cat > /etc/systemd/system/${SERVICE_NAME}.service << EOFSVC
[Unit]
Description=Lab Website Calculation Server (Flask)
After=network.target

[Service]
Type=notify
User=www-data
Group=www-data
WorkingDirectory=${APP_DIR}/backend
Environment="PATH=${APP_DIR}/backend/venv/bin"
Environment="DATA_DIR=${DATA_DIR}"
Environment="PYTHONUNBUFFERED=1"

# Gunicorn with optimized settings for data processing
ExecStart=${APP_DIR}/backend/venv/bin/gunicorn \\
    --workers 4 \\
    --worker-class sync \\
    --bind 0.0.0.0:5000 \\
    --timeout 300 \\
    --graceful-timeout 60 \\
    --max-requests 1000 \\
    --max-requests-jitter 50 \\
    --access-logfile - \\
    --error-logfile - \\
    --log-level info \\
    app:app

# Restart policy
Restart=always
RestartSec=10

# Resource limits (adjust based on your server)
MemoryLimit=8G
CPUQuota=200%

[Install]
WantedBy=multi-user.target
EOFSVC

echo -e "${GREEN}✓ Systemd service created${NC}"
echo ""

echo -e "${YELLOW}Step 5: Setting up data directory...${NC}"
if [ ! -d "$DATA_DIR" ]; then
    echo "Creating data directory: $DATA_DIR"
    mkdir -p "$DATA_DIR"
fi

# Set permissions
chown -R www-data:www-data "$APP_DIR"
chown -R www-data:www-data "$DATA_DIR"
chmod -R 755 "$APP_DIR"
chmod -R 755 "$DATA_DIR"

echo -e "${GREEN}✓ Data directory ready${NC}"
echo ""

echo -e "${YELLOW}Step 6: Configuring firewall (optional)...${NC}"
read -p "Configure firewall to only allow static server? (y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    # Enable UFW
    ufw --force enable
    
    # Allow SSH (important!)
    ufw allow 22/tcp
    
    # Allow port 5000 only from static server
    ufw allow from $ALLOWED_IP to any port 5000 proto tcp
    
    # If this IS the static server (same machine), allow localhost
    ufw allow from 127.0.0.1 to any port 5000 proto tcp
    
    ufw reload
    echo -e "${GREEN}✓ Firewall configured${NC}"
    ufw status
else
    echo "Skipping firewall configuration"
fi
echo ""

echo -e "${YELLOW}Step 7: Starting Flask service...${NC}"
systemctl daemon-reload
systemctl enable ${SERVICE_NAME}
systemctl restart ${SERVICE_NAME}

# Wait for service to start
sleep 3

# Check status
if systemctl is-active --quiet ${SERVICE_NAME}; then
    echo -e "${GREEN}✓ Flask service is running${NC}"
else
    echo -e "${RED}✗ Flask service failed to start${NC}"
    echo "Check logs with: sudo journalctl -u ${SERVICE_NAME} -n 50"
    exit 1
fi
echo ""

echo -e "${YELLOW}Step 8: Testing Flask application...${NC}"

# Test health endpoint
HTTP_CODE=$(curl -s -o /dev/null -w "%{http_code}" http://localhost:5000/api/articles || echo "000")
if [ "$HTTP_CODE" = "200" ]; then
    echo -e "${GREEN}✓ API responding (HTTP $HTTP_CODE)${NC}"
else
    echo -e "${YELLOW}⚠ API returned HTTP $HTTP_CODE${NC}"
    echo "This might be normal if no datasets are registered yet"
fi

# Test cache endpoint
HTTP_CODE=$(curl -s -o /dev/null -w "%{http_code}" http://localhost:5000/api/cache || echo "000")
if [ "$HTTP_CODE" = "200" ]; then
    echo -e "${GREEN}✓ Cache endpoint OK${NC}"
fi

echo ""

echo -e "${GREEN}=========================================="
echo "Calculation Server Deployment Complete!"
echo "==========================================${NC}"
echo ""
echo "Server Details:"
echo "  - Type: Flask Calculation Server"
echo "  - Binding: 0.0.0.0:5000"
echo "  - Data Directory: $DATA_DIR"
echo "  - Workers: 4"
echo "  - Timeout: 300s"
echo ""
echo "Service Management:"
echo "  - Check status: sudo systemctl status ${SERVICE_NAME}"
echo "  - Restart: sudo systemctl restart ${SERVICE_NAME}"
echo "  - View logs: sudo journalctl -u ${SERVICE_NAME} -f"
echo "  - Stop: sudo systemctl stop ${SERVICE_NAME}"
echo ""
echo "API Endpoints:"
echo "  - http://localhost:5000/api/articles"
echo "  - http://localhost:5000/api/datasets"
echo "  - http://localhost:5000/api/cache"
echo "  - http://localhost:5000/explorer"
echo ""
echo "Update Application:"
echo "  cd $APP_DIR"
echo "  git pull origin main"
echo "  sudo systemctl restart ${SERVICE_NAME}"
echo ""
echo "Register H5AD Datasets:"
echo "  cd $APP_DIR/backend"
echo "  source venv/bin/activate"
echo "  python register_datasets.py"
echo ""
echo -e "${YELLOW}IMPORTANT: Update static server Nginx config with this server's IP!${NC}"
echo ""
