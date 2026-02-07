#!/bin/bash

#############################################
# Static Server Deployment Script
# Server Role: Serve static content + Reverse proxy
# Technology: Nginx
# Domain: pinolilab.org
# Server IP: 175.41.197.121
#############################################

set -e

echo "=========================================="
echo "Deploying Static Server (pinolilab.org)"
echo "=========================================="
echo ""

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

# Check if running as root
if [ "$EUID" -ne 0 ]; then
    echo -e "${RED}Please run with sudo: sudo bash deploy_static_server.sh${NC}"
    exit 1
fi

# Configuration
DOMAIN="pinolilab.org"
APP_DIR="/var/www/lab-website"
STATIC_DIR="/var/www/lab-website/static-pages"
REPO_URL="https://github.com/wujiaxi/lab-website.git"

# Calculation server address
CALC_SERVER="192.168.25.2:5000"
# Examples:
# CALC_SERVER="192.168.1.100:5000"  # Private IP
# CALC_SERVER="calc.pinolilab.org:5000"  # Domain name

echo -e "${YELLOW}Configuration:${NC}"
echo "Domain: $DOMAIN"
echo "App Directory: $APP_DIR"
echo "Static Directory: $STATIC_DIR"
echo "Calculation Server: $CALC_SERVER"
echo ""
read -p "Press Enter to continue or Ctrl+C to abort..."

echo -e "${YELLOW}Step 1: Installing Nginx...${NC}"
apt update
apt install -y nginx git python3 python3-jinja2
echo -e "${GREEN}✓ Nginx installed${NC}"
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

echo -e "${YELLOW}Step 3: Extracting static HTML pages from Flask templates...${NC}"
mkdir -p "$STATIC_DIR"

# Create a simple Python script to render templates to static HTML
cat > /tmp/extract_static.py << 'EOFPYTHON'
#!/usr/bin/env python3
import os
import sys
from pathlib import Path
from jinja2 import Environment, FileSystemLoader

# Setup paths
template_dir = Path(sys.argv[1])  # backend/templates
static_dir = Path(sys.argv[2])    # /var/www/lab-website/static-pages

# Create Jinja environment
env = Environment(loader=FileSystemLoader(str(template_dir)))

# Templates to render as static pages
static_pages = {
    'index.html': 'index.html',
    'publications.html': 'publications/index.html',
    'research.html': 'research/index.html',
    'contact.html': 'contact/index.html',
}

for template_name, output_path in static_pages.items():
    print(f"Rendering {template_name} -> {output_path}")
    
    # Load and render template
    template = env.get_template(template_name)
    html_content = template.render()
    
    # Write to output
    output_file = static_dir / output_path
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text(html_content)
    print(f"  ✓ Created {output_file}")

print("\nStatic pages extracted successfully!")
EOFPYTHON

# Run extraction script
python3 /tmp/extract_static.py "$APP_DIR/backend/templates" "$STATIC_DIR"
rm /tmp/extract_static.py

# Copy static assets
cp -r "$APP_DIR/backend/static" "$APP_DIR/"
echo -e "${GREEN}✓ Static content extracted${NC}"
echo ""

echo -e "${YELLOW}Step 4: Configuring Nginx...${NC}"

# Backup existing config
if [ -f /etc/nginx/sites-available/pinolilab.org ]; then
    cp /etc/nginx/sites-available/pinolilab.org /etc/nginx/sites-available/pinolilab.org.backup.$(date +%Y%m%d%H%M%S)
fi

# Create Nginx configuration
cat > /etc/nginx/sites-available/pinolilab.org << EOFNGINX
# Upstream for calculation server
upstream calc_server {
    server ${CALC_SERVER};
    keepalive 32;
}

server {
    listen 80;
    listen [::]:80;
   
    server_name pinolilab.org www.pinolilab.org;
   
    # For Let's Encrypt SSL verification
    location ~ /.well-known/acme-challenge {
        allow all;
        root /var/www/pinolilab.org;
    }
   
    # Serve static pages directly
    location = / {
        root ${STATIC_DIR};
        try_files /index.html =404;
    }
    
    location = /publications {
        root ${STATIC_DIR};
        try_files /publications/index.html =404;
    }
    
    location = /research {
        root ${STATIC_DIR};
        try_files /research/index.html =404;
    }
    
    location = /contact {
        root ${STATIC_DIR};
        try_files /contact/index.html =404;
    }
   
    # Serve static assets (CSS, JS, images)
    location /static {
        alias ${APP_DIR}/backend/static;
        expires 30d;
        add_header Cache-Control "public, immutable";
    }
   
    # Proxy explorer page to calculation server
    location /explorer {
        proxy_pass http://calc_server;
        proxy_http_version 1.1;
        proxy_set_header Connection "";
        proxy_set_header Host \$host;
        proxy_set_header X-Real-IP \$remote_addr;
        proxy_set_header X-Forwarded-For \$proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto \$scheme;
        
        # Timeouts for large data transfers
        proxy_connect_timeout 120s;
        proxy_send_timeout 120s;
        proxy_read_timeout 120s;
    }
    
    # Proxy all API requests to calculation server
    location /api {
        proxy_pass http://calc_server;
        proxy_http_version 1.1;
        proxy_set_header Connection "";
        proxy_set_header Host \$host;
        proxy_set_header X-Real-IP \$remote_addr;
        proxy_set_header X-Forwarded-For \$proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto \$scheme;
        
        # Longer timeouts for data processing
        proxy_connect_timeout 180s;
        proxy_send_timeout 180s;
        proxy_read_timeout 180s;
        
        # Buffering settings for large responses
        proxy_buffering on;
        proxy_buffer_size 4k;
        proxy_buffers 8 4k;
        proxy_busy_buffers_size 8k;
    }
}
EOFNGINX

# Enable site
ln -sf /etc/nginx/sites-available/pinolilab.org /etc/nginx/sites-enabled/

# Remove default site
rm -f /etc/nginx/sites-enabled/default

# Test configuration
echo -e "${YELLOW}Testing Nginx configuration...${NC}"
nginx -t

# Reload Nginx
systemctl reload nginx
echo -e "${GREEN}✓ Nginx configured and reloaded${NC}"
echo ""

echo -e "${YELLOW}Step 5: Setting permissions...${NC}"
chown -R www-data:www-data "$APP_DIR"
chmod -R 755 "$APP_DIR"
echo -e "${GREEN}✓ Permissions set${NC}"
echo ""

echo -e "${YELLOW}Step 6: Testing static pages...${NC}"
for page in "/" "/publications" "/research" "/contact"; do
    HTTP_CODE=$(curl -s -o /dev/null -w "%{http_code}" "http://localhost${page}" || echo "000")
    if [ "$HTTP_CODE" = "200" ]; then
        echo -e "${GREEN}✓ $page - HTTP $HTTP_CODE${NC}"
    else
        echo -e "${RED}✗ $page - HTTP $HTTP_CODE${NC}"
    fi
done
echo ""

echo -e "${YELLOW}Step 7: Installing SSL certificate (optional)...${NC}"
read -p "Install SSL certificate with Let's Encrypt? (y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    apt install -y certbot python3-certbot-nginx
    certbot --nginx -d pinolilab.org -d www.pinolilab.org
    echo -e "${GREEN}✓ SSL installed${NC}"
fi
echo ""

echo -e "${GREEN}=========================================="
echo "Static Server Deployment Complete!"
echo "==========================================${NC}"
echo ""
echo "Server Details:"
echo "  - Type: Static + Reverse Proxy"
echo "  - Domain: https://$DOMAIN"
echo "  - Static Pages: $STATIC_DIR"
echo "  - Proxying to: $CALC_SERVER"
echo ""
echo "Service Management:"
echo "  - Test Nginx: sudo nginx -t"
echo "  - Reload Nginx: sudo systemctl reload nginx"
echo "  - View logs: sudo tail -f /var/log/nginx/access.log"
echo ""
echo "Update Static Content:"
echo "  cd $APP_DIR"
echo "  git pull origin main"
echo "  python3 scripts/extract_static.py"
echo "  sudo systemctl reload nginx"
echo ""
echo -e "${YELLOW}NEXT STEP: Deploy calculation server at $CALC_SERVER${NC}"
echo ""
