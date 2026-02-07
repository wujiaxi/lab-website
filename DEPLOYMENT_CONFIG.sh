#!/bin/bash
#############################################
# Deployment Configuration
# Edit this file before running deployment scripts
#############################################

# ============================================
# STATIC SERVER CONFIGURATION
# ============================================

# Domain name
export DOMAIN="pinolilab.org"

# Calculation server address (IP:PORT)
# Options:
#   - Same machine: "localhost:5000" or "127.0.0.1:5000"
#   - Private network: "192.168.1.100:5000"
#   - Public IP: "12.34.56.78:5000"
#   - Domain name: "calc.pinolilab.org:5000"
export CALC_SERVER="localhost:5000"  # ⚠️ CHANGE THIS

# ============================================
# CALCULATION SERVER CONFIGURATION
# ============================================

# Directory where H5AD files are stored
export DATA_DIR="/home/jiaxi/spatial_data"  # ⚠️ VERIFY THIS PATH

# IP address of static server (for firewall rules)
# Used to restrict access to Flask backend
export STATIC_SERVER_IP="175.41.197.121"  # ⚠️ CHANGE THIS IF NEEDED

# ============================================
# FLASK CONFIGURATION
# ============================================

# Number of Gunicorn workers
# Recommendation: (2 × CPU cores) + 1
# Example: 4-core CPU → 9 workers
export GUNICORN_WORKERS="4"

# Request timeout in seconds
# Increase for large datasets
export GUNICORN_TIMEOUT="300"

# Memory limit for Flask service
# Adjust based on dataset sizes
export MEMORY_LIMIT="8G"

# Cache expiry in hours
# How long to keep datasets in memory
export CACHE_EXPIRY_HOURS="6"

# ============================================
# OPTIONAL: SSL CONFIGURATION
# ============================================

# Email for Let's Encrypt SSL certificate
export SSL_EMAIL=""  # ⚠️ ADD YOUR EMAIL

# Enable automatic SSL setup? (yes/no)
export AUTO_SSL="yes"

# ============================================
# DEPLOYMENT OPTIONS
# ============================================

# Skip firewall configuration? (yes/no)
export SKIP_FIREWALL="no"

# Enable verbose logging? (yes/no)
export VERBOSE="yes"

# ============================================
# VALIDATION
# ============================================

echo "=========================================="
echo "Deployment Configuration Summary"
echo "=========================================="
echo ""
echo "Static Server:"
echo "  Domain: $DOMAIN"
echo "  Calculation Server: $CALC_SERVER"
echo "  SSL Email: ${SSL_EMAIL:-'Not set'}"
echo ""
echo "Calculation Server:"
echo "  Data Directory: $DATA_DIR"
echo "  Allowed IP: $STATIC_SERVER_IP"
echo "  Workers: $GUNICORN_WORKERS"
echo "  Timeout: ${GUNICORN_TIMEOUT}s"
echo "  Memory Limit: $MEMORY_LIMIT"
echo ""
echo "=========================================="
echo ""

# Validation warnings
if [ "$CALC_SERVER" = "localhost:5000" ]; then
    echo "⚠️  WARNING: CALC_SERVER is set to localhost"
    echo "   This only works if both servers are on the same machine"
    echo ""
fi

if [ ! -d "$DATA_DIR" ]; then
    echo "⚠️  WARNING: DATA_DIR does not exist: $DATA_DIR"
    echo "   Make sure to create this directory or update the path"
    echo ""
fi

if [ -z "$SSL_EMAIL" ]; then
    echo "⚠️  WARNING: SSL_EMAIL is not set"
    echo "   SSL setup will be skipped unless you provide an email"
    echo ""
fi

echo "To use this configuration, run:"
echo "  source DEPLOYMENT_CONFIG.sh"
echo "  sudo -E bash scripts/deploy_static_server.sh"
echo "  sudo -E bash scripts/deploy_calculation_server.sh"
echo ""
