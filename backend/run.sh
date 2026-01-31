#!/bin/bash
# Start the Flask backend server

# Configuration
HOST="0.0.0.0"
PORT=5000
WORKERS=4

echo "🚀 Starting Spatial Data Explorer Backend"
echo "   Host: $HOST"
echo "   Port: $PORT"
echo ""

# Check if running in production or development
if [ "$1" = "prod" ]; then
    echo "Running in PRODUCTION mode with gunicorn..."
    gunicorn -w $WORKERS -b $HOST:$PORT app:app
else
    echo "Running in DEVELOPMENT mode..."
    python3 app.py
fi
