#!/bin/bash

# fetch_ketcher.sh - Download and setup Ketcher build locally
# This script creates a minimal working Ketcher installation

set -e

# Configuration
KETCHER_DIR="public/static/ketcher"

echo "🧪 Setting up minimal Ketcher for local hosting..."

# Check if already present (idempotent)
if [[ -f "${KETCHER_DIR}/index.html" ]]; then
    echo "✅ Ketcher already present at ${KETCHER_DIR}/index.html"
    echo "🌐 Local URL: http://localhost:8011/static/ketcher/index.html"
    exit 0
fi

# Create directories
echo "📁 Creating directory structure..."
mkdir -p "${KETCHER_DIR}"

# Create a minimal Ketcher HTML wrapper that uses CDN but serves locally
echo "📄 Creating Ketcher HTML wrapper..."
cat > "${KETCHER_DIR}/index.html" << 'EOF'
<!DOCTYPE html>
<html>
<head>
    <title>Ketcher - Molecular Editor</title>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <link rel="stylesheet" href="https://unpkg.com/ketcher@2.19.0/dist/ketcher.css">
    <style>
        body { margin: 0; padding: 0; }
        #ketcher-root { width: 100%; height: 100vh; }
    </style>
</head>
<body>
    <div id="ketcher-root"></div>
    <script src="https://unpkg.com/ketcher@2.19.0/dist/ketcher.js"></script>
    <script>
        window.addEventListener('DOMContentLoaded', function() {
            window.ketcher = new window.Ketcher('ketcher-root', {
                // Configuration options
                buttons: ['file', 'edit', 'view', 'zoom', 'about'],
                staticResourcesUrl: 'https://unpkg.com/ketcher@2.19.0/dist/'
            });
        });
        
        // Bridge functions for parent window communication
        window.getMolfile = function() {
            return window.ketcher ? window.ketcher.getMolfile() : '';
        };
        
        window.setMolfile = function(molfile) {
            if (window.ketcher) {
                return window.ketcher.setMolfile(molfile);
            }
        };
        
        window.getSmiles = function() {
            return window.ketcher ? window.ketcher.getSmiles() : '';
        };
        
        window.setSmiles = function(smiles) {
            if (window.ketcher) {
                return window.ketcher.setSmiles(smiles);
            }
        };
        
        window.clear = function() {
            if (window.ketcher) {
                window.ketcher.setMolfile('');
            }
        };
    </script>
</body>
</html>
EOF

# Verify installation
if [[ -f "${KETCHER_DIR}/index.html" ]]; then
    echo "✅ Ketcher successfully installed!"
    echo "📁 Location: ${KETCHER_DIR}/"
    echo "🌐 Local URL: http://localhost:8011/static/ketcher/index.html"
    echo "💡 This uses CDN resources but serves locally to avoid CORS issues"
else
    echo "❌ Installation failed - index.html not found"
    exit 1
fi

echo "🎉 Ketcher local hosting setup complete!"
