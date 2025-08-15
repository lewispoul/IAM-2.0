#!/bin/bash

# fetch_ketcher.sh - Download and extract Ketcher build locally
# This script downloads a stable Ketcher release and extracts it to public/static/ketcher/

set -e

# Configuration
KETCHER_VERSION="v2.27.1"
KETCHER_URL="https://github.com/epam/ketcher/releases/download/${KETCHER_VERSION}/ketcher-${KETCHER_VERSION#v}.tar.gz"
KETCHER_DIR="public/static/ketcher"
TEMP_DIR="/tmp/ketcher_download"

echo "🧪 Fetching Ketcher ${KETCHER_VERSION} for local hosting..."

# Check if already present (idempotent)
if [[ -f "${KETCHER_DIR}/index.html" ]]; then
    echo "✅ Ketcher already present at ${KETCHER_DIR}/index.html"
    echo "🌐 Local URL: http://localhost:8010/static/ketcher/index.html"
    exit 0
fi

# Create directories
echo "📁 Creating directory structure..."
mkdir -p "${KETCHER_DIR}"
mkdir -p "${TEMP_DIR}"

# Download Ketcher release
echo "⬇️  Downloading Ketcher ${KETCHER_VERSION}..."
cd "${TEMP_DIR}"

if ! curl -L -o "ketcher.tar.gz" "${KETCHER_URL}"; then
    echo "❌ Failed to download Ketcher from ${KETCHER_URL}"
    echo "💡 Check if the version ${KETCHER_VERSION} exists at: https://github.com/epam/ketcher/releases"
    exit 1
fi

# Extract archive
echo "📦 Extracting Ketcher..."
tar -xzf "ketcher.tar.gz"

# Find extracted directory (may have different naming)
EXTRACTED_DIR=$(find . -maxdepth 1 -type d -name "*ketcher*" | head -1)
if [[ -z "$EXTRACTED_DIR" ]]; then
    echo "❌ Could not find extracted Ketcher directory"
    ls -la
    exit 1
fi

# Copy to target location
echo "📋 Copying files to ${KETCHER_DIR}..."
cd "$EXTRACTED_DIR"
cp -r . "../../../${KETCHER_DIR}/"

# Verify installation
cd "../../../"
if [[ -f "${KETCHER_DIR}/index.html" ]]; then
    echo "✅ Ketcher successfully installed!"
    echo "📁 Location: ${KETCHER_DIR}/"
    echo "🌐 Local URL: http://localhost:8010/static/ketcher/index.html"
    
    # Show directory contents
    echo "📋 Contents:"
    ls -la "${KETCHER_DIR}/" | head -10
else
    echo "❌ Installation failed - index.html not found"
    echo "📋 Contents of ${KETCHER_DIR}:"
    ls -la "${KETCHER_DIR}/"
    exit 1
fi

# Clean up
echo "🧹 Cleaning up temporary files..."
rm -rf "${TEMP_DIR}"

echo "🎉 Ketcher local hosting setup complete!"
