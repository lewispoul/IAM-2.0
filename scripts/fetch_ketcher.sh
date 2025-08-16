#!/usr/bin/env bash
set -euo pipefail
DEST="IAM_GUI/static/ketcher"
mkdir -p "$DEST"
echo "[*] Fetching ketcher-standalone via npm (requires node/npm)."
TMP="$(mktemp -d)"
pushd "$TMP" >/dev/null
  npm init -y >/dev/null 2>&1 || true
  npm i ketcher-standalone --no-audit --no-fund
  # Copy a simple distribution; adjust path if package layout changes
  # Prefer the built app if available; fallback to dist/
  if [ -d node_modules/ketcher-standalone/dist ]; then
    rsync -a node_modules/ketcher-standalone/dist/ "$DEST"/
  else
    rsync -a node_modules/ketcher-standalone/ "$DEST"/
  fi
popd >/dev/null
rm -rf "$TMP"
echo "[+] Ketcher files placed in $DEST"
echo "[i] If you prefer a fixed release, download a 'standalone' ZIP from EPAM releases and unzip to $DEST."
