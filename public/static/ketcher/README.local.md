# Local Ketcher Installation

This directory contains a local copy of Ketcher molecular editor for IAM-2.0.

## Setup

To fetch/update Ketcher locally:

```bash
bash scripts/fetch_ketcher.sh
```

This will:

- Download the latest stable Ketcher release from EPAM
- Extract to `public/static/ketcher/`
- Make it available at `http://localhost:8010/static/ketcher/index.html`

## Usage

1. Run the IAM-2.0 backend: `make run`
2. Open Ketcher at: `http://localhost:8010/static/ketcher/index.html`
3. Use the main UI at: `http://localhost:8010/ketcher.html`

## Notes

- Ketcher files are ignored by git (see `.gitignore`)
- The fetch script is idempotent - safe to re-run
- Uses stable release v2.27.1 by default
- No cross-site iframe issues since it's hosted locally
