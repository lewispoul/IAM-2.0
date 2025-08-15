## Quick Start: Environment Setup

> **Note:** The `iam2` conda environment is now set to auto-activate for all new terminal sessions via your `~/.bashrc` file. You do not need to manually activate it each time.

```sh
# Create and activate the IAM2 environment (recommended)
mamba create -n iam2 -f chem-env.yaml   # or: conda env create -f chem-env.yaml
conda activate iam2
make verify    # Check required packages
make run       # Start backend API
```

## Local Ketcher (no iframe blocking)

IAM-2.0 now includes local Ketcher hosting to avoid cross-site iframe issues:

```bash
bash scripts/fetch_ketcher.sh    # Download Ketcher locally
make run                         # Start backend
open http://localhost:8011/      # Main UI with local Ketcher
```

**Available interfaces:**

- **Main UI**: `http://localhost:8011/` - Full molecular analysis interface
- **Direct Ketcher**: `http://localhost:8011/static/ketcher/index.html` - Molecular editor only
- **API Docs**: `http://localhost:8011/docs` - Interactive API documentation

**Troubleshooting:**

- If "Firefox Can't Open This Page": Use Chrome/Chromium instead
- If Ketcher doesn't load: Check `/static/ketcher` is ignored by git
- Static files are served from `public/static/ketcher/` (vendor files ignored)


# IAM 2.0 — Intelligent Agent for Molecules
![IAM UI](docs/images/IAM-2.0UI.png)

![IAM at a Glance](docs/images/iam_at_a_glance.png)

**Model. Simulate. Predict. Visualize.**  
A clean, modular platform for computational chemistry with a focus on energetic materials — built for local devices and HPC.

- ⚙️ Engines: XTB (local), Psi4 (local), Gaussian-ready (HPC)
- 📈 Predictors: Kamlet–Jacobs, Keshavarz, ML (VoD/Pcj/ΔHdet)
- 🧪 Pipelines: SMILES→XYZ→Optimization→Properties→Prediction
- 🖥️ UI: Molecule viewer, Orbitals, Spectra, Performance
- 🧰 Dev: FastAPI, Pydantic, RDKit, 3Dmol.js, GitHub Actions

📘 Read the full whitepaper: **[`/docs/IAM_2.0_Whitepaper.md`](docs/IAM_2.0_Whitepaper.md)**

## Repository layout

IAM-2.0/
  api/
    __init__.py
    xtb_routes.py
    psi4_routes.py
    empirical_routes.py
    cj_routes.py
  docs/
    API.md
    README_UI.md
    dev_notes.md
  iam/
    __init__.py
    backend/
      __init__.py
      utils/
        __init__.py
        persistence.py
    core/
      __init__.py
      empirical/
        __init__.py
        predictors.py
    runners/
      __init__.py
      xtb.py
      psi4.py
      empirical.py
      cj.py
  installers/
    Miniconda3-latest-Linux-x86_64.sh
  scripts/
    __init__.py
    verify_env.py
    README.md
  tests/
    __init__.py
    api/
      test_routes.py
    core/
      test_empirical.py
    runners/
      test_cj.py
    utils/
      test_persistence.py
  REPORTS/
  .copilot-instructions.md
  .gitignore
  chem-env.yaml
  Makefile
  README.md
  pyproject.toml
  Dockerfile
  docker-compose.yml

---

# IAM-2.0 × NOX — Consolidated Integration Plan & Status (2025-08-15)

## Objectives
- IAM-2.0 is a standalone compute API; NOX consumes via OpenAPI v1.
- One normalized response shape everywhere: {"ok": bool, "data": object|null, "errors": []}.
- Conda/Mamba for local dev; Docker is optional and does not disrupt conda flow.
- Ketcher is frontend-only (CDN); all chemistry and persistence are backend.

## Current State
- FastAPI app wired: /compute/{xtb,psi4,empirical,cj}, /convert/molfile, /ketcher/{to-smiles,to-xyz,run}.
- Persistence utilities: success-only writes under IAM_Knowledge (env override: IAM_RESULTS_BASE).
- All tests pass for current routes; docs and progress log updated.
- OpenAPI v1 spec prepared and served at /v1/openapi.json.

## Source of Truth
- OpenAPI: openapi/iam2.v1.yaml is the contract NOX will ingest.
- Error envelope: Always normalized schema; functional errors → ok=false and HTTP 200. Auth/transport → 4xx/5xx.
- Persistence: Write JSON snapshots to Results/ and append benchmark_auto.csv only on ok=true.

## Environments
- Local dev: Conda env iam2 (working). Makefile targets intact.
- Docker (optional): Existing Dockerfiles untouched. New Docker assets go under /deploy/ and use micromamba. Makefile and conda docs unchanged.

## Ketcher Integration
- Minimal HTML page (or NOX tab) embeds Ketcher via CDN.
- Flow: get molfile → POST /v1/convert/molfile → /v1/compute/* (or /v1/ketcher/run).
- REST/curl examples provided under examples/requests/.

## CI/CD
- IAM: lint + typecheck + pytest → docker image (optional) → preview deploy → publish OpenAPI artifact.
- NOX: import OpenAPI → regenerate SDKs → contract tests vs IAM preview.

## Guardrails
- Work only inside IAM-2.0 repo.
- Never overwrite or delete existing Docker files; new go under /deploy/.
- No large binaries or vendor JS bundles (use CDN for Ketcher).
- Small, Conventional Commits; append summary to REPORTS/IAM2_min_stack.md.

## Next Steps
1. OpenAPI v1 polish and validation.
2. Envelope conformance and shape tests.
3. Persistence hooks and /export/zip safety.
4. Minimal Ketcher UI or REST examples.
5. Docs update: API.md, README, persistence notes, examples.
6. Optional Docker: new assets under /deploy/.
7. Finalize & report: pytest, operator commands, summary in REPORTS/IAM2_min_stack.md.

## Operator Commands
```
make run
# then open http://localhost:8010/ketcher.html (if UI present)
# or use examples/requests/*.http
```

---
