# IAM2.0 Complete Directory Tree Analysis
**Generated:** August 16, 2025  
**Purpose:** Post-refactor cleanup analysis to identify duplicates and validate structure

## Directory Structure Overview

```text
IAM-2.0/
├── .copilot-instructions.md                   # GitHub Copilot configuration
├── .git/                                      # Git repository (excluded from detailed listing)
├── .github/
│   └── workflows/
│       └── ci.yml                            # CI/CD pipeline configuration
├── .gitignore                                # Git ignore rules
├── .pre-commit-config.yaml                   # Pre-commit hooks configuration
├── .pytest_cache/                           # Pytest cache directory
│   ├── .gitignore
│   ├── CACHEDIR.TAG
│   ├── README.md
│   └── v/
│       └── cache/
│           ├── lastfailed
│           └── nodeids
├── .vscode/
│   └── settings.json                         # VS Code workspace settings
├── CHANGELOG.md                              # Project changelog
├── CONTRIBUTING.md                           # Contribution guidelines
├── Dockerfile                                # Main Docker container configuration
├── IAM_GUI/
│   └── templates/
│       └── iam_viewer_connected.html         # Frontend HTML template with 3Dmol.js integration
├── IAM_Knowledge/
│   ├── Exports/                              # (Empty directory)
│   ├── Results/                              # Chemistry calculation results (77+ JSON files)
│   │   ├── 20250815_174521_convert.json
│   │   ├── 20250815_174521_ketcher_run.json
│   │   ├── 20250815_174544_cj.json
│   │   ├── 20250815_174544_convert.json
│   │   ├── 20250815_174544_empirical.json
│   │   ├── 20250815_174544_ketcher_run.json
│   │   ├── [...77+ more result files...]
│   │   └── 20250815_191518_ketcher_run.json
│   └── benchmark_auto.csv                    # Automated benchmark data
├── LICENSE                                   # Project license
├── Makefile                                  # Build and development commands
├── PHASE_COMPLETE_KETCHER_INTEGRATION.md     # Ketcher integration completion report
├── README.md                                 # Main project documentation
├── REPORTS/
│   └── IAM2_min_stack.md                    # Minimal stack implementation report
├── SESSION_PROGRESS_2025-08-15.md           # Development session progress
├── api/                                      # 🔴 DUPLICATE: Legacy API structure
│   ├── __init__.py
│   ├── __pycache__/
│   │   ├── __init__.cpython-310.pyc
│   │   └── main.cpython-310.pyc
│   ├── cj_routes.py
│   ├── empirical_routes.py
│   ├── main.py
│   ├── psi4_routes.py
│   └── xtb_routes.py
├── backend/                                  # ✅ NEW: FastAPI backend structure
│   ├── api/
│   │   ├── calc.py                          # Chemistry calculation endpoints
│   │   ├── convert.py                       # Molecular conversion endpoints
│   │   └── predict.py                       # Performance prediction endpoints
│   ├── converters/
│   │   └── rdkit.py                         # RDKit molecular format converters
│   └── jobs/
│       ├── psi4_stub.py                     # Psi4 calculation stub
│       └── xtb_stub.py                      # XTB calculation stub
├── chem-env.yaml                            # Conda environment specification
├── docker-compose.dev.yml                   # Development Docker Compose
├── docker-compose.yml                       # Production Docker Compose
├── docker/
│   └── nginx/
│       └── dev.conf                         # Nginx development configuration
├── docs/                                    # Comprehensive documentation
│   ├── API.md                              # API documentation
│   ├── IAM_2.0_Whitepaper.md              # Technical whitepaper
│   ├── IAM_Master_Technical_Report_FINAL.md # Master technical report
│   ├── MIGRATION.md                        # Migration guide
│   ├── README.md                           # Documentation index
│   ├── images/                             # Documentation images (6 files)
│   │   ├── IAM-2.0UI.png
│   │   ├── iam_at_a_glance.png
│   │   ├── iam_case_homo_lumo.png
│   │   ├── iam_case_vod_output.png
│   │   ├── iam_data_pipeline.png
│   │   └── iam_dependency_timeline.png
│   ├── plans/                              # Planning documents
│   │   ├── IAM2.0_master_plan.md
│   │   ├── IAM2.0_risks.md
│   │   ├── IAM2.0_todo.md
│   │   ├── endpoint_migration_plan.md
│   │   └── ui_migration_plan.md
│   ├── reports/                            # Analysis reports
│   │   ├── IAM2.0_exec_overview.md
│   │   ├── IAM2.0_feature_diff.csv
│   │   ├── cicd_and_deploy.md
│   │   ├── env_requirements_unified.md
│   │   ├── file_io_map.md
│   │   ├── iam_conversion_endpoints.md
│   │   ├── iam_gui_blueprint.md
│   │   ├── iam_ketcher_integration_points.md
│   │   ├── ketcher_bundle_audit.md
│   │   ├── nox_api_inventory.md
│   │   ├── nox_chemistry_capabilities.md
│   │   ├── nox_data_and_jobs.md
│   │   ├── runbooks_and_scripts.md
│   │   └── xtb_psi4_pipeline_trace.md
│   ├── specs/
│   │   └── ketcher_postmessage_contract.md  # Ketcher integration specification
│   └── tests/
│       ├── TESTPLAN_IAM2.0.md              # Comprehensive test plan
│       ├── api/                            # (Empty)
│       ├── chem/                           # (Empty)
│       └── ui_contracts/                   # (Empty)
├── examples/                                # API usage examples
│   ├── README.md
│   └── requests/                           # HTTP request examples
│       ├── compute.http
│       ├── convert.http
│       ├── export.http
│       ├── ketcher_local.http
│       ├── ketcher_run.http
│       ├── ketcher_to_smiles.http
│       └── ketcher_to_xyz.http
├── frontend/                                # 🔴 DUPLICATE: Legacy frontend files
│   ├── ketcher_test.html
│   ├── standalone_test.html
│   └── test_success.html
├── iam/                                     # 🔴 DUPLICATE: Original Flask backend
│   ├── __init__.py
│   ├── __pycache__/
│   │   └── __init__.cpython-310.pyc
│   ├── backend/
│   │   ├── __init__.py
│   │   ├── __pycache__/
│   │   │   ├── __init__.cpython-310.pyc
│   │   │   └── app.cpython-310.pyc
│   │   ├── app.py                          # Flask application
│   │   ├── routes/                         # Flask routes
│   │   │   ├── __pycache__/
│   │   │   │   ├── cj.cpython-310.pyc
│   │   │   │   ├── convert.cpython-310.pyc
│   │   │   │   ├── empirical.cpython-310.pyc
│   │   │   │   ├── export.cpython-310.pyc
│   │   │   │   ├── ketcher.cpython-310.pyc
│   │   │   │   ├── psi4.cpython-310.pyc
│   │   │   │   └── xtb.cpython-310.pyc
│   │   │   ├── cj.py
│   │   │   ├── convert.py
│   │   │   ├── empirical.py
│   │   │   ├── export.py
│   │   │   ├── ketcher.py
│   │   │   ├── psi4.py
│   │   │   └── xtb.py
│   │   ├── schemas/
│   │   │   ├── __pycache__/
│   │   │   │   ├── common.cpython-310.pyc
│   │   │   │   └── error.cpython-310.pyc
│   │   │   ├── common.py
│   │   │   └── error.py
│   │   └── utils/
│   │       ├── __init__.py
│   │       ├── __pycache__/
│   │       │   ├── __init__.cpython-310.pyc
│   │       │   └── persistence.cpython-310.pyc
│   │       └── persistence.py
│   ├── core/
│   │   ├── __init__.py
│   │   ├── __pycache__/
│   │   │   └── __init__.cpython-310.pyc
│   │   └── empirical/
│   │       ├── __init__.py
│   │       ├── __pycache__/
│   │       │   ├── __init__.cpython-310.pyc
│   │       │   └── iam_empirical_predictor.cpython-310.pyc
│   │       ├── iam_empirical_predictor.py
│   │       └── predictors.py
│   └── runners/
│       ├── __init__.py
│       ├── __pycache__/
│       │   ├── __init__.cpython-310.pyc
│       │   └── cantera_cj.cpython-310.pyc
│       ├── cantera_cj.py
│       ├── cj.py
│       ├── empirical.py
│       ├── psi4.py
│       └── xtb.py
├── installers/
│   └── Miniconda3-latest-Linux-x86_64.sh    # Miniconda installer
├── integration_test_final.py                # Final integration test
├── openapi/
│   └── iam2.v1.yaml                        # OpenAPI specification
├── public/                                  # Static web assets
│   ├── ketcher.html
│   └── static/
│       └── ketcher/
│           ├── README.local.md
│           ├── bridge.js
│           └── index.html
├── pyproject.toml                          # Python project configuration
├── requirements.txt                        # Python dependencies
├── scripts/                                # Operational scripts
│   ├── __init__.py
│   ├── fetch_ketcher.sh                   # Ketcher asset fetching
│   ├── run_backend.sh                     # Backend startup script
│   ├── run_tests.sh                       # Test execution script
│   ├── validate_openapi.py                # OpenAPI validation
│   └── verify_env.py                      # Environment verification
├── src/
│   └── iam/
│       └── __init__.py                    # Package initialization
├── test_ui_integration.py                  # UI integration test
├── tests/                                  # ✅ Comprehensive test suite
│   ├── __init__.py
│   ├── __pycache__/
│   │   └── __init__.cpython-310.pyc
│   ├── api/                               # API tests
│   │   ├── __init__.py
│   │   ├── __pycache__/
│   │   │   ├── __init__.cpython-310.pyc
│   │   │   ├── test_conformance.cpython-310-pytest-8.4.1.pyc
│   │   │   ├── test_convert.cpython-310-pytest-8.4.1.pyc
│   │   │   ├── test_envelope_edge_cases.cpython-310-pytest-8.4.1.pyc
│   │   │   ├── test_export.cpython-310-pytest-8.4.1.pyc
│   │   │   ├── test_http_examples.cpython-310-pytest-8.4.1.pyc
│   │   │   ├── test_ketcher.cpython-310-pytest-8.4.1.pyc
│   │   │   ├── test_ketcher_smoke.cpython-310-pytest-8.4.1.pyc
│   │   │   ├── test_persistence_safety.cpython-310-pytest-8.4.1.pyc
│   │   │   └── test_routes.cpython-310-pytest-8.4.1.pyc
│   │   ├── test_conformance.py
│   │   ├── test_convert.py
│   │   ├── test_convert_and_ketcher.py
│   │   ├── test_envelope_edge_cases.py
│   │   ├── test_export.py
│   │   ├── test_export_and_persistence.py
│   │   ├── test_http_examples.py
│   │   ├── test_ketcher.py
│   │   ├── test_ketcher_run.py
│   │   ├── test_ketcher_smoke.py
│   │   ├── test_molfile_to_xyz.py
│   │   ├── test_persistence_safety.py
│   │   ├── test_psi4_stub.py
│   │   ├── test_routes.py
│   │   ├── test_smiles_to_xyz.py
│   │   ├── test_static_local_ketcher.py
│   │   └── test_xtb_stub.py
│   ├── chem/
│   │   └── test_regression_snapshots.py
│   ├── conftest.py                        # Pytest configuration
│   ├── core/
│   │   └── __init__.py
│   ├── integration/
│   │   └── test_ketcher_live.py
│   ├── runners/
│   │   ├── __init__.py
│   │   ├── __pycache__/
│   │   │   ├── __init__.cpython-310.pyc
│   │   │   └── test_cj.cpython-310-pytest-8.4.1.pyc
│   │   └── test_cj.py
│   ├── test_smoke.py
│   ├── ui_contracts/
│   │   ├── test_postmessage_mock.py
│   │   └── test_static_ketcher.py
│   ├── unit/
│   │   ├── __pycache__/
│   │   │   ├── test_empirical.cpython-310-pytest-8.4.1.pyc
│   │   │   └── test_health.cpython-310-pytest-8.4.1.pyc
│   │   ├── test_empirical.py
│   │   └── test_health.py
│   └── utils/
│       ├── __init__.py
│       ├── __pycache__/
│       │   ├── __init__.cpython-310.pyc
│       │   └── test_persistence.cpython-310-pytest-8.4.1.pyc
│       └── test_persistence.py
└── timer.dat                              # Performance timing data
```

## Summary Statistics

- **Total Directories**: 69
- **Total Files**: 300+ (including cached/generated files)
- **Main Python Files**: ~150
- **Test Files**: ~40
- **Documentation Files**: ~30
- **Configuration Files**: ~15
- **Generated/Cache Files**: ~100+

## Identified Duplicates and Cleanup Targets

### 🔴 Major Duplicates Requiring Cleanup

1. **API Structure Duplication**:
   - `api/` (legacy Flask structure) - **REMOVE**
   - `backend/` (new FastAPI structure) - **KEEP**

2. **Backend Duplication**:
   - `iam/backend/` (original Flask backend) - **REMOVE**
   - `backend/` (new FastAPI backend) - **KEEP**

3. **Frontend Duplication**:
   - `frontend/` (legacy test files) - **REVIEW/REMOVE**
   - `IAM_GUI/templates/` (new integrated frontend) - **KEEP**
   - `public/` (static assets) - **KEEP**

4. **Runner/Core Duplication**:
   - `iam/runners/` (original chemistry runners) - **MIGRATE THEN REMOVE**
   - `iam/core/empirical/` (empirical predictors) - **MIGRATE THEN REMOVE**
   - `backend/jobs/` (new stub implementations) - **KEEP AND EXPAND**

### 🟡 Minor Duplicates/Review Needed

1. **Test Structure**:
   - Some test files may have overlapping coverage
   - `docs/tests/` directories are empty - **REMOVE**

2. **Cache/Generated Files**:
   - `.pytest_cache/` - normal operation
   - `__pycache__/` directories - normal operation
   - `IAM_Knowledge/Results/` - 77+ result files - **REVIEW FOR RELEVANCE**

3. **Documentation**:
   - Multiple markdown files with potentially overlapping content
   - Some may be outdated

### 🟢 Clean Structure (Keep As-Is)

1. **New IAM2.0 Architecture**:
   - `backend/` - FastAPI implementation
   - `tests/` - comprehensive test suite
   - `docker/` - containerization
   - `.github/workflows/` - CI/CD

2. **Configuration**:
   - `pyproject.toml`, `requirements.txt` - dependency management
   - `docker-compose*.yml` - orchestration
   - `Makefile` - build automation

3. **Documentation**:
   - `docs/` - well-organized documentation structure

## Recommended Cleanup Actions

### Phase 1: Remove Legacy Duplicates

```bash
# Remove legacy API structure
rm -rf api/

# Remove legacy Flask backend 
rm -rf iam/

# Remove legacy frontend tests
rm -rf frontend/

# Remove empty documentation directories
rm -rf docs/tests/api/
rm -rf docs/tests/chem/
rm -rf docs/tests/ui_contracts/
```

### Phase 2: Migrate Valuable Components

1. Extract reusable logic from `iam/runners/` to `backend/jobs/`
2. Migrate empirical predictors from `iam/core/empirical/` to `backend/`
3. Consolidate any unique functionality

### Phase 3: Clean Results/Cache

1. Archive or remove old result files from `IAM_Knowledge/Results/`
2. Review and consolidate documentation

## Migration Status

- ✅ **FastAPI Backend**: Complete implementation in `backend/`
- ✅ **Test Suite**: Comprehensive coverage in `tests/`
- ✅ **CI/CD**: GitHub Actions in `.github/workflows/`
- ✅ **Docker**: Multi-service setup in `docker-compose.yml`
- ✅ **Frontend**: New template in `IAM_GUI/templates/`
- 🔄 **Legacy Cleanup**: Needs execution of cleanup actions above

---

*This document provides a complete inventory for systematic cleanup and consolidation of the IAM-2.0 codebase.*
