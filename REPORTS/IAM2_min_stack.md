## Step 0-3 Summary (2025-08-15)
- Added repo guardrails (.copilot-instructions.md), VS Code settings, and chem-env.yaml for conda/mamba environment.
- Updated README with install instructions and Makefile to default to python3 and run FastAPI app via uvicorn.
- Created backend app wiring (FastAPI) and common schema helpers (ok/fail) in iam/backend/app.py and iam/backend/schemas/common.py.

## Part B / Step 1 (2025-08-15)
- Added FastAPI app wiring and included routers for xtb, psi4, empirical, cj.
- Implemented shared schema helpers (ok/fail) and ComputeRequest model.
- Created compute route stubs for /compute/xtb, /compute/psi4, /compute/empirical, /compute/cj, all returning normalized schema.
- Updated Makefile to run uvicorn for FastAPI app.
- Added pytest smoke tests for all compute endpoints.
- Documented endpoints and schema in docs/API.md.

## Part B / Step 2 (2025-08-15)
- Added empirical predictor stubs for Kamlet–Jacobs and Keshavarz methods with deterministic outputs.
- Implemented predict_empirical interface to select method.
- Wired /compute/empirical route to use empirical predictor.
- Added unit tests for empirical predictor and route.
