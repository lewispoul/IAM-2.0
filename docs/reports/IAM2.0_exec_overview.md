# IAM2.0 Migration Executive Overview

## Scope

IAM2.0 is a major refactor of the IAM computational chemistry platform, integrating:
- A modular UI (Ketcher sketcher, 3Dmol.js viewer, job runner, results & performance tabs)
- Backend migration to NOX API (FastAPI), with robust REST endpoints
- Scientific computing integration (RDKit, XTB, Psi4, ML predictors, Cantera/CJ)

## Key Features Migrated

- **Ketcher Integration:** Static sketcher at `/static/ketcher/index.html`, secure postMessage round-trip, robust conversion of drawn molecules to MOL/SMILES.
- **3D Visualization:** Real-time rendering using 3Dmol.js; seamless conversion of SMILES/MOL to XYZ via `/smiles_to_xyz` and `/molfile_to_xyz`.
- **Job Execution:** Modular status panel for job queue/progress; endpoints for XTB and Psi4 calculation stubs, with parseable JSON output.
- **Results & Performance:** Tabbed interface for summary, output, log, and performance prediction (KJ, ML, CJ).
- **Scientific Integration:** Support for RDKit, XTB, Psi4, Cantera libraries; prediction endpoints for detonation and energetic materials.

## Validation & Testing

- Pytest validation plan covers UI contracts, API endpoints, chemistry outputs, and regression snapshots.
- Scaffolding for tests under `tests/api`, `tests/chem`, and `tests/ui_contracts`.

## Migration & Adapter Strategy

- Legacy Flask routes mapped to new NOX FastAPI endpoints, with standardized JSON request/response schema and adapters as needed.
- UI fetch patterns updated for new backend routes; each module isolated for maintainability.

## Risks & Priorities

- Immediate focus: UI, Ketcher, SMILES/MOL→XYZ flow (lowest effort, highest impact).
- Next: XTB/Psi4 job runner and output parsing.
- Then: Performance prediction (KJ, ML, CJ/Cantera).
- See attached TODO and risk register for actionable items and mitigation.

## Next Steps

- Complete endpoint migration and UI module wiring.
- Execute validation plan; unblock core molecule input→3D→job→results flow.
- Progressively integrate scientific backends and advanced predictors.

*IAM2.0 is well-positioned for a robust relaunch, with a modern, modular architecture and validated migration plan.*
