# IAM2.0 Migration Pytest Validation Plan

This document outlines the test plan and scaffolding for validating IAM2.0 key features after migration.

---

## 1. UI Static Assets & Contracts

- **Test:** Static Ketcher is served at `/static/ketcher/index.html`
- **Test:** postMessage round-trip for molecule data (mocked)

### Pytest Files:
- `tests/ui_contracts/test_static_ketcher.py`
- `tests/ui_contracts/test_postmessage_mock.py`

---

## 2. API Endpoints

- **Test:** `/smiles_to_xyz` endpoint returns valid XYZ coordinates for simple molecules (e.g., benzene, methane)
- **Test:** `/molfile_to_xyz` endpoint returns valid XYZ for MOL input
- **Test:** XTB and Psi4 stub endpoints return parseable JSON (structure, energy fields)

### Pytest Files:
- `tests/api/test_smiles_to_xyz.py`
- `tests/api/test_molfile_to_xyz.py`
- `tests/api/test_xtb_stub.py`
- `tests/api/test_psi4_stub.py`

---

## 3. Chemistry Outputs & Regression

- **Test:** Regression snapshots for JSON outputs from calculation endpoints (including results, logs)
- **Test:** Ensure XYZ and result fields remain consistent for standard inputs

### Pytest Files:
- `tests/chem/test_regression_snapshots.py`

---

## 4. Pytest Scaffolding

- All tests use `pytest` and may use `pytest-snapshot` for regression.
- Use `requests` or `httpx` for API calls.
- UI contracts use `requests` for static assets, and `unittest.mock` to simulate postMessage.

---

## 5. Example Test Matrix

| Test Name                        | Endpoint/Feature                    | Input                        | Expected Output                    |
|----------------------------------|-------------------------------------|------------------------------|------------------------------------|
| test_static_ketcher_served       | /static/ketcher/index.html          | GET                          | 200 OK, HTML with Ketcher app      |
| test_postmessage_roundtrip       | JS postMessage (mocked)             | MOL block                    | Receives same MOL in response      |
| test_smiles_to_xyz_benzene       | /smiles_to_xyz                      | SMILES: "C1=CC=CC=C1"        | XYZ block with 6 C, 6 H            |
| test_molfile_to_xyz_methane      | /molfile_to_xyz                     | MOL: methane                 | XYZ block with 1 C, 4 H            |
| test_xtb_stub_json               | /run_xtb (stubbed)                  | XYZ: methane                 | JSON with 'energy', 'dipole'       |
| test_psi4_stub_json              | /run_psi4 (stubbed)                 | XYZ: methane                 | JSON with 'energy', 'orbitals'     |
| test_regression_xtb_output       | /run_xtb                            | XYZ: methane                 | Output matches snapshot            |
| test_regression_psi4_output      | /run_psi4                           | XYZ: benzene                 | Output matches snapshot            |

---

## 6. Folder Structure

```
tests/
├── api/
│   ├── test_smiles_to_xyz.py
│   ├── test_molfile_to_xyz.py
│   ├── test_xtb_stub.py
│   ├── test_psi4_stub.py
├── chem/
│   ├── test_regression_snapshots.py
├── ui_contracts/
│   ├── test_static_ketcher.py
│   ├── test_postmessage_mock.py
```

---

## 7. Next Steps

- Scaffold test files with boilerplate and docstrings.
- Implement fixtures for API client and regression snapshots.
- Ensure coverage for all migrated endpoints.

---
