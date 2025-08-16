# IAM Legacy Endpoint Migration Plan for IAM2.0 (nox backend)

> **Note:** Only the first 10 search results are shown. For more legacy endpoints, [see full code search results in GitHub](https://github.com/lewispoul/IAM/search?q=endpoint+route+api+flask+handler+schema+request+response+migrate+adapter&type=code).

---

## Goal

Expose IAM legacy endpoints (Flask) on the new IAM2.0 nox backend.  
For each endpoint:  
- Map to new path, handler module, request/response schema  
- Specify adapter logic if existing code needs translation  
- Ensure compatibility for frontend/UI calls

---

## Migration Table

| Legacy Path        | Proposed IAM2.0 Path     | Handler Module            | Request Schema                        | Response Schema                      | Adapter/Notes                      |
|--------------------|-------------------------|--------------------------|---------------------------------------|--------------------------------------|-------------------------------------|
| `/smiles_to_xyz`   | `/api/convert/smiles-xyz` | `nox.handlers.convert`   | `{ smiles: str }` (JSON POST)         | `{ success: bool, xyz: str, error?: str }` | Minor: call RDKit, wrap error msg   |
| `/molfile_to_xyz`  | `/api/convert/mol-xyz`    | `nox.handlers.convert`   | `{ molfile: str }` (JSON POST)        | `{ success: bool, xyz: str, error?: str }` | Minor: call RDKit, wrap error msg   |
| `/run_xtb`         | `/api/calc/xtb`           | `nox.handlers.calc`      | `{ xyz: str, options?: dict }`        | `{ success: bool, results: dict, error?: str }` | Adapter: spawn xtb, parse output    |
| `/run_psi4`        | `/api/calc/psi4`          | `nox.handlers.calc`      | `{ xyz: str, method?: str, options?: dict }` | `{ success: bool, results: dict, error?: str }` | Adapter: spawn psi4, parse output   |
| `/physical_chemistry/electron_repulsion/calculate` | `/api/edu/electron-repulsion` | `nox.handlers.edu` | `{ xyz_content: str, method: str, charge?: int, multiplicity?: int }` | `{ success: bool, ...results }` | Direct, move logic to `nox.handlers.edu` |
| `/physical_chemistry/molecular_orbitals/calculate` | `/api/edu/molecular-orbitals` | `nox.handlers.edu` | `{ xyz_content?: str, smiles?: str, method?: str, ... }` | `{ success: bool, ...results }` | Direct/adapter as needed            |

---

## Example Handler Outlines

### 1. Module: `nox.handlers.convert`

```python
# /api/convert/smiles-xyz
def smiles_to_xyz(request):
    smiles = request.json['smiles']
    try:
        xyz = rdkit_smiles_to_xyz(smiles)
        return { "success": True, "xyz": xyz }
    except Exception as e:
        return { "success": False, "error": str(e) }
```

### 2. Adapter Example for XTB

```python
# /api/calc/xtb
def run_xtb(request):
    xyz = request.json['xyz']
    options = request.json.get('options', {})
    # Adapter: format input, call xtb binary, parse stdout
    results = xtb_runner.run(xyz, options)
    return { "success": results['ok'], "results": results, "error": results.get('error') }
```

---

## Request/Response Schemas

- **Request:** Always JSON, fields named as in legacy but snake_case.
- **Response:** Standardized `{ success: bool, ... }`, always includes error field if not success.

**Example:**
```json
// Request:
{ "smiles": "C1=CC=CC=C1" }
// Response:
{ "success": true, "xyz": "7\nenergy: ...\nC ..." }
```

---

## Adapter Logic

- Wrap legacy function calls with error handling
- Translate input field names if needed
- Parse stdout/stderr from CLI tools (xtb, psi4)
- Forward additional fields (e.g., charge, method) as needed

---

## Endpoint Naming Conventions

- All endpoints under `/api/`
- Use resource/action pattern: `/api/convert/...`, `/api/calc/...`, `/api/edu/...`
- Pluralize resource names only if multiple objects returned

---

## Migration Steps

1. Identify all Flask route functions and their schemas
2. Write nox handler modules (convert, calc, edu)
3. Implement adapters for CLI tools and legacy logic
4. Test frontend integration, update fetch paths in JS/UI
5. Document all new endpoint schemas and error conventions

---

**For full legacy endpoint details, refer to:**  
[GitHub code search results](https://github.com/lewispoul/IAM/search?q=endpoint+route+api+flask+handler+schema+request+response+migrate+adapter&type=code)

**Results may be incomplete due to search limits.**
