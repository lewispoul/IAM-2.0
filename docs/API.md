# API.md

> Placeholder for IAM 2.0 API specification.

> **Note:** The `iam2` conda environment is set to auto-activate for all new terminal sessions. All API development and testing should be performed in this environment for consistency.

## /compute/cj
- **POST /compute/cj**
- Request JSON example:
  ```json
  {
    "stoich": {"H2": 2, "O2": 1},
    "rho0": 1.6,
    "dhf": -50.0,
    "metals": {"Fe": 0.1},
    "species_db": "default"
  }
  ```
- Normalized response shape:
  ```json
  {
    "ok": true,
    "data": {
      "Pcj": 1.23,
      "Tcj": 3456.0,
      "VoD": null,
      "products": [],
      "artifacts": {"species_db": "default", "input_sum": 3.0, "used_metals": true}
    },
    "errors": []
  }
  ```
- On error:
  ```json
  {
    "ok": false,
    "data": {},
    "errors": ["error message"]
  }
  ```
- Note: VoD is None in stub; full CJ will be implemented later.
