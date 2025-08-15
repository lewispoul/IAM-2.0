# API.md

IAM 2.0 API Specification

## Error Envelope (Standard)

All error responses use a normalized envelope:

```json
{
  "code": 400,
  "message": "Invalid input",
  "details": {},
  "correlation_id": "uuid"
}
```

- `code`: HTTP status code (e.g., 400, 404, 429, 500)
- `message`: Human-readable error message
- `details`: Additional error context (e.g., field, engine, quota)
- `correlation_id`: Unique request ID for UI/debugging

All endpoints return this envelope for errors, with appropriate status codes.

## Endpoints

### /convert/molfile
- **POST /convert/molfile**
- Request: `{ "molfile": "..." }`
- Response: Normalized schema

### /ketcher/to-smiles
- **POST /ketcher/to-smiles**
- Request: `{ "molfile": "..." }`
- Response: Normalized schema

### /ketcher/to-xyz
- **POST /ketcher/to-xyz**
- Request: `{ "molfile": "..." }` or `{ "smiles": "..." }`
- Response: Normalized schema

### /ketcher/run
- **POST /ketcher/run**
- Request: `{ "molfile": "...", "smiles": "...", "method": "xtb|psi4|empirical|cj", "options": { ... } }`
- Response: Normalized schema

### /compute/xtb
- **POST /compute/xtb**
- Request: `{ "payload": { "smiles": "...", "options": { ... } } }`
- Response: Normalized schema

### /compute/psi4
- **POST /compute/psi4**
- Request: `{ "payload": { "smiles": "...", "options": { ... } } }`
- Response: Normalized schema

### /compute/empirical
- **POST /compute/empirical**
- Request: `{ "payload": { "formula": "...", "method": "KJ|Keshavarz" } }`
- Response: Normalized schema

### /compute/cj
- **POST /compute/cj**
- Request: `{ "stoich": { ... }, "rho0": ..., "dhf": ..., "metals": { ... }, "species_db": "..." }`
- Response: Normalized schema

### /export/zip
- **POST /export/zip**
- Request: `{ "artifacts": ["path1", "path2", ...] }`
- Response: Normalized schema

## Request Schemas

- **ConvertMolfileRequest**: `{ "molfile": "..." }`
- **KetcherToXYZRequest**: `{ "molfile": "..." }` or `{ "smiles": "..." }`
- **KetcherRunRequest**: `{ "molfile": "...", "smiles": "...", "method": "xtb|psi4|empirical|cj", "options": { ... } }`
- **ComputeXTBRequest**: `{ "payload": { "smiles": "...", "options": { ... } } }`
- **ComputePsi4Request**: `{ "payload": { "smiles": "...", "options": { ... } } }`
- **ComputeEmpiricalRequest**: `{ "payload": { "formula": "...", "method": "KJ|Keshavarz" } }`
- **ComputeCJRequest**: `{ "stoich": { ... }, "rho0": ..., "dhf": ..., "metals": { ... }, "species_db": "..." }`
- **ExportZipRequest**: `{ "artifacts": ["path1", "path2", ...] }`

## Error Handling Examples

### Invalid Input
```json
{
  "code": 400,
  "message": "Invalid SMILES",
  "details": {"field": "input_data"},
  "correlation_id": "uuid"
}
```

### File Not Found
```json
{
  "code": 404,
  "message": "File not found: Results/fake.txt",
  "details": {"field": "Results/fake.txt"},
  "correlation_id": "uuid"
}
```

### Quota Exceeded
```json
{
  "code": 429,
  "message": "Quota exceeded",
  "details": {"limit": "MAX_JOBS"},
  "correlation_id": "uuid"
}
```
