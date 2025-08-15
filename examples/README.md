# Examples: REST Client & Curl Usage

This folder contains comprehensive example requests for the IAM-2.0 backend endpoints, demonstrating both successful operations and error handling with proper error envelope responses. You can use these with the VS Code REST Client extension or via curl scripts.

## Usage

- Open any `.http` file in VS Code and click "Send Request" above the request block.
- Or, copy the request body and use curl:

```sh
curl -X POST http://localhost:8010/convert/molfile \
  -H "Content-Type: application/json" \
  -d '{"molfile": "..."}'
```

## Included Examples

### Core Endpoints
- **convert.http**: Convert molfile to SMILES (success/error cases)
- **ketcher_to_smiles.http**: Ketcher molfile to SMILES conversion
- **ketcher_to_xyz.http**: Ketcher molfile/smiles to XYZ coordinates  
- **ketcher_run.http**: Ketcher computation runs (empirical, cj methods)

### Compute Endpoints  
- **compute.http**: All compute endpoints (/compute/xtb, /compute/psi4, /compute/empirical, /compute/cj)
  - Includes success cases with proper payloads
  - Includes error cases showing envelope error responses

### Data Export
- **export.http**: Export/zip functionality with security testing
  - Path traversal protection demonstrations
  - Missing file handling
  - Empty artifact list handling

## Error Envelope Examples

All error cases in these examples demonstrate the normalized error envelope format:
```json
{
  "code": 400,
  "message": "Error description",  
  "details": {"field": "problematic_field"},
  "correlation_id": "uuid-for-tracing"
}
```

## Success Response Examples  

All success cases demonstrate the normalized success envelope format:
```json
{
  "ok": true,
  "data": {"result": "data"},
  "errors": []
}
```

## Testing Notes

- All examples use a corrected V2000 molfile for methane
- Error cases test validation, security, and edge conditions
- Examples cover all NOX integration endpoints and envelope conformance
- Use these examples to verify API behavior and error handling
