# Examples: REST Client & Curl Usage

This folder contains example requests for the IAM-2.0 backend endpoints. You can use these with the VS Code REST Client extension or via curl scripts.

## Usage

- Open any `.http` file in VS Code and click "Send Request" above the request block.
- Or, copy the request body and use curl:

```sh
curl -X POST http://localhost:8010/convert/molfile \
  -H "Content-Type: application/json" \
  -d '{"molfile": "..."}'
```

## Included Examples
- convert.http: Convert molfile to SMILES
- ketcher_to_smiles.http: Ketcher molfile to SMILES
- ketcher_to_xyz.http: Ketcher molfile/smiles to XYZ
- ketcher_run.http: Ketcher run (empirical, cj)

All examples use a small V2000 molfile for methane inline.
