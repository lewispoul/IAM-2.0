# IAM: Conversion and Calculation Endpoints

**Note:** Only the first 10 results are shown. More endpoints may exist.
[View more results in GitHub](https://github.com/lewispoul/IAM/search?q=SMILES+mol+xyz+XTB+Psi4+return+json+error+temp+parse+wrapper&type=code)

---

## 1. `/molfile_to_xyz` (Convert MOL to XYZ)

- **Route:** `/molfile_to_xyz`
- **Method:** POST
- **Parameters:** `molfile` (string, MOL block from Ketcher or input)
- **Return JSON Schema:**
  ```json
  {
    "success": true,
    "xyz": "<xyz string>"
  }
  ```
- **Error Handling:** Catches missing molecule, conversion error, or network issues. Returns specific error message.
- **Files Used:**
  - `rdkit.Chem` for MOL parsing and embedding
  - Parser: RDKit

---

## 2. `/run_xtb` (Run XTB Calculation)

- **Route:** `/run_xtb`
- **Method:** POST
- **Parameters:** `xyz` (string, XYZ coordinates)
- **Return JSON Schema:**
  ```json
  {
    "success": true,
    "results": {
      "energy": -123.456,
      "homo_lumo_gap": 0.234,
      "dipole": [0.1, 0.2, 0.3]
    }
  }
  ```
- **Error Handling:** Uses try/except around XTB runner, catches subprocess failures, output parse errors.
- **Files Used:**
  - Wrapper: `xtb_wrapper.py` (calls CLI, parses output)
  - Example input: `tools/xtb/xtbopt.xyz`

---

## 3. `/run_psi4` (Run Psi4 Calculation)

- **Route:** `/run_psi4`
- **Method:** POST
- **Parameters:** `xyz` (string), `method` (string, e.g., 'HF', 'DFT'), optional calculation parameters.
- **Return JSON Schema:**
  ```json
  {
    "success": true,
    "results": {
      "energy": -234.567,
      "orbitals": {...},
      "method": "B3LYP"
    }
  }
  ```
- **Error Handling:** Handles Psi4 import errors, calculation failures, input validation.
- **Files Used:**
  - Psi4 wrapper (not shown in top 10, but implied in env validation and performance modules)
  - Parser: likely custom in backend

---

## 4. JSON Results

- All endpoints above return JSON with either `success: true` and results, or `success: false` and error.
- Schemas are nested, e.g., energy, orbital gaps, dipole moments.
- Example: `xtbout.json` shows structure for XTB results.

---

## 5. Validation, Parsing, and Performance

- **Input Validation:** Optimized via functions like `validate_xyz_format_optimized` (see `performance_utils.py`) for XYZ content.
- **Parsers:** RDKit used for MOL/SMILES ↔ XYZ, custom parser for XTB output.
- **Wrappers:** `xtb_wrapper.py` for XTB; likely similar for Psi4.
- **Performance Tracking:** Decorators track requests, errors, cache stats (see `performance_utils.py`).

---

## 6. Error Handling

- All endpoints catch and report errors in JSON response.
- Frontend notifies users with toast messages for conversion/calculation errors.

---

## 7. Example JS Usage (`iam_fixed.js`)

```javascript
fetch('/molfile_to_xyz', {
  method: 'POST',
  body: JSON.stringify({ molfile: event.data.molfile })
})
.then(response => response.json())
.then(data => {
  if (data.success) {
    // Handle successful conversion
  }
})
.catch(error => {
  showToast(`Erreur réseau: ${error.message}`, 'error');
});
```

---

## 8. Example Python (XTB Wrapper)

```python
# name=xtb_wrapper.py
def run_xtb(input_file, output_file="xtb_output/xtb_output.log"):
    try:
        # Call XTB executable
        subprocess.run(['xtb', input_file], capture_output=True, check=True)
        return parse_xtb_output(output_file)
    except Exception as e:
        return {"success": False, "error": str(e)}
```

---

**For more endpoints and details:**  
[Search the repo for more code](https://github.com/lewispoul/IAM/search?q=SMILES+mol+xyz+XTB+Psi4+return+json+error+temp+parse+wrapper&type=code)

--- 

**Incomplete results:** Only top 10 code matches shown.
