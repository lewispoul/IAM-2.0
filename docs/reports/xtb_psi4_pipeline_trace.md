# IAM End-to-End Flow: "Run" Button to Results Display (XTB & Psi4)

**Note:** Only the first 10 search results are shown. For more details, [see full code search results in GitHub](https://github.com/lewispoul/IAM/search?q=run+xtb+psi4+fetch+response+json+output+parse+display+render&type=code).

---

## A. XTB Pipeline

### 1. **Front-End Call Site**
- **File:** `IAM_GUI/static/inputs.js`, `webui/templates/index.html`
- **Function:** User clicks "Run XTB" button/form.
- **Action:** JS (or HTML form) calls backend route via `fetch` or form POST.
  ```html
  <button onclick="runXTB()">Run XTB</button>
  <script>
  function runXTB() {
    const xyz = getXYZFromKetcher();
    fetch('/run_xtb', {
      method: 'POST',
      body: JSON.stringify({xyz: xyz})
    });
  }
  </script>
  ```

### 2. **Backend Route**
- **File:** `backend.py` or equivalent Flask app (not in top 10, inferred from system)
- **Route:** `/run_xtb` (POST)
- **Function:** Receives XYZ data, invokes XTB wrapper.

### 3. **Wrapper Invocation**
- **File:** `xtb_wrapper.py`
- **Function:** `run_xtb(input_file, output_file="xtb_output/xtb_output.log")`
- **Action:** Calls XTB executable with XYZ file, writes output log.

### 4. **Output Parsing**
- **File:** `xtb_wrapper.py`
- **Function:** `parse_xtb_output(output_file)`
- **Action:** Parses output log for key results (energy, HOMO-LUMO gap, dipole).

### 5. **Response JSON Fields**
- **Returned from backend:**
  ```json
  {
    "success": true,
    "energy": -123.456,
    "homo_lumo_gap": 0.234,
    "dipole": [0.1, 0.2, 0.3],
    "optimized_xyz": "7\n\nC 0.0 0.0 0.0\n..."
  }
  ```
  *(See also: `xtbout.json` for full schema)*

### 6. **UI Rendering Site**
- **File:** `webui/templates/index.html`, `IAM_GUI/static/outputs.js`
- **Function:** Renders results in list/table.
  ```html
  <div id="xtb-results">
    <table>
      <tr><td>Energy:</td><td id="energy"></td></tr>
      <tr><td>HOMO-LUMO Gap:</td><td id="gap"></td></tr>
    </table>
  </div>
  ```

---

## B. Psi4 Pipeline

### 1. **Front-End Call Site**
- **File:** `IAM_GUI/static/inputs.js`, (similar to XTB)
- **Function:** User clicks "Run Psi4" button/form.
- **Action:** JS sends `fetch('/run_psi4', ...)` with XYZ and method.

### 2. **Backend Route**
- **File:** `backend.py` or equivalent Flask app
- **Route:** `/run_psi4` (POST)
- **Function:** Receives XYZ and method, invokes Psi4 wrapper.

### 3. **Wrapper Invocation**
- **File:** (Psi4 wrapper not in top 10; typical name: `psi4_wrapper.py`)
- **Function:** `run_psi4(input_file, ...)`
- **Action:** Calls Psi4 executable, writes output files/logs.

### 4. **Output Parsing**
- **File:** (inferred, `psi4_wrapper.py`)
- **Function:** `parse_psi4_output(output_file)`
- **Action:** Parses output for quantum chemistry results (energy, orbitals).

### 5. **Response JSON Fields**
- **Returned from backend:**
  ```json
  {
    "success": true,
    "energy": -234.567,
    "method": "B3LYP",
    "orbitals": {
      "homo": -0.123,
      "lumo": 0.045
    },
    "cube_files": ["homo.cube", "lumo.cube"]
  }
  ```

### 6. **UI Rendering Site**
- **File:** `IAM_GUI/static/outputs.js`, `webui/templates/index.html`
- **Function:** Renders results (energy, orbitals) in web UI.

---

## File & Function Map

| Stage              | XTB                                  | Psi4                                 |
|--------------------|--------------------------------------|--------------------------------------|
| Frontend Call      | `inputs.js`, `index.html`            | `inputs.js`, `index.html`            |
| Backend Route      | `/run_xtb` (POST)                    | `/run_psi4` (POST)                   |
| Wrapper            | `xtb_wrapper.py:run_xtb`              | `psi4_wrapper.py:run_psi4` (inferred)|
| Output Parse       | `xtb_wrapper.py:parse_xtb_output`     | `psi4_wrapper.py:parse_psi4_output`  |
| Response Schema    | Energy, HOMO-LUMO gap, dipole       | Energy, orbitals, method, cube files |
| UI Rendering       | `outputs.js`, results table         | `outputs.js`, orbital visualization  |

**Note:** Results may be incomplete due to search limits.
