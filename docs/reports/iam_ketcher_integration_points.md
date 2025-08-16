# Ketcher Integration Points in IAM

_Results may be incomplete due to search limitations. Update with a full file audit as needed._

---

## 1. Ketcher References: iframe, postMessage, molecule set/get, export

| File                   | Line(s)       | Message Types / Actions              | Payload Shape                  | Receiving Code              | Backend Routes Called      |
|------------------------|---------------|--------------------------------------|-------------------------------|-----------------------------|----------------------------|
| `templates/ketcher.html`<br>`index.html` | (various)         | Ketcher iframe is embedded           | `<iframe src="...ketcher.html">` | JS parent in main app      | N/A (static load)          |
| `static/ketcher_bridge.js` | (various)         | `postMessage` to/from Ketcher iframe  | `{type: 'getMol', format: 'mol'}`, `{type: 'setMol', data: ...}`, `{type: 'exportSmiles'}` | Ketcher iframe JS, parent JS | `/load_mol`, `/save_mol`   |
| `static/inputs.js`     | (various)         | Sets/gets molecule from Ketcher       | SMILES string, MOL block       | Ketcher bridge, input forms | `/predict`, `/run_xtb`, `/run_psi4` |
| `static/outputs.js`    | (various)         | Exports molecule to output/download   | SMILES/MOL block               | Output pane JS              | N/A                        |
| `static/app.js`        | (various)         | Tab logic, triggers Ketcher actions   | N/A                            | Main controller             | Various                    |

---

## 2. Message Types and Payloads

- **`postMessage` (to Ketcher iframe):**
  - `{type: "setMol", data: <mol block|smiles>, format: "mol"|"smiles"}`
  - `{type: "exportSmiles"}`
- **`postMessage` (from Ketcher iframe):**
  - `{type: "molChanged", mol: <mol block>, smiles: <smiles string>}`
  - `{type: "exportSmiles", smiles: <smiles string>}`

**Payload shapes:**
- MOL block: multi-line string
- SMILES: string
- Object: `{type, data, format}`

---

## 3. Receiving Code

- **In main app JS** (`static/app.js`, `static/ketcher_bridge.js`):
  - Listens for `message` events from Ketcher iframe.
  - Calls backend with molecule payload when user triggers save/load/export.
- **Input/output panes**:
  - Use the Ketcher bridge to get/set molecule data for prediction and export.

---

## 4. Backend Routes Called

- **Molecule load/save:** `/load_mol`, `/save_mol`
- **Prediction endpoints:** `/predict`, `/run_xtb`, `/run_psi4`
- **May also interact with:** `/performance` (if molecule used for performance prediction)

---

## Example Snippet

```js
// Send molecule to Ketcher
ketcherFrame.contentWindow.postMessage({type: "setMol", data: molBlock, format: "mol"}, "*");

// Listen for molecule changes
window.addEventListener("message", function(event) {
  if (event.data.type === "molChanged") {
    // Handle molecule change
    const molData = event.data.mol;
    const smilesData = event.data.smiles;
    // Process molecule data...
  }
});
```

---

## Summary Table

| Integration Point   | JS Controller        | HTML Template        | Backend Route    |
|---------------------|---------------------|---------------------|------------------|
| Ketcher iframe      | ketcher_bridge.js   | ketcher.html        | N/A              |
| Molecule input      | inputs.js           | inputs.html         | /predict, /run_* |
| Molecule export     | outputs.js          | outputs.html        | N/A              |
| Tab navigation      | app.js              | index.html          | Various          |

---

**Results may be incomplete due to search limits.**
