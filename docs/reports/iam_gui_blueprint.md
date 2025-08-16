# IAM GUI Blueprint: Component Map

_Results may be incomplete due to search limitations. Update with full file audit as needed._

---

| Component           | File(s)                       | Responsibilities                                                     | Depends On                | Backend Route(s) Called                |
|---------------------|------------------------------|---------------------------------------------------------------------|---------------------------|----------------------------------------|
| Main Layout/Template| `templates/index.html`        | Overall app structure, tab panel, loads panes for viewer/Ketcher     | JS controllers, CSS       | `/`, `/predict`, `/run_xtb`, `/run_psi4`, `/performance` |
| Ketcher Editor Pane | `templates/ketcher.html`, `static/ketcher_bridge.js` | Embeds Ketcher (iframe), handles molecule drawing/editing           | Ketcher JS, postMessage   | `/load_mol`, `/save_mol`               |
| 3D Viewer Pane      | `templates/viewer.html`, `static/viewer.js`, `static/3dmol.js` | Renders molecule in 3D, display HOMO/LUMO, handles visualization    | 3Dmol.js, molecule state  | `/get_xyz`, `/get_cube`, `/get_orbitals`|
| Input Section       | `templates/inputs.html`, `static/inputs.js` | Handles SMILES/XYZ file input, parameter forms                       | HTML forms, fetch API     | `/predict`, `/run_xtb`, `/run_psi4`    |
| Output Section      | `templates/outputs.html`, `static/outputs.js` | Shows calculated results, tables/plots, downloads                    | Results state, Plotly.js  | `/predict`, `/performance`             |
| Performance Prediction | `templates/performance.html`, `static/performance.js` | Displays combustion/detonation predictions, CJ calculations         | JS logic, input/output    | `/performance`, `/run_cantera`         |
| Navigation/Tabs     | `static/app.js`, `templates/index.html` | Tab switching, layout management, state syncing                      | HTML structure, CSS       | Various (routes per tab)               |
| CSS/Styling         | `static/style.css`            | App layout, responsive design, tab/pane appearance                   | HTML structure            | N/A                                    |

---

## Component Details

**Main Layout/Template**
- Responsible for page structure, loading all panes via tabs.
- Depends on JS for controller logic and CSS for layout.
- Calls backend for main page load and all prediction endpoints.

**Ketcher Editor Pane**
- Loads Ketcher via iframe, bridges molecule data via postMessage.
- JS controller syncs editor state and sends molecular data to backend.

**3D Viewer Pane**
- Uses 3Dmol.js to render molecules, visualize orbitals (HOMO/LUMO), cube files.
- JS controller fetches geometry/cube data from backend.

**Input Section**
- Handles user inputs (SMILES, XYZ, parameters).
- JS fetches prediction and quantum chemistry endpoints.

**Output Section**
- Displays property prediction results, tables, plots.
- JS updates UI from backend responses.

**Performance Prediction**
- Pane for combustion/detonation predictions, CJ computations.
- JS fetches performance endpoints, renders results.

**Navigation/Tabs**
- Controls which pane is visible, manages state between tabs.
- Coordinates data flow between components.

---

## Tab Structure

```
Main App
├── Ketcher Tab (molecule editor)
├── 3D Viewer Tab (visualization)
├── Input Tab (parameters, file upload)
├── Output Tab (results display)
└── Performance Tab (predictions)
```

---

## Data Flow

1. **User draws molecule** in Ketcher → postMessage to parent
2. **Parent receives molecule** → sends to Input section
3. **Input section** → calls backend `/predict` or `/run_*`
4. **Backend responds** → Output section displays results
5. **Results include XYZ** → 3D Viewer renders molecule
6. **Performance prediction** → calls `/performance` endpoint

---

**Results may be incomplete due to search limits.**
