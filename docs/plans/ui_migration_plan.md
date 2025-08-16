# IAM2.0 UI Modularization & Route Wiring Plan

## Overview

Reassemble the IAM UI as modular, maintainable components in IAM2.0. Each module is wired to backend routes (IAM2.0 API or NOX API) and uses robust fetch patterns.

---

## 1. **Ketcher Panel**

- **Purpose:** Molecular sketcher; users draw or edit molecules, export MOL/SMILES.
- **Route wiring:** 
  - Static asset: `/static/ketcher/index.html`
  - PostMessage integration with parent UI; interface fetches MOL/SMILES output from the iframe.
  - Conversion endpoint: `POST /api/convert/mol-xyz` (NOX or IAM2.0)
- **Fetch pattern:** 
  - On "Export" or "Load to 3D" button:
    ```js
    const molfile = getMolfileFromKetcher(); // via postMessage
    const resp = await fetch('/api/convert/mol-xyz', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ molfile })
    });
    const { success, xyz } = await resp.json();
    if (success) showXYZIn3DViewer(xyz);
    ```

---

## 2. **3D Viewer**

- **Purpose:** Renders 3D molecular structure (XYZ or MOL); uses 3Dmol.js.
- **Route wiring:**
  - Receives XYZ from Ketcher panel, file upload, or result endpoints.
  - Optionally: `GET /api/molecule/:job_id/xyz`
- **Fetch pattern:**
  ```js
  // After successful conversion or job completion
  fetch(`/api/molecule/${job_id}/xyz`)
    .then(r => r.json())
    .then(data => ViewerManager.renderMolecule(data.xyz));
  ```

---

## 3. **Job Runner Status**

- **Purpose:** Tracks status of submitted jobs (XTB, Psi4, ML, CJ).
- **Route wiring:**
  - `GET /api/job/:job_id/status`
  - `GET /api/jobs/active`
- **Fetch pattern:**
  ```js
  setInterval(async () => {
    const resp = await fetch(`/api/job/${job_id}/status`);
    const { status, progress, error } = await resp.json();
    updateJobStatusUI(status, progress, error);
  }, 2000);
  ```

---

## 4. **Results Tabs (Summary, Output, Log)**

- **Purpose:** Tabbed view for output: Summary, raw Output, backend Log.
- **Route wiring:**
  - `GET /api/job/:job_id/summary`
  - `GET /api/job/:job_id/output`
  - `GET /api/job/:job_id/log`
- **Fetch pattern:**
  ```js
  async function loadTab(tabType) {
    const resp = await fetch(`/api/job/${job_id}/${tabType}`);
    const data = await resp.json();
    renderTabContent(tabType, data);
  }
  // Tabs: ["summary", "output", "log"]
  ```

---

## 5. **Performance Prediction Tab (KJ, ML, CJ/Cantera)**

- **Purpose:** Predict detonation performance (velocity, pressure, enthalpy) with ML, empirical, and CJ/Cantera models.
- **Route wiring:**
  - `POST /api/predict/performance` with `{ job_id, model: "ML"|"KJ"|"CJ" }`
  - `GET /api/job/:job_id/prediction?model=ML|KJ|CJ`
- **Fetch pattern:**
  ```js
  const resp = await fetch('/api/predict/performance', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ job_id, model: selectedModel })
  });
  const { prediction, graphs } = await resp.json();
  showPrediction(prediction, graphs);
  ```

---

## 6. **Component Structure**

```jsx
<App>
  <KetcherPanel />
  <ThreeDViewer />
  <JobRunnerStatus />
  <ResultsTabs>
    <Tab label="Summary" />
    <Tab label="Output" />
    <Tab label="Log" />
    <Tab label="Performance" />
  </ResultsTabs>
</App>
```

---

## 7. **Route Mapping Table**

| UI Module                 | Route(s)                                 | Data Flow              |
|---------------------------|------------------------------------------|------------------------|
| KetcherPanel              | /static/ketcher/index.html, /api/convert/mol-xyz | MOL → XYZ             |
| ThreeDViewer              | /api/molecule/:job_id/xyz                | XYZ → 3Dmol.js         |
| JobRunnerStatus           | /api/job/:job_id/status, /api/jobs/active| Job status             |
| ResultsTabs (Summary...)  | /api/job/:job_id/{summary,output,log}    | Output JSON/text/log   |
| PerformancePredictionTab  | /api/predict/performance, /api/job/:job_id/prediction | Metrics, graphs        |

---

## 8. **Fetch & State Patterns**

- Use async/await and robust error handling throughout.
- Prefer POST for conversions, GET for result retrieval.
- Use polling for job status.
- Use Redux/Context/Pinia for global state where applicable.

---

## 9. **Implementation Notes**

- Each module is isolated for testability and maintenance.
- All API endpoints accept/return JSON with `{ success, ... }`.
- Ketcher integration uses postMessage for secure communication and round-trip confirmation.
- 3Dmol.js is loaded as a static asset (`/static/3Dmol-min.js`).
- Performance tab includes ML and CJ/Cantera results, supported via backend libraries.

---

## 10. **Next Steps**

- Scaffold each module as a standalone JS/React/Vue component.
- Implement API compatibility layer to allow legacy and new routes.
- Add regression tests for round-trip and output formats.
- Document API contract for each module.

---

**For more details, see:**
- [IAM_GUI/static/script.js](https://github.com/lewispoul/IAM/blob/main/IAM_GUI/static/script.js)
- [IAM2.0 API FastAPI routes](https://github.com/lewispoul/nox)
- [IAM Front-End Refactor Summary](https://github.com/lewispoul/IAM/blob/main/FRONT_END_REFACTOR_SUMMARY.md)

**Results may be incomplete due to search limits.**
