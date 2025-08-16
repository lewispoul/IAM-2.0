# IAM2.0 Migration TODO List

| Priority | Task | Owner | Difficulty | Est. Effort | Notes |
|----------|------|-------|------------|-------------|-------|
| 1 | Finalize static Ketcher panel & postMessage round-trip | lewispoul | Medium | 4h | Unblocks molecule input; UI + backend integration |
| 2 | Wire `/smiles_to_xyz` and `/molfile_to_xyz` endpoints to NOX backend | lewispoul | Low | 2h | Adapter for RDKit, standardized schema |
| 3 | Refactor 3Dmol.js viewer for new XYZ fetch patterns | lewispoul | Low | 2h | UI update to support new backend |
| 4 | Implement job runner status module (polling, error handling) | lewispoul | Medium | 3h | Needed for calculation feedback |
| 5 | Stub XTB and Psi4 endpoints (parseable JSON, error handling) | lewispoul | Medium | 4h | Unblocks quantum calculation flow |
| 6 | Tabbed results panel (summary, output, log) | lewispoul | Low | 2h | UI wiring, fetch per tab |
| 7 | Implement regression snapshot pytest tests | lewispoul | Medium | 3h | Contract for outputs, avoids regressions |
| 8 | Add performance prediction tab (KJ, ML, CJ/Cantera) | lewispoul | High | 6h | Last major scientific component |
| 9 | Integrate CJ/Cantera backend and prediction endpoint | lewispoul | High | 8h | Scientific/engineering work |
| 10 | Document API contract and update developer guide | lewispoul | Low | 2h | Finalization/maintenance |

*Prioritization: unblock UI and molecule conversion flow first, then quantum jobs, then prediction tabs.*
