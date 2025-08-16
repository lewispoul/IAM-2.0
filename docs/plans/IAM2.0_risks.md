# IAM2.0 Migration Risk Register

| Risk | Impact | Likelihood | Mitigation | Owner | Status |
|------|--------|------------|------------|-------|--------|
| Ketcher postMessage/iframe fails in browser edge cases | High | Medium | Use robust JS handler, cross-browser QA | lewispoul | Ongoing |
| RDKit missing/fails in container | High | Low | Add validation, fallback error message | lewispoul | Mitigated |
| /smiles_to_xyz or /molfile_to_xyz returns invalid/empty XYZ | High | Medium | Add regression tests, error logging | lewispoul | Ongoing |
| XTB/Psi4 backend stub not returning correct JSON | Medium | Medium | Contract test, stub fallback, error schema | lewispoul | Planned |
| UI and backend route mismatch (fetch fails) | High | Medium | Update UI fetch patterns, adapter layer | lewispoul | Ongoing |
| Performance prediction tab delayed (CJ/Cantera integration) | Medium | Medium | Prioritize core predictors first, fallback to empirical | lewispoul | Acceptable |
| Snapshot regression tests not covering edge cases | Medium | Medium | Expand tests, add molecule diversity | lewispoul | Ongoing |
| Scientific library install/compatibility (Psi4, Cantera) | High | Low | Scripted env checks, containerize | lewispoul | Planned |
| API contract drift (schema mismatch) | Medium | Medium | Centralize schema, document endpoints | lewispoul | Planned |
| Unclear ownership of future modules | Low | Medium | Assign lead for each module, track PRs | lewispoul | Open |

*Focus mitigation on UI + molecule conversion flow, then scientific jobs, then prediction. Review and update risks monthly or per milestone.*
