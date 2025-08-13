# Contributing to IAM 2.0

## Dev Setup
- Python 3.10+
- `pip install -e .[dev]`
- `pre-commit install`

## PR Guidelines
- Small, focused PRs tied to an issue.
- Reference `docs/IAM_Master_Technical_Report_FINAL.md` §section.
- Include tests (unit/integration) and docs updates.

## Contracts
- Use the unified results schema & error envelopes defined in the FINAL doc.
- Keep API examples valid against `io/schema.py` (when added).
