# Contributing to IAM 2.0

## Getting Started
- Python 3.10+
- Create env with conda or venv
- `pip install -r requirements.txt` (or `pip install -e .[dev]` if using pyproject)
- `pre-commit install` (if configured)

## Branching & PRs
- Use feature branches: `feat/<area>`, `fix/<bug>`, `docs/<topic>`
- Reference spec sections by ID in issues/PRs (e.g., `[Spec-ID: ENG.XTB]`)
- Include tests and update docs when behavior changes

## Coding Standards
- PEP8, type hints where practical
- Lint with ruff, type-check with mypy
- Return the unified results schema and standard error envelope from all endpoints

## Adding a New Module (Quick Steps)
1. Open a feature request issue referencing the Whitepaper section.
2. Create module under `engines/`, `predictors/`, or `pipelines/` as appropriate.
3. Implement per Acceptance Criteria in the spec.
4. Add unit/integration tests.
5. Update API examples and docs if endpoints change.
6. Open a PR linking the issue and section IDs.
