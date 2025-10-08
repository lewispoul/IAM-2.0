#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="$(pwd)"
REPORT_ROOT="$PROJECT_ROOT/reports/Nox_Audit_$(date +%Y%m%d_%H%M%S)"
PY_MOD_DIRS=("nox" "iam" "pinox")
mkdir -p "$REPORT_ROOT"/{inventory,graphs,quality,security,tests,docs,raw}
echo "▶ Reports -> $REPORT_ROOT"
python3 -m pip install --user --upgrade -r requirements-dev.txt || true
cloc . --exclude-dir=.git,venv,__pycache__,node_modules --json > "$REPORT_ROOT/inventory/cloc.json" || true
git log --oneline --decorate --graph --all > "$REPORT_ROOT/inventory/git_history.txt" || true
git remote -v > "$REPORT_ROOT/inventory/git_remotes.txt" || true
{ echo "# Tree (depth 3)"; echo; find . -maxdepth 3 -type d -printf "%p\n" | sort; } > "$REPORT_ROOT/inventory/tree.txt"
python3 -m pip freeze > "$REPORT_ROOT/inventory/pip_freeze.txt" || true
pipdeptree --warn silence > "$REPORT_ROOT/inventory/pipdeptree.txt" || true
for m in "${PY_MOD_DIRS[@]}"; do
  if [ -d "$m" ]; then
    pydeps "$m" --max-bacon=4 --noshow --output "$REPORT_ROOT/graphs/${m}_imports.svg" || true
  fi
done
rg -n --no-heading -e "@app\.route|@.*\.get\(|@.*\.post\(|APIRouter\(" -g "**/*.py" > "$REPORT_ROOT/docs/routes_hits.txt" || true
rg -n --no-heading -e "FastAPI\(|Flask\(" -g "**/*.py" >> "$REPORT_ROOT/docs/routes_hits.txt" || true
rg -n --no-heading -e "openapi|swagger" -g "**/*.{yml,yaml,json,py}" > "$REPORT_ROOT/docs/openapi_hits.txt" || true
for f in $(rg -l -e "openapi|swagger" -g "**/*.y*ml" -g "**/*.json"); do
  openapi-spec-validator "$f" && echo "VALID: $f" >> "$REPORT_ROOT/docs/openapi_validation.txt" || echo "INVALID: $f" >> "$REPORT_ROOT/docs/openapi_validation.txt"
done
ruff check . > "$REPORT_ROOT/quality/ruff.txt" || true
mypy --ignore-missing-imports . > "$REPORT_ROOT/quality/mypy.txt" || true
radon cc . -s -a > "$REPORT_ROOT/quality/radon_cc.txt" || true
radon mi . -s > "$REPORT_ROOT/quality/radon_mi.txt" || true
bandit -r . -q -f txt -o "$REPORT_ROOT/security/bandit.txt" || true
safety check --full-report > "$REPORT_ROOT/security/safety.txt" || true
pytest -q --maxfail=1 --disable-warnings --cov=. --cov-report=term-missing || true
pytest --cov=. --cov-report=html:"$REPORT_ROOT/tests/htmlcov" || true
pytest --cov=. --cov-report=xml:"$REPORT_ROOT/tests/coverage.xml" || true
rg -n --no-heading -e "TODO|FIXME|HACK" -g "**/*.*" > "$REPORT_ROOT/raw/todos.txt" || true
cat > "$REPORT_ROOT/SUMMARY.md" <<'MD'
# Nox – Audit d’Architecture (auto-généré)
Dossier : inventory/ (cloc, git), graphs/ (imports), docs/ (routes, openapi), quality/ (ruff, mypy, radon), security/ (bandit, safety), tests/ (coverage), raw/todos.txt.
## Priorités
1) Corriger ruff/mypy/bandit critiques. 2) Monter MI radon. 3) Ajouter tests sur modules centraux. 4) Stabiliser OpenAPI. 5) MAJ deps vulnérables. 6) Brancher CI.
MD
echo "✅ Report: $REPORT_ROOT/SUMMARY.md"