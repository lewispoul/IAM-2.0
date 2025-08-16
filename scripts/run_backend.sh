#!/bin/bash
set -e

# IAM2.0 Backend Runner Script
# Starts the FastAPI backend with appropriate configuration

echo "🚀 Starting IAM2.0 Backend Server"
echo "================================="

# Check if we're in the right directory
if [[ ! -f "backend/api/convert.py" ]]; then
    echo "❌ Error: Must run from IAM-2.0 root directory"
    echo "Current directory: $(pwd)"
    exit 1
fi

# Set environment variables
export PYTHONPATH="$(pwd):$PYTHONPATH"
export IAM_RESULTS_BASE="${IAM_RESULTS_BASE:-$(pwd)/IAM_Knowledge}"
export LOG_LEVEL="${LOG_LEVEL:-INFO}"
export ENVIRONMENT="${ENVIRONMENT:-development}"

# Create required directories
mkdir -p "$IAM_RESULTS_BASE"
mkdir -p logs

# Check if conda environment is active
if [[ -n "$CONDA_DEFAULT_ENV" ]]; then
    echo "✅ Conda environment: $CONDA_DEFAULT_ENV"
else
    echo "⚠️  No conda environment detected"
    echo "Consider running: conda activate iam2"
fi

# Check dependencies
echo "📦 Checking dependencies..."
python -c "import fastapi; print(f'✅ FastAPI {fastapi.__version__}')" || {
    echo "❌ FastAPI not found. Run: pip install -r requirements.txt"
    exit 1
}

python -c "import uvicorn; print(f'✅ Uvicorn available')" || {
    echo "❌ Uvicorn not found. Run: pip install uvicorn"
    exit 1
}

# Check optional chemistry dependencies
echo "🧪 Checking chemistry dependencies..."
python -c "
try:
    import rdkit
    print('✅ RDKit available')
except ImportError:
    print('⚠️  RDKit not available - using stub mode')
"

python -c "
try:
    import psi4
    print('✅ Psi4 available')
except ImportError:
    print('⚠️  Psi4 not available - using stub mode')
"

# Show configuration
echo ""
echo "🔧 Configuration:"
echo "  PYTHONPATH: $PYTHONPATH"
echo "  IAM_RESULTS_BASE: $IAM_RESULTS_BASE"
echo "  LOG_LEVEL: $LOG_LEVEL"
echo "  ENVIRONMENT: $ENVIRONMENT"
echo ""

# Start the server
echo "🌟 Starting server on http://localhost:8000"
echo "📚 API docs will be available at http://localhost:8000/docs"
echo "🔄 Press Ctrl+C to stop"
echo ""

exec uvicorn iam.backend.app:app \
    --host 0.0.0.0 \
    --port 8000 \
    --reload \
    --reload-dir backend \
    --reload-dir iam \
    --log-level "${LOG_LEVEL,,}" \
    --access-log
