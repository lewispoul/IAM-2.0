#!/bin/bash
set -e

# IAM2.0 Test Runner Script
# Runs the complete test suite with proper environment setup

echo "🧪 IAM2.0 Test Suite Runner"
echo "============================"

# Check if we're in the right directory
if [[ ! -f "tests/conftest.py" ]]; then
    echo "❌ Error: Must run from IAM-2.0 root directory"
    echo "Current directory: $(pwd)"
    exit 1
fi

# Set environment variables
export PYTHONPATH="$(pwd):$PYTHONPATH"
export IAM_RESULTS_BASE="${IAM_RESULTS_BASE:-$(pwd)/IAM_Knowledge_test}"
export ENVIRONMENT="testing"

# Create test directories
mkdir -p "$IAM_RESULTS_BASE"
mkdir -p logs

# Clean previous test artifacts
echo "🧹 Cleaning previous test artifacts..."
find . -name "*.pyc" -delete 2>/dev/null || true
find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
rm -rf .pytest_cache htmlcov .coverage 2>/dev/null || true

# Check dependencies
echo "📦 Checking test dependencies..."
python -c "import pytest; print(f'✅ Pytest {pytest.__version__}')" || {
    echo "❌ Pytest not found. Run: pip install pytest"
    exit 1
}

# Parse command line arguments
TEST_TYPE="all"
COVERAGE=false
VERBOSE=false
MARKERS=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --smoke)
            TEST_TYPE="smoke"
            shift
            ;;
        --integration)
            TEST_TYPE="integration"
            shift
            ;;
        --api)
            TEST_TYPE="api"
            shift
            ;;
        --chem)
            TEST_TYPE="chem"
            shift
            ;;
        --ui)
            TEST_TYPE="ui"
            shift
            ;;
        --coverage)
            COVERAGE=true
            shift
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        --markers)
            MARKERS="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --smoke        Run smoke tests only"
            echo "  --integration  Run integration tests only"
            echo "  --api          Run API tests only"
            echo "  --chem         Run chemistry tests only"
            echo "  --ui           Run UI contract tests only"
            echo "  --coverage     Generate coverage report"
            echo "  -v, --verbose  Verbose output"
            echo "  --markers M    Run tests with specific markers"
            echo "  -h, --help     Show this help"
            exit 0
            ;;
        *)
            echo "❌ Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Build pytest command
PYTEST_ARGS="--tb=short"

if [[ "$VERBOSE" == "true" ]]; then
    PYTEST_ARGS="$PYTEST_ARGS -v"
fi

if [[ "$COVERAGE" == "true" ]]; then
    PYTEST_ARGS="$PYTEST_ARGS --cov=backend --cov-report=html --cov-report=term-missing"
fi

# Add test selection based on type
case $TEST_TYPE in
    "smoke")
        PYTEST_ARGS="$PYTEST_ARGS tests/test_smoke.py"
        echo "🔥 Running smoke tests..."
        ;;
    "integration")
        PYTEST_ARGS="$PYTEST_ARGS -m integration"
        echo "🔗 Running integration tests..."
        ;;
    "api")
        PYTEST_ARGS="$PYTEST_ARGS tests/api/"
        echo "🌐 Running API tests..."
        ;;
    "chem")
        PYTEST_ARGS="$PYTEST_ARGS tests/chem/"
        echo "🧪 Running chemistry tests..."
        ;;
    "ui")
        PYTEST_ARGS="$PYTEST_ARGS tests/ui_contracts/"
        echo "🖥️  Running UI contract tests..."
        ;;
    "all")
        PYTEST_ARGS="$PYTEST_ARGS tests/"
        echo "🎯 Running all tests..."
        ;;
esac

if [[ -n "$MARKERS" ]]; then
    PYTEST_ARGS="$PYTEST_ARGS -m '$MARKERS'"
    echo "🏷️  Using markers: $MARKERS"
fi

# Show configuration
echo ""
echo "🔧 Test Configuration:"
echo "  PYTHONPATH: $PYTHONPATH"
echo "  IAM_RESULTS_BASE: $IAM_RESULTS_BASE"
echo "  Test type: $TEST_TYPE"
echo "  Coverage: $COVERAGE"
echo "  Verbose: $VERBOSE"
echo ""

# Check for conda environment
if [[ -n "$CONDA_DEFAULT_ENV" ]]; then
    echo "✅ Conda environment: $CONDA_DEFAULT_ENV"
else
    echo "⚠️  No conda environment detected"
fi

# Run the tests
echo "🏃 Executing: pytest $PYTEST_ARGS"
echo ""

pytest $PYTEST_ARGS

# Show results summary
RESULT=$?
echo ""
if [[ $RESULT -eq 0 ]]; then
    echo "✅ All tests passed!"
    if [[ "$COVERAGE" == "true" ]]; then
        echo "📊 Coverage report generated in htmlcov/"
    fi
else
    echo "❌ Some tests failed (exit code: $RESULT)"
fi

# Clean up test artifacts if successful
if [[ $RESULT -eq 0 && "$COVERAGE" == "false" ]]; then
    echo "🧹 Cleaning up test artifacts..."
    rm -rf "$IAM_RESULTS_BASE" 2>/dev/null || true
fi

exit $RESULT
