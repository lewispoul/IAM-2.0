# tests/conftest.py
"""
Pytest configuration and fixtures for IAM2.0 testing.
"""
import os
import socket
import pytest
import http.client
from unittest.mock import Mock
from fastapi.testclient import TestClient

from backend.main import app


def _port_open(host: str, port: int) -> bool:
    """Check if a port is open with a quick timeout."""
    try:
        with socket.create_connection((host, port), timeout=0.3):
            return True
    except OSError:
        return False


def _health_ok() -> bool:
    """Check if external service on :5000 is healthy."""
    if not _port_open("127.0.0.1", 5000):
        return False
    try:
        conn = http.client.HTTPConnection("127.0.0.1", 5000, timeout=0.3)
        conn.request("GET", "/healthz")
        resp = conn.getresponse()
        return resp.status == 200
    except Exception:
        return False


def pytest_collection_modifyitems(config, items):
    """Mark external integration tests and skip if service not available."""
    require_external = (os.getenv("IAM_E2E") == "1") and _health_ok()
    
    # Known external integration test files
    external_patterns = [
        "test_molfile_to_xyz.py",
        "test_psi4_stub.py", 
        "test_smiles_to_xyz.py",
        "test_xtb_stub.py"
    ]
    
    for item in items:
        nid = item.nodeid
        if any(pattern in nid for pattern in external_patterns):
            item.add_marker(pytest.mark.external)
            if not require_external:
                item.add_marker(pytest.mark.skip(reason="External service not running (set IAM_E2E=1 and start service on :5000)"))


@pytest.fixture(scope="session")
def client():
    """FastAPI TestClient for API tests."""
    return TestClient(app)


@pytest.fixture()
def tmp_results_dir(monkeypatch, tmp_path):
    """Isolate persistence by pointing IAM_RESULTS_BASE to a temp dir."""
    base = tmp_path / "IAM_Knowledge"
    base.mkdir(parents=True, exist_ok=True)
    monkeypatch.setenv("IAM_RESULTS_BASE", str(base))
    yield base


@pytest.fixture
def flask_test_client():
    """Flask test client for legacy endpoint compatibility testing."""
    from flask import Flask
    app_flask = Flask(__name__)
    app_flask.config['TESTING'] = True
    
    @app_flask.route('/health')
    def health():
        return {'status': 'ok'}
    
    with app_flask.test_client() as client:
        yield client


@pytest.fixture  
def fastapi_test_client():
    """Enhanced FastAPI test client for IAM2.0 endpoints."""
    from fastapi import FastAPI
    
    test_app = FastAPI()
    
    @test_app.get("/health")
    def health():
        return {"status": "ok"}
    
    @test_app.post("/api/convert/smiles-to-xyz")
    def smiles_to_xyz(data: dict):
        # Mock response for testing
        return {
            "success": True,
            "xyz": "1\ntest\nC 0.0 0.0 0.0",
            "atoms": 1
        }
    
    @test_app.post("/api/calc/xtb")
    def xtb_calc(data: dict):
        # Mock response for testing
        return {
            "success": True,
            "results": {
                "energy": -100.0,
                "status": "completed"
            }
        }
    
    client = TestClient(test_app)
    yield client


@pytest.fixture
def mock_rdkit():
    """Mock RDKit functionality for testing without dependencies."""
    mock_mol = Mock()
    mock_mol.GetNumAtoms.return_value = 6
    
    mock_rdkit = Mock()
    mock_rdkit.Chem.MolFromSmiles.return_value = mock_mol
    mock_rdkit.Chem.MolToMolBlock.return_value = "mock mol block"
    
    return mock_rdkit


@pytest.fixture
def sample_molecules():
    """Sample molecular data for testing."""
    return {
        'benzene_smiles': 'C1=CC=CC=C1',
        'methane_mol': """
  Methane
  RDKit

  5 4  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0
    0.0000    0.0000    1.0890 H   0  0
    1.0267    0.0000   -0.3630 H   0  0
   -0.5133   -0.8892   -0.3630 H   0  0
   -0.5133    0.8892   -0.3630 H   0  0
  1 2 1  0  0  0  0
  1 3 1  0  0  0  0
  1 4 1  0  0  0  0
  1 5 1  0  0  0  0
M  END
""",
        'methane_xyz': """1
energy
C 0.0 0.0 0.0"""
    }


@pytest.fixture
def test_server_url():
    """Base URL for integration tests."""
    return "http://localhost:5000"


def has_local_ketcher() -> bool:
    """Return True if a local Ketcher build is present."""
    return os.path.exists("public/static/ketcher/index.html")


# Pytest markers
def pytest_configure(config):
    """Configure custom pytest markers."""
    config.addinivalue_line(
        "markers", "integration: mark test as integration test requiring running server"
    )
    config.addinivalue_line(
        "markers", "slow: mark test as slow running"
    )
    config.addinivalue_line(
        "markers", "chemistry: mark test as requiring chemistry dependencies"
    )
