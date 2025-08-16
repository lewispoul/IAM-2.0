import pytest
from fastapi.testclient import TestClient
from backend.main import app

client = TestClient(app)

@pytest.mark.parametrize("route, payload", [
    ("/convert/molfile", {"molfile": "\n  Methane\n  OpenAI2025\n\n  5  4  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6291    0.6291    0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6291   -0.6291    0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6291   -0.6291   -0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6291    0.6291   -0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  1  3  1  0\n  1  4  1  0\n  1  5  1  0\nM  END\n"}),
    ("/ketcher/to-smiles", {"molfile": "C"}),
    ("/ketcher/to-xyz", {"smiles": "C"}),
    ("/ketcher/run", {"smiles": "C", "method": "empirical"}),
    ("/export/zip", {"artifacts": ["Results/fake.txt"]}),
])
def test_normalized_response_shape(route, payload):
    resp = client.post(route, json=payload)
    # Define expected status codes for each route
    error_routes = {
        "/convert/molfile": 400,  # invalid molfile
        "/ketcher/to-smiles": 400,  # invalid molfile
        "/export/zip": 404,  # file not found
    }
    if route in error_routes:
        assert resp.status_code == error_routes[route]
        data = resp.json()
        # Error envelope shape
        assert set(data.keys()) == {"code", "message", "details", "correlation_id"}
        assert isinstance(data["code"], int)
        assert isinstance(data["message"], str)
        assert isinstance(data["details"], dict)
        assert isinstance(data["correlation_id"], str)
    else:
        assert resp.status_code == 200
        data = resp.json()
        assert set(data.keys()) == {"ok", "data", "errors"}
        assert isinstance(data["ok"], bool)
        assert isinstance(data["data"], dict)
        assert isinstance(data["errors"], list)

@pytest.mark.parametrize("route,payload,expected_status", [
    ("/compute/xtb", {"payload": {"smiles": "C", "options": {}}}, 200),
    ("/compute/xtb", {"payload": None}, 400),
    ("/compute/psi4", {"payload": {"smiles": "C", "options": {}}}, 200),
    ("/compute/psi4", {"payload": None}, 400),
    ("/compute/empirical", {"payload": {"formula": "H2O", "method": "KJ"}}, 200),
    ("/compute/empirical", {"payload": None}, 400),
    ("/compute/cj", {"stoich": {"H2": 2, "O2": 1}, "rho0": 1.6}, 200),
    ("/compute/cj", {"stoich": {}, "rho0": None}, 422),
])
def test_compute_endpoints(route, payload, expected_status):
    resp = client.post(route, json=payload)
    assert resp.status_code == expected_status
    if expected_status == 200:
        data = resp.json()
        assert set(data.keys()) == {"ok", "data", "errors"}
        assert isinstance(data["ok"], bool)
        assert isinstance(data["data"], dict)
        assert isinstance(data["errors"], list)
    else:
        data = resp.json()
        # Accept either error envelope or FastAPI validation error
        if set(data.keys()) == {"code", "message", "details", "correlation_id"}:
            assert isinstance(data["code"], int)
            assert isinstance(data["message"], str)
            assert isinstance(data["details"], dict)
            assert isinstance(data["correlation_id"], str)
        elif "detail" in data:
            assert isinstance(data["detail"], list) or isinstance(data["detail"], str)
