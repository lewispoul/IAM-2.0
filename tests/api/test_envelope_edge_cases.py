import pytest
from fastapi.testclient import TestClient
from backend.main import app

client = TestClient(app)

@pytest.mark.parametrize("route,payload,expected_status", [
    # Large payload tests
    ("/convert/molfile", {"molfile": "x" * 3_000_000}, 400),  # >2MB
    ("/ketcher/to-smiles", {"molfile": "x" * 3_000_000}, 400),  # >2MB

    # Empty/null payload tests
    ("/convert/molfile", {}, 422),  # Missing required field
    ("/ketcher/to-smiles", {}, 422),  # Missing required field
    ("/convert/molfile", {"molfile": ""}, 400),  # Empty string molfile

    # Invalid method tests
    ("/ketcher/run", {"smiles": "C", "method": "invalid"}, 501),

    # Invalid stoichiometry for CJ
    ("/compute/cj", {"stoich": "invalid", "rho0": 1.6}, 422),
    ("/compute/cj", {"stoich": {}, "rho0": "invalid"}, 422),
])
def test_edge_case_error_envelopes(route, payload, expected_status):
    resp = client.post(route, json=payload)
    assert resp.status_code == expected_status
    data = resp.json()
    
    # For 422 errors (FastAPI validation), accept either error envelope or FastAPI format
    if expected_status == 422:
        if "detail" in data:
            # FastAPI validation error format
            assert isinstance(data["detail"], (list, str))
        else:
            # Our error envelope format
            assert set(data.keys()) == {"code", "message", "details", "correlation_id"}
    else:
        # Our error envelope format for other errors
        assert set(data.keys()) == {"code", "message", "details", "correlation_id"}
        assert isinstance(data["code"], int)
        assert isinstance(data["message"], str)
        assert isinstance(data["details"], dict)
        assert isinstance(data["correlation_id"], str)
        assert data["code"] == expected_status

@pytest.mark.parametrize("route,valid_payload", [
    ("/convert/molfile", {"molfile": "\n  Methane\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6291    0.6291    0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6291   -0.6291    0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6291   -0.6291   -0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6291    0.6291   -0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1  3  1  0  0  0  0\n  1  4  1  0  0  0  0\n  1  5  1  0  0  0  0\nM  END\n"}),
    ("/ketcher/to-xyz", {"smiles": "C"}),
    ("/ketcher/run", {"smiles": "C", "method": "empirical"}),
    ("/compute/xtb", {"payload": {"smiles": "C", "options": {}}}),
    ("/compute/psi4", {"payload": {"smiles": "C", "options": {}}}),
    ("/compute/empirical", {"payload": {"formula": "H2O", "method": "KJ"}}),
    ("/compute/cj", {"stoich": {"H2": 2, "O2": 1}, "rho0": 1.6}),
])
def test_success_envelope_shape(route, valid_payload):
    resp = client.post(route, json=valid_payload)
    assert resp.status_code == 200
    
    data = resp.json()
    assert set(data.keys()) == {"ok", "data", "errors"}
    assert isinstance(data["ok"], bool)
    assert isinstance(data["data"], dict)
    assert isinstance(data["errors"], list)
    assert data["ok"] == True
    assert data["errors"] == []
