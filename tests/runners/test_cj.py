import pytest
from fastapi.testclient import TestClient
from backend.main import app

client = TestClient(app)

def test_cj_good_input():
    payload = {
        "stoich": {"H2": 2, "O2": 1},
        "rho0": 1.6
    }
    resp = client.post("/compute/cj", json=payload)
    assert resp.status_code == 200
    data = resp.json()
    assert data["ok"] is True
    assert set(data["data"].keys()) >= {"Pcj", "Tcj", "VoD", "products", "artifacts"}
    assert data["data"]["VoD"] is None
    assert isinstance(data["data"]["products"], list)
    assert isinstance(data["data"]["artifacts"], dict)

@pytest.mark.parametrize("bad_payload", [
    {},  # missing stoich, rho0
    {"stoich": {}},  # empty stoich
    {"stoich": {"": 1}, "rho0": 1.6},  # empty species key
    {"stoich": {"H2": -1}, "rho0": 1.6},  # non-positive fraction
    {"stoich": {"H2": 2, "O2": 1}, "rho0": 0},  # rho0 <= 0
    {"stoich": {"H2": 2, "O2": 1}, "rho0": "bad"},  # rho0 not a number
])
def test_cj_bad_inputs(bad_payload):
    resp = client.post("/compute/cj", json=bad_payload)
    # Invalid payloads should return 422 Unprocessable Entity
    assert resp.status_code == 422
    data = resp.json()
    # FastAPI validation errors return 'detail' key, not custom schema
    assert "detail" in data
