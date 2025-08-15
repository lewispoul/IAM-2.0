import pytest
from fastapi.testclient import TestClient
from iam.backend.app import app

client = TestClient(app)

@pytest.mark.parametrize("route", [
    "/compute/xtb",
    "/compute/psi4",
    "/compute/empirical",
    "/compute/cj",
])
def test_compute_routes(route):
    resp = client.post(route, json={"payload": {"test": True}})
    assert resp.status_code == 200
    data = resp.json()
    assert set(data.keys()) == {"ok", "data", "errors"}
    assert isinstance(data["ok"], bool)
    assert isinstance(data["data"], dict)
    assert isinstance(data["errors"], list)
