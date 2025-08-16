from fastapi.testclient import TestClient
import pytest
from backend.main import app

def test_health():
    client = TestClient(app)
    r = client.get("/api/v1/health")
    assert r.status_code == 200
    assert r.json().get("status") == "ok"
