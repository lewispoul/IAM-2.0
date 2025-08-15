import os
import tempfile
import zipfile
import pytest
from fastapi.testclient import TestClient
from iam.backend.app import app

client = TestClient(app)

@pytest.fixture
def temp_results_dir(monkeypatch):
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.setenv("IAM_RESULTS_BASE", tmpdir)
        yield tmpdir

def test_export_zip_success(temp_results_dir):
    # Create a fake artifact
    artifact_path = os.path.join(temp_results_dir, "Results", "test.txt")
    os.makedirs(os.path.dirname(artifact_path), exist_ok=True)
    with open(artifact_path, "w") as f:
        f.write("test content")
    rel_path = os.path.relpath(artifact_path, temp_results_dir)
    resp = client.post("/export/zip", json={"artifacts": [rel_path]})
    assert resp.status_code == 200
    data = resp.json()
    assert data["ok"] is True
    zip_rel_path = data["data"]["zip_path"]
    zip_abs_path = os.path.join(temp_results_dir, zip_rel_path)
    assert os.path.isfile(zip_abs_path)
    # Check zip contents
    with zipfile.ZipFile(zip_abs_path, "r") as zipf:
        assert rel_path in zipf.namelist()


def test_export_zip_traversal(temp_results_dir):
    # Attempt path traversal
    resp = client.post("/export/zip", json={"artifacts": ["../etc/passwd"]})
    assert resp.status_code == 400  # Proper error status
    data = resp.json()
    assert "correlation_id" in data
    assert "Invalid path" in data["message"]
