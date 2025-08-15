# tests/conftest.py
import os
import pytest
from fastapi.testclient import TestClient

from iam.backend.app import app


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


def has_local_ketcher() -> bool:
    """Return True if a local Ketcher build is present."""
    return os.path.exists("public/static/ketcher/index.html")
