import pytest
from fastapi.testclient import TestClient
from iam.backend.app import app

client = TestClient(app)

def test_ketcher_to_smiles_good():
    molfile = """
  MJ201100                      

  6  5  0  0  0  0            999 V2000
    1.2080   -0.6980    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2080   -0.6980    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2080   -2.1220    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -2.8200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2080   -2.1220    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
M  END
"""
    resp = client.post("/ketcher/to-smiles", json={"molfile": molfile})
    assert resp.status_code == 200
    data = resp.json()
    assert data["ok"] is True
    assert "smiles" in data["data"]

def test_ketcher_to_xyz_placeholder():
    resp = client.post("/ketcher/to-xyz", json={"smiles": "C"})
    assert resp.status_code == 200
    data = resp.json()
    assert data["ok"] is True
    assert "xyz" in data["data"]

def test_ketcher_run_empirical():
    resp = client.post("/ketcher/run", json={"smiles": "C", "method": "empirical"})
    assert resp.status_code == 200
    data = resp.json()
    assert data["ok"] is True
    assert "Pcj" in data["data"] or "Tcj" in data["data"]

def test_ketcher_run_bad_method():
    resp = client.post("/ketcher/run", json={"smiles": "C", "method": "unknown"})
    assert resp.status_code == 200
    data = resp.json()
    assert data["ok"] is False
    assert data["errors"]
