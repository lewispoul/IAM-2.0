import pytest
from fastapi.testclient import TestClient
from iam.backend.app import app

client = TestClient(app)

def test_convert_molfile_good():
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
    resp = client.post("/convert/molfile", json={"molfile": molfile})
    assert resp.status_code == 200
    data = resp.json()
    assert data["ok"] is True
    assert "smiles" in data["data"]
    assert "formula" in data["data"]

def test_convert_molfile_bad():
    resp = client.post("/convert/molfile", json={"molfile": ""})
    assert resp.status_code == 400  # Proper error status
    data = resp.json()
    assert "correlation_id" in data
    assert "molfile must be a non-empty string" in data["message"]
