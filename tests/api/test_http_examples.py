"""
Test that the example HTTP requests work correctly
"""
import pytest
from fastapi.testclient import TestClient
from backend.main import app

client = TestClient(app)

class TestHTTPExamples:
    """Test cases matching the examples in examples/requests/ folder."""
    
    def test_convert_molfile_success_example(self):
        """Test the success case from convert.http"""
        molfile = "\n  Methane\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6291    0.6291    0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6291   -0.6291    0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6291   -0.6291   -0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6291    0.6291   -0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1  3  1  0  0  0  0\n  1  4  1  0  0  0  0\n  1  5  1  0  0  0  0\nM  END\n"
        
        resp = client.post("/convert/molfile", json={"molfile": molfile})
        assert resp.status_code == 200
        data = resp.json()
        assert data["ok"] is True
        assert "smiles" in data["data"]
        assert data["data"]["smiles"] == "C"  # Methane
    
    def test_convert_molfile_error_example(self):
        """Test the error case from convert.http"""
        resp = client.post("/convert/molfile", json={"molfile": "invalid molfile structure"})
        assert resp.status_code == 400
        data = resp.json()
        assert "correlation_id" in data
        assert "Could not parse molfile" in data["message"]
    
    def test_ketcher_run_success_example(self):
        """Test the empirical success case from ketcher_run.http"""
        resp = client.post("/ketcher/run", json={
            "smiles": "C",
            "method": "empirical"
        })
        assert resp.status_code == 200
        data = resp.json()
        assert data["ok"] is True
        assert "Pcj" in data["data"]  # Empirical results
    
    def test_ketcher_run_error_example(self):
        """Test the invalid method error case from ketcher_run.http"""
        resp = client.post("/ketcher/run", json={
            "smiles": "C",
            "method": "invalid_method"
        })
        assert resp.status_code == 501
        data = resp.json()
        assert "correlation_id" in data
        assert "not implemented" in data["message"]
    
    def test_compute_xtb_success_example(self):
        """Test the XTB success case from compute.http"""
        resp = client.post("/compute/xtb", json={
            "payload": {
                "smiles": "C",
                "options": {
                    "method": "GFN2-xTB",
                    "opt": True
                }
            }
        })
        assert resp.status_code == 200
        data = resp.json()
        assert data["ok"] is True
        assert "data" in data
    
    def test_compute_cj_success_example(self):
        """Test the CJ success case from compute.http"""
        resp = client.post("/compute/cj", json={
            "stoich": {"H2": 2, "O2": 1},
            "rho0": 1.6
        })
        assert resp.status_code == 200
        data = resp.json()
        assert data["ok"] is True
        assert "Pcj" in data["data"]  # CJ results
    
    def test_compute_cj_error_example(self):
        """Test the CJ error case from compute.http"""
        resp = client.post("/compute/cj", json={
            "stoich": "invalid_string",
            "rho0": 1.6
        })
        assert resp.status_code == 422
        # This should trigger FastAPI validation error
    
    def test_export_zip_error_example(self):
        """Test the path traversal error case from export.http"""
        resp = client.post("/export/zip", json={
            "artifacts": ["../../../etc/passwd"]
        })
        assert resp.status_code == 400
        data = resp.json()
        assert "correlation_id" in data
        assert "Invalid path" in data["message"]
        
    def test_ketcher_to_xyz_success_example(self):
        """Test the XYZ success case from ketcher_to_xyz.http"""
        resp = client.post("/ketcher/to-xyz", json={"smiles": "C"})
        assert resp.status_code == 200
        data = resp.json()
        assert data["ok"] is True
        assert "xyz" in data["data"]
        
    def test_ketcher_to_smiles_success_example(self):
        """Test the SMILES success case from ketcher_to_smiles.http""" 
        molfile = "\n  Methane\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6291    0.6291    0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6291   -0.6291    0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6291   -0.6291   -0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6291    0.6291   -0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1  3  1  0  0  0  0\n  1  4  1  0  0  0  0\n  1  5  1  0  0  0  0\nM  END\n"
        
        resp = client.post("/ketcher/to-smiles", json={"molfile": molfile})
        assert resp.status_code == 200
        data = resp.json()
        assert data["ok"] is True
        assert "smiles" in data["data"]
        assert data["data"]["smiles"] == "C"
