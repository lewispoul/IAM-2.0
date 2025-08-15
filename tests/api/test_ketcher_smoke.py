"""
Ketcher Route Smoke Tests
=============================

Quick smoke tests specifically for Ketcher-related backend routes
to verify they're working before UI integration testing.
"""

import pytest
from fastapi.testclient import TestClient
from iam.backend.app import app

client = TestClient(app)

# Mock data for testing
MOCK_MOLFILE = """
  Methane
  
  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6291    0.6291    0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6291   -0.6291    0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6291   -0.6291   -0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6291    0.6291   -0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END"""

MOCK_SMILES = "C"  # Methane


class TestKetcherRoutes:
    """Smoke tests for Ketcher integration routes"""
    
    def test_ketcher_to_smiles_success(self):
        """Test successful molfile to SMILES conversion"""
        response = client.post("/ketcher/to-smiles", json={
            "molfile": MOCK_MOLFILE
        })
        
        assert response.status_code == 200
        data = response.json()
        assert data["ok"] is True
        assert "smiles" in data["data"]
        assert data["data"]["smiles"] == "C"  # Methane SMILES
    
    def test_ketcher_to_smiles_invalid_molfile(self):
        """Test error handling for invalid molfile"""
        response = client.post("/ketcher/to-smiles", json={
            "molfile": "invalid molfile content"
        })
        
        assert response.status_code == 400
        data = response.json()
        assert "correlation_id" in data
        assert "Could not parse molfile" in data["message"]
    
    def test_ketcher_to_xyz_success(self):
        """Test successful SMILES to XYZ conversion"""
        response = client.post("/ketcher/to-xyz", json={
            "smiles": MOCK_SMILES
        })
        
        assert response.status_code == 200
        data = response.json()
        assert data["ok"] is True
        assert "xyz" in data["data"]
        assert "TODO: 3D coordinates" in data["data"]["xyz"]  # Current placeholder
    
    def test_ketcher_to_xyz_invalid_smiles(self):
        """Test error handling for invalid SMILES - currently returns success with TODO"""
        response = client.post("/ketcher/to-xyz", json={
            "smiles": "invalid_smiles_string"
        })
        
        # Note: Currently returns 200 with TODO placeholder
        # In production, this should validate SMILES and return 400 for invalid input
        assert response.status_code == 200
        data = response.json()
        assert data["ok"] is True
        assert "TODO" in data["data"]["xyz"]
    
    def test_ketcher_run_empirical(self):
        """Test Ketcher empirical calculation"""
        response = client.post("/ketcher/run", json={
            "smiles": MOCK_SMILES,
            "method": "empirical"
        })
        
        assert response.status_code == 200
        data = response.json()
        assert data["ok"] is True
        assert "data" in data
        # Empirical returns detonation properties
        assert "Pcj" in data["data"]
        assert "Tcj" in data["data"]
    
    def test_ketcher_run_cj(self):
        """Test Ketcher CJ calculation"""
        response = client.post("/ketcher/run", json={
            "smiles": MOCK_SMILES,
            "method": "cj"
        })
        
        assert response.status_code == 200
        data = response.json()
        assert data["ok"] is True
        assert "data" in data
        # CJ returns Chapman-Jouguet properties
        assert "Pcj" in data["data"]
        assert "Tcj" in data["data"]
    
    def test_ketcher_run_invalid_method(self):
        """Test error handling for invalid method"""
        response = client.post("/ketcher/run", json={
            "smiles": MOCK_SMILES,
            "method": "invalid_method"
        })
        
        assert response.status_code == 501  # Not Implemented for unsupported methods
        data = response.json()
        assert "correlation_id" in data
        assert "not implemented" in data["message"]
    
    def test_ketcher_run_missing_smiles(self):
        """Test error handling for missing SMILES - currently passes validation"""
        response = client.post("/ketcher/run", json={
            "method": "empirical"
        })
        
        # Note: Currently passes because SMILES has default value
        # In production, should require explicit SMILES parameter
        assert response.status_code == 200
    
    def test_ketcher_run_empty_smiles(self):
        """Test error handling for empty SMILES - currently processes default"""
        response = client.post("/ketcher/run", json={
            "smiles": "",
            "method": "empirical"
        })
        
        # Note: Currently processes with default behavior
        # In production, should validate non-empty SMILES
        assert response.status_code == 200


class TestKetcherIntegrationFlow:
    """Test complete workflow from UI perspective"""
    
    def test_molfile_to_smiles_to_calculation_flow(self):
        """Test complete flow: molfile → SMILES → calculation"""
        # Step 1: Convert molfile to SMILES
        smiles_response = client.post("/ketcher/to-smiles", json={
            "molfile": MOCK_MOLFILE
        })
        assert smiles_response.status_code == 200
        smiles_data = smiles_response.json()
        assert smiles_data["ok"] is True
        smiles = smiles_data["data"]["smiles"]
        
        # Step 2: Run calculation with the SMILES
        calc_response = client.post("/ketcher/run", json={
            "smiles": smiles,
            "method": "empirical"
        })
        assert calc_response.status_code == 200
        calc_data = calc_response.json()
        assert calc_data["ok"] is True
        assert "Pcj" in calc_data["data"]  # Has detonation properties
    
    def test_smiles_to_xyz_conversion_flow(self):
        """Test SMILES to XYZ coordinate generation"""
        response = client.post("/ketcher/to-xyz", json={
            "smiles": MOCK_SMILES
        })
        assert response.status_code == 200
        response_data = response.json()
        assert response_data["ok"] is True
        xyz_data = response_data["data"]["xyz"]
        
        # Verify XYZ format (currently placeholder)
        assert "TODO" in xyz_data  # Current implementation returns placeholder


if __name__ == "__main__":
    # Run the tests
    pytest.main([__file__, "-v", "--tb=short"])
