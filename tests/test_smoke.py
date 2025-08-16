"""
Smoke tests for IAM2.0 API endpoints.
Basic functionality tests to ensure core endpoints are working.
"""
import pytest


def test_health_endpoint(fastapi_test_client):
    """Test health check endpoint returns 200."""
    response = fastapi_test_client.get("/health")
    assert response.status_code == 200
    assert response.json() == {"status": "ok"}


def test_smiles_to_xyz_endpoint(fastapi_test_client):
    """/api/convert/smiles-to-xyz returns 200 with valid structure."""
    payload = {"smiles": "C1=CC=CC=C1"}
    response = fastapi_test_client.post("/api/convert/smiles-to-xyz", json=payload)
    
    assert response.status_code == 200
    data = response.json()
    assert data["success"] is True
    assert "xyz" in data
    assert "atoms" in data
    assert isinstance(data["atoms"], int)


def test_xtb_calc_endpoint(fastapi_test_client):
    """/api/calc/xtb returns JSON with mock energy."""
    payload = {"xyz": "1\nenergy\nC 0.0 0.0 0.0"}
    response = fastapi_test_client.post("/api/calc/xtb", json=payload)
    
    assert response.status_code == 200
    data = response.json()
    assert data["success"] is True
    assert "results" in data
    assert isinstance(data["results"], dict)
    assert "energy" in data["results"]
    assert isinstance(data["results"]["energy"], (int, float))


@pytest.mark.integration
def test_integration_server_health(test_server_url):
    """Integration test: ensure test server is responding."""
    import requests
    try:
        response = requests.get(f"{test_server_url}/health", timeout=5)
        assert response.status_code == 200
    except requests.ConnectionError:
        pytest.skip("Integration test server not available")


@pytest.mark.chemistry
def test_chemistry_pipeline_smoke(sample_molecules, mock_rdkit):
    """Smoke test for chemistry pipeline without real dependencies."""
    benzene = sample_molecules['benzene_smiles']
    methane_mol = sample_molecules['methane_mol']
    
    # Test that sample data is valid
    assert benzene == "C1=CC=CC=C1"
    assert "Methane" in methane_mol
    assert "C   0  0" in methane_mol
    
    # Test mock RDKit responds correctly
    mol = mock_rdkit.Chem.MolFromSmiles(benzene)
    assert mol is not None
    assert mock_rdkit.Chem.MolFromSmiles.called
    assert mol.GetNumAtoms() == 6
