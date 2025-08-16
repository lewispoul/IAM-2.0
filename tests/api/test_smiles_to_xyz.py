import requests

def test_smiles_to_xyz_benzene():
    """/smiles_to_xyz returns valid XYZ for benzene"""
    smiles = "C1=CC=CC=C1"
    resp = requests.post("http://localhost:5000/smiles_to_xyz", json={"smiles": smiles})
    assert resp.status_code == 200
    data = resp.json()
    assert data.get("success")
    xyz = data.get("xyz")
    assert xyz and "C" in xyz and "H" in xyz
    assert xyz.count("\n") > 6  # Should be >6 lines for benzene
