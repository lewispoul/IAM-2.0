import requests

def test_molfile_to_xyz_methane():
    """/molfile_to_xyz returns valid XYZ for methane"""
    methane_mol = """
  Methane
  RDKit

  5 4  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0
    0.0000    0.0000    1.0890 H   0  0
    1.0267    0.0000   -0.3630 H   0  0
   -0.5133   -0.8892   -0.3630 H   0  0
   -0.5133    0.8892   -0.3630 H   0  0
  1 2 1  0  0  0  0
  1 3 1  0  0  0  0
  1 4 1  0  0  0  0
  1 5 1  0  0  0  0
M  END
"""
    resp = requests.post("http://localhost:5000/molfile_to_xyz", json={"molfile": methane_mol})
    assert resp.status_code == 200
    data = resp.json()
    assert data.get("success")
    xyz = data.get("xyz")
    assert xyz and "C" in xyz and "H" in xyz
    assert xyz.count("\n") >= 4
