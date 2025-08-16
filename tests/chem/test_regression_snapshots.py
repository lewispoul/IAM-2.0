import pytest
import requests

@pytest.mark.parametrize("endpoint, payload, snapshot_name", [
    ("smiles_to_xyz", {"smiles": "C1=CC=CC=C1"}, "xyz_benzene"),
    ("molfile_to_xyz", {"molfile": """
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
"""}, "xyz_methane"),
    ("run_xtb", {"xyz": "1\nenergy\nC 0.0 0.0 0.0"}, "xtb_stub_methane"),
    ("run_psi4", {"xyz": "1\nenergy\nC 0.0 0.0 0.0"}, "psi4_stub_methane"),
])
def test_regression_snapshots(snapshot, endpoint, payload, snapshot_name):
    """Regression snapshot for endpoint JSON output"""
    resp = requests.post(f"http://localhost:5000/{endpoint}", json=payload)
    assert resp.status_code == 200
    snapshot.assert_match(resp.json(), name=snapshot_name)
