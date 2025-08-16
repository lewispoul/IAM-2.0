# tests/api/test_convert_and_ketcher.py
from pathlib import Path

MOCK_MOLFILE = """\

     RDKit          3D

  5  4  0  0  0  0  0  0  0  0999 V2000
   -0.0221    0.0032    0.0165 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6665    0.8884   -0.1014 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3778   -0.8578   -0.5883 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0964   -0.3151    1.0638 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.9699    0.2812   -0.3906 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
M  END
"""
MOCK_SMILES = "C"  # methane


def assert_ok(resp):
    assert resp.status_code == 200
    body = resp.json()
    assert set(body.keys()) == {"ok", "data", "errors"}
    assert body["ok"] is True
    assert isinstance(body["data"], dict)
    assert body["errors"] == []
    return body


def assert_fail(resp):
    assert resp.status_code == 200
    body = resp.json()
    assert set(body.keys()) == {"ok", "data", "errors"}
    assert body["ok"] is False
    assert isinstance(body["errors"], list) and body["errors"]
    return body


def test_to_smiles_success(client, tmp_results_dir):
    r = client.post("/ketcher/to-smiles", json={"molfile": MOCK_MOLFILE})
    body = assert_ok(r)
    assert body["data"]["smiles"] == "C"
    assert "formula" in body["data"]
    assert "inchi" in body["data"]


def test_to_smiles_invalid_molfile(client, tmp_results_dir):
    r = client.post("/ketcher/to-smiles", json={"molfile": "garbage"})
    # Invalid molfile format should return 400 (matches original Flask behavior)
    assert r.status_code == 400
    body = r.json()
    assert "correlation_id" in body
    assert "Could not parse molfile" in body["message"]


def test_to_xyz_success_from_smiles(client, tmp_results_dir):
    r = client.post("/ketcher/to-xyz", json={"smiles": MOCK_SMILES})
    body = assert_ok(r)
    assert "xyz" in body["data"]
    assert isinstance(body["data"]["xyz"], str)


def test_to_xyz_invalid_smiles(client, tmp_results_dir):
    r = client.post("/ketcher/to-xyz", json={"smiles": "not_a_smiles"})
    body = assert_fail(r)
    assert any("smiles" in e.lower() for e in body["errors"]) or any(
        "parse" in e.lower() for e in body["errors"]
    )


def test_to_xyz_validation_missing(client):
    r = client.post("/ketcher/to-xyz", json={})
    assert r.status_code == 422
