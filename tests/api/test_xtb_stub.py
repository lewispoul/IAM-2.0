import requests
import json

def test_xtb_stub_json():
    """Stubbed /run_xtb returns parseable JSON with expected fields"""
    xyz = "1\nenergy\nC 0.0 0.0 0.0"
    resp = requests.post("http://localhost:5000/run_xtb", json={"xyz": xyz})
    assert resp.status_code == 200
    data = resp.json()
    assert data.get("success")
    results = data.get("results")
    assert isinstance(results, dict)
    assert "energy" in results
