import pytest
import requests

@pytest.mark.external
def test_static_ketcher_served():
    """Static Ketcher app should be served from /static/ketcher/index.html"""
    resp = requests.get("http://localhost:5000/static/ketcher/index.html")
    assert resp.status_code == 200
    assert "<title>Ketcher" in resp.text or "Ketcher" in resp.text
