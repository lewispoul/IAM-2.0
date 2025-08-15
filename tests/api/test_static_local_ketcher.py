# tests/api/test_static_local_ketcher.py
import pytest
from tests.conftest import has_local_ketcher


@pytest.mark.skipif(not has_local_ketcher(), reason="local Ketcher not fetched")
def test_local_ketcher_index_served(client):
    r = client.get("/static/ketcher/index.html")
    assert r.status_code == 200


def test_ketcher_shell_served(client):
    r = client.get("/ketcher.html")
    assert r.status_code == 200
