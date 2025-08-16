# tests/api/test_ketcher_run.py
import os
from pathlib import Path
from unittest.mock import patch


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


def test_run_basic_success(client, tmp_results_dir):
    req = {"xyz": "2\nMethane\nC 0.0 0.0 0.0\nH 0.0 1.0 0.0"}
    r = client.post("/run/", json=req)
    body = assert_ok(r)
    assert isinstance(body["data"], dict)


def test_run_missing_xyz(client):
    r = client.post("/run/", json={})
    assert r.status_code == 422


def test_run_empty_xyz(client, tmp_results_dir):
    r = client.post("/run/", json={"xyz": ""})
    body = assert_fail(r)
    assert any("empty" in e.lower() or "xyz" in e.lower() for e in body["errors"])


def test_run_invalid_xyz_format(client, tmp_results_dir):
    r = client.post("/run/", json={"xyz": "not_valid_xyz"})
    body = assert_fail(r)
    assert body["errors"]


def test_run_with_options(client, tmp_results_dir):
    req = {
        "xyz": "2\nMethane\nC 0.0 0.0 0.0\nH 0.0 1.0 0.0",
        "options": {"method": "B3LYP", "basis": "6-31G"}
    }
    r = client.post("/run/", json=req)
    body = assert_ok(r)
    assert "calculation_id" in body["data"]


@patch("backend.main.run_calc_task")
def test_run_internal_task_call(mock_run, client, tmp_results_dir):
    mock_run.return_value = {"calc_id": "test123", "status": "submitted"}
    req = {"xyz": "2\nMethane\nC 0.0 0.0 0.0\nH 0.0 1.0 0.0"}
    r = client.post("/run/", json=req)
    body = assert_ok(r)
    mock_run.assert_called_once()
    args, kwargs = mock_run.call_args
    assert "xyz" in kwargs or len(args) >= 1


def test_run_creates_files_in_results(client, tmp_results_dir):
    req = {"xyz": "2\nMethane\nC 0.0 0.0 0.0\nH 0.0 1.0 0.0"}
    r = client.post("/run/", json=req)
    body = assert_ok(r)
    
    # Check if any files were created in the results directory
    results_files = list(tmp_results_dir.rglob("*"))
    assert len(results_files) > 0, "Expected files to be created in results directory"
