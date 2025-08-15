# tests/api/test_export_and_persistence.py
import json
from pathlib import Path


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


def test_export_csv_success(client, tmp_results_dir):
    # Setup: create a test calculation result
    test_data = {
        "calculation_id": "test123",
        "molecule": "C",
        "energy": -40.123,
        "homo": -0.5,
        "lumo": 0.2
    }
    result_file = tmp_results_dir / "test123" / "result.json"
    result_file.parent.mkdir(exist_ok=True)
    result_file.write_text(json.dumps(test_data))
    
    # Test export
    r = client.post("/export/csv", json={"calculation_ids": ["test123"]})
    body = assert_ok(r)
    
    assert "csv_content" in body["data"]
    assert "filename" in body["data"]
    csv_content = body["data"]["csv_content"]
    assert "calculation_id" in csv_content
    assert "test123" in csv_content


def test_export_csv_nonexistent_calc(client, tmp_results_dir):
    r = client.post("/export/csv", json={"calculation_ids": ["nonexistent"]})
    body = assert_fail(r)
    assert any("not found" in e.lower() or "exist" in e.lower() for e in body["errors"])


def test_export_csv_empty_list(client):
    r = client.post("/export/csv", json={"calculation_ids": []})
    assert r.status_code == 422


def test_export_csv_missing_field(client):
    r = client.post("/export/csv", json={})
    assert r.status_code == 422


def test_export_json_success(client, tmp_results_dir):
    # Setup: create a test calculation result
    test_data = {
        "calculation_id": "test456",
        "molecule": "CCO",
        "energy": -154.789,
        "properties": {"dipole": 1.23}
    }
    result_file = tmp_results_dir / "test456" / "result.json"
    result_file.parent.mkdir(exist_ok=True)
    result_file.write_text(json.dumps(test_data))
    
    # Test export
    r = client.post("/export/json", json={"calculation_ids": ["test456"]})
    body = assert_ok(r)
    
    assert "json_content" in body["data"]
    assert "filename" in body["data"]
    json_content = json.loads(body["data"]["json_content"])
    assert isinstance(json_content, list)
    assert len(json_content) == 1
    assert json_content[0]["calculation_id"] == "test456"


def test_export_json_multiple_calculations(client, tmp_results_dir):
    # Setup: create multiple test calculations
    for calc_id in ["calc1", "calc2"]:
        test_data = {
            "calculation_id": calc_id,
            "molecule": "C",
            "energy": -40.0
        }
        result_file = tmp_results_dir / calc_id / "result.json"
        result_file.parent.mkdir(exist_ok=True)
        result_file.write_text(json.dumps(test_data))
    
    r = client.post("/export/json", json={"calculation_ids": ["calc1", "calc2"]})
    body = assert_ok(r)
    
    json_content = json.loads(body["data"]["json_content"])
    assert len(json_content) == 2
    calc_ids = {item["calculation_id"] for item in json_content}
    assert calc_ids == {"calc1", "calc2"}


def test_list_calculations_success(client, tmp_results_dir):
    # Setup: create test calculations
    for i, calc_id in enumerate(["list_test1", "list_test2"]):
        test_data = {
            "calculation_id": calc_id,
            "molecule": f"C{i}",
            "timestamp": f"2024-01-{i+10}T12:00:00"
        }
        result_file = tmp_results_dir / calc_id / "result.json"
        result_file.parent.mkdir(exist_ok=True)
        result_file.write_text(json.dumps(test_data))
    
    r = client.get("/calculations/")
    body = assert_ok(r)
    
    assert "calculations" in body["data"]
    calculations = body["data"]["calculations"]
    assert isinstance(calculations, list)
    assert len(calculations) >= 2
    
    # Check that our test calculations are included
    calc_ids = {calc["calculation_id"] for calc in calculations}
    assert "list_test1" in calc_ids
    assert "list_test2" in calc_ids


def test_list_calculations_empty(client, tmp_results_dir):
    # Ensure results directory is empty
    for item in tmp_results_dir.iterdir():
        if item.is_dir():
            for subitem in item.rglob("*"):
                if subitem.is_file():
                    subitem.unlink()
            item.rmdir()
    
    r = client.get("/calculations/")
    body = assert_ok(r)
    
    assert "calculations" in body["data"]
    assert body["data"]["calculations"] == []


def test_get_calculation_details_success(client, tmp_results_dir):
    # Setup: create a detailed test calculation
    test_data = {
        "calculation_id": "detail_test",
        "molecule": "C2H6",
        "energy": -79.123,
        "homo": -0.45,
        "lumo": 0.15,
        "properties": {
            "dipole_moment": 0.0,
            "total_charge": 0
        },
        "files": ["input.xyz", "output.log"]
    }
    result_file = tmp_results_dir / "detail_test" / "result.json"
    result_file.parent.mkdir(exist_ok=True)
    result_file.write_text(json.dumps(test_data))
    
    r = client.get("/calculations/detail_test")
    body = assert_ok(r)
    
    assert "calculation" in body["data"]
    calc_data = body["data"]["calculation"]
    assert calc_data["calculation_id"] == "detail_test"
    assert calc_data["molecule"] == "C2H6"
    assert "properties" in calc_data


def test_get_calculation_details_not_found(client, tmp_results_dir):
    r = client.get("/calculations/nonexistent_calc")
    body = assert_fail(r)
    assert any("not found" in e.lower() for e in body["errors"])


def test_calculation_persistence_across_requests(client, tmp_results_dir):
    # First, run a calculation
    req = {"xyz": "2\nTest molecule\nC 0.0 0.0 0.0\nH 0.0 1.0 0.0"}
    r1 = client.post("/run/", json=req)
    body1 = assert_ok(r1)
    calc_id = body1["data"].get("calculation_id")
    
    if calc_id:
        # Then, verify it persists in the list
        r2 = client.get("/calculations/")
        body2 = assert_ok(r2)
        calc_ids = {calc["calculation_id"] for calc in body2["data"]["calculations"]}
        assert calc_id in calc_ids
