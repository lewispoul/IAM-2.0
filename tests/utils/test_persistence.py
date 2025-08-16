import json
from backend.utils import persistence

def test_save_result_json(tmp_path):
    payload = {"ok": True, "data": {"foo": 1}, "errors": []}
    name = "test"
    # Patch RESULTS to tmp_path
    persistence.RESULTS = tmp_path / "Results"
    persistence.RESULTS.mkdir(parents=True, exist_ok=True)
    fp = persistence.save_result_json(name, payload)
    assert fp.endswith(f"_{name}.json")
    data = json.loads((tmp_path / "Results" / fp.split("/")[-1]).read_text())
    assert data == payload

def test_append_benchmark_row(tmp_path):
    persistence.BASE = tmp_path
    row = {"name": "test", "ok": True, "Pcj": 1.0, "Tcj": 2.0}
    fp = persistence.append_benchmark_row(row)
    assert fp.endswith("benchmark_auto.csv")
    lines = (tmp_path / "benchmark_auto.csv").read_text().splitlines()
    assert lines[0].startswith("Pcj,Tcj,name,ok") or lines[0].startswith("name,ok,Pcj,Tcj")
    assert any("test" in l for l in lines)
