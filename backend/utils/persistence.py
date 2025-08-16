from pathlib import Path
import json, csv, time
import os
from typing import Dict, List, Optional, Any

BASE = Path(os.environ.get("IAM_RESULTS_BASE", "IAM_Knowledge"))
RESULTS = BASE / "Results"
EXPORTS = BASE / "Exports"
for d in (RESULTS, EXPORTS):
    d.mkdir(parents=True, exist_ok=True)

def save_result_json(name: str, payload: dict) -> str:
    ts = time.strftime("%Y%m%d_%H%M%S")
    fp = RESULTS / f"{ts}_{name}.json"
    fp.write_text(json.dumps(payload, indent=2))
    return str(fp)

def append_benchmark_row(row: dict) -> str:
    BASE.mkdir(parents=True, exist_ok=True)
    fp = BASE / "benchmark_auto.csv"
    new = not fp.exists()
    with fp.open("a", newline="") as f:
        fieldnames = sorted(row.keys())
        w = csv.DictWriter(f, fieldnames=fieldnames)
        if new: w.writeheader()
        w.writerow(row)
    return str(fp)

def save_result(calc_id: str, result_data: Dict[str, Any]) -> str:
    """Save calculation result by calculation ID."""
    calc_dir = RESULTS / calc_id
    calc_dir.mkdir(exist_ok=True)
    result_file = calc_dir / "result.json"
    result_file.write_text(json.dumps(result_data, indent=2))
    return str(result_file)

def get_result(calc_id: str) -> Optional[Dict[str, Any]]:
    """Get calculation result by calculation ID."""
    calc_dir = RESULTS / calc_id
    result_file = calc_dir / "result.json"
    
    if result_file.exists():
        try:
            return json.loads(result_file.read_text())
        except (json.JSONDecodeError, IOError):
            return None
    return None

def list_results() -> List[Dict[str, Any]]:
    """List all available calculation results."""
    results = []
    
    if not RESULTS.exists():
        return results
    
    for calc_dir in RESULTS.iterdir():
        if calc_dir.is_dir():
            result_file = calc_dir / "result.json"
            if result_file.exists():
                try:
                    result_data = json.loads(result_file.read_text())
                    results.append(result_data)
                except (json.JSONDecodeError, IOError):
                    continue
    
    return results
