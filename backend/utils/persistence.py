"""Persistence layer for IAM-2.0."""

import json
import csv
import time
import re
import os
from pathlib import Path
from typing import Dict, Any, List, Optional

# Backward compatibility for tests
BASE = Path.cwd() / "IAM_Knowledge"
RESULTS = BASE / "Results"

# Dynamic path resolution to respect environment changes
def get_results_base() -> Path:
    """Get the results base directory from environment or default."""
    # Use the global BASE if it's been modified (for tests), otherwise use environment
    global BASE
    if BASE != Path.cwd() / "IAM_Knowledge":
        return BASE
    base_path = os.getenv("IAM_RESULTS_BASE", str(Path.cwd() / "IAM_Knowledge"))
    return Path(base_path)

def get_results_dir() -> Path:
    """Get the directory where calc results are stored (directly in base, not a subdirectory)."""
    # Use the global RESULTS if it's been modified (for tests), otherwise compute
    global RESULTS
    if RESULTS != BASE / "Results":
        return RESULTS
    return get_results_base()

def get_exports_dir() -> Path:
    """Get the Exports subdirectory."""
    return get_results_base() / "Exports"

def ensure_base() -> Path:
    """Create results directories if missing. Return the results base path."""
    results_dir = get_results_dir()
    exports_dir = get_exports_dir()
    results_dir.mkdir(parents=True, exist_ok=True)
    exports_dir.mkdir(parents=True, exist_ok=True)
    return get_results_base()

def calc_dir(calc_id: str) -> Path:
    """Return Results / calc_id with path traversal protection."""
    # Validate calc_id to prevent path traversal
    if not re.match(r'^[A-Za-z0-9._-]+$', calc_id):
        raise ValueError(f"Invalid calc_id: {calc_id}")
    if '..' in calc_id or '/' in calc_id or '\\' in calc_id:
        raise ValueError(f"Path traversal attempt in calc_id: {calc_id}")
    
    return get_results_dir() / calc_id

def get_results_bulk(calc_ids: List[str]) -> List[Dict[str, Any]]:
    """Return list of {calc_id, data} records for those found; ignore missing."""
    results = []
    
    for calc_id in calc_ids:
        data = get_calc(calc_id)
        if data is not None:
            results.append({"calc_id": calc_id, "data": data})
    
    return results

def save_result_json(name: str, payload: dict) -> str:
    """Save a JSON result with timestamp."""
    results_dir = get_results_dir()
    ensure_base()
    ts = time.strftime("%Y%m%d_%H%M%S")
    fp = results_dir / f"{ts}_{name}.json"
    fp.write_text(json.dumps(payload, indent=2))
    return str(fp)

def append_benchmark_row(row: dict) -> str:
    """Append a row to the benchmark CSV file."""
    base = get_results_base()
    base.mkdir(parents=True, exist_ok=True)
    fp = base / "benchmark_auto.csv"
    new = not fp.exists()
    with fp.open("a", newline="") as f:
        fieldnames = sorted(row.keys())
        w = csv.DictWriter(f, fieldnames=fieldnames)
        if new: w.writeheader()
        w.writerow(row)
    return str(fp)

def save_calc(calc_id: str, payload: dict) -> Path:
    """Save calculation result by calculation ID."""
    calc_directory = calc_dir(calc_id)
    calc_directory.mkdir(parents=True, exist_ok=True)
    result_file = calc_directory / "result.json"
    result_file.write_text(json.dumps(payload, indent=2))
    return result_file

def list_calcs() -> List[str]:
    """Return sorted calc_id list of directories that contain result.json."""
    ensure_base()
    calc_ids = []
    
    results_dir = get_results_dir()
    if not results_dir.exists():
        return calc_ids
    
    for calc_directory in results_dir.iterdir():
        if calc_directory.is_dir():
            result_file = calc_directory / "result.json"
            if result_file.exists():
                calc_ids.append(calc_directory.name)
    
    return sorted(calc_ids)

def get_calc(calc_id: str) -> Optional[Dict[str, Any]]:
    """Get calculation result by calculation ID."""
    try:
        calc_directory = calc_dir(calc_id)
        result_file = calc_directory / "result.json"
        
        if result_file.exists():
            try:
                return json.loads(result_file.read_text())
            except (json.JSONDecodeError, IOError):
                return None
    except ValueError:
        # Invalid calc_id
        return None
    return None

# Legacy compatibility wrappers
def save_result(calc_id: str, payload: dict) -> Path:
    """Legacy wrapper for save_calc."""
    return save_calc(calc_id, payload)

def get_result(calc_id: str) -> Optional[Dict[str, Any]]:
    """Legacy wrapper for get_calc."""
    return get_calc(calc_id)

def list_results() -> List[str]:
    """Legacy wrapper for list_calcs."""
    return list_calcs()
