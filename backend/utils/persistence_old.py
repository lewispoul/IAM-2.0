from pathlib import Path
import json, csv, time
import os
import re
from typing import Dict, List, Optional, Any

# Honor the IAM_RESULTS_BASE environment variable used by tests
"""Persistence layer for IAM-2.0."""

import json
import csv
import time
import re
import os
from pathlib import Path
from typing import Dict, Any, List, Optional

# Dynamic path resolution to respect environment changes
def get_results_base() -> Path:
    """Get the results base directory from environment or default."""
    base_path = os.getenv("IAM_RESULTS_BASE", str(Path.cwd() / "IAM_Knowledge"))
    return Path(base_path)

def get_results_dir() -> Path:
    """Get the Results subdirectory."""
    return get_results_base() / "Results"

def get_exports_dir() -> Path:
    """Get the Exports subdirectory."""
    return get_results_base() / "Exports"
BASE = RESULTS_BASE
RESULTS = RESULTS_BASE / "Results"
EXPORTS = RESULTS_BASE / "Exports"

def ensure_base() -> Path:
    """Create RESULTS_BASE if missing. Return the Path."""
    RESULTS.mkdir(parents=True, exist_ok=True)
    EXPORTS.mkdir(parents=True, exist_ok=True)
    return RESULTS_BASE

def calc_dir(calc_id: str) -> Path:
    """Return RESULTS_BASE / calc_id with path traversal protection."""
    # Validate calc_id to prevent path traversal
    if not re.match(r'^[A-Za-z0-9._-]+$', calc_id):
        raise ValueError(f"Invalid calc_id: {calc_id}")
    if '..' in calc_id or '/' in calc_id or '\\' in calc_id:
        raise ValueError(f"Path traversal attempt in calc_id: {calc_id}")
    
    ensure_base()
    return RESULTS / calc_id

# Legacy compatibility functions
def save_result(calc_id: str, result_data: Dict[str, Any]) -> str:
    """Save calculation result by calculation ID (legacy compatibility)."""
    result_file = save_calc(calc_id, result_data)
    return str(result_file)

def get_result(calc_id: str) -> Optional[Dict[str, Any]]:
    """Get calculation result by calculation ID (legacy compatibility)."""
    return get_calc(calc_id)

def list_results() -> List[Dict[str, Any]]:
    """List all available calculation results (legacy compatibility)."""
    results = []
    calc_ids = list_calcs()
    
    for calc_id in calc_ids:
        data = get_calc(calc_id)
        if data is not None:
            results.append(data)
    
    return results

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
    
    if not RESULTS.exists():
        return calc_ids
    
    for calc_directory in RESULTS.iterdir():
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

def get_results_bulk(calc_ids: List[str]) -> List[Dict[str, Any]]:
    """Return list of {calc_id, data} records for those found; ignore missing."""
    results = []
    
    for calc_id in calc_ids:
        data = get_calc(calc_id)
        if data is not None:
            results.append({"calc_id": calc_id, "data": data})
    
    return results
