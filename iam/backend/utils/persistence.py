from pathlib import Path
import json, csv, time
import os
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
