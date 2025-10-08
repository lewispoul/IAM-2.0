#!/usr/bin/env python3
import os
import pathlib
import datetime
ROOT = pathlib.Path(os.getcwd())
modules = ["nox", "iam", "pinox"]
edges = set()


def scan(pkg):
    p = ROOT/pkg
    if not p.is_dir():
        return
    for f in p.rglob("*.py"):
        text = f.read_text(encoding="utf-8", errors="ignore")
        for line in text.splitlines():
            s = line.strip()
            if s.startswith(("from ", "import ")):
                for m in modules:
                    if f"{m}." in s or s.startswith(f"from {m} ") or s.startswith(f"import {m}"):
                        if pkg != m:
                            edges.add((pkg, m))


for m in modules:
    scan(m)
out = ["flowchart LR"]+[f'  {m}["{m.upper()}"]' for m in modules] + \
    [f"  {a} --> {b}" for a, b in sorted(edges)]
report_dir = ROOT/"reports" / \
    f"Arch_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
report_dir.mkdir(parents=True, exist_ok=True)
(report_dir/"architecture.mmd").write_text("\n".join(out), encoding="utf-8")
print(f"Generated: {report_dir/'architecture.mmd'}")
