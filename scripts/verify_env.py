#!/usr/bin/env python
import importlib
import sys

pkgs = [
    ("rdkit", "rdkit.__version__"),
    ("psi4", "psi4.__version__"),
    ("cantera", "cantera.__version__"),
    # xtb might be via python wrapper or subprocess; try common options
    ("xtb", "__version__"),
]

def get_version(mod_name: str, attr_expr: str):
    try:
        mod = importlib.import_module(mod_name)
        # Evaluate attr path safely
        obj = mod
        for part in attr_expr.split(".")[1:]:
            obj = getattr(obj, part)
        return str(obj)
    except Exception as e:
        return f"unavailable ({e.__class__.__name__}: {e})"

def main():
    print("=== IAM 2.0 environment check ===")
    for name, attr in pkgs:
        ver = get_version(name, attr)
        print(f"{name:8s} : {ver}")

    # Also show Python and platform
    print(f"Python   : {sys.version.split()[0]}")
    print(f"Platform : {sys.platform}")
    print("Note: If rdkit/psi4/cantera are 'unavailable', activate chem-env or install requirements.")

if __name__ == "__main__":
    main()
