#!/usr/bin/env python
import importlib
import sys
import subprocess
import shutil
from pathlib import Path

# Python packages to check
pkgs = [
    ("rdkit", "rdkit.__version__"),
    ("psi4", "psi4.__version__"),
    ("cantera", "cantera.__version__"),
    ("fastapi", "fastapi.__version__"),
    ("uvicorn", "uvicorn.__version__"), 
    ("pydantic", "pydantic.__version__"),
    ("numpy", "numpy.__version__"),
    ("pandas", "pandas.__version__"),
    ("pytest", "pytest.__version__"),
    # xtb might be via python wrapper or subprocess; try common options
    ("xtb", "__version__"),
]

# CLI binaries to check
cli_binaries = [
    "xtb",
    "psi4", 
    "obabel",
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

def check_cli_binary(binary_name: str):
    """Check if CLI binary is available in PATH"""
    try:
        if shutil.which(binary_name):
            # Try to get version
            result = subprocess.run([binary_name, "--version"], 
                                  capture_output=True, text=True, timeout=5)
            if result.returncode == 0:
                # Extract first line which usually contains version
                version_line = result.stdout.strip().split('\n')[0]
                return f"available ({version_line})"
            else:
                return "available (version check failed)"
        else:
            return "not in PATH"
    except (subprocess.TimeoutExpired, subprocess.SubprocessError, FileNotFoundError) as e:
        return f"unavailable ({e.__class__.__name__})"

def main():
    print("=== IAM 2.0 Environment Check ===")
    print()
    
    # Check Python packages
    print("Python Packages:")
    print("-" * 40)
    all_available = True
    for name, attr in pkgs:
        ver = get_version(name, attr)
        status = "✅" if "unavailable" not in ver else "❌"
        print(f"{status} {name:12s} : {ver}")
        if "unavailable" in ver:
            all_available = False
    
    print()
    
    # Check CLI binaries  
    print("CLI Binaries:")
    print("-" * 40)
    for binary in cli_binaries:
        status_info = check_cli_binary(binary)
        status = "✅" if "available" in status_info else "❌"
        print(f"{status} {binary:12s} : {status_info}")
        if "available" not in status_info:
            all_available = False
    
    print()
    print("System Information:")
    print("-" * 40)
    print(f"✅ Python       : {sys.version.split()[0]}")
    print(f"✅ Platform     : {sys.platform}")
    
    # Check for Ketcher availability (via CDN - always available)
    print(f"✅ Ketcher      : available (via CDN)")
    
    print()
    if all_available:
        print("🎉 All dependencies are available!")
        print("Ready for IAM-2.0 development and testing.")
    else:
        print("⚠️  Some dependencies are missing.")
        print("Run 'conda activate chem-env' or install missing packages.")
        print("See requirements.txt or chem-env.yaml for installation guide.")
    
    return 0 if all_available else 1

if __name__ == "__main__":
    main()
