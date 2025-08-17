"""Environment configuration utilities for IAM-2.0."""

import os
from pathlib import Path


def flag(name: str) -> bool:
    """Check if an environment variable is set to a truthy value.
    
    Args:
        name: Environment variable name
        
    Returns:
        True if the variable is set to "1", "true", "TRUE", "yes", or "YES"
    """
    return os.getenv(name, "0") in {"1", "true", "TRUE", "yes", "YES"}


def results_base() -> Path:
    """Get the results base directory from environment or default.
    
    Returns:
        Path to the results base directory (IAM_Knowledge/Results by default)
    """
    base_path = os.getenv("IAM_RESULTS_BASE", "IAM_Knowledge/Results")
    return Path(base_path)


def exports_base() -> Path:
    """Get the exports base directory.
    
    Returns:
        Path to the exports directory (IAM_Knowledge/Exports by default)
    """
    return Path("IAM_Knowledge/Exports")
