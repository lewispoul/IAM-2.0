"""
XTB calculation stub for IAM2.0.
Provides deterministic mock results for development and testing.
"""
import asyncio
import logging
from typing import Dict, Any, Optional
import hashlib

logger = logging.getLogger(__name__)


async def run_xtb_calculation(
    xyz: str,
    method: str = "gfn2",
    charge: int = 0,
    multiplicity: int = 1
) -> Dict[str, Any]:
    """
    Stub implementation of XTB calculation.
    
    Args:
        xyz: XYZ coordinates string
        method: XTB method (gfn1, gfn2, gfn-ff)
        charge: Molecular charge
        multiplicity: Spin multiplicity
        
    Returns:
        Deterministic mock results based on input hash
    """
    logger.info(f"Running XTB stub calculation: {method}")
    
    # Simulate calculation time
    await asyncio.sleep(0.1)
    
    # Generate deterministic results based on input hash
    input_hash = hashlib.md5(f"{xyz}{method}{charge}{multiplicity}".encode()).hexdigest()
    hash_int = int(input_hash[:8], 16)
    
    # Count atoms for realistic scaling
    lines = xyz.strip().split('\n')
    atom_count = int(lines[0]) if lines and lines[0].isdigit() else 1
    
    # Deterministic "results" based on hash
    base_energy = -50.0 * atom_count  # Rough scaling with atom count
    energy_variation = (hash_int % 1000) / 1000.0 * 10.0  # 0-10 variation
    final_energy = base_energy - energy_variation
    
    return {
        "energy": final_energy,
        "energy_units": "hartree",
        "method": method,
        "converged": True,
        "optimization_steps": hash_int % 50 + 10,
        "dipole_moment": {
            "magnitude": (hash_int % 100) / 100.0 * 5.0,  # 0-5 Debye
            "x": ((hash_int % 37) - 18) / 18.0,
            "y": (((hash_int >> 8) % 37) - 18) / 18.0,
            "z": (((hash_int >> 16) % 37) - 18) / 18.0
        },
        "homo_lumo_gap": (hash_int % 200) / 1000.0 + 0.1,  # 0.1-0.3 hartree
        "total_charge": float(charge),
        "multiplicity": multiplicity,
        "atom_count": atom_count,
        "calculation_time": 0.1,
        "stub": True,
        "hash_seed": input_hash[:8]
    }


def get_available_methods() -> list:
    """Return list of available XTB methods."""
    return ["gfn1", "gfn2", "gfn-ff"]


def validate_method(method: str) -> bool:
    """Validate XTB method name."""
    return method in get_available_methods()


async def run_xtb_optimization(
    xyz: str,
    method: str = "gfn2",
    charge: int = 0,
    multiplicity: int = 1
) -> Dict[str, Any]:
    """
    Stub implementation of XTB geometry optimization.
    
    Returns optimized geometry along with energy.
    """
    logger.info(f"Running XTB optimization stub: {method}")
    
    # Run base calculation
    results = await run_xtb_calculation(xyz, method, charge, multiplicity)
    
    # Add optimization-specific data
    lines = xyz.strip().split('\n')
    if len(lines) > 2:
        # Return "optimized" geometry (slightly perturbed original)
        optimized_xyz = lines[0] + '\n' + "Optimized by XTB stub\n"
        hash_int = int(results["hash_seed"], 16)
        
        for i, line in enumerate(lines[2:], 2):
            parts = line.split()
            if len(parts) >= 4:
                symbol = parts[0]
                x = float(parts[1]) + ((hash_int >> (i*2)) % 21 - 10) / 1000.0  # ±0.01 perturbation
                y = float(parts[2]) + ((hash_int >> (i*2+8)) % 21 - 10) / 1000.0
                z = float(parts[3]) + ((hash_int >> (i*2+16)) % 21 - 10) / 1000.0
                optimized_xyz += f"{symbol} {x:.6f} {y:.6f} {z:.6f}\n"
    else:
        optimized_xyz = xyz
    
    results.update({
        "optimized_xyz": optimized_xyz,
        "optimization_converged": True,
        "max_force": (hash_int % 100) / 100000.0,  # Small force
        "rms_force": (hash_int % 50) / 200000.0,
    })
    
    return results
