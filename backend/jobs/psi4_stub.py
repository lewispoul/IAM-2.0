"""
Psi4 calculation stub for IAM2.0.
Provides deterministic mock results for quantum chemistry calculations.
"""
import asyncio
import logging
from typing import Dict, Any, Optional, List
import hashlib
import math

logger = logging.getLogger(__name__)


async def run_psi4_calculation(
    xyz: str,
    method: str = "HF",
    charge: int = 0,
    multiplicity: int = 1,
    basis: str = "6-31G"
) -> Dict[str, Any]:
    """
    Stub implementation of Psi4 calculation.
    
    Args:
        xyz: XYZ coordinates string
        method: Quantum chemistry method
        charge: Molecular charge
        multiplicity: Spin multiplicity
        basis: Basis set name
        
    Returns:
        Deterministic mock results based on input hash
    """
    logger.info(f"Running Psi4 stub calculation: {method}/{basis}")
    
    # Simulate longer calculation time for QM methods
    await asyncio.sleep(0.5)
    
    # Generate deterministic results
    input_hash = hashlib.md5(f"{xyz}{method}{charge}{multiplicity}{basis}".encode()).hexdigest()
    hash_int = int(input_hash[:8], 16)
    
    # Count atoms for realistic scaling
    lines = xyz.strip().split('\n')
    atom_count = int(lines[0]) if lines and lines[0].isdigit() else 1
    
    # Method-dependent energy scaling
    method_factors = {
        "HF": 1.0,
        "B3LYP": 1.1,
        "MP2": 1.05,
        "CCSD": 1.02,
        "CCSD(T)": 1.01
    }
    method_factor = method_factors.get(method, 1.0)
    
    # Basis set correction
    basis_corrections = {
        "STO-3G": 0.95,
        "6-31G": 1.0,
        "6-31G*": 1.02,
        "6-31+G*": 1.03,
        "cc-pVDZ": 1.04,
        "cc-pVTZ": 1.05
    }
    basis_correction = basis_corrections.get(basis, 1.0)
    
    # Calculate mock energy
    base_energy = -75.0 * atom_count * method_factor * basis_correction
    energy_variation = (hash_int % 10000) / 10000.0 * 5.0
    final_energy = base_energy - energy_variation
    
    # Generate orbital energies (mock HOMO/LUMO)
    num_orbitals = atom_count * 7  # Rough estimate
    homo_energy = -0.3 - (hash_int % 100) / 1000.0
    lumo_energy = 0.1 + (hash_int % 150) / 1000.0
    
    results = {
        "energy": final_energy,
        "energy_units": "hartree",
        "method": method,
        "basis": basis,
        "converged": True,
        "scf_iterations": (hash_int % 30) + 10,
        "homo_energy": homo_energy,
        "lumo_energy": lumo_energy,
        "homo_lumo_gap": lumo_energy - homo_energy,
        "total_charge": float(charge),
        "multiplicity": multiplicity,
        "atom_count": atom_count,
        "dipole_moment": {
            "magnitude": (hash_int % 200) / 100.0,  # 0-2 Debye
            "x": ((hash_int % 73) - 36) / 36.0,
            "y": (((hash_int >> 8) % 73) - 36) / 36.0,
            "z": (((hash_int >> 16) % 73) - 36) / 36.0,
            "units": "debye"
        },
        "calculation_time": 0.5,
        "stub": True,
        "hash_seed": input_hash[:8]
    }
    
    # Add method-specific data
    if method in ["MP2", "CCSD", "CCSD(T)"]:
        results["correlation_energy"] = -(hash_int % 1000) / 10000.0  # Small negative
        results["reference_energy"] = final_energy - results["correlation_energy"]
    
    if method in ["B3LYP", "PBE", "M06"]:
        results["exchange_energy"] = -(hash_int % 500) / 1000.0
        results["correlation_energy"] = -(hash_int % 300) / 2000.0
    
    return results


def get_available_methods() -> List[str]:
    """Return list of available Psi4 methods."""
    return [
        "HF", "B3LYP", "PBE", "M06", "wB97X-D",
        "MP2", "CCSD", "CCSD(T)",
        "CASSCF", "CASPT2"
    ]


def get_available_basis_sets() -> List[str]:
    """Return list of available basis sets."""
    return [
        "STO-3G", "6-31G", "6-31G*", "6-31+G*", "6-311G*",
        "cc-pVDZ", "cc-pVTZ", "cc-pVQZ",
        "aug-cc-pVDZ", "aug-cc-pVTZ"
    ]


def validate_method(method: str) -> bool:
    """Validate Psi4 method name."""
    return method in get_available_methods()


def validate_basis(basis: str) -> bool:
    """Validate basis set name."""
    return basis in get_available_basis_sets()


async def run_psi4_optimization(
    xyz: str,
    method: str = "HF",
    charge: int = 0,
    multiplicity: int = 1,
    basis: str = "6-31G"
) -> Dict[str, Any]:
    """
    Stub implementation of Psi4 geometry optimization.
    
    Returns optimized geometry along with final energy.
    """
    logger.info(f"Running Psi4 optimization stub: {method}/{basis}")
    
    # Run base calculation
    results = await run_psi4_calculation(xyz, method, charge, multiplicity, basis)
    
    # Add optimization-specific data
    lines = xyz.strip().split('\n')
    if len(lines) > 2:
        # Return "optimized" geometry
        optimized_xyz = lines[0] + '\n' + f"Optimized by Psi4 stub ({method}/{basis})\n"
        hash_int = int(results["hash_seed"], 16)
        
        for i, line in enumerate(lines[2:], 2):
            parts = line.split()
            if len(parts) >= 4:
                symbol = parts[0]
                # Smaller perturbations for QM optimization
                x = float(parts[1]) + ((hash_int >> (i*3)) % 11 - 5) / 2000.0  # ±0.0025 perturbation
                y = float(parts[2]) + ((hash_int >> (i*3+8)) % 11 - 5) / 2000.0
                z = float(parts[3]) + ((hash_int >> (i*3+16)) % 11 - 5) / 2000.0
                optimized_xyz += f"{symbol} {x:.6f} {y:.6f} {z:.6f}\n"
    else:
        optimized_xyz = xyz
    
    results.update({
        "optimized_xyz": optimized_xyz,
        "optimization_converged": True,
        "optimization_steps": (hash_int % 20) + 5,
        "max_gradient": (hash_int % 50) / 500000.0,  # Very small gradient
        "rms_gradient": (hash_int % 30) / 1000000.0,
        "final_energy": results["energy"] - 0.001  # Slightly lower after optimization
    })
    
    return results


async def run_psi4_frequency(
    xyz: str,
    method: str = "HF",
    charge: int = 0,
    multiplicity: int = 1,
    basis: str = "6-31G"
) -> Dict[str, Any]:
    """
    Stub implementation of Psi4 frequency calculation.
    
    Returns vibrational frequencies and thermodynamic data.
    """
    logger.info(f"Running Psi4 frequency stub: {method}/{basis}")
    
    # Run base calculation
    results = await run_psi4_calculation(xyz, method, charge, multiplicity, basis)
    
    # Count atoms for proper frequency count
    lines = xyz.strip().split('\n')
    atom_count = int(lines[0]) if lines and lines[0].isdigit() else 1
    num_frequencies = 3 * atom_count - 6  # 3N-6 for non-linear molecules
    
    if num_frequencies < 0:
        num_frequencies = 0
    
    hash_int = int(results["hash_seed"], 16)
    
    # Generate mock frequencies (cm⁻¹)
    frequencies = []
    for i in range(num_frequencies):
        # Typical vibrational range: 200-4000 cm⁻¹
        freq = 200 + ((hash_int >> (i % 24)) % 3800)
        frequencies.append(freq)
    
    frequencies.sort()
    
    # Mock thermodynamic data (298.15 K)
    results.update({
        "frequencies": frequencies,
        "frequency_units": "cm^-1",
        "zero_point_energy": sum(f * 0.00012 for f in frequencies[:10]),  # Rough ZPE
        "thermal_energy": results["energy"] + 0.002,
        "enthalpy": results["energy"] + 0.0024,
        "entropy": (hash_int % 100) / 1000.0 + 0.05,  # cal/mol/K
        "gibbs_free_energy": results["energy"] - 0.005,
        "temperature": 298.15,
        "pressure": 1.0,
        "units": {
            "energy": "hartree",
            "entropy": "cal/mol/K",
            "temperature": "K",
            "pressure": "atm"
        }
    })
    
    return results
