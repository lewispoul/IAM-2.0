"""
Enhanced Psi4 integration for IAM-2.0.
Supports real Psi4 quantum chemistry calculations when enabled via environment flag.
"""
import asyncio
import logging
import os
import tempfile
import subprocess
import json
from typing import Dict, Any, Optional
from pathlib import Path

from ..utils.env import flag
from .psi4_stub import run_psi4_calculation as stub_run_psi4_calculation

logger = logging.getLogger(__name__)


class Psi4Calculator:
    """Enhanced Psi4 calculator with environment flag support."""
    
    def __init__(self):
        self.enabled = flag("IAM_ENABLE_PSI4")
        self.psi4_executable = self._find_psi4_executable()
        self.python_psi4_available = self._check_python_psi4()
        
        if self.enabled and not (self.psi4_executable or self.python_psi4_available):
            logger.warning("IAM_ENABLE_PSI4=true but Psi4 not found, using stub")
            self.enabled = False
    
    def _find_psi4_executable(self) -> Optional[str]:
        """Find Psi4 executable in PATH or specified location."""
        # Check common locations
        for possible_name in ['psi4', 'psi4.exe']:
            try:
                result = subprocess.run(
                    ['which', possible_name], 
                    capture_output=True, 
                    text=True, 
                    check=False
                )
                if result.returncode == 0:
                    return possible_name
            except Exception:
                pass
        
        # Check custom path if provided
        custom_path = os.getenv("PSI4_EXECUTABLE")
        if custom_path and Path(custom_path).exists():
            return custom_path
            
        return None
    
    def _check_python_psi4(self) -> bool:
        """Check if Psi4 is available as a Python module."""
        try:
            import psi4
            return True
        except ImportError:
            return False
    
    async def run_calculation(
        self,
        xyz: str,
        method: str = "HF",
        charge: int = 0,
        multiplicity: int = 1,
        basis: str = "6-31G",
        optimization: bool = False
    ) -> Dict[str, Any]:
        """
        Run Psi4 calculation with environment flag control.
        
        Args:
            xyz: XYZ coordinates string
            method: Quantum chemistry method (HF, B3LYP, MP2, etc.)
            charge: Molecular charge
            multiplicity: Spin multiplicity
            basis: Basis set name
            optimization: Whether to perform geometry optimization
            
        Returns:
            Calculation results dict
        """
        if not self.enabled:
            logger.info("Psi4 calculations disabled (IAM_ENABLE_PSI4=false), using stub")
            return await stub_run_psi4_calculation(xyz, method, charge, multiplicity, basis)
        
        logger.info(f"Running real Psi4 calculation: {method}/{basis} (enabled via IAM_ENABLE_PSI4)")
        
        try:
            if self.python_psi4_available:
                return await self._run_python_psi4(xyz, method, charge, multiplicity, basis, optimization)
            elif self.psi4_executable:
                return await self._run_executable_psi4(xyz, method, charge, multiplicity, basis, optimization)
            else:
                raise RuntimeError("No Psi4 implementation available")
        except Exception as e:
            logger.error(f"Real Psi4 calculation failed: {e}, falling back to stub")
            return await stub_run_psi4_calculation(xyz, method, charge, multiplicity, basis)
    
    async def _run_python_psi4(
        self,
        xyz: str,
        method: str,
        charge: int,
        multiplicity: int,
        basis: str,
        optimization: bool
    ) -> Dict[str, Any]:
        """Run Psi4 calculation using Python API."""
        try:
            import psi4
        except ImportError:
            raise RuntimeError("Psi4 Python module not available")
        
        # Set up Psi4
        psi4.core.set_output_file("psi4_output.dat", False)
        psi4.set_memory('500 MB')
        
        # Parse XYZ and create molecule
        lines = xyz.strip().split('\n')
        atom_count = int(lines[0])
        
        # Build molecule string for Psi4
        mol_string = f"{charge} {multiplicity}\n"
        for line in lines[2:2+atom_count]:
            parts = line.strip().split()
            if len(parts) >= 4:
                symbol, x, y, z = parts[0], parts[1], parts[2], parts[3]
                mol_string += f"{symbol} {x} {y} {z}\n"
        
        molecule = psi4.geometry(mol_string)
        
        # Set basis set
        psi4.set_options({'basis': basis})
        
        # Run calculation
        if optimization:
            energy = psi4.optimize(method, molecule=molecule)
            # Get optimized geometry
            opt_xyz = self._extract_optimized_xyz(molecule, atom_count)
        else:
            energy = psi4.energy(method, molecule=molecule)
            opt_xyz = None
        
        # Extract properties
        wavefunction = psi4.core.get_global_option("WFN")
        
        results = {
            "energy": float(energy),
            "energy_units": "hartree",
            "method": method,
            "basis": basis,
            "charge": float(charge),
            "multiplicity": multiplicity,
            "converged": True,
            "atom_count": atom_count,
            "stub": False
        }
        
        if opt_xyz:
            results.update({
                "optimized_xyz": opt_xyz,
                "optimization_converged": True
            })
        
        return results
    
    def _extract_optimized_xyz(self, molecule, atom_count: int) -> str:
        """Extract optimized geometry from Psi4 molecule object."""
        try:
            # Get geometry from molecule
            geom = molecule.geometry()
            xyz_lines = [str(atom_count), "Optimized by Psi4"]
            
            for i in range(atom_count):
                symbol = molecule.symbol(i)
                x = geom.get(i, 0) * 0.52917721067  # Convert bohr to angstrom
                y = geom.get(i, 1) * 0.52917721067
                z = geom.get(i, 2) * 0.52917721067
                xyz_lines.append(f"{symbol} {x:.6f} {y:.6f} {z:.6f}")
            
            return "\n".join(xyz_lines)
        except Exception:
            return None
    
    async def _run_executable_psi4(
        self,
        xyz: str,
        method: str,
        charge: int,
        multiplicity: int,
        basis: str,
        optimization: bool
    ) -> Dict[str, Any]:
        """Run Psi4 calculation using external executable."""
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create Psi4 input file
            input_file = temp_path / "input.dat"
            psi4_input = self._generate_psi4_input(xyz, method, charge, multiplicity, basis, optimization)
            input_file.write_text(psi4_input)
            
            # Run Psi4 calculation
            cmd = [self.psi4_executable, str(input_file), "-o", "output.dat"]
            
            logger.info(f"Executing: {' '.join(cmd)}")
            
            process = await asyncio.create_subprocess_exec(
                *cmd,
                cwd=temp_dir,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )
            
            stdout, stderr = await process.communicate()
            
            if process.returncode != 0:
                error_msg = stderr.decode() if stderr else "Unknown Psi4 error"
                raise RuntimeError(f"Psi4 failed with return code {process.returncode}: {error_msg}")
            
            # Parse Psi4 output
            output_file = temp_path / "output.dat"
            if output_file.exists():
                output_content = output_file.read_text()
                return self._parse_psi4_output(output_content, method, basis, charge, multiplicity)
            else:
                raise RuntimeError("Psi4 output file not found")
    
    def _generate_psi4_input(
        self,
        xyz: str,
        method: str,
        charge: int,
        multiplicity: int,
        basis: str,
        optimization: bool
    ) -> str:
        """Generate Psi4 input file content."""
        
        lines = xyz.strip().split('\n')
        atom_count = int(lines[0])
        
        # Build molecule block
        mol_block = f"""
molecule {{
    {charge} {multiplicity}
"""
        
        for line in lines[2:2+atom_count]:
            parts = line.strip().split()
            if len(parts) >= 4:
                symbol, x, y, z = parts[0], parts[1], parts[2], parts[3]
                mol_block += f"    {symbol} {x} {y} {z}\n"
        
        mol_block += "}\n"
        
        # Build full input
        input_content = f"""
memory 500 MB

set basis {basis}
set reference rhf

{mol_block}

"""
        
        if optimization:
            input_content += f"optimize('{method}')\n"
        else:
            input_content += f"energy('{method}')\n"
        
        return input_content
    
    def _parse_psi4_output(
        self,
        output: str,
        method: str,
        basis: str,
        charge: int,
        multiplicity: int
    ) -> Dict[str, Any]:
        """Parse Psi4 output and extract relevant data."""
        
        results = {
            "method": method,
            "basis": basis,
            "charge": float(charge),
            "multiplicity": multiplicity,
            "converged": True,
            "stub": False
        }
        
        # Look for energy in output
        for line in output.split('\n'):
            if 'Final Energy:' in line or 'Total Energy =' in line:
                parts = line.split()
                for i, part in enumerate(parts):
                    try:
                        energy = float(part)
                        results["energy"] = energy
                        results["energy_units"] = "hartree"
                        break
                    except ValueError:
                        continue
        
        # Look for HOMO/LUMO
        homo_found = False
        for line in output.split('\n'):
            if 'HOMO' in line and 'LUMO' in line:
                try:
                    # Try to extract HOMO/LUMO energies
                    parts = line.split()
                    if len(parts) >= 6:
                        results["homo_energy"] = float(parts[2])
                        results["lumo_energy"] = float(parts[5])
                        results["homo_lumo_gap"] = results["lumo_energy"] - results["homo_energy"]
                        homo_found = True
                except (ValueError, IndexError):
                    pass
        
        # Look for dipole moment
        for line in output.split('\n'):
            if 'Dipole Moment:' in line or 'Total Dipole' in line:
                try:
                    # Try to extract dipole components
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if 'X=' in part:
                            results["dipole_moment"] = {
                                "x": float(parts[i+1]),
                                "y": float(parts[i+3]),
                                "z": float(parts[i+5]),
                                "units": "debye"
                            }
                            break
                except (ValueError, IndexError):
                    pass
        
        return results


# Global calculator instance
calculator = Psi4Calculator()


async def run_psi4_calculation_enhanced(
    xyz: str,
    method: str = "HF",
    charge: int = 0,
    multiplicity: int = 1,
    basis: str = "6-31G"
) -> Dict[str, Any]:
    """
    Enhanced Psi4 calculation with environment flag support.
    
    This function respects the IAM_ENABLE_PSI4 environment variable:
    - When True: Attempts to run real Psi4 calculations
    - When False (default): Uses deterministic stub implementation
    """
    return await calculator.run_calculation(xyz, method, charge, multiplicity, basis)


async def run_psi4_optimization_enhanced(
    xyz: str,
    method: str = "HF",
    charge: int = 0,
    multiplicity: int = 1,
    basis: str = "6-31G"
) -> Dict[str, Any]:
    """
    Enhanced Psi4 geometry optimization with environment flag support.
    """
    return await calculator.run_calculation(xyz, method, charge, multiplicity, basis, optimization=True)


def get_psi4_status() -> Dict[str, Any]:
    """Get current Psi4 integration status."""
    return {
        "enabled": calculator.enabled,
        "executable_found": calculator.psi4_executable is not None,
        "executable_path": calculator.psi4_executable,
        "python_module_available": calculator.python_psi4_available,
        "environment_flag": flag("IAM_ENABLE_PSI4"),
        "mode": "real" if calculator.enabled else "stub",
        "implementation": "python" if calculator.python_psi4_available else "executable" if calculator.psi4_executable else "none"
    }
