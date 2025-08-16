"""
Enhanced XTB integration for IAM-2.0.
Supports real XTB calculations when enabled via environment flag.
"""
import asyncio
import logging
import os
import tempfile
import subprocess
from typing import Dict, Any, Optional
from pathlib import Path

from ..utils.env import flag
from .xtb_stub import run_xtb_calculation as stub_run_xtb_calculation

logger = logging.getLogger(__name__)


class XTBCalculator:
    """Enhanced XTB calculator with environment flag support."""
    
    def __init__(self):
        self.enabled = flag("IAM_ENABLE_XTB")
        self.xtb_executable = self._find_xtb_executable()
        
        if self.enabled and not self.xtb_executable:
            logger.warning("IAM_ENABLE_XTB=true but XTB executable not found, using stub")
            self.enabled = False
    
    def _find_xtb_executable(self) -> Optional[str]:
        """Find XTB executable in PATH or specified location."""
        # Check common locations
        for possible_name in ['xtb', 'xtb.exe']:
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
        custom_path = os.getenv("XTB_EXECUTABLE")
        if custom_path and Path(custom_path).exists():
            return custom_path
            
        return None
    
    async def run_calculation(
        self,
        xyz: str,
        method: str = "gfn2",
        charge: int = 0,
        multiplicity: int = 1,
        optimization: bool = False
    ) -> Dict[str, Any]:
        """
        Run XTB calculation with environment flag control.
        
        Args:
            xyz: XYZ coordinates string
            method: XTB method (gfn1, gfn2, gfn-ff)
            charge: Molecular charge
            multiplicity: Spin multiplicity
            optimization: Whether to perform geometry optimization
            
        Returns:
            Calculation results dict
        """
        if not self.enabled:
            logger.info("XTB calculations disabled (IAM_ENABLE_XTB=false), using stub")
            return await stub_run_xtb_calculation(xyz, method, charge, multiplicity)
        
        logger.info(f"Running real XTB calculation: {method} (enabled via IAM_ENABLE_XTB)")
        
        try:
            return await self._run_real_xtb(xyz, method, charge, multiplicity, optimization)
        except Exception as e:
            logger.error(f"Real XTB calculation failed: {e}, falling back to stub")
            return await stub_run_xtb_calculation(xyz, method, charge, multiplicity)
    
    async def _run_real_xtb(
        self,
        xyz: str,
        method: str,
        charge: int,
        multiplicity: int,
        optimization: bool
    ) -> Dict[str, Any]:
        """Run actual XTB calculation using external executable."""
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Write XYZ file
            xyz_file = temp_path / "input.xyz"
            xyz_file.write_text(xyz)
            
            # Prepare XTB command
            cmd = [
                self.xtb_executable,
                str(xyz_file),
                f"--{method}",
                f"--chrg", str(charge),
                f"--uhf", str(multiplicity - 1),
                "--json"  # Request JSON output if supported
            ]
            
            if optimization:
                cmd.append("--opt")
            
            # Run XTB calculation
            logger.info(f"Executing: {' '.join(cmd)}")
            
            process = await asyncio.create_subprocess_exec(
                *cmd,
                cwd=temp_dir,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )
            
            stdout, stderr = await process.communicate()
            
            if process.returncode != 0:
                error_msg = stderr.decode() if stderr else "Unknown XTB error"
                raise RuntimeError(f"XTB failed with return code {process.returncode}: {error_msg}")
            
            # Parse XTB output
            return self._parse_xtb_output(
                stdout.decode(),
                stderr.decode(),
                temp_path,
                method,
                charge,
                multiplicity
            )
    
    def _parse_xtb_output(
        self,
        stdout: str,
        stderr: str,
        work_dir: Path,
        method: str,
        charge: int,
        multiplicity: int
    ) -> Dict[str, Any]:
        """Parse XTB output and extract relevant data."""
        
        results = {
            "method": method,
            "charge": float(charge),
            "multiplicity": multiplicity,
            "converged": True,
            "stub": False
        }
        
        # Look for energy in output
        for line in stdout.split('\n'):
            if 'TOTAL ENERGY' in line:
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        results["energy"] = float(parts[2])
                        results["energy_units"] = "hartree"
                    except (ValueError, IndexError):
                        pass
        
        # Look for optimized geometry if available
        opt_xyz_file = work_dir / "xtbopt.xyz"
        if opt_xyz_file.exists():
            try:
                results["optimized_xyz"] = opt_xyz_file.read_text()
                results["optimization_converged"] = True
            except Exception:
                pass
        
        # Look for additional properties in stderr/stdout
        self._extract_properties(stdout, stderr, results)
        
        return results
    
    def _extract_properties(self, stdout: str, stderr: str, results: Dict[str, Any]):
        """Extract additional molecular properties from XTB output."""
        
        # Extract dipole moment
        for line in stdout.split('\n'):
            if 'molecular dipole' in line.lower():
                # Try to parse dipole moment
                parts = line.split()
                try:
                    if len(parts) >= 4:
                        dipole_mag = float(parts[-1])
                        results["dipole_moment"] = {"magnitude": dipole_mag}
                except (ValueError, IndexError):
                    pass
        
        # Extract HOMO-LUMO gap
        for line in stdout.split('\n'):
            if 'HOMO-LUMO GAP' in line or 'HL-Gap' in line:
                parts = line.split()
                try:
                    gap_value = float(parts[-2])  # Usually second to last
                    results["homo_lumo_gap"] = gap_value
                except (ValueError, IndexError):
                    pass
        
        # Count atoms from original calculation
        if "atom_count" not in results:
            atom_lines = [l for l in stdout.split('\n') if 'atoms' in l.lower()]
            if atom_lines:
                try:
                    parts = atom_lines[0].split()
                    results["atom_count"] = int(parts[0])
                except (ValueError, IndexError):
                    results["atom_count"] = 1


# Global calculator instance
calculator = XTBCalculator()


async def run_xtb_calculation_enhanced(
    xyz: str,
    method: str = "gfn2",
    charge: int = 0,
    multiplicity: int = 1
) -> Dict[str, Any]:
    """
    Enhanced XTB calculation with environment flag support.
    
    This function respects the IAM_ENABLE_XTB environment variable:
    - When True: Attempts to run real XTB calculations
    - When False (default): Uses deterministic stub implementation
    """
    return await calculator.run_calculation(xyz, method, charge, multiplicity)


async def run_xtb_optimization_enhanced(
    xyz: str,
    method: str = "gfn2",
    charge: int = 0,
    multiplicity: int = 1
) -> Dict[str, Any]:
    """
    Enhanced XTB geometry optimization with environment flag support.
    """
    return await calculator.run_calculation(xyz, method, charge, multiplicity, optimization=True)


def get_xtb_status() -> Dict[str, Any]:
    """Get current XTB integration status."""
    return {
        "enabled": calculator.enabled,
        "executable_found": calculator.xtb_executable is not None,
        "executable_path": calculator.xtb_executable,
        "environment_flag": flag("IAM_ENABLE_XTB"),
        "mode": "real" if calculator.enabled else "stub"
    }
