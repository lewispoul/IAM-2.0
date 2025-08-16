"""
RDKit-based molecular converters for IAM2.0.
Provides SMILES→XYZ, MOL→XYZ, and SMILES→PDB conversion functionality.
"""
import logging
from typing import Dict, Optional, Any
from dataclasses import dataclass, field

try:
    from rdkit import Chem  # type: ignore
    from rdkit.Chem import AllChem, rdMolDescriptors  # type: ignore
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available - using stub implementations")
    Chem = None  # type: ignore
    AllChem = None  # type: ignore
    rdMolDescriptors = None  # type: ignore

logger = logging.getLogger(__name__)


@dataclass
class ConversionResult:
    """Result of molecular conversion operation."""
    success: bool
    xyz: Optional[str] = None
    atoms: Optional[int] = None
    error: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)


class RDKitConverter:
    """RDKit-based molecular format converter."""
    
    def __init__(self):
        self.available = RDKIT_AVAILABLE
        if not self.available:
            logger.warning("RDKit not available - converter will return stub results")
    
    def smiles_to_xyz(self, smiles: str, optimize: bool = True) -> ConversionResult:
        """
        Convert SMILES string to XYZ coordinates.
        
        Args:
            smiles: SMILES representation of molecule
            optimize: Whether to optimize geometry with force field
            
        Returns:
            ConversionResult with XYZ coordinates or error
        """
        if not self.available:
            return self._stub_smiles_to_xyz(smiles)
            
        try:
            # Parse SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return ConversionResult(
                    success=False,
                    error=f"Invalid SMILES: {smiles}"
                )
            
            # Add hydrogens
            mol = Chem.AddHs(mol)
            
            # Generate 3D coordinates
            if AllChem.EmbedMolecule(mol) != 0:
                return ConversionResult(
                    success=False,
                    error="Failed to generate 3D coordinates"
                )
            
            # Optional geometry optimization
            if optimize:
                try:
                    AllChem.MMFFOptimizeMolecule(mol)
                except Exception as e:
                    logger.warning(f"Optimization failed: {e}")
            
            # Convert to XYZ format
            xyz_content = self._mol_to_xyz(mol)
            atom_count = mol.GetNumAtoms()
            
            return ConversionResult(
                success=True,
                xyz=xyz_content,
                atoms=atom_count,
                metadata={
                    "smiles": smiles,
                    "optimized": optimize,
                    "molecular_weight": rdMolDescriptors.CalcExactMolWt(mol)
                }
            )
            
        except Exception as e:
            logger.error(f"SMILES conversion failed: {e}")
            return ConversionResult(
                success=False,
                error=str(e)
            )
    
    def molfile_to_xyz(self, molfile: str) -> ConversionResult:
        """
        Convert MOL file content to XYZ coordinates.
        
        Args:
            molfile: MOL file content as string
            
        Returns:
            ConversionResult with XYZ coordinates or error
        """
        if not self.available:
            return self._stub_molfile_to_xyz(molfile)
            
        try:
            # Parse MOL file
            mol = Chem.MolFromMolBlock(molfile)
            if mol is None:
                return ConversionResult(
                    success=False,
                    error="Invalid MOL file format"
                )
            
            # Convert to XYZ format
            xyz_content = self._mol_to_xyz(mol)
            atom_count = mol.GetNumAtoms()
            
            return ConversionResult(
                success=True,
                xyz=xyz_content,
                atoms=atom_count,
                metadata={
                    "source": "molfile",
                    "molecular_weight": rdMolDescriptors.CalcExactMolWt(mol) if mol else None
                }
            )
            
        except Exception as e:
            logger.error(f"MOL file conversion failed: {e}")
            return ConversionResult(
                success=False,
                error=str(e)
            )
    
    def _mol_to_xyz(self, mol) -> str:
        """
        Convert RDKit molecule to XYZ format string.
        
        Args:
            mol: RDKit molecule object
            
        Returns:
            XYZ format string
        """
        if not self.available:
            return "1\nstub\nC 0.0 0.0 0.0"
            
        conformer = mol.GetConformer()
        atom_count = mol.GetNumAtoms()
        
        lines = [str(atom_count)]
        lines.append("Generated by IAM2.0 RDKit converter")
        
        for i in range(atom_count):
            atom = mol.GetAtomWithIdx(i)
            pos = conformer.GetAtomPosition(i)
            symbol = atom.GetSymbol()
            
            lines.append(f"{symbol} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")
        
        return "\n".join(lines)
    
    def _mol_to_pdb(self, mol, title: str = "IAM2.0 Generated") -> str:
        """
        Convert RDKit molecule to PDB format string.
        
        Args:
            mol: RDKit molecule object
            title: PDB title/header
            
        Returns:
            PDB format string
        """
        if not self.available:
            return "HEADER    STUB PDB\nATOM      1  C   MOL A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n"
            
        conformer = mol.GetConformer()
        lines = []
        
        # PDB header
        lines.append(f"HEADER    {title}")
        lines.append("REMARK   Generated by IAM2.0 RDKit converter")
        
        # ATOM records
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            pos = conformer.GetAtomPosition(i)
            symbol = atom.GetSymbol()
            
            # PDB ATOM record format
            lines.append(
                f"ATOM  {i+1:5d}  {symbol:<3s} MOL A   1    "
                f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}  1.00  0.00           {symbol:>2s}"
            )
        
        lines.append("END")
        return "\n".join(lines)
    
    def smiles_to_pdb(self, smiles: str, optimize: bool = True, title: str = "IAM2.0 Generated") -> ConversionResult:
        """
        Convert SMILES string to PDB format.
        
        Args:
            smiles: SMILES representation of molecule
            optimize: Whether to optimize geometry with force field
            title: PDB header title
            
        Returns:
            ConversionResult with PDB content or error
        """
        if not self.available:
            return self._stub_smiles_to_pdb(smiles, title)
            
        try:
            # Parse SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return ConversionResult(
                    success=False,
                    error=f"Invalid SMILES: {smiles}"
                )
            
            # Add hydrogens
            mol = Chem.AddHs(mol)
            
            # Generate 3D coordinates
            if AllChem.EmbedMolecule(mol) != 0:
                return ConversionResult(
                    success=False,
                    error="Failed to generate 3D coordinates"
                )
            
            # Optional geometry optimization
            if optimize:
                try:
                    AllChem.MMFFOptimizeMolecule(mol)
                except Exception as e:
                    logger.warning(f"Optimization failed: {e}")
            
            # Convert to PDB format
            pdb_content = self._mol_to_pdb(mol, title)
            atom_count = mol.GetNumAtoms()
            
            return ConversionResult(
                success=True,
                xyz=pdb_content,  # Using xyz field for PDB content
                atoms=atom_count,
                metadata={
                    "format": "pdb",
                    "smiles": smiles,
                    "optimized": optimize,
                    "title": title,
                    "molecular_weight": rdMolDescriptors.CalcExactMolWt(mol)
                }
            )
            
        except Exception as e:
            logger.error(f"SMILES to PDB conversion failed: {e}")
            return ConversionResult(
                success=False,
                error=str(e)
            )
    
    def _stub_smiles_to_pdb(self, smiles: str, title: str = "IAM2.0 Generated") -> ConversionResult:
        """Stub implementation for PDB conversion when RDKit is not available."""
        atom_count = len(smiles.replace("=", "").replace("#", "")) // 2 + 1
        pdb_content = f"HEADER    {title}\nREMARK   Stub conversion of {smiles}\n"
        
        for i in range(atom_count):
            pdb_content += f"ATOM  {i+1:5d}  C   MOL A   1    {0.0:8.3f}{0.0:8.3f}{i*1.5:8.3f}  1.00  0.00           C\n"
        
        pdb_content += "END\n"
        
        return ConversionResult(
            success=True,
            xyz=pdb_content,
            atoms=atom_count,
            metadata={"format": "pdb", "stub": True, "smiles": smiles, "title": title}
        )
    
    def _stub_smiles_to_xyz(self, smiles: str) -> ConversionResult:
        """Stub implementation when RDKit is not available."""
        atom_count = len(smiles.replace("=", "").replace("#", "")) // 2 + 1  # Rough estimate
        xyz_content = f"{atom_count}\nStub conversion of {smiles}\n"
        xyz_content += "\n".join([f"C 0.0 0.0 {i*1.5:.1f}" for i in range(atom_count)])
        
        return ConversionResult(
            success=True,
            xyz=xyz_content,
            atoms=atom_count,
            metadata={"stub": True, "smiles": smiles}
        )
    
    def _stub_molfile_to_xyz(self, molfile: str) -> ConversionResult:
        """Stub implementation when RDKit is not available."""
        # Count atoms from MOL file header
        lines = molfile.strip().split('\n')
        atom_count = 1
        
        for line in lines:
            if line.strip() and not line.startswith('M  END'):
                parts = line.split()
                if len(parts) >= 4 and parts[3].isalpha():
                    atom_count += 1
                    
        xyz_content = f"{atom_count}\nStub conversion from MOL file\n"
        xyz_content += "\n".join([f"C 0.0 0.0 {i*1.5:.1f}" for i in range(atom_count)])
        
        return ConversionResult(
            success=True,
            xyz=xyz_content,
            atoms=atom_count,
            metadata={"stub": True, "source": "molfile"}
        )


# Global converter instance
converter = RDKitConverter()


def smiles_to_xyz(smiles: str, optimize: bool = True) -> ConversionResult:
    """Convert SMILES to XYZ coordinates."""
    return converter.smiles_to_xyz(smiles, optimize)


def molfile_to_xyz(molfile: str) -> ConversionResult:
    """Convert MOL file to XYZ coordinates."""
    return converter.molfile_to_xyz(molfile)


def smiles_to_pdb(smiles: str, optimize: bool = True, title: str = "IAM2.0 Generated") -> ConversionResult:
    """Convert SMILES to PDB format."""
    return converter.smiles_to_pdb(smiles, optimize, title)
