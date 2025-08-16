"""
FastAPI conversion endpoints for IAM2.0.
Provides modern REST API for molecular format conversion.
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import Dict, Any, Optional

from ..converters.rdkit import smiles_to_xyz, molfile_to_xyz, smiles_to_pdb

router = APIRouter()


class ConvertSMILESRequest(BaseModel):
    """Request model for SMILES conversion."""
    smiles: str
    optimize: bool = True


class ConvertMOLRequest(BaseModel):
    """Request model for MOL file conversion."""
    molfile: str


class ConvertSMILESToPDBRequest(BaseModel):
    """Request model for SMILES to PDB conversion."""
    smiles: str
    optimize: bool = True
    title: Optional[str] = "IAM2.0 Generated"


class ConvertResponse(BaseModel):
    """Response model for conversion operations."""
    success: bool
    xyz: Optional[str] = None
    atoms: Optional[int] = None
    error: Optional[str] = None
    metadata: Dict[str, Any] = {}


@router.post("/smiles-to-xyz", response_model=ConvertResponse)
async def convert_smiles_to_xyz(request: ConvertSMILESRequest):
    """
    Convert SMILES string to XYZ coordinates.
    
    Args:
        request: SMILES conversion request
        
    Returns:
        ConvertResponse with XYZ coordinates or error
    """
    try:
        result = smiles_to_xyz(request.smiles, request.optimize)
        
        if not result.success:
            raise HTTPException(status_code=400, detail=result.error)
        
        return ConvertResponse(
            success=result.success,
            xyz=result.xyz,
            atoms=result.atoms,
            metadata=result.metadata
        )
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/mol-to-xyz", response_model=ConvertResponse)
async def convert_mol_to_xyz(request: ConvertMOLRequest):
    """
    Convert MOL file content to XYZ coordinates.
    
    Args:
        request: MOL file conversion request
        
    Returns:
        ConvertResponse with XYZ coordinates or error
    """
    try:
        result = molfile_to_xyz(request.molfile)
        
        if not result.success:
            raise HTTPException(status_code=400, detail=result.error)
        
        return ConvertResponse(
            success=result.success,
            xyz=result.xyz,
            atoms=result.atoms,
            metadata=result.metadata
        )
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/smiles-to-pdb", response_model=ConvertResponse)
async def convert_smiles_to_pdb(request: ConvertSMILESToPDBRequest):
    """
    Convert SMILES string to PDB format.
    
    Args:
        request: SMILES to PDB conversion request
        
    Returns:
        ConvertResponse with PDB content or error
    """
    try:
        result = smiles_to_pdb(request.smiles, request.optimize, request.title or "IAM2.0 Generated")
        
        if not result.success:
            raise HTTPException(status_code=400, detail=result.error)
        
        return ConvertResponse(
            success=result.success,
            xyz=result.xyz,  # Contains PDB content
            atoms=result.atoms,
            metadata=result.metadata
        )
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
