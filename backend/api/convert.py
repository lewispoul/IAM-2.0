"""
Conversion API endpoints for IAM2.0.
Handles molecular format conversions (SMILES, MOL → XYZ).
"""
from fastapi import APIRouter, HTTPException, status
from pydantic import BaseModel, Field
from typing import Optional
import logging

from ..converters.rdkit import smiles_to_xyz, molfile_to_xyz

logger = logging.getLogger(__name__)
router = APIRouter(prefix="/api/convert", tags=["conversion"])


class SmilesRequest(BaseModel):
    """Request model for SMILES to XYZ conversion."""
    smiles: str = Field(..., description="SMILES string to convert")
    optimize: bool = Field(True, description="Apply force field optimization")


class MolfileRequest(BaseModel):
    """Request model for MOL file to XYZ conversion."""
    molfile: str = Field(..., description="MOL file content")


class ConversionResponse(BaseModel):
    """Response model for conversion operations."""
    success: bool
    xyz: Optional[str] = None
    atoms: Optional[int] = None
    error: Optional[str] = None
    metadata: dict = {}


@router.post("/smiles-to-xyz", response_model=ConversionResponse)
async def convert_smiles_to_xyz(request: SmilesRequest):
    """
    Convert SMILES string to XYZ coordinates.
    
    - **smiles**: SMILES representation of the molecule
    - **optimize**: Whether to optimize geometry (default: true)
    
    Returns XYZ coordinates with atom count and metadata.
    """
    logger.info(f"Converting SMILES to XYZ: {request.smiles[:50]}...")
    
    try:
        result = smiles_to_xyz(request.smiles, request.optimize)
        
        return ConversionResponse(
            success=result.success,
            xyz=result.xyz,
            atoms=result.atoms,
            error=result.error,
            metadata=result.metadata or {}
        )
        
    except Exception as e:
        logger.error(f"SMILES conversion failed: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Conversion failed: {str(e)}"
        )


@router.post("/molfile-to-xyz", response_model=ConversionResponse)
async def convert_molfile_to_xyz(request: MolfileRequest):
    """
    Convert MOL file content to XYZ coordinates.
    
    - **molfile**: Complete MOL file content as string
    
    Returns XYZ coordinates with atom count and metadata.
    """
    logger.info("Converting MOL file to XYZ...")
    
    try:
        result = molfile_to_xyz(request.molfile)
        
        return ConversionResponse(
            success=result.success,
            xyz=result.xyz,
            atoms=result.atoms,
            error=result.error,
            metadata=result.metadata or {}
        )
        
    except Exception as e:
        logger.error(f"MOL file conversion failed: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Conversion failed: {str(e)}"
        )


@router.get("/health")
async def conversion_health():
    """Health check for conversion services."""
    from ..converters.rdkit import converter
    
    return {
        "status": "ok",
        "rdkit_available": converter.available,
        "endpoints": [
            "/api/convert/smiles-to-xyz",
            "/api/convert/molfile-to-xyz"
        ]
    }
