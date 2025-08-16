"""
Calculation API endpoints for IAM2.0.
Handles chemistry calculations (XTB, Psi4, etc.).
"""
from fastapi import APIRouter, HTTPException, status, BackgroundTasks
from pydantic import BaseModel, Field
from typing import Optional, Dict, Any
import logging
import asyncio
from uuid import uuid4

from ..jobs.xtb_stub import run_xtb_calculation
from ..jobs.psi4_stub import run_psi4_calculation

logger = logging.getLogger(__name__)
router = APIRouter(prefix="/api/calc", tags=["calculation"])


class XYZRequest(BaseModel):
    """Request model for XYZ-based calculations."""
    xyz: str = Field(..., description="XYZ coordinate data")
    method: Optional[str] = Field("gfn2", description="Calculation method/level")
    charge: int = Field(0, description="Molecular charge")
    multiplicity: int = Field(1, description="Spin multiplicity")


class CalculationResponse(BaseModel):
    """Response model for calculation operations."""
    success: bool
    job_id: Optional[str] = None
    results: Optional[Dict[str, Any]] = None
    status: str = "completed"
    error: Optional[str] = None
    metadata: Dict[str, Any] = {}


class JobStatusResponse(BaseModel):
    """Response model for job status queries."""
    job_id: str
    status: str  # pending, running, completed, failed
    progress: Optional[float] = None
    results: Optional[Dict[str, Any]] = None
    error: Optional[str] = None


# In-memory job tracking (replace with Redis/database in production)
active_jobs: Dict[str, Dict[str, Any]] = {}


@router.post("/xtb", response_model=CalculationResponse)
async def run_xtb(request: XYZRequest, background_tasks: BackgroundTasks):
    """
    Run XTB (extended tight-binding) calculation.
    
    - **xyz**: Molecule coordinates in XYZ format
    - **method**: XTB method (gfn1, gfn2, gfn-ff)
    - **charge**: Molecular charge
    - **multiplicity**: Spin multiplicity
    
    Returns energy, orbitals, and other molecular properties.
    """
    logger.info(f"Starting XTB calculation with method: {request.method}")
    
    try:
        job_id = str(uuid4())
        
        # For now, run synchronously (stub implementation)
        result = await run_xtb_calculation(
            xyz=request.xyz,
            method=request.method,
            charge=request.charge,
            multiplicity=request.multiplicity
        )
        
        return CalculationResponse(
            success=True,
            job_id=job_id,
            results=result,
            status="completed",
            metadata={
                "method": request.method,
                "charge": request.charge,
                "multiplicity": request.multiplicity
            }
        )
        
    except Exception as e:
        logger.error(f"XTB calculation failed: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Calculation failed: {str(e)}"
        )


@router.post("/psi4", response_model=CalculationResponse)
async def run_psi4(request: XYZRequest, background_tasks: BackgroundTasks):
    """
    Run Psi4 quantum chemistry calculation.
    
    - **xyz**: Molecule coordinates in XYZ format
    - **method**: Quantum chemistry method (HF, B3LYP, etc.)
    - **charge**: Molecular charge
    - **multiplicity**: Spin multiplicity
    
    Returns energy, orbitals, and quantum properties.
    """
    logger.info(f"Starting Psi4 calculation with method: {request.method}")
    
    try:
        job_id = str(uuid4())
        
        # For now, run synchronously (stub implementation)
        result = await run_psi4_calculation(
            xyz=request.xyz,
            method=request.method or "HF",
            charge=request.charge,
            multiplicity=request.multiplicity
        )
        
        return CalculationResponse(
            success=True,
            job_id=job_id,
            results=result,
            status="completed",
            metadata={
                "method": request.method or "HF",
                "charge": request.charge,
                "multiplicity": request.multiplicity
            }
        )
        
    except Exception as e:
        logger.error(f"Psi4 calculation failed: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Calculation failed: {str(e)}"
        )


@router.get("/job/{job_id}", response_model=JobStatusResponse)
async def get_job_status(job_id: str):
    """
    Get status of a calculation job.
    
    Returns job status, progress, and results if completed.
    """
    if job_id not in active_jobs:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Job {job_id} not found"
        )
    
    job = active_jobs[job_id]
    
    return JobStatusResponse(
        job_id=job_id,
        status=job.get("status", "unknown"),
        progress=job.get("progress"),
        results=job.get("results"),
        error=job.get("error")
    )


@router.get("/health")
async def calculation_health():
    """Health check for calculation services."""
    return {
        "status": "ok",
        "active_jobs": len(active_jobs),
        "endpoints": [
            "/api/calc/xtb",
            "/api/calc/psi4",
            "/api/calc/job/{job_id}"
        ],
        "methods": {
            "xtb": ["gfn1", "gfn2", "gfn-ff"],
            "psi4": ["HF", "B3LYP", "MP2", "CCSD"]
        }
    }
