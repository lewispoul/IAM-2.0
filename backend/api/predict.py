"""
Performance prediction API endpoints for IAM2.0.
Handles CJ detonation velocity and pressure calculations.
"""
from fastapi import APIRouter, HTTPException, status
from pydantic import BaseModel, Field, validator
from typing import Optional, Dict, Any, List
import logging
import hashlib
import math

logger = logging.getLogger(__name__)
router = APIRouter(prefix="/api/predict", tags=["prediction"])


class CJRequest(BaseModel):
    """Request model for CJ detonation calculations."""
    formula: str = Field(..., description="Chemical formula (e.g., C2H4N4O4)")
    density: float = Field(..., description="Density in g/cm³", gt=0, le=10.0)
    temperature: float = Field(298.15, description="Temperature in Kelvin", gt=0, le=5000)
    pressure: float = Field(1.0, description="Pressure in atm", gt=0, le=1000)
    
    @validator('formula')
    def validate_formula(cls, v):
        """Basic validation of chemical formula format."""
        if not v or len(v.strip()) == 0:
            raise ValueError("Formula cannot be empty")
        # Allow letters, numbers, parentheses, brackets
        import re
        if not re.match(r'^[A-Za-z0-9()[\]\.]+$', v.strip()):
            raise ValueError("Invalid characters in formula")
        return v.strip()


class CJResponse(BaseModel):
    """Response model for CJ detonation calculations."""
    success: bool
    formula: str
    density: float
    pcj: Optional[float] = Field(None, description="CJ pressure in GPa")
    vod: Optional[float] = Field(None, description="Detonation velocity in m/s")
    temperature_cj: Optional[float] = Field(None, description="CJ temperature in K")
    gamma: Optional[float] = Field(None, description="Heat capacity ratio")
    error: Optional[str] = None
    metadata: Dict[str, Any] = {}


def parse_formula(formula: str) -> Dict[str, int]:
    """
    Parse chemical formula into element counts.
    
    Simplified parser for basic formulas like C2H4N4O4.
    """
    import re
    elements = {}
    
    # Simple regex pattern for element followed by optional number
    pattern = r'([A-Z][a-z]?)(\d*)'
    matches = re.findall(pattern, formula)
    
    for element, count in matches:
        count = int(count) if count else 1
        elements[element] = elements.get(element, 0) + count
    
    return elements


def calculate_molecular_weight(elements: Dict[str, int]) -> float:
    """Calculate molecular weight from element composition."""
    # Atomic weights (simplified set)
    atomic_weights = {
        'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.811,
        'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
        'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.086, 'P': 30.974,
        'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078
    }
    
    molecular_weight = 0.0
    for element, count in elements.items():
        if element in atomic_weights:
            molecular_weight += atomic_weights[element] * count
        else:
            # Default weight for unknown elements
            molecular_weight += 50.0 * count
    
    return molecular_weight


async def calculate_cj_stub(
    formula: str,
    density: float,
    temperature: float = 298.15,
    pressure: float = 1.0
) -> Dict[str, Any]:
    """
    Stub implementation of CJ detonation calculation.
    
    Uses empirical correlations and deterministic results for development.
    In production, this would interface with Cantera or other thermochemistry libraries.
    """
    
    # Parse formula for composition
    try:
        elements = parse_formula(formula)
        molecular_weight = calculate_molecular_weight(elements)
    except Exception as e:
        raise ValueError(f"Invalid formula: {e}")
    
    # Generate deterministic results based on input hash
    input_hash = hashlib.md5(f"{formula}{density}{temperature}{pressure}".encode()).hexdigest()
    hash_int = int(input_hash[:8], 16)
    
    # Element-based scaling factors
    c_atoms = elements.get('C', 0)
    h_atoms = elements.get('H', 0) 
    n_atoms = elements.get('N', 0)
    o_atoms = elements.get('O', 0)
    
    # Oxygen balance calculation (simplified)
    oxygen_balance = (o_atoms - 2*c_atoms - 0.5*h_atoms) * 15.999 / molecular_weight * 100
    
    # Empirical correlations (simplified Kamlet-Jacobs equations)
    # These are rough approximations for demonstration
    
    # Detonation velocity (m/s)
    # VoD ≈ A * (N*M*Q)^0.5 * (ρ/ρ0)^β
    
    # Nitrogen content factor
    nitrogen_factor = (n_atoms * 14.007) / molecular_weight
    
    # Base velocity from composition
    base_vod = 1800 + nitrogen_factor * 3000  # Base velocity
    base_vod += abs(oxygen_balance) * 5  # Oxygen balance contribution
    
    # Density scaling (typical β ≈ 0.5)
    density_factor = (density / 1.0) ** 0.5
    vod = base_vod * density_factor
    
    # Add deterministic variation
    vod += (hash_int % 800 - 400)  # ±400 m/s variation
    vod = max(1000, min(10000, vod))  # Reasonable bounds
    
    # CJ Pressure (GPa) - empirical correlation
    # Pcj ≈ ρ * VoD^2 / (γ+1) / 10^9
    gamma = 2.5 + (hash_int % 100) / 1000.0  # 2.5-2.6 typical range
    pcj = density * vod**2 / (gamma + 1) / 1e9
    
    # CJ Temperature (K) - rough estimate
    # Higher for more energetic compositions
    temp_cj = 2500 + nitrogen_factor * 1000 + (hash_int % 500)
    temp_cj = max(2000, min(5000, temp_cj))
    
    return {
        "pcj": round(pcj, 2),
        "vod": round(vod, 0),
        "temperature_cj": round(temp_cj, 0),
        "gamma": round(gamma, 3),
        "composition": elements,
        "molecular_weight": round(molecular_weight, 3),
        "oxygen_balance": round(oxygen_balance, 1),
        "nitrogen_content": round(nitrogen_factor * 100, 1),
        "density": density,
        "stub": True,
        "hash_seed": input_hash[:8]
    }


@router.post("/cj/v1", response_model=CJResponse)
async def predict_cj_detonation(request: CJRequest):
    """
    Calculate Chapman-Jouguet (CJ) detonation properties.
    
    Uses thermochemical calculations to predict:
    - CJ pressure (GPa)
    - Detonation velocity (m/s) 
    - CJ temperature (K)
    - Heat capacity ratio
    
    Currently uses stub implementation with empirical correlations.
    Production version will integrate with Cantera thermochemistry solver.
    """
    logger.info(f"CJ calculation for {request.formula} at {request.density} g/cm³")
    
    try:
        results = await calculate_cj_stub(
            request.formula,
            request.density,
            request.temperature,
            request.pressure
        )
        
        return CJResponse(
            success=True,
            formula=request.formula,
            density=request.density,
            pcj=results["pcj"],
            vod=results["vod"],
            temperature_cj=results["temperature_cj"],
            gamma=results["gamma"],
            metadata=results
        )
        
    except ValueError as e:
        logger.warning(f"Invalid input for CJ calculation: {e}")
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )
    except Exception as e:
        logger.error(f"CJ calculation failed: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Calculation failed: {str(e)}"
        )


@router.get("/methods")
async def get_prediction_methods():
    """Get available prediction methods and their descriptions."""
    return {
        "methods": {
            "cj/v1": {
                "name": "Chapman-Jouguet Detonation",
                "description": "Predicts detonation velocity and pressure using thermochemical calculations",
                "inputs": ["formula", "density", "temperature", "pressure"],
                "outputs": ["pcj", "vod", "temperature_cj", "gamma"],
                "status": "stub",
                "version": "1.0.0"
            }
        },
        "planned_methods": [
            "equilibrium_products",
            "ignition_delay",
            "flame_speed",
            "shock_hugoniot"
        ]
    }


@router.get("/health")
async def prediction_health():
    """Health check for prediction services."""
    return {
        "status": "ok",
        "methods": ["cj/v1"],
        "backend": "stub",
        "cantera_available": False,  # Will be True when integrated
        "endpoints": [
            "/api/predict/cj/v1",
            "/api/predict/methods"
        ]
    }
