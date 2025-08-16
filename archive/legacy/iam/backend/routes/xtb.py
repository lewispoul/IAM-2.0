from fastapi import APIRouter, Request
from iam.backend.schemas.common import ok, fail
from pydantic import BaseModel

router = APIRouter()

class ComputeRequest(BaseModel):
    payload: dict | None = None

@router.post("/compute/xtb")
async def compute_xtb(req: ComputeRequest):
    from fastapi.responses import JSONResponse
    try:
        # Minimal input check
        if req.payload is None:
            return JSONResponse(fail(["Missing payload"], code=400), status_code=400)
        return JSONResponse(ok({"result": "xtb stub", "input": req.payload}))
    except Exception as e:
        return JSONResponse(fail([str(e)], code=500), status_code=500)
