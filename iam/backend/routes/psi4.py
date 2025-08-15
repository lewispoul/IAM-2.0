from fastapi import APIRouter, Request
from iam.backend.schemas.common import ok, fail
from pydantic import BaseModel

router = APIRouter()

class ComputeRequest(BaseModel):
    payload: dict | None = None

@router.post("/compute/psi4")
async def compute_psi4(req: ComputeRequest):
    try:
        if req.payload is None:
            return fail(["Missing payload"])
        return ok({"result": "psi4 stub", "input": req.payload})
    except Exception as e:
        return fail([str(e)])
