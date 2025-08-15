from fastapi import APIRouter, Request
from iam.backend.schemas.common import ok, fail
from pydantic import BaseModel

router = APIRouter()

class ComputeRequest(BaseModel):
    payload: dict | None = None

@router.post("/compute/cj")
async def compute_cj(req: ComputeRequest):
    try:
        if req.payload is None:
            return fail(["Missing payload"])
        return ok({"result": "cj stub", "input": req.payload})
    except Exception as e:
        return fail([str(e)])
