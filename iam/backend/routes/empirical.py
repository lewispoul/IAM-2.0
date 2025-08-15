from fastapi import APIRouter, Request
from iam.backend.schemas.common import ok, fail
from pydantic import BaseModel

router = APIRouter()

class ComputeRequest(BaseModel):
    payload: dict | None = None

@router.post("/compute/empirical")
async def compute_empirical(req: ComputeRequest):
    try:
        if req.payload is None:
            return fail(["Missing payload"])
        return ok({"result": "empirical stub", "input": req.payload})
    except Exception as e:
        return fail([str(e)])
