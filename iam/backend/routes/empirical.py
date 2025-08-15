from fastapi import APIRouter, Request
from iam.backend.schemas.common import ok, fail
from pydantic import BaseModel

router = APIRouter()

class ComputeRequest(BaseModel):
    payload: dict | None = None

from iam.core.empirical.iam_empirical_predictor import predict_empirical

@router.post("/compute/empirical")
async def compute_empirical(req: ComputeRequest):
    try:
        if req.payload is None:
            return fail(["Missing payload"])
        result = predict_empirical(req.payload)
        return ok(result)
    except Exception as e:
        return fail([str(e)])
