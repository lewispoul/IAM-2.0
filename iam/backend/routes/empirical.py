from fastapi import APIRouter, Request
from iam.backend.schemas.common import ok, fail
from pydantic import BaseModel

router = APIRouter()

class ComputeRequest(BaseModel):
    payload: dict | None = None


from iam.core.empirical.iam_empirical_predictor import predict_empirical
from iam.backend.utils.persistence import save_result_json, append_benchmark_row

@router.post("/compute/empirical")
async def compute_empirical(req: ComputeRequest):
    from fastapi.responses import JSONResponse
    try:
        if req.payload is None:
            return JSONResponse(fail(["Missing payload"], code=400), status_code=400)
        result = predict_empirical(req.payload)
        resp = ok(result)
        # Persistence only if ok
        if resp["ok"]:
            save_result_json("empirical", resp)
            bench = {"name": "empirical", "ok": True}
            for k in ["method", "Pcj", "Tcj", "VoD"]:
                if k in result:
                    bench[k] = result[k]
            append_benchmark_row(bench)
        return JSONResponse(resp)
    except Exception as e:
        return JSONResponse(fail([str(e)], code=500), status_code=500)
