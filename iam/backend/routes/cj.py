from fastapi import APIRouter
from iam.backend.schemas.common import ok, fail
from pydantic import BaseModel

from iam.runners.cantera_cj import predict_cj
from iam.backend.utils.persistence import save_result_json, append_benchmark_row

router = APIRouter()

class CJRequest(BaseModel):
    stoich: dict
    rho0: float
    dhf: float | None = None
    metals: dict | None = None
    species_db: str | None = None

@router.post("/compute/cj")
async def compute_cj(req: CJRequest):
    try:
        result = predict_cj(
            stoich=req.stoich,
            rho0=req.rho0,
            dhf=req.dhf,
            metals=req.metals,
            species_db=req.species_db,
        )
        resp = ok(result)
        if resp["ok"]:
            save_result_json("cj", resp)
            bench = {"name": "cj", "ok": True}
            for k in ["Pcj", "Tcj", "VoD"]:
                if k in result:
                    bench[k] = result[k]
            append_benchmark_row(bench)
        return resp
    except Exception as e:
        return fail([str(e)])
