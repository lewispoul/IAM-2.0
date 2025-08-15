from fastapi import APIRouter
from pydantic import BaseModel
from typing import Optional, Dict, Any
from iam.backend.schemas.common import ok, fail
try:
    from iam.backend.utils.persistence import save_result_json, append_benchmark_row
except ImportError:
    save_result_json = None
    append_benchmark_row = None
from rdkit import Chem

class MolfileRequest(BaseModel):
    molfile: str

class ToXYZRequest(BaseModel):
    molfile: Optional[str] = None
    smiles: Optional[str] = None

class RunRequest(BaseModel):
    molfile: Optional[str] = None
    smiles: Optional[str] = None
    method: str
    options: Optional[Dict[str, Any]] = None

router = APIRouter()

@router.post("/ketcher/to-smiles")
async def ketcher_to_smiles(req: MolfileRequest):
    from iam.backend.routes.convert import convert_molfile
    return await convert_molfile(req)

@router.post("/ketcher/to-xyz")
async def ketcher_to_xyz(req: ToXYZRequest):
    from fastapi.responses import JSONResponse
    # TODO: Proper 3D embedding, for now placeholder
    xyz = "TODO: 3D coordinates from RDKit"
    result = {"xyz": xyz}
    return JSONResponse(ok(result))

@router.post("/ketcher/run")
async def ketcher_run(req: RunRequest):
    from fastapi.responses import JSONResponse
    smiles = req.smiles
    if req.molfile and not smiles:
        mol = Chem.MolFromMolBlock(req.molfile, sanitize=True)
        if mol is None:
            return JSONResponse(fail(["Could not parse molfile"], code=400, details={"field": "molfile"}), status_code=400)
        smiles = Chem.MolToSmiles(mol)
    payload = {"smiles": smiles, "options": req.options or {}}
    method = req.method.lower()
    if method == "empirical":
        from iam.core.empirical.iam_empirical_predictor import predict_empirical
        data = predict_empirical(payload)
        resp = ok(data)
        if save_result_json and resp["ok"]:
            save_result_json("ketcher_run", resp)
            bench = {"name": "ketcher", "method": method, "ok": True}
            for k in ["Pcj", "Tcj", "VoD"]:
                if k in data:
                    bench[k] = data[k]
            if append_benchmark_row:
                append_benchmark_row(bench)
        return JSONResponse(resp)
    elif method == "cj":
        from iam.runners.cantera_cj import predict_cj
        data = predict_cj(stoich={"H2":2,"O2":1}, rho0=1.6)
        resp = ok(data)
        if save_result_json and resp["ok"]:
            save_result_json("ketcher_run", resp)
            bench = {"name": "ketcher", "method": method, "ok": True}
            for k in ["Pcj", "Tcj", "VoD"]:
                if k in data:
                    bench[k] = data[k]
            if append_benchmark_row:
                append_benchmark_row(bench)
        return JSONResponse(resp)
    else:
        return JSONResponse(fail([f"Method '{method}' not implemented in stub"], code=501), status_code=501)
