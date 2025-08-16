from fastapi import APIRouter
from pydantic import BaseModel
from iam.backend.schemas.common import ok, fail
try:
    from iam.backend.utils.persistence import save_result_json
except ImportError:
    save_result_json = None
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

class MolfileRequest(BaseModel):
    molfile: str

router = APIRouter()

@router.post("/convert/molfile")
async def convert_molfile(req: MolfileRequest):
    from fastapi.responses import JSONResponse
    molfile = req.molfile
    if not isinstance(molfile, str) or not (1 < len(molfile) < 2_000_000):
        return JSONResponse(fail(["molfile must be a non-empty string <2MB"]), status_code=400)
    try:
        mol = Chem.MolFromMolBlock(molfile, sanitize=True)
        if mol is None:
            return JSONResponse(fail(["Could not parse molfile"], code=400, details={"field": "molfile"}), status_code=400)
        smiles = Chem.MolToSmiles(mol)
        try:
            from rdkit.Chem import inchi
            inchi_str = inchi.MolToInchi(mol)
        except Exception:
            inchi_str = None
        formula = rdMolDescriptors.CalcMolFormula(mol)
        result = {"smiles": smiles, "inchi": inchi_str, "formula": formula}
        resp = ok(result)
        if save_result_json and resp["ok"]:
            save_result_json("convert", result)
            from iam.backend.utils.persistence import append_benchmark_row
            append_benchmark_row({"name": "convert", "ok": True})
        return JSONResponse(resp)
    except Exception as e:
        return JSONResponse(fail([str(e)], code=500), status_code=500)
