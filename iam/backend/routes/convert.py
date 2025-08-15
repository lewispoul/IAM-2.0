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
    molfile = req.molfile
    if not isinstance(molfile, str) or not (1 < len(molfile) < 2_000_000):
        return fail(["molfile must be a non-empty string <2MB"])
    try:
        mol = Chem.MolFromMolBlock(molfile, sanitize=True)
        if mol is None:
            return fail(["Could not parse molfile"])
        smiles = Chem.MolToSmiles(mol)
        try:
            from rdkit.Chem import inchi
            inchi_str = inchi.MolToInchi(mol)
        except Exception:
            inchi_str = None
        formula = rdMolDescriptors.CalcMolFormula(mol)
        result = {"smiles": smiles, "inchi": inchi_str, "formula": formula}
        if save_result_json:
            save_result_json("convert", result)
        return ok(result)
    except Exception as e:
        return fail([str(e)])
