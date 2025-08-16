from fastapi import APIRouter
from pydantic import BaseModel
from typing import List
import os
import zipfile
import time
from iam.backend.schemas.common import ok, fail
from iam.backend.schemas.common import ok, fail

router = APIRouter()

class ExportZipRequest(BaseModel):
    artifacts: List[str]

@router.post("/export/zip")
async def export_zip(req: ExportZipRequest):
    from fastapi.responses import JSONResponse
    base_dir = os.environ.get("IAM_RESULTS_BASE", "IAM_Knowledge")
    exports_dir = os.path.join(base_dir, "Exports")
    os.makedirs(exports_dir, exist_ok=True)
    timestamp = int(time.time())
    zip_name = f"{timestamp}_export.zip"
    zip_path = os.path.join(exports_dir, zip_name)
    rel_zip_path = os.path.relpath(zip_path, base_dir)
    files_to_zip = []
    for rel_path in req.artifacts:
        abs_path = os.path.abspath(os.path.join(base_dir, rel_path))
        # Validate path is under base_dir
        if not abs_path.startswith(os.path.abspath(base_dir)):
            return JSONResponse(fail([f"Invalid path: {rel_path}"], code=400, details={"field": rel_path}), status_code=400)
        if not os.path.isfile(abs_path):
            return JSONResponse(fail([f"File not found: {rel_path}"], code=404, details={"field": rel_path}), status_code=404)
        files_to_zip.append((abs_path, rel_path))
    # Create zip
    with zipfile.ZipFile(zip_path, "w") as zipf:
        for abs_path, rel_path in files_to_zip:
            zipf.write(abs_path, arcname=rel_path)
    return JSONResponse(ok({"zip_path": rel_zip_path}))
