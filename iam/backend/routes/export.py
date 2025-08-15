from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import FileResponse
from pydantic import BaseModel
from typing import List
import os
import zipfile
import time

router = APIRouter()

class ExportZipRequest(BaseModel):
    artifacts: List[str]

@router.post("/export/zip")
async def export_zip(req: ExportZipRequest, request: Request):
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
            return {"ok": False, "data": {}, "errors": [f"Invalid path: {rel_path}"]}
        if not os.path.isfile(abs_path):
            return {"ok": False, "data": {}, "errors": [f"File not found: {rel_path}"]}
        files_to_zip.append((abs_path, rel_path))
    # Create zip
    with zipfile.ZipFile(zip_path, "w") as zipf:
        for abs_path, rel_path in files_to_zip:
            zipf.write(abs_path, arcname=rel_path)
    return {"ok": True, "data": {"zip_path": rel_zip_path}, "errors": []}
