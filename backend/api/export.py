"""Export API endpoints for IAM-2.0."""

from fastapi import APIRouter, Body, Query
from fastapi.responses import StreamingResponse, JSONResponse
from typing import List, Optional
from io import StringIO, BytesIO
import csv, json, zipfile
import re
from pathlib import Path

from backend.api.envelope import ok, err
from backend.utils.persistence import list_calcs, get_calc, get_results_bulk, get_results_base

router = APIRouter()

@router.get("/export/list")
async def export_list():
    """List all available calculation IDs."""
    calc_ids = list_calcs()
    return ok({"calc_ids": calc_ids})

@router.get("/export/details")
async def export_details(id: str = Query(...)):
    """Get details for a specific calculation."""
    result = get_calc(id)
    if result is None:
        return ok({"calc_id": id, "result": None})
    return ok({"calc_id": id, "result": result})

@router.get("/calculations/")
async def list_calculations():
    """List all available calculations."""
    calc_ids = list_calcs()
    calculations = []
    
    for calc_id in calc_ids:
        calc_data = get_calc(calc_id)
        if calc_data:
            # Add calc_id to the data if not present
            if "calculation_id" not in calc_data:
                calc_data["calculation_id"] = calc_id
            calculations.append(calc_data)
    
    return ok({"calculations": calculations})

@router.get("/calculations/{calc_id}")
async def get_calculation_details(calc_id: str):
    """Get details for a specific calculation."""
    result = get_calc(calc_id)
    if result is None:
        from backend.api.envelope import fail
        return JSONResponse(content=fail(["Calculation not found"]), status_code=200)
    return ok({"calculation": result})

@router.post("/export/csv")
async def export_csv(payload: dict = Body(...)):
    """Export calculation results as CSV."""
    calc_ids = payload.get("calculation_ids", [])
    
    if not calc_ids:
        from backend.api.envelope import fail
        return JSONResponse(content=fail(["Missing or empty 'calculation_ids' field"]), status_code=422)
    
    results = get_results_bulk(calc_ids)
    
    if not results:
        missing_ids = [cid for cid in calc_ids if get_calc(cid) is None]
        from backend.api.envelope import fail
        return JSONResponse(content=fail([f"No calculations found: {missing_ids}"]), status_code=200)
    
    # Generate CSV content
    csv_buffer = StringIO()
    
    if results:
        # Get all unique keys from all results for headers
        all_keys = set(["calc_id"])
        for result in results:
            if isinstance(result.get("data"), dict):
                # Add scalar keys directly, nested objects as JSON strings
                for key, value in result["data"].items():
                    if isinstance(value, (str, int, float, bool, type(None))):
                        all_keys.add(key)
                    else:
                        all_keys.add(key)  # Will be JSON dumped
        
        headers = sorted(all_keys)
        writer = csv.DictWriter(csv_buffer, fieldnames=headers)
        writer.writeheader()
        
        for result in results:
            row = {"calc_id": result["calc_id"]}
            if isinstance(result.get("data"), dict):
                for key, value in result["data"].items():
                    if isinstance(value, (str, int, float, bool, type(None))):
                        row[key] = value
                    else:
                        row[key] = json.dumps(value)
            writer.writerow(row)
    
    csv_content = csv_buffer.getvalue()
    
    # Return CSV content in JSON response (legacy format)
    return ok({
        "csv_content": csv_content,
        "filename": "calculations.csv"
    })

@router.post("/export/json")
async def export_json(payload: dict = Body(...)):
    """Export calculation results as JSON."""
    calc_ids = payload.get("calculation_ids", [])
    
    if not calc_ids:
        from backend.api.envelope import fail
        return JSONResponse(content=fail(["Missing or empty 'calculation_ids' field"]), status_code=200)
    
    results = get_results_bulk(calc_ids)
    
    if not results:
        missing_ids = [cid for cid in calc_ids if get_calc(cid) is None]
        from backend.api.envelope import fail
        return JSONResponse(content=fail([f"No calculations found: {missing_ids}"]), status_code=200)
    
    # Transform results for export - return the actual calculation data
    items = []
    for result in results:
        # Return the raw data from the file, which should include calculation_id
        items.append(result["data"])
    
    # Return JSON content in envelope (legacy format)
    return ok({
        "json_content": json.dumps(items, indent=2),
        "filename": "calculations.json"
    })

@router.post("/export/zip")
async def export_zip(payload: dict = Body(...)):
    """Export artifacts as ZIP file."""
    artifacts = payload.get("artifacts", [])
    
    if not artifacts:
        return err("Missing or empty 'artifacts' field", 422, {"field": "artifacts"})
    
    # Validate artifact paths for security
    missing_files = []
    valid_files = []
    
    for artifact in artifacts:
        # Path traversal protection
        if '..' in artifact or artifact.startswith('/') or '\\' in artifact:
            return err("Path traversal blocked", 400, {"invalid_path": artifact})
        
        # Clean and resolve path
        results_base = get_results_base()
        artifact_path = results_base / artifact
        try:
            # Ensure the resolved path is still under results_base
            artifact_path = artifact_path.resolve()
            if not str(artifact_path).startswith(str(results_base.resolve())):
                return err("Path traversal blocked", 400, {"invalid_path": artifact})
            
            if artifact_path.exists():
                valid_files.append((artifact, artifact_path))
            else:
                missing_files.append(artifact)
        except (OSError, ValueError):
            missing_files.append(artifact)
    
    # If any files are missing, return error
    if missing_files:
        return err("Missing artifacts", 400, {"missing": missing_files})
    
    # Create ZIP file in memory
    zip_buffer = BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        for artifact_name, artifact_path in valid_files:
            # Use safe archive name (relative path)
            safe_name = artifact_name.replace('\\', '/').lstrip('/')
            if artifact_path.is_file():
                zip_file.write(artifact_path, safe_name)
            elif artifact_path.is_dir():
                # Add directory and its contents
                for file_path in artifact_path.rglob('*'):
                    if file_path.is_file():
                        rel_path = file_path.relative_to(artifact_path.parent)
                        zip_file.write(file_path, str(rel_path).replace('\\', '/'))
    
    zip_buffer.seek(0)
    
    return StreamingResponse(
        BytesIO(zip_buffer.read()),
        media_type="application/zip",
        headers={"Content-Disposition": "attachment; filename=exports.zip"}
    )
