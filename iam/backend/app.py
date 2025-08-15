from fastapi import FastAPI
from fastapi.responses import FileResponse, JSONResponse
import yaml
import os
from iam.backend.routes import xtb, psi4, empirical, cj, convert, ketcher, export

app = FastAPI()

@app.get("/v1/openapi.json")
async def get_openapi_json():
	yaml_path = os.path.join(os.path.dirname(__file__), "../../openapi/iam2.v1.yaml")
	with open(yaml_path, "r") as f:
		spec = yaml.safe_load(f)
	return JSONResponse(spec)

app.include_router(xtb.router)
app.include_router(psi4.router)
app.include_router(empirical.router)
app.include_router(cj.router)
app.include_router(convert.router)
app.include_router(ketcher.router)
app.include_router(export.router)

from fastapi import FastAPI
from iam.backend.routes import xtb, psi4, empirical, cj, convert, ketcher, export

app = FastAPI()
@app.get("/v1/openapi.json")
async def get_openapi_json():
	yaml_path = os.path.join(os.path.dirname(__file__), "../../openapi/iam2.v1.yaml")
	with open(yaml_path, "r") as f:
		spec = yaml.safe_load(f)
	return JSONResponse(spec)
app.include_router(xtb.router)
app.include_router(psi4.router)
app.include_router(empirical.router)
app.include_router(cj.router)
app.include_router(convert.router)
app.include_router(ketcher.router)
app.include_router(export.router)
