from fastapi import FastAPI
from fastapi.responses import FileResponse, JSONResponse, HTMLResponse
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
import yaml
import os
from iam.backend.routes import xtb, psi4, empirical, cj, convert, ketcher, export

app = FastAPI(
    title="IAM-2.0 API",
    description="Integrated Analysis of Materials 2.0 - Molecular Analysis API",
    version="1.0.0"
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, specify actual origins
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount static files for local Ketcher hosting
static_dir = os.path.join(os.path.dirname(__file__), "../../public/static")
if os.path.exists(static_dir):
    app.mount("/static", StaticFiles(directory=static_dir), name="static")

@app.get("/")
async def serve_ui():
    """Serve the Ketcher UI for easy access"""
    ui_file = os.path.join(os.path.dirname(__file__), "../../frontend/ketcher_test.html")
    if os.path.exists(ui_file):
        return FileResponse(ui_file, media_type="text/html")
    else:
        return HTMLResponse("""
        <html>
            <body>
                <h1>🧪 IAM-2.0 API Server</h1>
                <p>API is running! Access the UI at:</p>
                <p><code>file:///home/lppoulin/IAM-2.0/frontend/ketcher_test.html</code></p>
                <p>Or view API docs at: <a href="/docs">/docs</a></p>
            </body>
        </html>
        """)

@app.get("/health")
async def health_check():
    """Health check endpoint"""
    return {"status": "healthy", "service": "iam-2.0-api"}

@app.get("/v1/openapi.json")
async def get_openapi_json():
    yaml_path = os.path.join(os.path.dirname(__file__), "../../openapi/iam2.v1.yaml")
    with open(yaml_path, "r") as f:
        spec = yaml.safe_load(f)
    return JSONResponse(spec)

# Include all route modules
app.include_router(xtb.router)
app.include_router(psi4.router)
app.include_router(empirical.router)
app.include_router(cj.router)
app.include_router(convert.router)
app.include_router(ketcher.router)
app.include_router(export.router)
