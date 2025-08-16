"""
FastAPI main application for IAM-2.0 backend.
Replaces the legacy Flask app from iam.backend.app.
"""

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

from backend.api import convert, calc, predict
from backend.api.legacy import router as legacy_router
from backend.api.export import router as export_router

# Create FastAPI application instance
app = FastAPI(
    title="IAM-2.0 API",
    description="Integrated Application for Molecular energetics - Version 2.0",
    version="2.0.0"
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Configure appropriately for production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount static files
app.mount("/static", StaticFiles(directory="IAM_GUI/static"), name="static")

# Include API routers
app.include_router(convert.router, prefix="/api/convert", tags=["conversion"])
app.include_router(calc.router, prefix="/api/calc", tags=["calculation"])
app.include_router(predict.router, prefix="/api/predict", tags=["prediction"])
app.include_router(export_router, prefix="/api/export", tags=["export"])

# Include legacy compatibility router (no prefix for exact path matching)
app.include_router(legacy_router, tags=["legacy"])

# Health check endpoint
@app.get("/healthz")
async def health_check():
    """Health check endpoint for monitoring."""
    return {"ok": True}

# Root endpoint
@app.get("/")
async def root():
    """Root endpoint with API information."""
    return {
        "service": "IAM-2.0 Backend",
        "version": "2.0.0",
        "status": "operational",
        "docs": "/docs"
    }

# Ketcher shell endpoint for tests
@app.get("/ketcher.html")
async def ketcher_shell():
    """Serve ketcher shell HTML directly."""
    from fastapi.responses import FileResponse
    return FileResponse("IAM_GUI/templates/ketcher_debug.html")
