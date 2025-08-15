
from fastapi import FastAPI
from iam.backend.routes import xtb, psi4, empirical, cj, convert, ketcher, export

app = FastAPI()
app.include_router(xtb.router)
app.include_router(psi4.router)
app.include_router(empirical.router)
app.include_router(cj.router)
app.include_router(convert.router)
app.include_router(ketcher.router)
app.include_router(export.router)
