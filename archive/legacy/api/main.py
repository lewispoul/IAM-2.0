from fastapi import FastAPI

app = FastAPI(title="IAM 2.0 API", version="0.0.1")

@app.get("/api/v1/health")
def health():
    return {"status": "ok"}
