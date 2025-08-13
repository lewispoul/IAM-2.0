# syntax=docker/dockerfile:1
FROM python:3.11-slim
WORKDIR /app
COPY pyproject.toml README.md /app/
COPY api /app/api
COPY src /app/src
RUN pip install --no-cache-dir --upgrade pip && pip install fastapi uvicorn pydantic
EXPOSE 8000
CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8000"]
