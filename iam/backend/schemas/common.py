def ok(data: dict) -> dict:
    return {"ok": True, "data": data, "errors": []}

from iam.backend.schemas.error import error_envelope
from typing import Optional
def fail(errors: list[str], code: int = 400, details: Optional[dict] = None) -> dict:
    # Use first error as message, rest as details
    message = errors[0] if errors else "Unknown error"
    return error_envelope(code, message, details or {"errors": errors})
