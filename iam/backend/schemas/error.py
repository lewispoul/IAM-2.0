from typing import Any, Optional
from uuid import uuid4

def error_envelope(code: int, message: str, details: Optional[Any] = None, correlation_id: Optional[str] = None) -> dict:
    """
    Returns a normalized error envelope for API responses.
    """
    if correlation_id is None:
        correlation_id = str(uuid4())
    return {
        "code": code,
        "message": message,
        "details": details or {},
        "correlation_id": correlation_id
    }

# Example usage:
# return error_envelope(400, "Invalid input", {"field": "molfile"})
