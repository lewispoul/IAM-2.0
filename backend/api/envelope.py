"""
Response envelope helpers for IAM2.0 API compatibility.
Provides standardized success/error response formats expected by tests.
"""
import uuid
from typing import Any, Dict, List, Optional, Union
from fastapi import HTTPException
from fastapi.responses import JSONResponse


def generate_correlation_id() -> str:
    """Generate a unique correlation ID for request tracking."""
    return str(uuid.uuid4())


def ok(data: Dict[str, Any], errors: Optional[List[str]] = None) -> Dict[str, Any]:
    """
    Create a standardized success envelope.
    
    Args:
        data: The response data
        errors: Optional list of non-critical errors/warnings
        
    Returns:
        Success envelope with format: {"ok": True, "data": {...}, "errors": [...]}
    """
    return {
        "ok": True,
        "data": data,
        "errors": errors or []
    }


def fail(errors: List[str], data: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Create a standardized failure envelope with 200 status.
    
    Args:
        errors: List of error messages
        data: Optional partial data
        
    Returns:
        Failure envelope with format: {"ok": False, "data": {...}, "errors": [...]}
    """
    return {
        "ok": False,
        "data": data or {},
        "errors": errors
    }


def err(message: str, code: int, details: Optional[Dict[str, Any]] = None) -> JSONResponse:
    """
    Create a standardized error response with proper HTTP status.
    
    Args:
        message: Main error message
        code: HTTP status code
        details: Optional additional error details
        
    Returns:
        JSONResponse with error envelope format: {"code": int, "message": str, "details": {...}, "correlation_id": str}
    """
    response_data = {
        "code": code,
        "message": message,
        "details": details or {},
        "correlation_id": generate_correlation_id()
    }
    
    return JSONResponse(
        status_code=code,
        content=response_data
    )


def validation_error(message: str, field: str, value: Any = None) -> JSONResponse:
    """
    Create a 422 validation error response.
    
    Args:
        message: Validation error message
        field: Field that failed validation
        value: Optional invalid value
        
    Returns:
        422 JSONResponse with validation error envelope
    """
    details = {"field": field}
    if value is not None:
        details["invalid_value"] = str(value)
    
    return err(f"Validation error: {message}", 422, details)


def not_found_error(resource: str, identifier: Optional[str] = None) -> JSONResponse:
    """
    Create a 404 not found error response.
    
    Args:
        resource: Type of resource that was not found
        identifier: Optional identifier of the missing resource
        
    Returns:
        404 JSONResponse with not found error envelope
    """
    message = f"{resource} not found"
    if identifier:
        message += f": {identifier}"
        
    details = {"resource": resource}
    if identifier:
        details["identifier"] = identifier
    
    return err(message, 404, details)


def bad_request_error(message: str, details: Optional[Dict[str, Any]] = None) -> JSONResponse:
    """
    Create a 400 bad request error response.
    
    Args:
        message: Error message describing what's wrong with the request
        details: Optional additional details about the error
        
    Returns:
        400 JSONResponse with bad request error envelope
    """
    return err(f"Bad request: {message}", 400, details)


def not_implemented_error(feature: str) -> JSONResponse:
    """
    Create a 501 not implemented error response.
    
    Args:
        feature: Feature or method that is not implemented
        
    Returns:
        501 JSONResponse with not implemented error envelope
    """
    return err("not implemented", 501, {"feature": feature})


def internal_error(message: str, error_id: Optional[str] = None) -> JSONResponse:
    """
    Create a 500 internal server error response.
    
    Args:
        message: Error message
        error_id: Optional error tracking ID
        
    Returns:
        500 JSONResponse with internal error envelope
    """
    details = {"error_id": error_id or generate_correlation_id()}
    return err(f"Internal server error: {message}", 500, details)


def path_traversal_error() -> JSONResponse:
    """
    Create a 400 error for path traversal attempts.
    
    Returns:
        400 JSONResponse for security violation
    """
    return bad_request_error(
        "Invalid path: path traversal not allowed",
        {"security_violation": "path_traversal"}
    )
