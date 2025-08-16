from typing import Dict, Any

# Deterministic placeholder values
KJ_PLACEHOLDER = {"Pcj": 25.0, "Tcj": 3500.0, "VoD": 7.5}
KESH_PLACEHOLDER = {"Pcj": 22.0, "Tcj": 3300.0, "VoD": 6.8}


def predict_kamlet_jacobs(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Stub for Kamlet–Jacobs empirical prediction."""
    return {**KJ_PLACEHOLDER, "input": payload}


def predict_keshavarz(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Stub for Keshavarz empirical prediction."""
    return {**KESH_PLACEHOLDER, "input": payload}


def predict_empirical(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Selects method based on 'method' field in payload."""
    method = payload.get("method", "kj").lower()
    if method == "keshavarz":
        return predict_keshavarz(payload)
    return predict_kamlet_jacobs(payload)
