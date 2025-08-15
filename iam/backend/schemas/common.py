def ok(data: dict) -> dict:
    return {"ok": True, "data": data, "errors": []}

def fail(errors: list[str]) -> dict:
    return {"ok": False, "data": {}, "errors": errors}
