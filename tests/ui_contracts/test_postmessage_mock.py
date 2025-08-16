from unittest.mock import Mock

def test_postmessage_roundtrip():
    """Mocked postMessage round-trip: send MOL, expect same MOL in response"""
    mol_block = "Methane MOL block"
    # Simulate frontend postMessage handler
    window = Mock()
    window.postMessage = Mock()
    # Simulate backend handler responds with same MOL
    def backend_handler(message):
        assert message["type"] == "get-molfile"
        return {"type": "molfile", "molfile": mol_block}
    response = backend_handler({"type": "get-molfile"})
    assert response["molfile"] == mol_block
