import pytest
from backend.empirical.iam_empirical_predictor import predict_kamlet_jacobs, predict_keshavarz, predict_empirical

@pytest.mark.parametrize("func,payload,expected_keys", [
    (predict_kamlet_jacobs, {"foo": 1}, {"Pcj", "Tcj", "VoD", "input"}),
    (predict_keshavarz, {"bar": 2}, {"Pcj", "Tcj", "VoD", "input"}),
    (predict_empirical, {"method": "kj"}, {"Pcj", "Tcj", "VoD", "input"}),
    (predict_empirical, {"method": "keshavarz"}, {"Pcj", "Tcj", "VoD", "input"}),
])
def test_empirical_predictors(func, payload, expected_keys):
    result = func(payload)
    assert set(result.keys()) == expected_keys
    assert isinstance(result["Pcj"], float)
    assert isinstance(result["Tcj"], float)
    assert isinstance(result["VoD"], float)
    assert isinstance(result["input"], dict)
