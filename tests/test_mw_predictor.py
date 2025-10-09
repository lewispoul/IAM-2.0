"""
Unit tests for molecular weight calculation in empirical predictor.
Tests both RDKit-based and fallback formula parsing implementations.
"""
import pytest
from backend.jobs.empirical_predictor import EmpiricalPredictor


class TestMolecularWeight:
    """Test molecular weight calculations."""

    def setup_method(self):
        """Setup test environment."""
        self.predictor = EmpiricalPredictor()

    def test_caffeine_mw(self):
        """Test caffeine molecular weight calculation."""
        smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # Caffeine
        result = self.predictor.predict_properties(
            smiles, ["molecular_weight"])

        assert result['success'] is True
        mw = result['properties'].get('molecular_weight')
        assert mw is not None
        assert abs(mw - 194.19) < 0.05  # ±0.05 g/mol tolerance

    def test_water_mw(self):
        """Test water molecular weight calculation."""
        smiles = "O"  # Water
        result = self.predictor.predict_properties(
            smiles, ["molecular_weight"])

        assert result['success'] is True
        mw = result['properties'].get('molecular_weight')
        assert mw is not None
        assert abs(mw - 18.015) < 0.05  # ±0.05 g/mol tolerance

    def test_co2_mw(self):
        """Test CO2 molecular weight calculation."""
        smiles = "O=C=O"  # Carbon dioxide
        result = self.predictor.predict_properties(
            smiles, ["molecular_weight"])

        assert result['success'] is True
        mw = result['properties'].get('molecular_weight')
        assert mw is not None
        assert abs(mw - 44.009) < 0.05  # ±0.05 g/mol tolerance

    def test_methane_mw(self):
        """Test methane molecular weight calculation."""
        smiles = "C"  # Methane
        result = self.predictor.predict_properties(
            smiles, ["molecular_weight"])

        assert result['success'] is True
        mw = result['properties'].get('molecular_weight')
        assert mw is not None
        assert abs(mw - 16.043) < 0.05  # ±0.05 g/mol tolerance


if __name__ == "__main__":
    pytest.main([__file__])
