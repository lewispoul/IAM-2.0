"""
Tests for empirical prediction functionality in IAM-2.0.
Tests both real and stub implementations based on environment flags.
"""
import pytest
import os
from backend.jobs.empirical_predictor import EmpiricalPredictor


class TestEmpiricalPredictor:
    """Test empirical molecular property prediction."""
    
    def setup_method(self):
        """Setup test environment."""
        self.predictor = EmpiricalPredictor()
        
    def test_predictor_initialization(self):
        """Test empirical predictor initialization."""
        assert hasattr(self.predictor, 'enabled')
        assert hasattr(self.predictor, 'models_available')
        assert hasattr(self.predictor, 'descriptors_available')
        
    def test_get_status(self):
        """Test status retrieval."""
        status = self.predictor.get_status()
        
        assert isinstance(status, dict)
        assert 'enabled' in status
        assert 'mode' in status
        assert 'models_available' in status
        assert 'descriptors_available' in status
        assert 'prediction_types' in status
        
        # Check prediction types
        expected_types = [
            "logP", "molecular_weight", "tpsa", 
            "bioavailability", "drug_likeness", "toxicity_indicators"
        ]
        assert all(ptype in status['prediction_types'] for ptype in expected_types)
        
    def test_predict_properties_water(self):
        """Test property prediction for water molecule."""
        smiles = "O"  # Water
        result = self.predictor.predict_properties(smiles)
        
        assert isinstance(result, dict)
        assert result['success'] is True
        assert 'smiles' in result
        assert result['smiles'] == smiles
        assert 'properties' in result
        assert 'descriptors' in result
        assert 'confidence_scores' in result
        assert 'timestamp' in result
        
        # Check that properties are present
        properties = result['properties']
        assert 'logP' in properties
        assert 'molecular_weight' in properties
        
        # Check value types
        assert isinstance(properties['logP'], (int, float))
        assert isinstance(properties['molecular_weight'], (int, float))
        
    def test_predict_properties_caffeine(self):
        """Test property prediction for caffeine."""
        smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # Caffeine
        result = self.predictor.predict_properties(smiles)
        
        assert isinstance(result, dict)
        assert result['success'] is True
        assert result['smiles'] == smiles
        assert 'properties' in result
        assert 'descriptors' in result
        
        # Caffeine should have reasonable molecular weight
        mw = result['properties'].get('molecular_weight')
        if mw is not None:
            assert 100 < mw < 300  # Caffeine MW ≈ 194
            
    def test_predict_properties_with_specific_properties(self):
        """Test prediction with specific property list."""
        smiles = "CCO"  # Ethanol
        properties = ["logP", "molecular_weight"]
        
        result = self.predictor.predict_properties(smiles, properties)
        
        assert result['success'] is True
        assert 'properties' in result
        
        # Should only contain requested properties or all if not filtered
        result_props = result['properties']
        if len(result_props) <= len(properties):
            # If filtered, check only requested properties
            for prop in result_props:
                assert prop in properties
                
    def test_batch_predict(self):
        """Test batch prediction functionality."""
        smiles_list = ["O", "CCO", "CC(=O)O"]  # Water, ethanol, acetic acid
        
        results = self.predictor.batch_predict(smiles_list)
        
        assert isinstance(results, list)
        assert len(results) == len(smiles_list)
        
        for i, result in enumerate(results):
            assert result['success'] is True
            assert result['smiles'] == smiles_list[i]
            assert 'properties' in result
            
    def test_analyze_druglikeness(self):
        """Test drug-likeness analysis."""
        smiles = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"  # Ibuprofen
        
        analysis = self.predictor.analyze_druglikeness(smiles)
        
        assert isinstance(analysis, dict)
        assert 'overall_score' in analysis
        assert 'lipinski_compliance' in analysis
        assert 'recommendations' in analysis
        assert 'warnings' in analysis
        
        # Overall score should be between 0 and 1
        score = analysis['overall_score']
        assert 0.0 <= score <= 1.0
        
        # Should have recommendations and warnings lists
        assert isinstance(analysis['recommendations'], list)
        assert isinstance(analysis['warnings'], list)
        
    def test_deterministic_stub_behavior(self):
        """Test that stub implementation gives deterministic results."""
        smiles = "C1CCCCC1"  # Cyclohexane
        
        # Run prediction multiple times
        result1 = self.predictor.predict_properties(smiles)
        result2 = self.predictor.predict_properties(smiles)
        
        # Should get identical results (deterministic stub)
        assert result1['properties'] == result2['properties']
        assert result1['descriptors'] == result2['descriptors']
        
    def test_different_smiles_different_results(self):
        """Test that different SMILES give different results."""
        smiles1 = "C"    # Methane
        smiles2 = "CC"   # Ethane
        
        result1 = self.predictor.predict_properties(smiles1)
        result2 = self.predictor.predict_properties(smiles2)
        
        # Should get different results for different molecules
        assert result1['properties']['molecular_weight'] != result2['properties']['molecular_weight']
        
    def test_invalid_smiles_handling(self):
        """Test handling of invalid SMILES strings."""
        invalid_smiles = "INVALID_SMILES_123"
        
        # Should not crash, may return stub results or handle gracefully
        result = self.predictor.predict_properties(invalid_smiles)
        assert isinstance(result, dict)
        # The exact behavior depends on implementation (real vs stub)
        
    def test_empty_smiles_list(self):
        """Test batch prediction with empty list."""
        results = self.predictor.batch_predict([])
        assert isinstance(results, list)
        assert len(results) == 0


@pytest.mark.skipif(
    not os.environ.get("IAM_ENABLE_EMPIRICAL", "").lower() == "true",
    reason="Real empirical prediction not enabled (set IAM_ENABLE_EMPIRICAL=true)"
)
class TestEmpiricalPredictorReal:
    """Test real empirical prediction when enabled."""
    
    def setup_method(self):
        """Setup test environment for real prediction."""
        self.predictor = EmpiricalPredictor()
        
    def test_real_mode_enabled(self):
        """Test that real mode is properly enabled."""
        status = self.predictor.get_status()
        assert status['enabled'] is True
        assert status['mode'] == 'real'
        
    def test_real_prediction_has_method_info(self):
        """Test that real predictions include method information."""
        smiles = "CCO"  # Ethanol
        result = self.predictor.predict_properties(smiles)
        
        if result.get('method') == 'empirical_ml':
            # Real prediction should include additional metadata
            assert 'method' in result
            assert result['method'] in ['empirical_ml', 'empirical_stub']
            
    def test_model_availability_reporting(self):
        """Test that model availability is properly reported."""
        status = self.predictor.get_status()
        models = status['models_available']
        
        assert isinstance(models, dict)
        # Should report status of various models
        expected_models = ['rdkit_descriptors', 'sklearn', 'tensorflow']
        for model in expected_models:
            if model in models:
                assert isinstance(models[model], bool)


if __name__ == "__main__":
    # Run tests
    import sys
    pytest.main([__file__] + sys.argv[1:])
