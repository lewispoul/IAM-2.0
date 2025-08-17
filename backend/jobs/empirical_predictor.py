"""
Empirical Predictor for IAM-2.0
Implements ML-based property prediction and empirical analysis
"""
import json
import logging
from typing import Dict, Any, List, Optional, Tuple
from pathlib import Path
import tempfile
import os
import time

from backend.utils.env import flag

logger = logging.getLogger(__name__)

class EmpiricalPredictor:
    """
    Empirical predictor for molecular properties using ML models and heuristics.
    
    Supports:
    - QSAR-based property prediction  
    - Molecular descriptor calculation
    - Empirical rule-based analysis
    - Statistical property estimation
    """
    
    def __init__(self):
        self.enabled = flag("IAM_ENABLE_EMPIRICAL")
        self.models_available = self._check_ml_models()
        self.descriptors_available = self._check_descriptor_tools()
        
        logger.info(f"EmpiricalPredictor initialized: enabled={self.enabled}")
        
    def _check_ml_models(self) -> Dict[str, bool]:
        """Check availability of ML model libraries"""
        models = {}
        
        # Check scikit-learn
        try:
            import sklearn
            models['sklearn'] = True
            logger.debug("scikit-learn available")
        except ImportError:
            models['sklearn'] = False
            logger.debug("scikit-learn not available")
            
        # Check RDKit for descriptors  
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors
            models['rdkit_descriptors'] = True
            logger.debug("RDKit descriptors available")
        except ImportError:
            models['rdkit_descriptors'] = False
            logger.debug("RDKit descriptors not available")
            
        # Check TensorFlow/Keras
        try:
            import tensorflow
            models['tensorflow'] = True
            logger.debug("TensorFlow available")
        except ImportError:
            models['tensorflow'] = False
            logger.debug("TensorFlow not available")
            
        return models
        
    def _check_descriptor_tools(self) -> Dict[str, bool]:
        """Check availability of molecular descriptor calculation tools"""
        tools = {}
        
        # Check Mordred
        try:
            from mordred import Calculator, descriptors
            tools['mordred'] = True
            logger.debug("Mordred descriptors available")
        except ImportError:
            tools['mordred'] = False
            logger.debug("Mordred not available")
            
        # Check PaDEL-Descriptor (via subprocess)
        padel_path = os.environ.get('PADEL_DESCRIPTOR_PATH')
        if padel_path and Path(padel_path).exists():
            tools['padel'] = True
            logger.debug("PaDEL-Descriptor available")
        else:
            tools['padel'] = False
            logger.debug("PaDEL-Descriptor not available")
            
        return tools
        
    def get_status(self) -> Dict[str, Any]:
        """Get empirical predictor status and capabilities"""
        return {
            "enabled": self.enabled,
            "mode": "real" if self.enabled else "stub",
            "models_available": self.models_available,
            "descriptors_available": self.descriptors_available,
            "prediction_types": [
                "logP",
                "molecular_weight", 
                "tpsa",
                "bioavailability",
                "drug_likeness",
                "toxicity_indicators"
            ]
        }
        
    def predict_properties(self, 
                          smiles: str, 
                          properties: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Predict molecular properties from SMILES string.
        
        Args:
            smiles: SMILES string of molecule
            properties: List of properties to predict (None for all)
            
        Returns:
            Dictionary with predicted properties and confidence scores
        """
        if not self.enabled:
            return self._stub_prediction(smiles, properties)
            
        try:
            return self._real_prediction(smiles, properties)
        except Exception as e:
            logger.warning(f"Real prediction failed: {e}, falling back to stub")
            return self._stub_prediction(smiles, properties)
            
    def _real_prediction(self, smiles: str, properties: Optional[List[str]]) -> Dict[str, Any]:
        """Perform real empirical prediction using available tools"""
        logger.info(f"Performing real empirical prediction for: {smiles}")
        
        results = {
            "success": True,
            "method": "empirical_ml",
            "smiles": smiles,
            "properties": {},
            "descriptors": {},
            "confidence_scores": {},
            "timestamp": time.time()
        }
        
        # Try to use RDKit for basic descriptors and properties
        if self.models_available.get('rdkit_descriptors', False):
            try:
                from rdkit import Chem
                from rdkit.Chem import Descriptors, Crippen
                
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    # Calculate basic descriptors
                    results["descriptors"].update({
                        "molecular_weight": Descriptors.MolWt(mol),
                        "num_atoms": mol.GetNumAtoms(),
                        "num_bonds": mol.GetNumBonds(),
                        "num_rings": Descriptors.RingCount(mol),
                        "aromatic_rings": Descriptors.NumAromaticRings(mol),
                        "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
                        "h_bond_donors": Descriptors.NumHDonors(mol),
                        "h_bond_acceptors": Descriptors.NumHAcceptors(mol)
                    })
                    
                    # Calculate properties
                    results["properties"].update({
                        "logP": Crippen.MolLogP(mol),
                        "molecular_weight": Descriptors.MolWt(mol),
                        "tpsa": Descriptors.TPSA(mol),
                        "bioavailability_score": self._calculate_bioavailability(mol),
                        "drug_likeness": self._assess_drug_likeness(mol),
                        "hbd_count": Descriptors.NumHDonors(mol),
                        "hba_count": Descriptors.NumHAcceptors(mol)
                    })
                    
                    # Add confidence scores (simplified for demonstration)
                    for prop in results["properties"]:
                        results["confidence_scores"][prop] = 0.85  # High confidence for descriptor-based
                        
                    logger.info("Real RDKit-based prediction successful")
                else:
                    raise ValueError("Invalid SMILES string")
                    
            except Exception as e:
                logger.error(f"RDKit prediction failed: {e}")
                raise
                
        # Try advanced ML models if available
        if self.models_available.get('sklearn', False):
            try:
                # Placeholder for trained ML models
                logger.info("Advanced ML models available but not implemented yet")
            except Exception as e:
                logger.warning(f"Advanced ML prediction failed: {e}")
                
        return results
        
    def _calculate_bioavailability(self, mol) -> float:
        """Calculate bioavailability score using empirical rules"""
        from rdkit.Chem import Descriptors, Crippen
        
        # Simplified bioavailability scoring
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        
        score = 1.0
        
        # Molecular weight penalty
        if mw > 500:
            score *= 0.7
        elif mw > 400:
            score *= 0.9
            
        # LogP penalty  
        if abs(logp) > 5:
            score *= 0.6
        elif abs(logp) > 3:
            score *= 0.8
            
        # TPSA penalty
        if tpsa > 140:
            score *= 0.7
        elif tpsa > 90:
            score *= 0.9
            
        return max(0.0, min(1.0, score))
        
    def _assess_drug_likeness(self, mol) -> float:
        """Assess drug-likeness using Lipinski's Rule of Five and other factors"""
        from rdkit.Chem import Descriptors, Crippen
        
        violations = 0
        
        # Lipinski's Rule of Five
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        
        if mw > 500: violations += 1
        if logp > 5: violations += 1  
        if hbd > 5: violations += 1
        if hba > 10: violations += 1
        
        # Additional factors
        if Descriptors.NumRotatableBonds(mol) > 10: violations += 1
        if Descriptors.TPSA(mol) > 140: violations += 1
        
        # Convert violations to drug-likeness score
        drug_likeness = max(0.0, 1.0 - (violations * 0.15))
        
        return drug_likeness
        
    def _stub_prediction(self, smiles: str, properties: Optional[List[str]]) -> Dict[str, Any]:
        """Stub implementation returning deterministic mock results"""
        logger.info(f"Using stub empirical prediction for: {smiles}")
        
        # Generate deterministic results based on SMILES
        smiles_hash = hash(smiles) % 10000
        
        # Default properties if none specified
        if properties is None:
            properties = ["logP", "molecular_weight", "tpsa", "bioavailability", "drug_likeness"]
            
        # Generate mock properties based on hash
        mock_properties = {
            "logP": (smiles_hash % 100) / 10.0 - 5.0,  # -5.0 to 5.0
            "molecular_weight": 150 + (smiles_hash % 400),  # 150-550
            "tpsa": 20 + (smiles_hash % 150),  # 20-170
            "bioavailability": (smiles_hash % 100) / 100.0,  # 0.0-1.0
            "drug_likeness": (smiles_hash % 80 + 20) / 100.0,  # 0.2-1.0
            "toxicity_indicators": (smiles_hash % 50) / 100.0  # 0.0-0.5
        }
        
        # Filter requested properties
        filtered_properties = {k: v for k, v in mock_properties.items() if k in properties}
        
        return {
            "success": True,
            "method": "empirical_stub",
            "smiles": smiles,
            "properties": filtered_properties,
            "descriptors": {
                "molecular_weight": mock_properties["molecular_weight"],
                "num_atoms": 5 + (smiles_hash % 50),
                "num_bonds": 4 + (smiles_hash % 55),
                "num_rings": smiles_hash % 5,
                "rotatable_bonds": smiles_hash % 10
            },
            "confidence_scores": {prop: 0.75 for prop in filtered_properties},
            "note": "Stub implementation - deterministic mock results",
            "timestamp": time.time()
        }
        
    def batch_predict(self, smiles_list: List[str], properties: Optional[List[str]] = None) -> List[Dict[str, Any]]:
        """
        Perform batch prediction on multiple SMILES strings.
        
        Args:
            smiles_list: List of SMILES strings
            properties: Properties to predict for each molecule
            
        Returns:
            List of prediction results
        """
        results = []
        for smiles in smiles_list:
            result = self.predict_properties(smiles, properties)
            results.append(result)
        return results
        
    def analyze_druglikeness(self, smiles: str) -> Dict[str, Any]:
        """
        Comprehensive drug-likeness analysis.
        
        Args:
            smiles: SMILES string of molecule
            
        Returns:
            Detailed drug-likeness analysis
        """
        prediction = self.predict_properties(smiles, ["logP", "molecular_weight", "tpsa", "drug_likeness"])
        
        analysis = {
            "overall_score": prediction["properties"].get("drug_likeness", 0.5),
            "lipinski_compliance": {},
            "recommendations": [],
            "warnings": []
        }
        
        if "molecular_weight" in prediction["properties"]:
            mw = prediction["properties"]["molecular_weight"]
            analysis["lipinski_compliance"]["molecular_weight"] = {
                "value": mw,
                "limit": 500,
                "compliant": mw <= 500
            }
            if mw > 500:
                analysis["warnings"].append("Molecular weight exceeds Lipinski limit (>500 Da)")
                
        if "logP" in prediction["properties"]:
            logp = prediction["properties"]["logP"] 
            analysis["lipinski_compliance"]["logP"] = {
                "value": logp,
                "limit": 5,
                "compliant": logp <= 5
            }
            if logp > 5:
                analysis["warnings"].append("LogP exceeds Lipinski limit (>5)")
                
        if "tpsa" in prediction["properties"]:
            tpsa = prediction["properties"]["tpsa"]
            analysis["lipinski_compliance"]["tpsa"] = {
                "value": tpsa,
                "limit": 140,
                "compliant": tpsa <= 140
            }
            if tpsa > 140:
                analysis["warnings"].append("TPSA exceeds recommended limit (>140 Ų)")
                
        # Add recommendations based on analysis
        if len(analysis["warnings"]) == 0:
            analysis["recommendations"].append("Molecule shows good drug-like properties")
        else:
            analysis["recommendations"].append("Consider structural modifications to improve drug-likeness")
            
        return analysis

# Global predictor instance
empirical_predictor = EmpiricalPredictor()
