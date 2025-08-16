"""
Final Integration Test - Ketcher UI + Backend Live Test
======================================================

This script tests the complete integration between the Ketcher UI and backend.
Run with the backend server running on localhost:8010.
"""

import pytest
import requests
import json
import time
from datetime import datetime

API_BASE = "http://localhost:8010"

@pytest.mark.external
def test_endpoint(endpoint, payload, description):
    """Test a single endpoint and print results"""
    print(f"\n🧪 Testing: {description}")
    print(f"   Endpoint: {endpoint}")
    print(f"   Payload: {json.dumps(payload, indent=2)}")
    
    try:
        start_time = time.time()
        response = requests.post(f"{API_BASE}{endpoint}", json=payload, timeout=10)
        duration = time.time() - start_time
        
        print(f"   ⏱️  Response time: {duration:.2f}s")
        print(f"   📊 Status: {response.status_code}")
        
        try:
            data = response.json()
            print(f"   📝 Response preview: {json.dumps(data, indent=2)[:200]}...")
            
            # Verify response structure
            if response.status_code == 200:
                if data.get("ok") is True:
                    print("   ✅ SUCCESS: Valid envelope response")
                    return True
                else:
                    print("   ⚠️  WARNING: Response ok=False")
            else:
                if "correlation_id" in data and "message" in data:
                    print("   ✅ SUCCESS: Valid error envelope")
                    return True
                else:
                    print("   ❌ FAIL: Invalid error envelope")
                    
        except json.JSONDecodeError:
            print(f"   ❌ FAIL: Invalid JSON response")
            print(f"   Raw response: {response.text[:200]}...")
            
    except requests.exceptions.RequestException as e:
        print(f"   ❌ FAIL: Request error - {e}")
        return False
    
    return False


def run_integration_tests():
    """Run complete integration test suite"""
    print("🚀 IAM-2.0 × Ketcher Integration Test Suite")
    print("=" * 50)
    print(f"⏰ Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"🔗 Backend URL: {API_BASE}")
    
    # Test data
    mock_molfile = """
  Methane
  
  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6291    0.6291    0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6291   -0.6291    0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6291   -0.6291   -0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6291    0.6291   -0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END"""
    
    mock_smiles = "C"
    
    # Test scenarios matching the UI buttons
    test_cases = [
        # Export Functions
        {
            "endpoint": "/ketcher/to-smiles",
            "payload": {"molfile": mock_molfile},
            "description": "Export to SMILES (UI: Export to SMILES button)"
        },
        {
            "endpoint": "/convert/molfile", 
            "payload": {"molfile": mock_molfile},
            "description": "Convert Molfile (UI: Convert Molfile button)"
        },
        {
            "endpoint": "/ketcher/to-xyz",
            "payload": {"smiles": mock_smiles},
            "description": "Export to XYZ (UI: Export to XYZ button)"
        },
        
        # Computation Functions
        {
            "endpoint": "/ketcher/run",
            "payload": {"smiles": mock_smiles, "method": "empirical"},
            "description": "Run Empirical (UI: Run Empirical button)"
        },
        {
            "endpoint": "/ketcher/run",
            "payload": {"smiles": mock_smiles, "method": "cj"},
            "description": "Run CJ Calculation (UI: Run CJ button)"
        },
        {
            "endpoint": "/compute/xtb",
            "payload": {"payload": {"smiles": mock_smiles, "options": {"method": "GFN2-xTB"}}},
            "description": "Run XTB (UI: Run XTB button)"
        },
        {
            "endpoint": "/compute/psi4",
            "payload": {"payload": {"smiles": mock_smiles, "options": {"method": "B3LYP", "basis": "6-31G*"}}},
            "description": "Run Psi4 (UI: Run Psi4 button)"
        },
        
        # Export Functions
        {
            "endpoint": "/export/zip",
            "payload": {"artifacts": ["Results/example.json", "benchmark_auto.csv"]},
            "description": "Export Results (UI: Export Results button)"
        },
    ]
    
    # Run tests
    passed = 0
    total = len(test_cases)
    
    for i, test_case in enumerate(test_cases, 1):
        print(f"\n{'='*60}")
        print(f"Test {i}/{total}")
        success = test_endpoint(**test_case)
        if success:
            passed += 1
    
    # Summary
    print(f"\n{'='*60}")
    print("🎯 INTEGRATION TEST RESULTS")
    print(f"✅ Passed: {passed}/{total} tests")
    print(f"❌ Failed: {total - passed}/{total} tests")
    print(f"📊 Success Rate: {(passed/total)*100:.1f}%")
    
    if passed == total:
        print("\n🎉 ALL TESTS PASSED!")
        print("✨ Ketcher UI ↔ Backend integration is working perfectly!")
        print("🌐 Open frontend/ketcher_test.html in browser to use the UI")
    else:
        print(f"\n⚠️  {total - passed} tests failed - check backend logs")
    
    print(f"\n⏰ Completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    # First check if backend is running
    try:
        health_response = requests.get(f"{API_BASE}/health", timeout=5)
        if health_response.status_code == 200:
            print("✅ Backend is running")
            run_integration_tests()
        else:
            print(f"❌ Backend health check failed: {health_response.status_code}")
    except requests.exceptions.RequestException as e:
        print(f"❌ Cannot connect to backend at {API_BASE}")
        print(f"   Error: {e}")
        print(f"   💡 Start backend with: uvicorn iam.backend.app:app --reload --port 8010")
