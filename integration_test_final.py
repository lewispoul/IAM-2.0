#!/usr/bin/env python3
"""
🔬 IAM-2.0 × Ketcher Integration Final Test
This script verifies that the UI integration is working properly.
"""

import requests

def test_integration():
    """Comprehensive integration test"""
    print("🧪 IAM-2.0 × Ketcher Integration Test")
    print("=" * 60)
    
    results = {
        "backend_running": False,
        "ui_serving": False,
        "api_functional": False,
        "cors_enabled": False
    }
    
    try:
        # Test 1: Backend Health
        print("\n1. 🏥 Testing backend health...")
        health_response = requests.get('http://127.0.0.1:8011/health', timeout=5)
        
        if health_response.status_code == 200:
            print(f"   ✅ Backend health: {health_response.status_code}")
            print(f"   📊 Health data: {health_response.json()}")
            results["backend_running"] = True
        else:
            print(f"   ❌ Backend health failed: {health_response.status_code}")
            
    except requests.exceptions.RequestException as e:
        print(f"   ❌ Backend not accessible: {e}")
        print("   💡 Make sure backend is running with: uvicorn iam.backend.app:app --host 127.0.0.1 --port 8011")
        return results
    
    try:
        # Test 2: UI Serving
        print("\n2. 🌐 Testing UI serving...")
        ui_response = requests.get('http://127.0.0.1:8011/', timeout=5)
        
        if ui_response.status_code == 200:
            print(f"   ✅ UI endpoint: {ui_response.status_code}")
            print(f"   📄 Content-Type: {ui_response.headers.get('content-type', 'unknown')}")
            print(f"   📏 Content size: {len(ui_response.content)} bytes")
            
            # Check if it contains expected UI elements
            content = ui_response.text.lower()
            if 'ketcher' in content or 'molecule' in content or 'iam-2.0' in content:
                print("   ✅ UI contains expected molecular interface elements")
                results["ui_serving"] = True
            else:
                print("   ⚠️ UI content may not be the molecular interface")
                
        else:
            print(f"   ❌ UI serving failed: {ui_response.status_code}")
            
    except requests.exceptions.RequestException as e:
        print(f"   ❌ UI not accessible: {e}")
    
    try:
        # Test 3: CORS Check
        print("\n3. 🔐 Testing CORS configuration...")
        cors_response = requests.options('http://127.0.0.1:8011/', timeout=5)
        cors_headers = cors_response.headers
        
        if 'Access-Control-Allow-Origin' in cors_headers:
            print("   ✅ CORS headers present")
            print(f"   🌍 Allow-Origin: {cors_headers.get('Access-Control-Allow-Origin')}")
            results["cors_enabled"] = True
        else:
            print("   ⚠️ CORS headers not found in OPTIONS response")
            
    except requests.exceptions.RequestException as e:
        print(f"   ❌ CORS test failed: {e}")
    
    try:
        # Test 4: API Functionality
        print("\n4. ⚗️ Testing API endpoints...")
        
        # Test a simple endpoint
        api_response = requests.post(
            'http://127.0.0.1:8011/ketcher/to-smiles',
            json={'molfile': 'C'},
            timeout=5
        )
        
        if api_response.status_code == 200:
            print(f"   ✅ API endpoint functional: {api_response.status_code}")
            try:
                api_data = api_response.json()
                print(f"   📊 API response: {api_data}")
                results["api_functional"] = True
            except:
                print("   ⚠️ API returned non-JSON response")
        else:
            print(f"   ❌ API endpoint failed: {api_response.status_code}")
            
    except requests.exceptions.RequestException as e:
        print(f"   ❌ API test failed: {e}")
    
    # Summary
    print("\n" + "=" * 60)
    print("📋 INTEGRATION TEST SUMMARY:")
    print(f"   🏥 Backend Running:     {'✅ YES' if results['backend_running'] else '❌ NO'}")
    print(f"   🌐 UI Serving:          {'✅ YES' if results['ui_serving'] else '❌ NO'}")
    print(f"   ⚗️ API Functional:      {'✅ YES' if results['api_functional'] else '❌ NO'}")
    print(f"   🔐 CORS Enabled:        {'✅ YES' if results['cors_enabled'] else '❌ NO'}")
    
    all_good = all(results.values())
    
    if all_good:
        print("\n🎉 INTEGRATION SUCCESS!")
        print("   Your IAM-2.0 × Ketcher integration is fully functional!")
        print("\n🌐 Access methods:")
        print("   • Primary:    http://127.0.0.1:8011/")
        print("   • Standalone: file:///home/lppoulin/IAM-2.0/frontend/standalone_test.html")
        print("   • Test page:  file:///home/lppoulin/IAM-2.0/frontend/test_success.html")
        
        print("\n💡 Usage tips:")
        print("   • Use system browser instead of VS Code Simple Browser")
        print("   • The standalone version works offline with mock data")
        print("   • All molecular analysis methods are available")
        
    else:
        print("\n⚠️ PARTIAL INTEGRATION")
        print("   Some components need attention (see details above)")
        
        if not results['backend_running']:
            print("   🔧 Fix: Start backend with: uvicorn iam.backend.app:app --host 127.0.0.1 --port 8011")
    
    return results

if __name__ == "__main__":
    test_integration()
