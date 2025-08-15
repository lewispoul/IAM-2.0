"""
Test static file serving for local Ketcher hosting.

These tests verify that static files are served correctly,
but skip if Ketcher is not installed locally.
"""

import pytest
import requests
import os
from pathlib import Path


BASE_URL = "http://127.0.0.1:8011"
KETCHER_INDEX = Path("public/static/ketcher/index.html")


def is_ketcher_installed():
    """Check if Ketcher is installed locally."""
    return KETCHER_INDEX.exists()


def is_server_running():
    """Check if the backend server is running."""
    try:
        response = requests.get(f"{BASE_URL}/health", timeout=2)
        return response.status_code == 200
    except requests.exceptions.RequestException:
        return False


pytestmark = pytest.mark.skipif(
    not is_server_running(),
    reason="Backend server not running on 8011"
)


@pytest.mark.skipif(not is_ketcher_installed(), reason="Ketcher not installed locally")
def test_ketcher_static_index():
    """Test that Ketcher index.html is served correctly."""
    response = requests.get(f"{BASE_URL}/static/ketcher/index.html")
    
    assert response.status_code == 200
    assert "text/html" in response.headers.get("content-type", "")
    assert "ketcher" in response.text.lower()


@pytest.mark.skipif(not is_ketcher_installed(), reason="Ketcher not installed locally")  
def test_ketcher_bridge_js():
    """Test that bridge.js is served correctly."""
    response = requests.get(f"{BASE_URL}/static/ketcher/bridge.js")
    
    assert response.status_code == 200
    assert "javascript" in response.headers.get("content-type", "") or "text/plain" in response.headers.get("content-type", "")
    assert "KetcherBridge" in response.text


def test_main_ui_serves():
    """Test that the main UI serves correctly."""
    response = requests.get(f"{BASE_URL}/")
    
    assert response.status_code == 200
    assert "text/html" in response.headers.get("content-type", "")
    
    # Should contain either the full UI or fallback message
    content = response.text.lower()
    assert "iam-2.0" in content
    

@pytest.mark.skipif(not is_ketcher_installed(), reason="Ketcher not installed locally")
def test_ketcher_integration_components():
    """Test that all integration components are present."""
    # Main UI
    response = requests.get(f"{BASE_URL}/")
    assert response.status_code == 200
    
    # Check for iframe pointing to local Ketcher
    if is_ketcher_installed():
        assert "/static/ketcher/index.html" in response.text
    
    # Static Ketcher
    response = requests.get(f"{BASE_URL}/static/ketcher/index.html")
    assert response.status_code == 200
    
    # Bridge JS
    response = requests.get(f"{BASE_URL}/static/ketcher/bridge.js")
    assert response.status_code == 200
    assert "getMolfile" in response.text
    assert "setSmiles" in response.text


def test_static_mount_available():
    """Test that static file mounting is working."""
    # Even if Ketcher is not installed, the static mount should be available
    response = requests.get(f"{BASE_URL}/static/")
    
    # Should return either directory listing or 404, but not 500
    assert response.status_code in [200, 403, 404]


@pytest.mark.skipif(is_ketcher_installed(), reason="Ketcher is installed")
def test_fallback_when_ketcher_missing():
    """Test behavior when Ketcher is not installed."""
    # Main UI should still work with fallback message
    response = requests.get(f"{BASE_URL}/")
    
    assert response.status_code == 200
    content = response.text.lower()
    assert "iam-2.0" in content
    
    # Should mention static/ketcher as alternative
    if "/static/ketcher" not in response.text:
        # Fallback UI should be shown
        assert "api docs" in content
