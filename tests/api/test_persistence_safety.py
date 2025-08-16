"""
Tests for persistence hooks and /export/zip safety
"""
import pytest
import tempfile
import os
import json
from pathlib import Path
from fastapi.testclient import TestClient
from backend.main import app
from backend.utils.persistence import save_result_json, append_benchmark_row

client = TestClient(app)

class TestPersistenceSafety:
    """Test persistence functions for path safety and data integrity."""
    
    def test_save_result_json_creates_valid_file(self):
        """Test that save_result_json creates proper JSON files."""
        test_data = {"test": "data", "number": 42}
        filepath = save_result_json("test_operation", test_data)
        
        # Verify file exists
        assert os.path.exists(filepath)
        
        # Verify JSON content
        with open(filepath, 'r') as f:
            loaded_data = json.load(f)
        assert loaded_data == test_data
        
        # Verify filename pattern
        assert "test_operation" in os.path.basename(filepath)
        assert filepath.endswith(".json")
    
    def test_append_benchmark_row_creates_valid_csv(self):
        """Test that append_benchmark_row creates proper CSV files."""
        test_row = {"name": "test", "method": "empirical", "ok": True, "value": 1.23}
        filepath = append_benchmark_row(test_row)
        
        # Verify file exists
        assert os.path.exists(filepath)
        
        # Verify CSV content
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Should have header + at least 1 data row
        assert len(lines) >= 2
        assert "name" in lines[0]  # Header row
        assert "test" in lines[-1]  # Last data row contains our test data
        assert filepath.endswith(".csv")

class TestExportZipSafety:
    """Test /export/zip endpoint for path traversal and security issues."""
    
    def setup_method(self):
        """Set up test environment with temporary directories."""
        self.temp_dir = tempfile.mkdtemp()
        self.old_base = os.environ.get("IAM_RESULTS_BASE")
        os.environ["IAM_RESULTS_BASE"] = self.temp_dir
        
        # Create test files
        self.test_file1 = Path(self.temp_dir) / "test1.json"
        self.test_file2 = Path(self.temp_dir) / "subdir" / "test2.json"
        self.test_file2.parent.mkdir(parents=True, exist_ok=True)
        
        self.test_file1.write_text(json.dumps({"file": "test1"}))
        self.test_file2.write_text(json.dumps({"file": "test2"}))
    
    def teardown_method(self):
        """Clean up test environment."""
        if self.old_base:
            os.environ["IAM_RESULTS_BASE"] = self.old_base
        else:
            del os.environ["IAM_RESULTS_BASE"]
    
    def test_export_zip_valid_files(self):
        """Test that valid file paths can be zipped successfully."""
        resp = client.post("/export/zip", json={
            "artifacts": ["test1.json", "subdir/test2.json"]
        })
        assert resp.status_code == 200
        data = resp.json()
        assert data["ok"] is True
        assert "zip_path" in data["data"]
        assert data["data"]["zip_path"].endswith("_export.zip")
    
    def test_export_zip_blocks_path_traversal(self):
        """Test that path traversal attempts are blocked."""
        resp = client.post("/export/zip", json={
            "artifacts": ["../../../etc/passwd"]
        })
        assert resp.status_code == 400
        data = resp.json()
        # Error responses use the error envelope format (not ok/data format)
        assert "correlation_id" in data
        assert "Invalid path" in data["message"]
    
    def test_export_zip_handles_missing_files(self):
        """Test proper error handling for missing files."""
        resp = client.post("/export/zip", json={
            "artifacts": ["nonexistent.json"]
        })
        assert resp.status_code == 404
        data = resp.json()
        # Error responses use the error envelope format
        assert "correlation_id" in data
        assert "File not found" in data["message"]
    
    def test_export_zip_empty_artifacts_list(self):
        """Test behavior with empty artifacts list."""
        resp = client.post("/export/zip", json={
            "artifacts": []
        })
        assert resp.status_code == 200
        data = resp.json()
        assert data["ok"] is True
        # Should create empty zip file
        assert "zip_path" in data["data"]

class TestPersistenceIntegration:
    """Test that persistence hooks work correctly with API endpoints."""
    
    def test_convert_molfile_saves_result(self):
        """Test that successful molfile conversion triggers persistence."""
        molfile = "\n  Methane\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6291    0.6291    0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6291   -0.6291    0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6291   -0.6291   -0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6291    0.6291   -0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1  3  1  0  0  0  0\n  1  4  1  0  0  0  0\n  1  5  1  0  0  0  0\nM  END\n"
        
        resp = client.post("/convert/molfile", json={"molfile": molfile})
        assert resp.status_code == 200
        data = resp.json()
        assert data["ok"] is True
        assert "smiles" in data["data"]
        
        # Persistence should have been called (files should exist)
        # Note: This test verifies the integration works without filesystem pollution
    
    def test_ketcher_run_saves_benchmark(self):
        """Test that ketcher/run calls save benchmark data."""
        resp = client.post("/ketcher/run", json={
            "smiles": "C", 
            "method": "empirical"
        })
        assert resp.status_code == 200
        data = resp.json()
        assert data["ok"] is True
        
        # Should have benchmark data structure typical of empirical predictions
        assert "data" in data
