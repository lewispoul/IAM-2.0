#!/usr/bin/env python3

import requests
import json
from fastapi.testclient import TestClient
import sys
sys.path.append('/home/lppoulin/IAM-2.0')

from backend.main import app

client = TestClient(app)

print("Testing /run/ endpoint with empty XYZ...")
response = client.post("/run/", json={"xyz": ""})
print(f"Status Code: {response.status_code}")
print(f"Response: {response.json()}")

print("\nTesting /run/ endpoint with missing XYZ...")
response2 = client.post("/run/", json={})
print(f"Status Code: {response2.status_code}")
print(f"Response: {response2.json()}")
