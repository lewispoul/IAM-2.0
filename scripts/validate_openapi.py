import sys
import yaml
from openapi_spec_validator import validate_spec

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "openapi/iam2.v1.yaml"
    with open(path, "r") as f:
        spec = yaml.safe_load(f)
    validate_spec(spec)
    print(f"OpenAPI spec {path} is valid.")

if __name__ == "__main__":
    main()
