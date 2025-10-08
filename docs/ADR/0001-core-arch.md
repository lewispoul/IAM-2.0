# ADR 0001: Core Architecture Decision

## Status
Accepted

## Context
The Nox project requires a modular architecture to handle IAM (Identity and Access Management) functionalities across multiple packages: nox, iam, and pinox. The system needs to support WSL/Ubuntu environments and potentially Raspberry Pi deployments.

## Decision
We adopt a multi-package Python structure with:
- **nox**: Core orchestration and main logic
- **iam**: Identity and access management components
- **pinox**: Peripheral integrations and extensions

Framework choices:
- Flask/FastAPI for web APIs based on module needs
- WSL/Ubuntu as primary development environment
- Raspberry Pi support for edge deployments

## Consequences
- Modular development allows independent evolution of components
- Framework flexibility enables optimal tool selection per module
- Cross-platform compatibility ensures broad deployment options