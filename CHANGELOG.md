# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2025-08-13
### Added
- New repository scaffold for IAM 2.0
- Full documentation set, including **IAM_2.0_Whitepaper.md** with roadmap and legacy appendix
- Pre-generated diagram PNGs in `/docs/images/` (Mermaid source embedded in the whitepaper)

### Changed
- Rebuilt architecture around unified results schema and FastAPI

### Removed
- Legacy IAM 1.x (Nox) fragile templates and ad-hoc outputs

## [Unreleased]
- Planned: Psi4 frequency jobs & spectra viewers, ML predictor registry, Gaussian HPC connector

## [2.0.0-alpha] - 2025-08-16
- FastAPI backend + legacy compatibility layer
- Export/persistence (CSV/JSON/ZIP) with path-safety
- Empirical predictors integrated
- Test harness: 113 passed, 2 skipped (external)
- Ketcher static served locally
