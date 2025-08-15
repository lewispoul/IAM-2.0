# MIGRATION.md

> Placeholder for IAM 1.x → 2.0 controlled rebuild plan.

## Environment Activation (2025-08-15)
- The `iam2` conda environment is now set to auto-activate for all new terminal sessions via `~/.bashrc`.
- All tests and development should be run in the `iam2` environment for consistent dependencies and results.
- To disable auto-activation, remove or comment out `conda activate iam2` in your `~/.bashrc`.
