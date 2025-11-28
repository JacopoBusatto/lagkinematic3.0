# AGENTS.md — Guidelines for Codex agents

This repository contains the Python rewrite of the LAGKINEMATIC 1.0 model.

## Project Layout
- `python/lagkinematic/` → main Python package
- `legacy_cpp/` → read-only reference C++ implementation (do NOT edit)
- `config/` → configuration files
- `tests/` → unit tests

## Codex Rules
1. **NEVER modify files in `legacy_cpp/`. They are read-only reference.**
2. For each new module translated from C++, create:
   - a Python file under `python/lagkinematic/...`
   - a related test under `tests/`
3. Prefer small, modular functions over monolithic classes.
4. Use:
   - type hints (`mypy`-friendly)
   - docstrings following numpy-style
   - Pydantic for validation (if needed)
5. Before proposing modifications:
   - run: `pytest -q`
   - run: `ruff check .`

## High-level goal for Codex
Gradually port the functionality from `legacy_cpp/` to Python, keeping:
- reproducibility
- testability
- clarity
- numerical equivalence to the C++ reference when possible

Codex: always explain what you plan to change, then provide a diff.
