# delta-phylo — Project Roadmap

This document tracks planned work for the **delta-phylo** package, organised by
milestone. Each item links to the corresponding GitHub issue.

---

## v0.1.x — Bug fixes (critical correctness)

These issues affect the correctness of results and should be resolved before
any new feature work.

| # | Issue | Description |
|---|-------|-------------|
| [#2](https://github.com/nathanufpb/morfology_analysis/issues/2) | **BUG** | Incorrect regex for ordered character types (`OM`) in DELTA parser |
| [#3](https://github.com/nathanufpb/morfology_analysis/issues/3) | **BUG** | Autapomorphy detection logic is incorrect in `MatrixCleaner` |
| [#4](https://github.com/nathanufpb/morfology_analysis/issues/4) | **BUG** | `mk_likelihood` ignores tree topology — not a real ML implementation |

**Definition of done for v0.1.x:**  
All three bugs fixed, each accompanied by a failing-before / passing-after
regression test.

---

## v0.2.0 — Test coverage (reliability foundation)

Before adding features, the project needs a comprehensive test suite. All
issues in this milestone should be resolved together (they are mostly
independent and can be worked on in parallel).

| # | Issue | File to create |
|---|-------|----------------|
| [#7](https://github.com/nathanufpb/morfology_analysis/issues/7) | **TESTS** | `tests/test_matrix_cleaning.py` — `MatrixCleaner` |
| [#8](https://github.com/nathanufpb/morfology_analysis/issues/8) | **TESTS** | `tests/test_encoding.py` — `MatrixEncoder` |
| [#9](https://github.com/nathanufpb/morfology_analysis/issues/9) | **TESTS** | `tests/test_likelihood.py` — `LikelihoodAnalysis` / `mk_likelihood` |
| [#10](https://github.com/nathanufpb/morfology_analysis/issues/10) | **TESTS** | `tests/test_validation.py` — validation utilities |
| [#11](https://github.com/nathanufpb/morfology_analysis/issues/11) | **TESTS** | `tests/test_tree_utils.py` — `TreeUtils` |
| [#12](https://github.com/nathanufpb/morfology_analysis/issues/12) | **TESTS** | `tests/test_cli.py` — CLI commands (Click `CliRunner`) |
| [#13](https://github.com/nathanufpb/morfology_analysis/issues/13) | **TESTS** | `tests/test_integration.py` — end-to-end with bundled example dataset |

**Target coverage goal:** ≥ 85% line coverage measured by `pytest --cov`.

**Definition of done for v0.2.0:**  
All test files created and passing; CI (`.github/workflows/ci.yml`) green on
Python 3.9 – 3.12.

---

## v0.3.0 — Feature improvements

| # | Issue | Description |
|---|-------|-------------|
| [#6](https://github.com/nathanufpb/morfology_analysis/issues/6) | **FEATURE** | NEXUS interleaved format (`interleave=True`) in `NexusWriter` |
| [#15](https://github.com/nathanufpb/morfology_analysis/issues/15) | **FEATURE** | Polymorphic strategy `'majority'` in `MatrixBuilder` |
| [#5](https://github.com/nathanufpb/morfology_analysis/issues/5) | **FEATURE** | Support numeric/continuous characters (`CharacterType.NUMERIC`) + Gower distance |

**Definition of done for v0.3.0:**  
All features implemented, each with corresponding unit tests; changelog updated.

---

## v0.4.0 — Refactoring & developer experience

| # | Issue | Description |
|---|-------|-------------|
| [#14](https://github.com/nathanufpb/morfology_analysis/issues/14) | **REFACTOR** | Add `__all__` exports to all `__init__.py` modules for a clean public API |

Additional housekeeping targeted for this release:
- Add `CONTRIBUTING.md` with contribution guidelines and branch strategy.
- Add GitHub issue/PR templates (`.github/ISSUE_TEMPLATE/`, `.github/PULL_REQUEST_TEMPLATE.md`).
- Configure `mypy` strict mode in `setup.cfg` / `pyproject.toml`.
- Pin development dependencies in `requirements-dev.txt`.
- Publish to PyPI as `delta-phylo 0.4.0`.

---

## Backlog (unscheduled)

Items identified during review that do not yet have an issue or a milestone:

- **Gamma rate variation (+G)** for the Mk likelihood model.
- **SPR (Subtree Pruning and Regrafting)** tree search in addition to NNI.
- **Bootstrap support values** on parsimony/distance trees.
- **Interactive visualisation** using Plotly (replace static matplotlib for tree plots).
- **Bayesian-ready MrBayes block** auto-tuning (ngen, samplefreq) based on matrix size.
- **Support for `conch` / `Intkey` DELTA directives** beyond the current subset.

---

## Issue label legend

| Label | Meaning |
|-------|---------|
| `bug` | Incorrect behaviour that must be fixed |
| `enhancement` | New capability |
| `testing` | Test coverage work |
| `refactor` | Code quality / architecture improvements |
| `documentation` | Docs only |

---

*Last updated: 2026-03-05*
