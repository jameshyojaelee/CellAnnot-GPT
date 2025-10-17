# GPT Cell Annotator 0.2.0 – Installation Matrix QA

Date: 2025-10-17  
Engineer: Codex CLI assistant (GPT-5)

## Environment
- Host Python: `/usr/bin/python3` (3.10.12) with project-managed 3.11.14 via Poetry
- Wheel under test: `dist/gpt_cell_annotator-0.2.0-py3-none-any.whl`
- Extras exercised: base, `[api]`, `[ui]`, `[scanpy]`, `[full]`
- Assets cached to default `~/.cache/gpt-cell-annotator` unless noted

## Summary

| Scenario | Install Command | Smoke Checks | Result |
| --- | --- | --- | --- |
| Base | `pip install dist/...whl` | `gca --version`; `gca annotate data/demo/pbmc_markers.csv --offline` | ✅ Pass *(required fix: added `pyyaml` dependency)* |
| API | `pip install 'dist/...whl[api]'` | `gca api --help` | ✅ Pass |
| Scanpy | `pip install 'dist/...whl[scanpy]'` | `gca scanpy annotate /tmp/gca-scanpy-demo.h5ad --cluster-key leiden --offline` | ✅ Pass *(uses synthetic AnnData stub)* |
| UI | `pip install 'dist/...whl[ui]'` | `python -c "import streamlit, altair, fpdf"` | ✅ Pass |
| Full | `pip install 'dist/...whl[full]'` | `gca --version`; `gca api --help` | ✅ Pass |

## Notable Findings
- Initial base install failed: `ModuleNotFoundError: No module named 'yaml'`. Resolved by adding `pyyaml = "^6.0.2"` to core dependencies (`pyproject.toml` + `poetry.lock` regenerated).
- `gca scanpy annotate` expects an AnnData input; repo lacks a bundled `.h5ad`. Generated a synthetic `/tmp/gca-scanpy-demo.h5ad` for smoke testing.
- Global pytest run with `poetry run poe test` tripped over external plugins. Running tests with `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 poetry run pytest` succeeds (65 passing).
- `poetry check` passes with deprecated-field warnings (migration to PEP 621 sections still pending).

## Artifact Integrity
- `poetry build` regenerated both sdist and wheel post-dependency update.
- `poetry run twine check dist/*` ⇒ PASSED.
- SHA256 hashes (post-fix):
  - `bfd84e9837d9cc61ff88661bee5ad66af7f63523d6eba6b55c6dd2118af2dc8d  dist/gpt_cell_annotator-0.2.0-py3-none-any.whl`
  - `f0b6531b06312e2f0cb49d718d8e7074ad0919b8efc5a46666be0ca55278161e  dist/gpt_cell_annotator-0.2.0.tar.gz`

## Next Actions
1. Capture pytest invocation guidance in `RELEASING.md` (note the `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1` requirement).
2. Consider shipping a lightweight AnnData demo or update docs with synthetic dataset instructions for scanpy smoke tests.
3. Proceed to publish via `twine upload dist/*` once release notes/changelog are final.
