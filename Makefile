.PHONY: install format lint test typecheck ruff ruff-fix build-marker-db api ui clean

install:
	poetry install --extras "dev" --extras "scanpy" --extras "api" --extras "ui"

format:
	poetry run black backend frontend evaluation scripts config tests

ruff:
	poetry run ruff check backend frontend evaluation scripts config tests

ruff-fix:
	poetry run ruff check --fix backend frontend evaluation scripts config tests

lint: ruff typecheck

typecheck:
	poetry run mypy backend

test:
	PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 poetry run pytest

release:
	poetry build
	poetry run twine check dist/*

build-marker-db:
	poetry run python scripts/build_marker_db.py

api:
	poetry run uvicorn backend.api.main:app --reload

ui:
	poetry run streamlit run frontend/streamlit_app.py

clean:
	rm -rf __pycache__ .pytest_cache .ruff_cache .mypy_cache
	rm -f marker_db.sqlite
