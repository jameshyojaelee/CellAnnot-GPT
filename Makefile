.PHONY: install format lint test typecheck ruff ruff-fix build-marker-db api ui clean

install:
	poetry install --extras "dev,scanpy,api,ui"

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
	poetry run pytest

build-marker-db:
	poetry run python scripts/build_marker_db.py

api:
	poetry run uvicorn backend.api.main:app --reload

ui:
	poetry run streamlit run frontend/streamlit_app.py

clean:
	rm -rf __pycache__ .pytest_cache .ruff_cache .mypy_cache
	rm -f marker_db.sqlite
