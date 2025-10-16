# Releasing the gptcellannotator R package

## Versioning

- Follow semantic versioning (major.minor.patch) aligned with the Python backend API contract.
- Bump the version in `clients/r/gptcellannotator/DESCRIPTION` and update `NEWS.md` (or release notes).
- Ensure the Python backend schema has not changed. If it has, coordinate a synchronized release with the backend team.

## Pre-flight checklist

1. Pull the latest `main` branch and ensure the working tree is clean.
2. Run the automation helpers from the repository root:
   ```bash
   make -C clients/r doc
   make -C clients/r test
   make -C clients/r check
   make -C clients/r pkgdown
   ```
3. Confirm `R CMD check --as-cran` passes on macOS and Linux (CI already covers this).
4. Inspect test coverage via `covr::report()` and ensure it stays above 80\% for exported functions.
5. Regenerate the vignette and README if examples changed.
6. Review the pkgdown site locally (`clients/r/gptcellannotator/docs/index.html`).
7. Update `clients/r/gptcellannotator/NEWS.md` (create if missing) with noteworthy changes.

## Dependency management

- Avoid pinning CRAN package versions unless the backend requires it; prefer minimum version comments in `DESCRIPTION` if needed.
- Changes to Python CLI behaviours must be reflected in tests that mock `gca` responses.
- When adding dependencies, update `Imports`/`Suggests` and ensure CI installs them.

## Publishing

### GitHub release

1. Create a git tag matching the version (e.g., `r-v0.1.0`).
2. Push the tag to GitHub; CI will build pkgdown and attach check artefacts.
3. Draft a GitHub release summarising changes, validation steps, and backend compatibility.

### CRAN (optional)

1. Follow the pre-flight checklist above and run `devtools::check_rhub()` and `devtools::check_win_devel()`.
2. Ensure reverse dependencies (if any) are notified.
3. Submit via `devtools::release()` and monitor CRAN feedback.
4. Once accepted, update installation instructions and tag the git commit.

## Post-release

- Increment the development version (append `.9000`) in `DESCRIPTION`.
- Communicate schema or contract changes with frontend/backend teams.
- Keep the pkgdown site updated with the release notes.
