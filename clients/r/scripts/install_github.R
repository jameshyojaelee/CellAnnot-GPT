#!/usr/bin/env Rscript
# Install gptcellannotator from the GitHub repo using pak (preferred) or remotes.

args <- commandArgs(trailingOnly = TRUE)
ref <- if (length(args) >= 1) args[[1]] else Sys.getenv("GPTCA_R_GITHUB_REF", "main")

repo_spec <- sprintf(
  "github::jameshyojaelee/CellAnnot-GPT@%s?subdir=clients/r/gptcellannotator",
  ref
)

message("Installing gptcellannotator from ", repo_spec)

ensure_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("pak", quietly = TRUE)) {
  message("pak not found; attempting to install from CRAN...")
  ensure_package("pak")
}

installed <- FALSE

if (requireNamespace("pak", quietly = TRUE)) {
  pak::pkg_install(repo_spec, ask = FALSE)
  installed <- TRUE
}

if (!installed) {
  message("Falling back to remotes::install_github()")
  ensure_package("remotes")
  remotes::install_github(
    "jameshyojaelee/CellAnnot-GPT",
    subdir = "clients/r/gptcellannotator",
    ref = ref,
    build = TRUE,
    dependencies = TRUE
  )
}

message("gptcellannotator installation complete.")
