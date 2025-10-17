#!/usr/bin/env Rscript
# Build and install gptcellannotator from a local source tarball.

pkg_root <- normalizePath(file.path("clients", "r", "gptcellannotator"))
dest_dir <- normalizePath(ifelse(nchar(Sys.getenv("GPTCA_R_TARBALL_DIR")) > 0,
  Sys.getenv("GPTCA_R_TARBALL_DIR"),
  tempdir()
))

message("Building package from ", pkg_root)

ensure_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

ensure_package("pkgbuild")
ensure_package("pak")

tarball_path <- pkgbuild::build(
  pkg = pkg_root,
  dest_path = dest_dir,
  binary = FALSE,
  quiet = TRUE
)

message("Built tarball: ", tarball_path)

pak::pkg_install(tarball_path, ask = FALSE)

message("gptcellannotator installed from local tarball.")
