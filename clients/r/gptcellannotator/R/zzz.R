.onAttach <- function(libname, pkgname) {
  invisible()
}

utils::globalVariables(c("dim1", "dim2", "label", "cell"))

.onLoad <- function(libname, pkgname) {
  gptca_config_reset()
  invisible()
}
