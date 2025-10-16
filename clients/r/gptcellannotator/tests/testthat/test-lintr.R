test_that("package style lints pass", {
  skip_if_not_installed("lintr")
  skip_on_cran()
  pkg_path <- system.file(package = "gptcellannotator")
  linters <- lintr::linters_with_defaults(
    line_length_linter = lintr::line_length_linter(160),
    object_length_linter = NULL
  )
  lintr::expect_lint_free(path = pkg_path, linters = linters, exclusions = list("tests/testthat"))
})
