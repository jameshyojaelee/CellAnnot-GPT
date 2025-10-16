test_that("package style lints pass", {
  skip_if_not_installed("lintr")
  skip_on_cran()
  lintr::expect_lint_free()
})
