test_that("HTTP transport can use httptest2 fixtures", {
  skip_if_not_installed("httptest2")
  cfg <- gptca_config(base_url = "https://mocked.api", offline = FALSE, cli_path = NA_character_)
  result <- httptest2::with_mock_dir("fixtures/batch", {
    gptca_annotate_markers(list(`0` = c("MS4A1", "CD79A")), config = cfg, fallback = FALSE)
  })
  expect_equal(result$clusters$primary_label, "B cell")
  expect_true(result$validated)
})
