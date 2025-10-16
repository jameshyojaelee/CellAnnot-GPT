test_that("HTTP transport can use httptest2 fixtures", {
  skip_if_not_installed("httptest2")
  cfg <- gptca_config(base_url = "https://mocked.api", offline = FALSE, cli_path = NA_character_)
  marker_list <- setNames(list(c("MS4A1", "CD79A")), "0")
  result <- httptest2::with_mock_dir("fixtures/http", {
    gptca_annotate_markers(marker_list, config = cfg, fallback = FALSE, return_validated = FALSE)
  })
  expect_equal(result$clusters$primary_label, "B cell")
  expect_false(result$validated)
})
