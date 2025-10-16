fake_validated_response <- function() {
  list(
    result = list(
      summary = list(
        total_clusters = 1L,
        supported_clusters = 1L,
        flagged_clusters = 0L,
        unknown_clusters = list()
      ),
      clusters = list(
        list(
          cluster_id = "0",
          status = "supported",
          confidence = "High",
          warnings = list(),
          annotation = list(
            primary_label = "B cell",
            ontology_id = "CL:0000236",
            confidence = "High",
            rationale = "Markers support B cell identity",
            markers = list("MS4A1", "CD79A")
          ),
          validation = list(
            cluster_id = "0",
            primary_label = "B cell",
            is_supported = TRUE
          )
        )
      )
    )
  )
}

test_that("gptca_prepare_clusters handles named lists", {
  clusters <- gptca_prepare_clusters(list(`0` = c("MS4A1", "CD79A")))
  expect_equal(length(clusters), 1L)
  expect_equal(clusters[[1]]$cluster_id, "0")
  expect_equal(clusters[[1]]$markers, c("MS4A1", "CD79A"))
})

test_that("gptca_prepare_clusters handles data frames", {
  df <- data.frame(
    cluster = c(0, 0, 1, 1),
    gene = c("MS4A1", "CD79A", "CD3E", "LCK"),
    avg_log2FC = c(2, 1.5, 1.8, 1.2)
  )
  clusters <- gptca_prepare_clusters(df, top_n = 1)
  expect_equal(length(clusters), 2L)
  expect_equal(clusters[[1]]$markers, "MS4A1")
})

test_that("gptca_annotate_markers parses validated response", {
  cfg <- gptca_config(base_url = "https://mock.api", offline = FALSE, cli_path = NA_character_)
  result <- with_mocked_bindings(
    gptca_annotate_markers(list(`0` = c("MS4A1", "CD79A")), config = cfg, fallback = FALSE),
    gptca_http_post = function(path, body, config) fake_validated_response()
  )
  expect_s3_class(result, "gptca_annotation")
  expect_true(result$validated)
  expect_equal(result$clusters$primary_label, "B cell")
  expect_equal(result$summary$total_clusters, 1L)
})

test_that("gptca_annotate_markers falls back to CLI", {
  cli_script <- withr::local_tempfile()
  writeLines(
    c(
      "#!/usr/bin/env bash",
      "while [[ $# -gt 0 ]]; do",
      "  if [[ $1 == \"--out-json\" ]]; then",
      "    shift",
      "    out=$1",
      "  fi",
      "  shift",
      "done",
      "cat <<'JSON' > \"$out\"",
      "{",
      "  \"summary\": {\"total_clusters\": 1, \"supported_clusters\": 0, \"flagged_clusters\": 1, \"unknown_clusters\": []},",
      "  \"clusters\": [",
      "    {",
      "      \"cluster_id\": \"0\",",
      "      \"status\": \"flagged\",",
      "      \"confidence\": \"Low\",",
      "      \"warnings\": [\"Mocked\"],",
      "      \"annotation\": {\"primary_label\": \"Unknown\", \"ontology_id\": null, \"confidence\": \"Low\", \"rationale\": \"\", \"markers\": []}",
      "    }",
      "  ]",
      "}",
      "JSON"
    ),
    cli_script
  )
  Sys.chmod(cli_script, mode = "0755")

  cfg <- gptca_config(base_url = "https://mock.api", cli_path = cli_script, offline = FALSE)

  warnings <- list()
  result <- withCallingHandlers(
    with_mocked_bindings(
      gptca_annotate_markers(list(`0` = c("X")), config = cfg, fallback = TRUE),
      gptca_http_post = function(path, body, config) stop("boom")
    ),
    warning = function(w) {
      warnings <<- c(warnings, list(w))
    }
  )
  expect_true(any(grepl("HTTP request failed", vapply(warnings, conditionMessage, character(1)), fixed = TRUE)))
  expect_true(result$validated)
  expect_equal(result$clusters$status, "flagged")
  expect_identical(result$clusters$warnings[[1]], list("Mocked"))
})
