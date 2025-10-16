test_that("PBMC integration flow works with mocked backend", {
  markers_path <- system.file("extdata", "pbmc_markers.csv", package = "gptcellannotator")
  pbmc <- utils::read.csv(markers_path, stringsAsFactors = FALSE)
  fake_response <- list(
    result = list(
      summary = list(
        total_clusters = 4L,
        supported_clusters = 3L,
        flagged_clusters = 1L,
        unknown_clusters = list("3")
      ),
      clusters = lapply(
        pbmc$cluster_id,
        function(cid) {
          list(
            cluster_id = as.character(cid),
            status = if (cid == 3) "unknown" else "supported",
            confidence = if (cid == 3) "Low" else "High",
            warnings = if (cid == 3) list("Insufficient markers") else list(),
            annotation = list(
              primary_label = switch(
                as.character(cid),
                "0" = "B cell",
                "1" = "T cell",
                "2" = "NK cell",
                "3" = "Unknown"
              ),
              ontology_id = switch(
                as.character(cid),
                "0" = "CL:0000236",
                "1" = "CL:0000625",
                "2" = "CL:0000629",
                NULL
              ),
              confidence = if (cid == 3) "Low" else "High",
              rationale = "Mock response",
              markers = jsonlite::fromJSON(pbmc$markers[pbmc$cluster_id == cid])
            )
          )
        }
      )
    )
  )

  cfg <- gptca_config(base_url = "https://mock.api", offline = FALSE, cli_path = NA_character_)
  annotations <- with_mocked_bindings(
    gptca_annotate_markers(pbmc, config = cfg, fallback = FALSE),
    gptca_http_post = function(path, body, config) fake_response
  )

  expect_equal(nrow(annotations$clusters), nrow(pbmc))
  expect_equal(annotations$summary$total_clusters, 4L)
  expect_equal(
    annotations$clusters$primary_label,
    c("B cell", "T cell", "NK cell", "Unknown")
  )
})
