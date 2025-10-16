test_that("PBMC integration flow works with mocked backend", {
  markers_path <- system.file("extdata", "pbmc_markers.csv", package = "gptcellannotator")
  pbmc <- utils::read.csv(markers_path, stringsAsFactors = FALSE)
  pbmc$marker_list <- lapply(pbmc$markers, function(x) {
    cleaned <- gsub("\\\\", "", x)
    cleaned <- gsub("\\[|\\]", "", cleaned)
    strsplit(cleaned, ",")[[1]] |> trimws() |> toupper()
  })
  label_map <- c(
    "0" = "B cell",
    "1" = "T cell",
    "2" = "NK cell",
    "3" = "Unknown",
    "4" = "Monocyte",
    "5" = "Cycling",
    "6" = "Platelet",
    "7" = "Erythrocyte"
  )
  ontology_map <- c(
    "0" = "CL:0000236",
    "1" = "CL:0000625",
    "2" = "CL:0000629",
    "4" = "CL:0000576",
    "6" = "CL:0000545"
  )
  fake_response <- list(
    result = list(
      summary = list(
        total_clusters = length(pbmc$cluster_id),
        supported_clusters = length(pbmc$cluster_id) - 1L,
        flagged_clusters = 1L,
        unknown_clusters = list("3")
      ),
      clusters = lapply(
        pbmc$cluster_id,
        function(cid) {
          cid_chr <- as.character(cid)
          list(
            cluster_id = cid_chr,
            status = if (cid == 3) "unknown" else "supported",
            confidence = if (cid == 3) "Low" else "High",
            warnings = if (cid == 3) list("Insufficient markers") else list(),
            annotation = list(
              primary_label = label_map[[cid_chr]],
              ontology_id = if (cid_chr %in% names(ontology_map)) ontology_map[[cid_chr]] else NULL,
              confidence = if (cid == 3) "Low" else "High",
              rationale = "Mock response",
              markers = pbmc$marker_list[[match(cid, pbmc$cluster_id)]]
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
  expect_equal(annotations$summary$total_clusters, length(pbmc$cluster_id))
  expect_equal(annotations$clusters$primary_label, unname(label_map))
})
