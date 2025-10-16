#' Internal: call the `gca` CLI for offline annotations
#'
#' @keywords internal
gptca_cli_available <- function(config = gptca_config_get()) {
  inherits(config, "GptcaConfig") && !is.null(config$cli_path) && file.exists(config$cli_path)
}

#' @keywords internal
gptca_cli_annotate <- function(clusters, dataset_context, config, offline = TRUE, marker_limit = NULL) {
  if (!gptca_cli_available(config)) {
    cli::cli_abort(
      c(
        "CLI fallback requested but {.code gca} binary not found.",
        "i" = "Install the Python package with `pip install gpt-cell-annotator` and ensure `gca` is on PATH."
      )
    )
  }

  if (length(clusters) == 0) {
    cli::cli_abort("No clusters supplied for CLI annotation.")
  }

  payload_df <- gptca_clusters_to_frame(clusters, marker_limit = marker_limit)
  input_path <- tempfile(fileext = ".csv")
  output_path <- tempfile(fileext = ".json")
  utils::write.csv(payload_df, input_path, row.names = FALSE, quote = TRUE)

  args <- c(
    "annotate",
    input_path,
    "--out-json",
    output_path
  )
  if (!is.null(dataset_context$species)) {
    args <- c(args, "--species", dataset_context$species)
  }
  if (!is.null(dataset_context$tissue)) {
    args <- c(args, "--tissue", dataset_context$tissue)
  }
  if (isTRUE(offline)) {
    args <- c(args, "--offline")
  }

  res <- system2(config$cli_path, args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(res, "status") %||% 0L
  if (!is.null(status) && status != 0) {
    cli::cli_abort(
      c(
        "CLI invocation failed with status {status}.",
        "x" = res
      )
    )
  }

  if (!file.exists(output_path)) {
    cli::cli_abort(
      c(
        "CLI invocation succeeded but JSON output missing.",
        "i" = "Ensure `gca annotate` supports `--out-json`."
      )
    )
  }

  jsonlite::fromJSON(output_path, simplifyVector = FALSE)
}

gptca_clusters_to_frame <- function(clusters, marker_limit = NULL) {
  tibble::tibble(
    cluster_id = vapply(clusters, function(cluster) cluster$cluster_id, "", USE.NAMES = FALSE),
    markers = vapply(
      clusters,
      function(cluster) {
        markers <- cluster$markers %||% character()
        if (!is.null(marker_limit)) {
          markers <- head(markers, marker_limit)
        }
        jsonlite::toJSON(markers, auto_unbox = TRUE)
      },
      "",
      USE.NAMES = FALSE
    )
  )
}
