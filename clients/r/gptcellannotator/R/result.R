#' Internal: normalise backend responses
#'
#' @keywords internal
gptca_parse_result <- function(result, validated = TRUE) {
  if (is.null(result)) {
    cli::cli_abort("Empty response received from GPT Cell Annotator.")
  }

  clusters_raw <- NULL
  summary <- NULL
  metrics <- NULL

  if (validated) {
    summary <- result$summary %||% NULL
    metrics <- result$metrics %||% NULL
    clusters_raw <- result$clusters %||% list()
  } else {
    clusters_raw <- list(result)
  }

  clusters_tbl <- gptca_tibble_from_clusters(clusters_raw)
  structure(
    list(
      clusters = clusters_tbl,
      summary = summary,
      metrics = metrics,
      validated = validated,
      raw = result
    ),
    class = c("gptca_annotation", "list")
  )
}

gptca_tibble_from_clusters <- function(clusters) {
  if (length(clusters) == 0) {
    return(tibble::tibble(
      cluster_id = character(),
      primary_label = character(),
      ontology_id = character(),
      confidence = character(),
      status = character(),
      rationale = character(),
      warnings = list(),
      validation = list(),
      markers = list(),
      alternatives = list(),
      raw = list()
    ))
  }

  tibble::tibble(
    cluster_id = vapply(
      clusters,
      function(cluster) {
        annotation <- cluster$annotation %||% cluster
        as.character(cluster$cluster_id %||% annotation$cluster_id %||% NA_character_)
      },
      "",
      USE.NAMES = FALSE
    ),
    primary_label = vapply(
      clusters,
      function(cluster) {
        annotation <- cluster$annotation %||% cluster
        annotation$primary_label %||% cluster$primary_label %||% NA_character_
      },
      "",
      USE.NAMES = FALSE
    ),
    ontology_id = vapply(
      clusters,
      function(cluster) {
        annotation <- cluster$annotation %||% cluster
        annotation$ontology_id %||% cluster$ontology_id %||% NA_character_
      },
      "",
      USE.NAMES = FALSE
    ),
    confidence = vapply(
      clusters,
      function(cluster) {
        annotation <- cluster$annotation %||% cluster
        cluster$confidence %||% annotation$confidence %||% NA_character_
      },
      "",
      USE.NAMES = FALSE
    ),
    status = vapply(
      clusters,
      function(cluster) {
        cluster$status %||% cluster$validation$status %||% NA_character_
      },
      "",
      USE.NAMES = FALSE
    ),
    rationale = vapply(
      clusters,
      function(cluster) {
        annotation <- cluster$annotation %||% cluster
        annotation$rationale %||% cluster$rationale %||% NA_character_
      },
      "",
      USE.NAMES = FALSE
    ),
    warnings = lapply(
      clusters,
      function(cluster) {
        metadata <- cluster[["metadata"]]
        warnings <- cluster$warnings %||% metadata$warnings %||% list()
        as_character_list(warnings)
      }
    ),
    validation = lapply(
      clusters,
      function(cluster) {
        metadata <- cluster[["metadata"]]
        cluster$validation %||% metadata$validation %||% NULL
      }
    ),
    markers = lapply(
      clusters,
      function(cluster) {
        annotation <- cluster$annotation %||% cluster
        annotation$markers %||% cluster$markers %||% character()
      }
    ),
    alternatives = lapply(
      clusters,
      function(cluster) {
        annotation <- cluster$annotation %||% cluster
        annotation$alternatives %||% list()
      }
    ),
    raw = lapply(clusters, identity)
  )
}

#' @export
print.gptca_annotation <- function(x, ...) {
  validated <- if (isTRUE(x$validated)) "validated" else "raw"
  cli::cli_text("{.strong GPT Cell Annotator annotations} ({validated})")
  clusters <- x$clusters
  if (nrow(clusters) == 0) {
    cli::cli_text("No clusters annotated.")
    return(invisible(x))
  }
  cli::cli_ul(vapply(
    seq_len(nrow(clusters)),
    function(i) {
      label <- clusters$primary_label[i] %||% "Unknown"
      cid <- clusters$cluster_id[i]
      confidence <- clusters$confidence[i] %||% ""
      status <- clusters$status[i] %||% ""
      sprintf("Cluster %s: %s (confidence: %s, status: %s)", cid, label, confidence, status)
    },
    "",
    USE.NAMES = FALSE
  ))
  invisible(x)
}

as_character_list <- function(x) {
  if (is.null(x)) {
    return(list())
  }
  if (is.character(x)) {
    return(as.list(x))
  }
  if (is.list(x)) {
    return(lapply(x, as.character))
  }
  as.list(as.character(x))
}
