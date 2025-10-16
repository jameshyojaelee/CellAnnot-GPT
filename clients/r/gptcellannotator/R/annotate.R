#' Annotate marker lists with GPT Cell Annotator
#'
#' @param markers Marker definitions. Accepts a named list of gene symbols,
#'   a data frame with `cluster`/`gene` columns (e.g. `FindAllMarkers`
#'   output), or a data frame with `cluster_id` and list-column `markers`.
#' @param species Optional species hint passed to the backend.
#' @param tissue Optional tissue hint passed to the backend.
#' @param return_validated When `TRUE` (default) returns the validated report
#'   produced by `/annotate_batch` or `return_validated = TRUE` on the cluster
#'   endpoint. Set to `FALSE` to access the raw LLM output for a single cluster.
#' @param marker_limit Maximum number of markers per cluster sent to the API.
#' @param offline Force the CLI fallback instead of HTTP calls.
#' @param fallback Attempt CLI fallback automatically if HTTP requests fail.
#' @param config An explicit [gptca_config()] object.
#' @param cluster_column Column containing cluster identifiers when `markers`
#'   is a data frame.
#' @param gene_column Column containing gene symbols when `markers` is a data
#'   frame of marker statistics.
#' @param markers_column Optional column with list-columns of markers per
#'   cluster.
#' @param score_column Optional column controlling ordering of markers (highest
#'   values are kept first). If omitted an appropriate column is guessed.
#' @param top_n Optional limit applied before `marker_limit` to keep the top-N
#'   genes per cluster.
#' @return A `gptca_annotation` object with tidy cluster annotations, summary,
#'   and raw JSON payload.
#' @export
gptca_annotate_markers <- function(
  markers,
  species = NULL,
  tissue = NULL,
  return_validated = TRUE,
  marker_limit = 50,
  offline = NULL,
  fallback = TRUE,
  config = gptca_config_get(),
  cluster_column = NULL,
  gene_column = NULL,
  markers_column = NULL,
  score_column = NULL,
  top_n = NULL
) {
  config <- config %||% gptca_config_get()
  if (!inherits(config, "GptcaConfig")) {
    cli::cli_abort("{.arg config} must be a {.cls GptcaConfig}.")
  }
  if (is.null(offline)) {
    offline <- config$offline
  }

  clusters <- gptca_prepare_clusters(
    markers = markers,
    cluster_column = cluster_column,
    gene_column = gene_column,
    markers_column = markers_column,
    score_column = score_column,
    top_n = top_n,
    marker_limit = marker_limit
  )

  if (length(clusters) == 0) {
    cli::cli_abort("No marker definitions detected.")
  }

  dataset_context <- gptca_dataset_context(species = species, tissue = tissue)

  if (!return_validated && length(clusters) > 1) {
    cli::cli_abort("`return_validated = FALSE` is only supported for a single cluster.")
  }

  if (!isTRUE(offline)) {
    result <- try(
      gptca_request_backend(clusters, dataset_context, config, return_validated),
      silent = TRUE
    )
    if (!inherits(result, "try-error")) {
      return(result)
    }
    if (isTRUE(fallback) && gptca_cli_available(config)) {
      cli::cli_warn(
        c(
          "HTTP request failed, attempting CLI fallback.",
          "i" = conditionMessage(attr(result, "condition"))
        )
      )
    } else {
      stop(conditionMessage(attr(result, "condition")), call. = FALSE)
    }
  }

  cli_result <- gptca_cli_annotate(
    clusters = clusters,
    dataset_context = dataset_context,
    config = config,
    offline = TRUE,
    marker_limit = marker_limit
  )
  gptca_parse_result(cli_result, validated = TRUE)
}

gptca_request_backend <- function(clusters, dataset_context, config, return_validated) {
  if (length(clusters) == 1L) {
    body <- list(
      cluster = clusters[[1]],
      dataset_context = dataset_context,
      return_validated = isTRUE(return_validated)
    )
    response <- gptca_http_post("/annotate_cluster", body = body, config = config)
    payload <- response$result %||% response
    gptca_parse_result(payload, validated = isTRUE(return_validated))
  } else {
    body <- list(
      clusters = clusters,
      dataset_context = dataset_context
    )
    response <- gptca_http_post("/annotate_batch", body = body, config = config)
    payload <- response$result %||% response
    gptca_parse_result(payload, validated = TRUE)
  }
}

gptca_dataset_context <- function(species = NULL, tissue = NULL) {
  context <- list()
  if (!is.null(species) && nzchar(species)) {
    context$species <- species
  }
  if (!is.null(tissue) && nzchar(tissue)) {
    context$tissue <- tissue
  }
  if (length(context) == 0) {
    NULL
  } else {
    context
  }
}

gptca_prepare_clusters <- function(
  markers,
  cluster_column = NULL,
  gene_column = NULL,
  markers_column = NULL,
  score_column = NULL,
  top_n = NULL,
  marker_limit = NULL
) {
  if (is.null(markers)) {
    cli::cli_abort("Markers input is required.")
  }

  if (is.list(markers) && !is.data.frame(markers)) {
    return(gptca_prepare_from_list(markers, marker_limit = marker_limit))
  }

  if (!is.data.frame(markers)) {
    cli::cli_abort("Unsupported marker input of class {.cls {class(markers)}}.")
  }

  if (!is.null(markers_column) && markers_column %in% names(markers)) {
    return(gptca_prepare_from_list_column(
      markers,
      markers_column = markers_column,
      cluster_column = cluster_column,
      marker_limit = marker_limit
    ))
  }

  cluster_column <- cluster_column %||% gptca_guess_cluster_column(names(markers))
  gene_column <- gene_column %||% gptca_guess_gene_column(names(markers))
  if (is.null(cluster_column) || is.null(gene_column)) {
    cli::cli_abort(
      "Unable to guess {.field cluster} or {.field gene} columns. Supply `cluster_column` and `gene_column`."
    )
  }

  score_column <- score_column %||% gptca_guess_score_column(names(markers))
  gptca_prepare_from_long(
    markers,
    cluster_column = cluster_column,
    gene_column = gene_column,
    score_column = score_column,
    top_n = top_n,
    marker_limit = marker_limit
  )
}

gptca_prepare_from_list <- function(markers, marker_limit = NULL) {
  ids <- names(markers)
  if (is.null(ids)) {
    ids <- as.character(seq_along(markers))
  }
  purrr::imap(
    markers,
    function(gene_list, cluster_id) {
      genes <- gptca_normalise_genes(gene_list, marker_limit = marker_limit)
      list(cluster_id = as.character(cluster_id), markers = genes)
    }
  )
}

gptca_prepare_from_list_column <- function(markers, markers_column, cluster_column, marker_limit = NULL) {
  if (is.null(cluster_column)) {
    cluster_column <- gptca_guess_cluster_column(names(markers))
  }
  if (is.null(cluster_column) || !cluster_column %in% names(markers)) {
    cli::cli_abort("Column with cluster identifiers not found.")
  }
  purrr::map2(
    markers[[cluster_column]],
    markers[[markers_column]],
    function(cid, genes) {
      gene_vec <- gptca_normalise_genes(genes, marker_limit = marker_limit)
      list(cluster_id = as.character(cid), markers = gene_vec)
    }
  )
}

gptca_prepare_from_long <- function(markers, cluster_column, gene_column, score_column, top_n, marker_limit) {
  split_markers <- split(markers, markers[[cluster_column]])
  purrr::imap(
    split_markers,
    function(cluster_df, cluster_id) {
      if (!is.null(score_column) && score_column %in% names(cluster_df)) {
        ord <- order(cluster_df[[score_column]], decreasing = TRUE, na.last = NA)
        cluster_df <- cluster_df[ord, , drop = FALSE]
      }
      genes <- cluster_df[[gene_column]]
      genes <- gptca_normalise_genes(genes, marker_limit = marker_limit, top_n = top_n)
      list(cluster_id = as.character(cluster_id), markers = genes)
    }
  )
}

gptca_guess_cluster_column <- function(columns) {
  candidates <- c("cluster_id", "cluster", "ident", "group", "clusterID")
  intersect(candidates, columns)[1]
}

gptca_guess_gene_column <- function(columns) {
  candidates <- c("gene", "feature", "symbol", "marker")
  intersect(candidates, columns)[1]
}

gptca_guess_score_column <- function(columns) {
  candidates <- c("avg_log2FC", "avg_logFC", "avg_log2fc", "avg_logfc", "avg_logFC")
  intersect(candidates, columns)[1]
}

gptca_normalise_genes <- function(genes, marker_limit = NULL, top_n = NULL) {
  if (is.list(genes) && !is.vector(genes)) {
    genes <- unlist(genes, recursive = TRUE, use.names = FALSE)
  }
  genes <- as.character(genes)
  genes <- genes[nzchar(genes)]
  genes <- toupper(trimws(genes))
  genes <- unique(genes)
  if (!is.null(top_n)) {
    genes <- head(genes, top_n)
  }
  if (!is.null(marker_limit)) {
    genes <- head(genes, marker_limit)
  }
  genes
}
