#' Annotate Seurat clusters via GPT Cell Annotator
#'
#' @param object A `Seurat` object.
#' @param markers Optional data frame from `Seurat::FindAllMarkers`. When
#'   omitted the function will compute positive markers using
#'   `Seurat::FindAllMarkers(object, only.pos = TRUE)`.
#' @param cluster_column Column name describing cluster IDs in `markers`.
#' @param gene_column Gene symbol column in `markers`.
#' @param score_column Score column used to rank markers.
#' @param top_n Keep the top-N genes per cluster before sending to the API.
#' @inheritParams gptca_annotate_markers
#'
#' @return A `gptca_annotation` object.
#' @export
gptca_annotate_seurat <- function(
  object,
  markers = NULL,
  species = NULL,
  tissue = NULL,
  return_validated = TRUE,
  marker_limit = 50,
  cluster_column = "cluster",
  gene_column = "gene",
  score_column = NULL,
  top_n = 25,
  offline = NULL,
  fallback = TRUE,
  config = gptca_config_get()
) {
  gptca_require_seurat()

  if (is.null(markers)) {
    cli::cli_inform("Computing markers via `Seurat::FindAllMarkers()`.")
    markers <- Seurat::FindAllMarkers(object, only.pos = TRUE, logfc.threshold = 0.25, verbose = FALSE)
  }

  gptca_annotate_markers(
    markers = markers,
    species = species,
    tissue = tissue,
    return_validated = return_validated,
    marker_limit = marker_limit,
    offline = offline,
    fallback = fallback,
    config = config,
    cluster_column = cluster_column,
    gene_column = gene_column,
    score_column = score_column,
    top_n = top_n
  )
}

#' Merge GPT Cell Annotator results into Seurat metadata
#'
#' @param object A `Seurat` object.
#' @param annotation A `gptca_annotation` object returned by
#'   [gptca_annotate_markers()] or [gptca_annotate_seurat()].
#' @param cluster_column Metadata column in `object` matching cluster IDs.
#' @param label_col Column name to store the primary label.
#' @param confidence_col Column for model confidence.
#' @param status_col Column for validation status (supported/flagged/...).
#' @param warnings_col Column storing concatenated warnings.
#'
#' @return The modified `Seurat` object.
#' @export
gptca_add_metadata <- function(
  object,
  annotation,
  cluster_column = "seurat_clusters",
  label_col = "gptca_label",
  confidence_col = "gptca_confidence",
  status_col = "gptca_status",
  warnings_col = "gptca_warnings"
) {
  gptca_require_seurat()
  if (!inherits(annotation, "gptca_annotation")) {
    cli::cli_abort("{.arg annotation} must be a {.cls gptca_annotation}.")
  }
  if (!cluster_column %in% colnames(object@meta.data)) {
    cli::cli_abort("Column {.field {cluster_column}} not found in Seurat metadata.")
  }

  clusters_tbl <- annotation$clusters
  cluster_ids <- clusters_tbl$cluster_id
  labels <- clusters_tbl$primary_label
  confidence <- clusters_tbl$confidence
  status <- clusters_tbl$status
  warnings <- vapply(
    clusters_tbl$warnings,
    function(x) if (length(x)) paste(unique(unlist(x)), collapse = "; ") else "",
    "",
    USE.NAMES = FALSE
  )

  meta <- object@meta.data
  key <- as.character(meta[[cluster_column]])
  idx <- match(key, cluster_ids)
  label_vec <- labels[idx]
  confidence_vec <- confidence[idx]
  status_vec <- status[idx]
  warnings_vec <- warnings[idx]

  new_meta <- data.frame(
    label_vec,
    confidence_vec,
    status_vec,
    warnings_vec,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  colnames(new_meta) <- c(label_col, confidence_col, status_col, warnings_col)
  rownames(new_meta) <- rownames(meta)
  object <- SeuratObject::AddMetaData(object, metadata = new_meta)
  object
}

gptca_require_seurat <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("SeuratObject", quietly = TRUE)) {
    cli::cli_abort(
      c(
        "Seurat integration requires the {.pkg Seurat} and {.pkg SeuratObject} packages.",
        "i" = "Install with `install.packages('Seurat')`."
      )
    )
  }
}
