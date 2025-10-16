#' Plot GPT Cell Annotator labels on UMAP
#'
#' @param object A `Seurat` object containing UMAP embeddings.
#' @param label_col Metadata column with GPT Cell Annotator labels.
#' @param status_col Metadata column with validation status.
#' @param reduction Dimensionality reduction slot (default `"umap"`).
#' @param flagged_status Status values highlighted with outlines.
#' @param point_size Point size for scatter plot.
#' @param alpha Point alpha.
#' @param palette Optional vector of colours to override the default palette.
#' @return A `ggplot2` object.
#' @export
gptca_plot_umap <- function(
  object,
  label_col = "gptca_label",
  status_col = "gptca_status",
  reduction = "umap",
  flagged_status = c("flagged", "unknown"),
  point_size = 0.6,
  alpha = 0.9,
  palette = NULL
) {
  gptca_require_seurat()
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("`gptca_plot_umap()` requires {.pkg ggplot2}. Install it to enable plotting.")
  }

  meta <- object@meta.data
  if (!label_col %in% colnames(meta)) {
    cli::cli_abort("Label column {.field {label_col}} not found in Seurat metadata.")
  }

  embeddings <- Seurat::Embeddings(object, reduction = reduction)
  if (ncol(embeddings) < 2) {
    cli::cli_abort("Reduction {.field {reduction}} does not provide 2-D embeddings.")
  }

  plot_df <- data.frame(
    dim1 = embeddings[, 1],
    dim2 = embeddings[, 2],
    label = meta[[label_col]],
    status = if (status_col %in% colnames(meta)) meta[[status_col]] else NA_character_,
    cell = rownames(meta),
    stringsAsFactors = FALSE
  )

  flagged_cells <- plot_df$cell[plot_df$status %in% flagged_status]
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = dim1, y = dim2, colour = label)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::labs(
      title = "GPT Cell Annotator labels",
      subtitle = if (length(flagged_cells)) "Flagged clusters highlighted" else NULL,
      colour = label_col,
      x = paste0(toupper(reduction), " 1"),
      y = paste0(toupper(reduction), " 2")
    ) +
    ggplot2::theme_minimal()

  if (!is.null(palette)) {
    p <- p + ggplot2::scale_colour_manual(values = palette)
  }

  if (length(flagged_cells)) {
    flagged_df <- subset(plot_df, cell %in% flagged_cells)
    p <- p + ggplot2::geom_point(
      data = flagged_df,
      colour = "black",
      size = point_size + 0.6,
      shape = 21,
      stroke = 0.5,
      inherit.aes = FALSE,
      ggplot2::aes(x = dim1, y = dim2)
    )
  }

  p
}
