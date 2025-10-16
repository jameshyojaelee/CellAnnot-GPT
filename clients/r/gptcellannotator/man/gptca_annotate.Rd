% Generated manually until roxygen can be run
\name{gptca_annotate_markers}
\alias{gptca_annotate_markers}
\alias{gptca_annotate_seurat}
\alias{gptca_add_metadata}
\alias{gptca_plot_umap}
\title{Annotate clusters with GPT Cell Annotator}
\description{
High-level helpers to send marker tables to GPT Cell Annotator, merge validated labels into Seurat metadata, and visualise annotations on UMAP embeddings.
}
\usage{
gptca_annotate_markers(
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
)

gptca_annotate_seurat(
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
)

gptca_add_metadata(
  object,
  annotation,
  cluster_column = "seurat_clusters",
  label_col = "gptca_label",
  confidence_col = "gptca_confidence",
  status_col = "gptca_status",
  warnings_col = "gptca_warnings"
)

gptca_plot_umap(
  object,
  label_col = "gptca_label",
  status_col = "gptca_status",
  reduction = "umap",
  flagged_status = c("flagged", "unknown"),
  point_size = 0.6,
  alpha = 0.9,
  palette = NULL
)
}
\arguments{
\item{markers}{Marker definitions: named list of gene symbols, a data frame with `cluster` and `gene` columns (e.g. \code{FindAllMarkers} output), or a data frame with `cluster_id` and list-column `markers`.}

\item{species}{Optional species hint included in the dataset context.}

\item{tissue}{Optional tissue hint.}

\item{return_validated}{When \code{TRUE} (default) request the validated report. Set to \code{FALSE} to retrieve raw LLM output for a single cluster.}

\item{marker_limit}{Maximum number of markers per cluster.}

\item{offline}{Logical flag forcing CLI fallback. Defaults to the active configuration.}

\item{fallback}{Attempt CLI fallback automatically when HTTP requests fail.}

\item{config}{A \code{GptcaConfig} object.}

\item{cluster_column}{Name of the column containing cluster IDs in marker data frames or Seurat metadata.}

\item{gene_column}{Name of the gene column in marker data frames.}

\item{markers_column}{Name of list-column storing marker vectors.}

\item{score_column}{Column used to rank markers (highest values first).}

\item{top_n}{Optional number of top markers to keep per cluster before applying \code{marker_limit}.}

\item{object}{A Seurat object.}

\item{label_col}{Metadata column used to store primary labels.}

\item{confidence_col}{Metadata column for confidence values.}

\item{status_col}{Metadata column storing validation status.}

\item{warnings_col}{Metadata column storing concatenated warnings.}

\item{annotation}{A \code{gptca_annotation} object returned by \code{gptca_annotate_markers()} or \code{gptca_annotate_seurat()}.}

\item{reduction}{Name of the reduction containing UMAP embeddings.}

\item{flagged_status}{Status values highlighted on the plot.}

\item{point_size}{Point size for scatter plot.}

\item{alpha}{Point alpha level.}

\item{palette}{Optional colour palette overriding defaults.}
}
\value{
\code{gptca_annotate_markers()} and \code{gptca_annotate_seurat()} return objects of class \code{gptca_annotation}. \code{gptca_add_metadata()} returns the modified Seurat object. \code{gptca_plot_umap()} returns a \pkg{ggplot2} object.
}
\examples{
## Not run: 
# markers <- read.csv(system.file("extdata", "pbmc_markers.csv", package = "gptcellannotator"))
# annotations <- gptca_annotate_markers(markers, species = "Homo sapiens")
# if (requireNamespace("Seurat", quietly = TRUE)) {
#   pbmc <- readRDS("path/to/pbmc.rds")
#   pbmc <- gptca_add_metadata(pbmc, annotations)
#   gptca_plot_umap(pbmc)
# }
## End(Not run)
}
