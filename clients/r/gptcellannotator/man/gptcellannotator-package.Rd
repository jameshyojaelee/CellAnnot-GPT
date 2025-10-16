% Generated manually until roxygen can be run
\name{gptcellannotator-package}
\alias{gptcellannotator}
\alias{gptcellannotator-package}
\docType{package}
\title{gptcellannotator: Connect Seurat workflows to GPT Cell Annotator}
\description{
An R interface to the GPT Cell Annotator backend. Provides helpers for sending marker tables, working with Seurat objects, parsing validation reports, and plotting annotations.
}
\details{
The package prefers HTTPS communication via the FastAPI service but can fall back to the `gca` CLI when offline. Set configuration defaults with \code{gptca_config()} and see \code{vignette("annotate-seurat")} for a complete workflow.
}
\author{GPT Cell Annotator Team <support@gpt-cell-annotator.org>}
\seealso{\code{gptca_config()}, \code{gptca_annotate_markers()}, \code{gptca_plot_umap()}}
