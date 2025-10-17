#!/usr/bin/env Rscript
# QA harness: exercise gptca_annotate_seurat against the REST backend.

suppressPackageStartupMessages({
  library(gptcellannotator)
  library(Seurat)
  library(jsonlite)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else file.path("docs", "reports", "seurat_integration", "artifacts")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

api_key <- Sys.getenv("GPTCA_API_KEY", unset = "")
if (!nzchar(api_key)) {
  stop("Set GPTCA_API_KEY before running the REST QA harness.", call. = FALSE)
}

base_url <- Sys.getenv("GPTCA_BASE_URL", unset = "https://api.gpt-cell-annotator.org")

cfg <- gptca_config(
  base_url = base_url,
  api_key = api_key,
  offline = FALSE,
  timeout = 120,
  retry_max = 3
)
gptca_config_set(cfg)

pbmc <- Seurat::pbmc_small
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, npcs = 10, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE)

markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.diff.pct = 0.1
)

annotations <- gptca_annotate_seurat(
  pbmc,
  markers = markers,
  species = "Homo sapiens",
  tissue = "Peripheral blood",
  return_validated = TRUE,
  top_n = 15
)

if (nrow(annotations$clusters) == 0) {
  stop("No clusters returned from REST annotation.", call. = FALSE)
}

clusters_path <- file.path(out_dir, "rest_clusters.csv")
utils::write.csv(annotations$clusters, clusters_path, row.names = FALSE)

summary_path <- file.path(out_dir, "rest_summary.json")
write_json(annotations$summary, summary_path, auto_unbox = TRUE, pretty = TRUE)

raw_path <- file.path(out_dir, "rest_annotations.json")
write_json(annotations$raw, raw_path, auto_unbox = TRUE, pretty = TRUE)

plot_path <- file.path(out_dir, "rest_umap.png")
umap_plot <- gptca_plot_umap(pbmc)
ggsave(plot_path, umap_plot, width = 6.5, height = 5, dpi = 180)

message("REST QA completed. Artefacts saved in ", normalizePath(out_dir))
