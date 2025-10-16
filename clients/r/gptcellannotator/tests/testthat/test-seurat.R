test_that("gptca_add_metadata merges annotations into Seurat object", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  counts <- matrix(rpois(200, lambda = 5), nrow = 20, dimnames = list(paste0("gene", 1:20), paste0("cell", 1:10)))
  obj <- Seurat::CreateSeuratObject(counts)
  obj$seurat_clusters <- rep(c("0", "1"), length.out = ncol(obj))

  annotation <- structure(
    list(
      clusters = tibble::tibble(
        cluster_id = c("0", "1"),
        primary_label = c("B cell", "T cell"),
        ontology_id = c("CL:0000236", "CL:0000625"),
        confidence = c("High", "Medium"),
        status = c("supported", "flagged"),
        rationale = c("Mock rationale", "Mock rationale"),
        warnings = list(character(), list("Potential doublet")),
        validation = list(
          list(is_supported = TRUE),
          list(is_supported = FALSE)
        ),
        markers = list(c("MS4A1", "CD79A"), c("CD3E")),
        alternatives = list(list(), list())
      ),
      summary = list(),
      metrics = NULL,
      validated = TRUE,
      raw = list()
    ),
    class = c("gptca_annotation", "list")
  )

  updated <- gptca_add_metadata(obj, annotation)
  meta_cols <- c("gptca_label", "gptca_confidence", "gptca_status", "gptca_warnings")
  expect_true(all(meta_cols %in% colnames(updated@meta.data)))
  expect_equal(unique(updated$gptca_label), c("B cell", "T cell"))
})

test_that("gptca_plot_umap returns ggplot", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("ggplot2")
  counts <- matrix(rpois(200, lambda = 5), nrow = 20, dimnames = list(paste0("gene", 1:20), paste0("cell", 1:10)))
  obj <- Seurat::CreateSeuratObject(counts)
  obj$seurat_clusters <- rep(c("0", "1"), length.out = ncol(obj))

  annotation <- structure(
    list(
      clusters = tibble::tibble(
        cluster_id = c("0", "1"),
        primary_label = c("B cell", "T cell"),
        ontology_id = c("CL:0000236", "CL:0000625"),
        confidence = c("High", "Medium"),
        status = c("supported", "flagged"),
        rationale = c("Mock rationale", "Mock rationale"),
        warnings = list(character(), list("Potential doublet")),
        validation = list(
          list(is_supported = TRUE),
          list(is_supported = FALSE)
        ),
        markers = list(c("MS4A1", "CD79A"), c("CD3E")),
        alternatives = list(list(), list())
      ),
      summary = list(),
      metrics = NULL,
      validated = TRUE,
      raw = list()
    ),
    class = c("gptca_annotation", "list")
  )

  obj <- gptca_add_metadata(obj, annotation)
  embeddings <- matrix(runif(20), ncol = 2, dimnames = list(colnames(obj), c("UMAP_1", "UMAP_2")))
  obj@reductions$umap <- Seurat::CreateDimReducObject(
    embeddings = embeddings,
    key = "UMAP_",
    assay = Seurat::DefaultAssay(obj)
  )

  plot <- gptca_plot_umap(obj)
  expect_s3_class(plot, "ggplot")
})
