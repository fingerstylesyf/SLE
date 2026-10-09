# Run from the repository root. Input: a QC-passed Seurat RNA object.
# Required metadata: Sample (donor identifier).
input_file <- "inputs/CLE/qc.rds"
output_file <- "results/CLE/integrated.rds"

library(Seurat)
set.seed(1234)

cle <- readRDS(input_file)
DefaultAssay(cle) <- "RNA"
cle_samples <- SplitObject(cle, split.by = "Sample")
cle_samples <- lapply(cle_samples, function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", verbose = FALSE)
  FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000,
                       verbose = FALSE)
})

cle_features <- SelectIntegrationFeatures(cle_samples, nfeatures = 2000)
cle_samples <- lapply(cle_samples, function(x) {
  x <- ScaleData(x, features = cle_features, verbose = FALSE)
  RunPCA(x, features = cle_features, npcs = 35, verbose = FALSE)
})
cle_anchors <- FindIntegrationAnchors(
  object.list = cle_samples, anchor.features = cle_features,
  reduction = "rpca", dims = 1:35, k.anchor = 5)
cle <- IntegrateData(cle_anchors, dims = 1:35)

DefaultAssay(cle) <- "integrated"
cle <- ScaleData(cle, verbose = FALSE)
cle <- RunPCA(cle, npcs = 35, verbose = FALSE)
cle <- RunUMAP(cle, reduction = "pca", dims = 1:20)
cle <- FindNeighbors(cle, reduction = "pca", dims = 1:20)
cle <- FindClusters(cle, resolution = 0.4)
DefaultAssay(cle) <- "RNA"

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(cle, output_file)
