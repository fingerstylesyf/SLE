# Run from the repository root. Input: a QC-passed Seurat RNA object.
# PBMC preprocessing includes standard DoubletFinder filtering.
# Required metadata: Name, percent.mt and percent.rb.
input_file <- "inputs/PBMC/qc.rds"
output_file <- "results/PBMC/integrated.rds"

library(Seurat)

pbmc <- readRDS(input_file)
DefaultAssay(pbmc) <- "RNA"
pbmc_samples <- SplitObject(pbmc, split.by = "Name")
pbmc_samples <- lapply(pbmc_samples, function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", verbose = FALSE)
  FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000,
                       verbose = FALSE)
})

pbmc_features <- SelectIntegrationFeatures(pbmc_samples, nfeatures = 2000)
pbmc_anchors <- FindIntegrationAnchors(
  object.list = pbmc_samples, anchor.features = pbmc_features,
  reduction = "cca", dims = 1:30)
pbmc <- IntegrateData(pbmc_anchors, dims = 1:30)

DefaultAssay(pbmc) <- "integrated"
pbmc <- ScaleData(pbmc, features = rownames(pbmc),
                  vars.to.regress = c("percent.mt", "percent.rb"), verbose = FALSE)
pbmc <- RunPCA(pbmc, npcs = 30, verbose = FALSE)
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:15)
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 1.3)
DefaultAssay(pbmc) <- "RNA"

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(pbmc, output_file)
