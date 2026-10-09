# Run from the repository root. Input: a QC-passed Seurat object with a joined RNA assay.
# Required metadata: Sample (donor identifier).
input_file <- "inputs/LN/qc.rds"
output_file <- "results/LN/integrated.rds"

library(Seurat)
library(harmony)
set.seed(1234)

ln <- readRDS(input_file)
DefaultAssay(ln) <- "RNA"
ln <- NormalizeData(ln, normalization.method = "LogNormalize", verbose = FALSE)
ln <- FindVariableFeatures(ln, selection.method = "vst", nfeatures = 2000,
                          verbose = FALSE)
ln <- ScaleData(ln, features = VariableFeatures(ln), verbose = FALSE)
ln <- RunPCA(ln, features = VariableFeatures(ln), npcs = 35, verbose = FALSE)

ln <- RunHarmony(ln, group.by.vars = "Sample", reduction.use = "pca",
                 dims.use = 1:35, reduction.save = "harmony")
ln <- RunUMAP(ln, reduction = "harmony", dims = 1:35, seed.use = 1234)
ln <- FindNeighbors(ln, reduction = "harmony", dims = 1:35)
ln <- FindClusters(ln, resolution = 0.4)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(ln, output_file)
