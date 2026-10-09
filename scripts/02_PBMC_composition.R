# PBMC annotation display and pooled composition.
library(Seurat)
pbmc <- readRDS("inputs/PBMC/TNK_annotated.rds")
DefaultAssay(pbmc) <- "RNA"
dir.create("results/PBMC/composition", recursive = TRUE, showWarnings = FALSE)

DimPlot(pbmc, reduction = "umap", group.by = "celltype", label = TRUE)
FeaturePlot(pbmc, features = "GZMK", reduction = "umap")
counts <- table(subset = pbmc$celltype, group = pbmc$Type)
percent <- 100 * prop.table(counts, margin = 2)
write.csv(as.data.frame.matrix(counts), "results/PBMC/composition/subset_counts.csv")
write.csv(as.data.frame.matrix(percent), "results/PBMC/composition/subset_percent.csv")
