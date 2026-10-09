# CLE skin annotation display and pooled composition.
library(Seurat)
cle <- readRDS("inputs/CLE/TNK_annotated.rds")
DefaultAssay(cle) <- "RNA"
dir.create("results/CLE/composition", recursive = TRUE, showWarnings = FALSE)

DimPlot(cle, reduction = "umap", group.by = "celltype", label = TRUE)
FeaturePlot(cle, features = "GZMK", reduction = "umap")
counts <- table(subset = cle$celltype, group = cle$Group)
percent <- 100 * prop.table(counts, margin = 2)
write.csv(as.data.frame.matrix(counts), "results/CLE/composition/subset_counts.csv")
write.csv(as.data.frame.matrix(percent), "results/CLE/composition/subset_percent.csv")
