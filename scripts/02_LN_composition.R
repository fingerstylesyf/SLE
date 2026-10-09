# LN kidney annotation display and pooled composition.
library(Seurat)
ln <- readRDS("inputs/LN/TNK_annotated.rds")
DefaultAssay(ln) <- "RNA"
dir.create("results/LN/composition", recursive = TRUE, showWarnings = FALSE)

DimPlot(ln, reduction = "umap", group.by = "celltype", label = TRUE)
FeaturePlot(ln, features = "GZMK", reduction = "umap")
counts <- table(subset = ln$celltype, group = ln$Group)
percent <- 100 * prop.table(counts, margin = 2)
write.csv(as.data.frame.matrix(counts), "results/LN/composition/subset_counts.csv")
write.csv(as.data.frame.matrix(percent), "results/LN/composition/subset_percent.csv")
