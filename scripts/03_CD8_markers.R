# CD8 subset markers: Wilcoxon tests followed by Bonferroni correction.
library(Seurat)
cd8 <- readRDS("inputs/cd8/T.rds")
DefaultAssay(cd8) <- "RNA"
Idents(cd8) <- "celltype"
sampled <- subset(cd8, downsample = 200)
dir.create("results/markers", recursive = TRUE, showWarnings = FALSE)

naive <- FindMarkers(sampled, ident.1 = "Naive", ident.2 = c("Tem", "CTL"),
                      assay = "RNA", slot = "data", test.use = "wilcox", only.pos = FALSE,
                      min.pct = 0.1, logfc.threshold = 0)
tem <- FindMarkers(sampled, ident.1 = "Tem", ident.2 = c("Naive", "CTL"),
                    assay = "RNA", slot = "data", test.use = "wilcox", only.pos = FALSE,
                    min.pct = 0.1, logfc.threshold = 0)
ctl <- FindMarkers(sampled, ident.1 = "CTL", ident.2 = c("Naive", "Tem"),
                    assay = "RNA", slot = "data", test.use = "wilcox", only.pos = FALSE,
                    min.pct = 0.1, logfc.threshold = 0)

# Bonferroni correction within each comparison.
naive$p_val_adj <- p.adjust(naive$p_val, method = "bonferroni")
tem$p_val_adj <- p.adjust(tem$p_val, method = "bonferroni")
ctl$p_val_adj <- p.adjust(ctl$p_val, method = "bonferroni")
write.csv(naive, "results/markers/Naive_markers.csv")
write.csv(tem, "results/markers/Tem_markers.csv")
write.csv(ctl, "results/markers/CTL_markers.csv")
