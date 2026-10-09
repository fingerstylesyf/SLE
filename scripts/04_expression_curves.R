# Expression curves from the original naive-region root.
suppressPackageStartupMessages(library(monocle3))
out <- "results/expression_curves"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cds <- readRDS("inputs/cd8/monocle3/cds.rds")
cds <- order_cells(cds, root_pr_nodes = "Y_2")
genes <- c("GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1")
pdf(file.path(out, "expression_pseudotime.pdf"), width = 10, height = 9)
print(plot_genes_in_pseudotime(cds[genes, ], color_cells_by = "celltype",
                               min_expr = 0.5, ncol = 2))
dev.off()
