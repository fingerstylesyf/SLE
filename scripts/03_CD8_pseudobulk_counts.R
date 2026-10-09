# Sum RNA UMI counts within donor and CD8 subset.
library(Seurat)
library(Matrix)
cd8 <- readRDS("inputs/cd8/T.rds")
counts <- GetAssayData(cd8, assay = "RNA", layer = "counts")
md <- cd8[[]]
stopifnot(identical(colnames(counts), rownames(md)),
          all(counts@x >= 0), all(counts@x == round(counts@x)))

groups <- paste(md$Name, md$celltype, sep = "__")
profiles <- unique(groups)
assignment <- sparseMatrix(i = seq_along(groups), j = match(groups, profiles), x = 1,
                            dims = c(length(groups), length(profiles)))
pb_counts <- counts %*% assignment
colnames(pb_counts) <- profiles
samples <- md[match(profiles, groups), c("Name", "celltype"), drop = FALSE]
samples$pb_id <- profiles
samples$n_cells <- as.integer(table(factor(groups, levels = profiles)))
rownames(samples) <- profiles
dir.create("results/pseudobulk", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(counts = pb_counts, samples = samples),
        "results/pseudobulk/CD8_pseudobulk.rds")
