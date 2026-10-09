# edgeR 4.4.2: donor-paired subset contrasts.
library(edgeR)
pb <- readRDS("results/pseudobulk/CD8_pseudobulk.rds")
samples <- pb$samples
samples <- samples[samples$celltype %in% c("Naive", "Tem", "CTL") &
                     samples$n_cells >= 20, ]
paired_donors <- names(which(table(samples$Name) == 3))
samples <- samples[samples$Name %in% paired_donors, ]
samples$Name <- factor(samples$Name)
samples$celltype <- factor(samples$celltype, levels = c("Naive", "Tem", "CTL"))
counts <- pb$counts[, match(samples$pb_id, colnames(pb$counts)), drop = FALSE]
stopifnot(all(table(samples$Name, samples$celltype) == 1),
          identical(colnames(counts), samples$pb_id))

design <- model.matrix(~ 0 + celltype + Name, samples)
stopifnot(qr(design)$rank == ncol(design))
y <- DGEList(counts)
y <- y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)

naive_contrast <- c(1, -0.5, -0.5, rep(0, ncol(design) - 3))
tem_contrast <- c(-0.5, 1, -0.5, rep(0, ncol(design) - 3))
ctl_contrast <- c(-0.5, -0.5, 1, rep(0, ncol(design) - 3))
naive_test <- glmQLFTest(fit, contrast = naive_contrast)
tem_test <- glmQLFTest(fit, contrast = tem_contrast)
ctl_test <- glmQLFTest(fit, contrast = ctl_contrast)
naive <- topTags(naive_test, n = Inf, adjust.method = "BH")$table
tem <- topTags(tem_test, n = Inf, adjust.method = "BH")$table
ctl <- topTags(ctl_test, n = Inf, adjust.method = "BH")$table

dir.create("results/pseudobulk", recursive = TRUE, showWarnings = FALSE)
write.csv(naive, "results/pseudobulk/Naive_DE.csv")
write.csv(tem, "results/pseudobulk/Tem_DE.csv")
write.csv(ctl, "results/pseudobulk/CTL_DE.csv")
