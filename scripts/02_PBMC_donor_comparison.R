# Exact donor-level comparisons from full-PBMC annotations.
md <- read.csv("inputs/PBMC/full_metadata.csv", row.names = 1, check.names = FALSE)
dir.create("results/PBMC/donors", recursive = TRUE, showWarnings = FALSE)
tnk_labels <- c("CD4_Naive", "CD4_Tcm", "Treg", "CD8_Naive", "CD8_Tem",
                "CD8_CTL", "MAIT", "gdT", "CD16_NK", "CD56_NK")

counts <- table(donor = md$Name, celltype = md$celltype)
proportions <- data.frame(
  donor = rownames(counts),
  group = c(Con = "HC", SLE = "SLE")[as.character(md$Type[match(rownames(counts), md$Name)])],
  total_PBMC_cells = rowSums(counts),
  TNK_cells = rowSums(counts[, tnk_labels, drop = FALSE]),
  Tem_cells = counts[, "CD8_Tem"])
proportions$Tem_pct_of_TNK <- 100 * proportions$Tem_cells / proportions$TNK_cells
proportions$Tem_pct_of_PBMC <- 100 * proportions$Tem_cells / proportions$total_PBMC_cells
hc <- proportions[proportions$group == "HC", ]
sle <- proportions[proportions$group == "SLE", ]
stopifnot(nrow(hc) == 6, nrow(sle) == 6,
          !anyDuplicated(proportions$Tem_pct_of_TNK),
          !anyDuplicated(proportions$Tem_pct_of_PBMC))

tnk_test <- wilcox.test(hc$Tem_pct_of_TNK, sle$Tem_pct_of_TNK,
                        alternative = "two.sided", exact = TRUE, correct = FALSE)
pbmc_test <- wilcox.test(hc$Tem_pct_of_PBMC, sle$Tem_pct_of_PBMC,
                         alternative = "two.sided", exact = TRUE, correct = FALSE)
comparisons <- data.frame(
  denominator = c("T/NK", "PBMC"),
  HC_mean_percent = c(mean(hc$Tem_pct_of_TNK), mean(hc$Tem_pct_of_PBMC)),
  SLE_mean_percent = c(mean(sle$Tem_pct_of_TNK), mean(sle$Tem_pct_of_PBMC)),
  U_HC = c(unname(tnk_test$statistic), unname(pbmc_test$statistic)),
  P_two_sided = c(tnk_test$p.value, pbmc_test$p.value))
write.csv(proportions, "results/PBMC/donors/Tem_proportions.csv", row.names = FALSE)
write.csv(comparisons, "results/PBMC/donors/Tem_group_comparisons.csv", row.names = FALSE)
