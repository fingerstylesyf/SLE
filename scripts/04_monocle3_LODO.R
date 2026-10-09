# Leave-one-donor-out graph refitting on the original UMAP.
suppressPackageStartupMessages(library(monocle3))
out <- "results/monocle3_lodo"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cds <- readRDS("inputs/cd8/monocle3/cds.rds")
scores <- read.csv("results/diffusion/cell_metadata_root_scores.csv", row.names = "cell_id")
colData(cds)$naive_marker_score <- scores[colnames(cds), "naive_marker_score"]
main <- read.csv("results/monocle3/monocle_cells_marker_root.csv")
reference <- setNames(main$pseudotime, main$cell_id)
donors <- as.character(colData(cds)$Name)
results <- list()

for (excluded in sort(unique(donors))) {
  set.seed(20260919)
  x <- cds[, donors != excluded]
  x <- cluster_cells(x, reduction_method = "UMAP", random_seed = 20260919)
  x <- learn_graph(x, use_partition = TRUE, close_loop = TRUE,
                   learn_graph_control = list(minimal_branch_len = 10))
  meta <- as.data.frame(colData(x))
  naive <- meta$celltype == "Naive"
  cutoff <- quantile(meta$naive_marker_score[naive], 0.75)
  candidates <- rownames(meta)[naive & meta$naive_marker_score >= cutoff]
  near <- x@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  node <- as.character(near[candidates, 1])
  donor <- as.character(meta[candidates, "Name"])
  weight <- 1 / as.numeric(table(donor)[donor])
  selected <- names(sort(tapply(weight, node, sum), decreasing = TRUE))[1]
  root <- paste0("Y_", sub("^Y_", "", selected))
  x <- order_cells(x, root_pr_nodes = root)

  pt <- pseudotime(x)
  common <- intersect(names(pt), names(reference))
  keep <- is.finite(pt[common]) & is.finite(reference[common])
  rho <- cor(pt[common][keep], reference[common][keep], method = "spearman")
  cells <- data.frame(cell_id = colnames(x), donor = as.character(meta$Name),
                      subset = as.character(meta$celltype),
                      pseudotime = unname(pt[colnames(x)]))
  cells$finite <- is.finite(cells$pseudotime)
  medians <- aggregate(pseudotime ~ donor + subset, cells[cells$finite, ], median)
  names(medians)[3] <- "median_pseudotime"
  prefix <- paste0("leave_out_", make.names(excluded))
  write.csv(cells, file.path(out, paste0(prefix, "_cells.csv")), row.names = FALSE)
  write.csv(medians, file.path(out, paste0(prefix, "_donor_medians.csv")), row.names = FALSE)
  results[[excluded]] <- data.frame(excluded_donor = excluded, n_cells = ncol(x),
      root_node = root, n_root_candidates = length(candidates),
      finite_fraction = mean(cells$finite), n_finite_pairs = sum(keep), spearman_vs_main = rho)
}
write.csv(do.call(rbind, results), file.path(out, "monocle_lodo_summary.csv"), row.names = FALSE)
