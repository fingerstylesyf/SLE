# Monocle3 marker-rooted trajectory.
suppressPackageStartupMessages(library(monocle3))
set.seed(20260919)
out <- "results/monocle3"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

cds <- readRDS("inputs/cd8/monocle3/cds.rds")
scores <- read.csv("results/diffusion/cell_metadata_root_scores.csv", row.names = "cell_id")
colData(cds)$naive_marker_score <- scores[colnames(cds), "naive_marker_score"]
meta <- as.data.frame(colData(cds))
naive <- meta$celltype == "Naive"
cutoff <- quantile(meta$naive_marker_score[naive], 0.75)
candidates <- rownames(meta)[naive & meta$naive_marker_score >= cutoff]
near <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
node <- as.character(near[candidates, 1])
donor <- as.character(meta[candidates, "Name"])
weight <- 1 / as.numeric(table(donor)[donor])
selected <- names(sort(tapply(weight, node, sum), decreasing = TRUE))[1]
root <- paste0("Y_", sub("^Y_", "", selected))
cds <- order_cells(cds, root_pr_nodes = root)

pt <- pseudotime(cds)
cells <- data.frame(cell_id = colnames(cds), donor = as.character(meta$Name),
                    subset = as.character(meta$celltype),
                    pseudotime = unname(pt[colnames(cds)]))
cells$finite <- is.finite(cells$pseudotime)
medians <- aggregate(pseudotime ~ donor + subset, cells[cells$finite, ], median)
names(medians)[3] <- "median_pseudotime"
roots <- data.frame(cell_id = candidates, donor = donor, closest_node = node,
                    weight = weight, selected_root = root)
write.csv(cells, file.path(out, "monocle_cells_marker_root.csv"), row.names = FALSE)
write.csv(medians, file.path(out, "monocle_donor_medians_marker_root.csv"), row.names = FALSE)
write.csv(roots, file.path(out, "monocle_roots_marker_root.csv"), row.names = FALSE)
write.csv(data.frame(root_node = root, n_root_candidates = length(candidates),
                    finite_fraction = mean(cells$finite)),
          file.path(out, "monocle_summary.csv"), row.names = FALSE)
