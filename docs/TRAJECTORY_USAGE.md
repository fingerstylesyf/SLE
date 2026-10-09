# Trajectory analysis

Inputs:

- `inputs/cd8/T.h5ad`: full CD8 object with sparse log-normalized RNA in `X`, gene symbols, integrated `X_pca`, and `Name`/`celltype` metadata.
- `inputs/cd8/monocle3/cds.rds`: fitted 3,000-cell Monocle3 graph with the imported UMAP and `Name`/`celltype` metadata.

Run from the repository root:

```sh
python scripts/05_diffusion_pseudotime.py
Rscript scripts/04_monocle3.R
python scripts/05_DPT_comparison.py
```

The DPT script first calculates marker scores across the full CD8 object and saves them in `results/diffusion/cell_metadata_root_scores.csv`. The Monocle3 script reads these scores, selects the root on the fitted graph, and writes cell pseudotime and donor medians. The comparison script aligns cell identifiers and annotations before calculating Spearman correlation.

Run the two additional analyses independently:

```sh
Rscript scripts/04_monocle3_LODO.R
Rscript scripts/04_expression_curves.R
```

LODO refits clustering, the graph and root selection once per excluded donor. Expression curves use the supplied original graph and its original naive-region root.
