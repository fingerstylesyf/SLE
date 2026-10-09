"""Main dynamical RNA velocity from matched spliced/unspliced counts."""
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt

out = Path("results/velocity")
out.mkdir(parents=True, exist_ok=True)
src = ad.read_h5ad("inputs/cd8/velocity_raw_matched.h5ad")
data = ad.AnnData(X=src.layers["spliced"].copy(),
                  obs=src.obs[["Name", "barcode", "ID", "celltype"]].copy(),
                  var=pd.DataFrame(index=src.var_names.copy()))
data.layers["spliced"] = src.layers["spliced"].copy()
data.layers["unspliced"] = src.layers["unspliced"].copy()
data.obsm["X_umap"] = src.obsm["X_umap"].copy()

np.random.seed(20260919)
scv.pp.filter_and_normalize(data, min_shared_counts=20, n_top_genes=2000)
sc.pp.pca(data, n_comps=30, random_state=20260919)
sc.pp.neighbors(data, n_neighbors=30, n_pcs=30,
                use_rep="X_pca", random_state=20260919)
scv.pp.moments(data, n_neighbors=None, n_pcs=None)
data.write_h5ad(out / "velocity_prepared.h5ad")

np.random.seed(20260919)
scv.tl.recover_dynamics(data, var_names=list(data.var_names), max_iter=30, n_jobs=8)
assert "fit_alpha" in data.var and data.var["fit_alpha"].notna().any()
data.write_h5ad(out / "fitted_iter30.h5ad")
scv.tl.velocity(data, mode="dynamical", min_likelihood=0.1, min_r2=None)
scv.tl.velocity_graph(data, n_jobs=8)
scv.tl.velocity_embedding(data, basis="umap")
data.var[["fit_likelihood", "fit_scaling", "velocity_genes"]].to_csv(out / "main_genes.csv")
data.write_h5ad(out / "main_velocity.h5ad")
ax = scv.pl.velocity_embedding_stream(data, basis="umap", color="celltype",
                                      density=1.4, arrow_size=1.2, show=False)
ax.figure.savefig(out / "main_velocity.pdf", bbox_inches="tight")
plt.close(ax.figure)
