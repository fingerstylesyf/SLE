"""Four likelihood thresholds applied to the same 30-iteration kinetic fit."""
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scvelo as scv

out = Path("results/velocity/thresholds")
out.mkdir(parents=True, exist_ok=True)
fitted = ad.read_h5ad("results/velocity/fitted_iter30.h5ad")
assert "fit_alpha" in fitted.var and fitted.var["fit_alpha"].notna().any()

data = fitted.copy()
np.random.seed(20260919)
scv.tl.velocity(data, mode="dynamical", min_likelihood=0.05, min_r2=None)
scv.tl.velocity_graph(data, n_jobs=8)
scv.tl.velocity_embedding(data, basis="umap")
n005 = int(data.var.velocity_genes.sum())
data.write_h5ad(out / "likelihood_005.h5ad")

data = fitted.copy()
np.random.seed(20260919)
scv.tl.velocity(data, mode="dynamical", min_likelihood=0.10, min_r2=None)
scv.tl.velocity_graph(data, n_jobs=8)
scv.tl.velocity_embedding(data, basis="umap")
n010 = int(data.var.velocity_genes.sum())
data.write_h5ad(out / "likelihood_010.h5ad")

data = fitted.copy()
np.random.seed(20260919)
scv.tl.velocity(data, mode="dynamical", min_likelihood=0.20, min_r2=None)
scv.tl.velocity_graph(data, n_jobs=8)
scv.tl.velocity_embedding(data, basis="umap")
n020 = int(data.var.velocity_genes.sum())
data.write_h5ad(out / "likelihood_020.h5ad")

data = fitted.copy()
np.random.seed(20260919)
scv.tl.velocity(data, mode="dynamical", min_likelihood=0.30, min_r2=None)
scv.tl.velocity_graph(data, n_jobs=8)
scv.tl.velocity_embedding(data, basis="umap")
n030 = int(data.var.velocity_genes.sum())
data.write_h5ad(out / "likelihood_030.h5ad")

pd.DataFrame({"min_likelihood": [0.05, 0.10, 0.20, 0.30],
              "n_velocity_genes": [n005, n010, n020, n030]}).to_csv(
    out / "selected_gene_counts.csv", index=False)
