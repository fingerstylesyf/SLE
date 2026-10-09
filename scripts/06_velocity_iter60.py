"""Independent 60-iteration sensitivity fit from the prepared counts."""
from pathlib import Path
import anndata as ad
import numpy as np
import scvelo as scv
import matplotlib.pyplot as plt

out = Path("results/velocity/iteration60")
out.mkdir(parents=True, exist_ok=True)
data = ad.read_h5ad("results/velocity/velocity_prepared.h5ad")
np.random.seed(20260919)
scv.tl.recover_dynamics(data, var_names=list(data.var_names),
                        max_iter=60, load_pars=False, n_jobs=8)
assert "fit_alpha" in data.var and data.var["fit_alpha"].notna().any()
scv.tl.velocity(data, mode="dynamical", min_likelihood=0.1, min_r2=None)
scv.tl.velocity_graph(data, n_jobs=8)
scv.tl.velocity_embedding(data, basis="umap")
data.var[["fit_likelihood", "fit_scaling", "velocity_genes"]].to_csv(out / "genes.csv")
data.write_h5ad(out / "velocity_iter60.h5ad")
ax = scv.pl.velocity_embedding_stream(data, basis="umap", color="celltype",
                                      density=1.4, arrow_size=1.2, show=False)
ax.figure.savefig(out / "velocity_iter60.pdf", bbox_inches="tight")
plt.close(ax.figure)
