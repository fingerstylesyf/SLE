"""Root-marker scores and diffusion pseudotime."""
from pathlib import Path
import json
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

out = Path("results/diffusion")
out.mkdir(parents=True, exist_ok=True)
source = sc.read_h5ad("inputs/cd8/T.h5ad")
donor = source.obs["Name"].astype(str).to_numpy()
subset = source.obs["celltype"].astype(str)
markers = ["CCR7", "TCF7", "LEF1", "SELL"]
expression = pd.DataFrame(source[:, markers].X.toarray(), index=source.obs_names, columns=markers)
scores = expression.rank(pct=True).mean(axis=1).to_numpy()
pd.DataFrame({"cell_id": source.obs_names, "Name": donor, "celltype": subset.to_numpy(),
              "naive_marker_score": scores}).to_csv(out / "cell_metadata_root_scores.csv", index=False)
means = expression.assign(celltype=subset).groupby("celltype", observed=True).mean()
means.to_csv(out / "root_marker_mean_expression.csv")
((means - means.mean()) / means.std(ddof=0)).to_csv(out / "root_marker_scaled_expression.csv")

pca = np.asarray(source.obsm["X_pca"][:, :30])
naive = subset.to_numpy() == "Naive"
candidates = np.flatnonzero(naive & (scores >= np.quantile(scores[naive], 0.75)))
candidate_pcs = pd.DataFrame(pca[candidates])
candidate_pcs["donor"] = donor[candidates]
donor_centroids = candidate_pcs.groupby("donor", sort=True).mean().to_numpy()
center = donor_centroids.mean(axis=0)
root = candidates[np.argmin(np.sum((pca[candidates] - center) ** 2, axis=1))]

data = ad.AnnData(obs=pd.DataFrame({"donor": donor, "subset": subset.to_numpy(),
                    "naive_marker_score": scores}, index=source.obs_names.copy()))
data.obsm["X_test"] = pca.copy()
data.uns["iroot"] = int(root)
sc.pp.neighbors(data, n_neighbors=30, use_rep="X_test", random_state=20260919)
sc.tl.diffmap(data, n_comps=15, random_state=20260919)
sc.tl.dpt(data, n_dcs=10, n_branchings=0)

cells = data.obs.copy()
cells.insert(0, "cell_id", data.obs_names)
cells["finite"] = np.isfinite(cells["dpt_pseudotime"])
cells.to_csv(out / "diffusion_cell_pseudotime.csv", index=False)
cells.groupby(["donor", "subset"], observed=True)["dpt_pseudotime"].median().to_csv(
    out / "diffusion_donor_medians.csv")
settings = {"n_pcs": 30, "n_neighbors": 30, "n_diffusion_components": 15, "n_dcs": 10,
            "seed": 20260919, "root_cell": str(data.obs_names[root]),
            "root_donor": str(donor[root]), "root_candidates": len(candidates),
            "n_cells": len(cells), "n_finite": int(cells["finite"].sum())}
(out / "diffusion_run_settings.json").write_text(json.dumps(settings, indent=2) + "\n", encoding="utf-8")
