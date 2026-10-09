"""Compare matched Monocle3 and diffusion pseudotimes."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

out = Path("results/diffusion")
monocle = pd.read_csv("results/monocle3/monocle_cells_marker_root.csv", dtype={"cell_id": str})
diffusion = pd.read_csv(out / "diffusion_cell_pseudotime.csv", dtype={"cell_id": str})
shared = monocle.merge(diffusion, on=["cell_id", "donor", "subset"], validate="one_to_one")
assert len(shared) == len(monocle), "All Monocle3 cells must match the DPT input."
finite = np.isfinite(shared["pseudotime"]) & np.isfinite(shared["dpt_pseudotime"])
rho = float(spearmanr(shared.loc[finite, "pseudotime"],
                      shared.loc[finite, "dpt_pseudotime"]).statistic)
shared["finite_pair"] = finite
shared.to_csv(out / "monocle3_diffusion_shared_cells.csv", index=False)
comparison = {"n_shared": len(shared), "n_finite_pairs": int(finite.sum()), "spearman_rho": rho}
(out / "monocle3_diffusion_comparison.json").write_text(
    json.dumps(comparison, indent=2) + "\n", encoding="utf-8")
