# Analysis and figure index

Figure numbers follow the revised manuscript.

| Figure | Script or input | Analysis |
|---|---|---|
| Figure 2A–C | Separate `01_PBMC`, `01_CLE`, `01_LN` scripts | CCA, RPCA and Harmony, respectively |
| Figure 2D–I | Separate `02_PBMC`, `02_CLE`, `02_LN` composition scripts | Annotated UMAP, GZMK expression and composition |
| Supplementary S1A–B | `02_PBMC_donor_comparison.R` | Donor-level exact comparisons |
| Supplementary S1C–E | `03_CD8_pseudobulk_counts.R`, `03_CD8_pseudobulk_edgeR.R` | Donor-paired subset contrasts |
| Figure 3A | `04_monocle3.R` | Marker-rooted pseudotime |
| Figure 3B | `06_rna_velocity.py` | Dynamical velocity |
| Figure 3C | `03_CD8_markers.R` | Cell-level subset markers |
| Figure 3D | `04_expression_curves.R` | Expression along pseudotime |
| Figure 3E | Annotated CD8 object | Granzyme and PRF1 expression on UMAP |
| Supplementary S2A | `05_diffusion_pseudotime.py` | Root-marker expression |
| Supplementary S2B | `04_monocle3.R` | Donor pseudotime medians |
| Supplementary S2C | `04_monocle3_LODO.R` | Leave-one-donor-out analysis |
| Supplementary S2D | `05_DPT_comparison.py` | DPT–Monocle3 comparison |
| Supplementary S2E | `tables/Figure_S2.xlsx` | Spliced/unspliced counts per donor |
| Velocity sensitivity | `06_velocity_thresholds.py`, `06_velocity_iter60.py` | Thresholds and independent iteration fit |
