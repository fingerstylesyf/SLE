# Inputs and running order

Run scripts from the repository root. Set input locations at the start of each script. Results are written under `results/`.

## Tissue datasets

| Dataset | Integration input | Composition input |
|---|---|---|
| PBMC | `inputs/PBMC/qc.rds`: QC-passed Seurat RNA object with `Name`, `percent.mt`, `percent.rb` | `inputs/PBMC/TNK_annotated.rds`: annotated object with `celltype`, `Type`, RNA data and UMAP |
| CLE | `inputs/CLE/qc.rds`: QC-passed Seurat RNA object with `Sample` | `inputs/CLE/TNK_annotated.rds`: annotated object with `celltype`, `Group`, RNA data and UMAP |
| LN | `inputs/LN/qc.rds`: QC-passed Seurat object with joined RNA assay and `Sample` | `inputs/LN/TNK_annotated.rds`: annotated object with `celltype`, `Group`, RNA data and UMAP |

Run the corresponding `01_*_integration.R` and `02_*_composition.R` scripts separately. Composition uses the annotated input objects. Annotation follows canonical markers and cluster-level expression.

PBMC donor comparisons use `inputs/PBMC/full_metadata.csv`, with one row per cell and columns `Name`, `Type` (`Con`/`SLE`) and `celltype`:

```sh
Rscript scripts/02_PBMC_donor_comparison.R
```

## CD8 expression

`inputs/cd8/T.rds` is the annotated CD8 Seurat object with a joined raw RNA counts layer, log-normalized RNA data, and `Name`/`celltype` metadata. Subsets are `Naive`, `Tem` and `CTL`.

```sh
Rscript scripts/03_CD8_pseudobulk_counts.R
Rscript scripts/03_CD8_pseudobulk_edgeR.R
```

The first script sums UMI counts within donor and subset. The second reads `results/pseudobulk/CD8_pseudobulk.rds` and fits the donor-paired edgeR model. An existing pseudobulk RDS can also be placed at this path; it must contain `counts` and `samples`, with sample columns `pb_id`, `Name`, `celltype`, and `n_cells`.

Run `Rscript scripts/03_CD8_markers.R` separately for cell-level markers. It downsamples each subset to at most 200 cells and applies Wilcoxon tests followed by Bonferroni correction.

## Trajectory

[TRAJECTORY_USAGE.md](docs/TRAJECTORY_USAGE.md) gives inputs and the sequence for DPT, Monocle3 and their comparison. LODO and expression curves have separate scripts.

## Velocity

`inputs/cd8/velocity_raw_matched.h5ad` contains barcode-matched raw `spliced` and `unspliced` counts, original `X_umap`, and `Name`, `barcode`, `ID`, `celltype` metadata.

```sh
python scripts/06_rna_velocity.py
```

The main script saves the prepared counts, the 30-iteration fit, velocities and a streamline plot. Run sensitivity analyses separately:

```sh
python scripts/06_velocity_thresholds.py
python scripts/06_velocity_iter60.py
```

Thresholds share the saved 30-iteration fit. The 60-iteration analysis starts from the prepared counts with fresh initialization.
