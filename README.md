# SLE single-cell analysis

Core analysis steps and source tables for Figures 2–3 and Supplementary Figures S1–S2.

PBMC, CLE skin and LN kidney have separate scripts. Each script follows one analysis from its input to its output. Integration starts from quality-controlled objects; PBMC doublets were removed using standard DoubletFinder.

| Analysis | Scripts |
|---|---|
| PBMC | `01_PBMC_integration.R`, `02_PBMC_composition.R`, `02_PBMC_donor_comparison.R` |
| CLE skin | `01_CLE_integration.R`, `02_CLE_composition.R` |
| LN kidney | `01_LN_integration.R`, `02_LN_composition.R` |
| CD8 markers | `03_CD8_markers.R` |
| Paired pseudobulk | `03_CD8_pseudobulk_counts.R`, `03_CD8_pseudobulk_edgeR.R` |
| Monocle3 | `04_monocle3.R`, `04_monocle3_LODO.R`, `04_expression_curves.R` |
| DPT | `05_diffusion_pseudotime.py`, `05_DPT_comparison.py` |
| RNA velocity | `06_rna_velocity.py`, `06_velocity_thresholds.py`, `06_velocity_iter60.py` |

All scripts are in `scripts/`. [INPUTS.md](INPUTS.md) lists the processed inputs and running order. Large matrices and raw sequencing files are supplied separately.

Four Excel workbooks contain Figure 2 cell counts and percentages, Figure 3 cell-level differential-expression results, Supplementary Figure S1 counts and paired pseudobulk results, and Supplementary Figure S2 cell-count summaries. Pseudobulk contrasts compare each CD8 subset with the equally weighted mean of the other two on the log-expression scale, accounting for donor pairing.

See the [figure index](docs/ANALYSIS_INDEX.md), [table list](docs/TABLE_INDEX.md), [software versions](docs/SOFTWARE.md) and [validation record](docs/VALIDATION.md).
