# Validation

The manuscript methods were checked against the source analysis scripts and the separate release scripts.

| Component | Check |
|---|---|
| Tissue integration | PBMC CCA, CLE RPCA and LN Harmony checked separately against their source workflows |
| Donor composition | Exact tests recalculated from full-PBMC annotations |
| Paired edgeR | All three contrasts rerun; all compared source statistics matched |
| Cell markers | Each subset sampled before comparison; Wilcoxon/Bonferroni and counts aggregation tested |
| Monocle3 main run | Pseudotimes compared cell by cell with saved results |
| DPT | Full analysis rerun and aligned to Monocle3 cells |
| LODO | Root selection/refitting checked against the original procedure and saved results |
| Velocity | Preprocessing, kinetic fit, selection thresholds and independent iteration fit checked separately |
| Tables | All four Excel exports checked cell by cell against their source values; percentage formulas, differential-expression adjustments and file checksums verified |

Sequencing preprocessing, all tissue integrations, all LODO fits and study-scale velocity dynamics were not rerun during packaging. Source fingerprints are in `SOURCE_PROVENANCE.csv`; release checksums are in the repository-root `FILE_MANIFEST.csv`.
