# Tables

| Workbook | Sheets | Content |
|---|---|---|
| [Figure_2.xlsx](../tables/Figure_2.xlsx) | PBMC, CLE, LN, Group_composition | Cell-type counts and percentages by sample, with pooled group composition |
| [Figure_3.xlsx](../tables/Figure_3.xlsx) | Naive_DE, Tem_DE, CTL_DE | Cell-level differential-expression results exported by `03_CD8_markers.R` |
| [Figure_S1.xlsx](../tables/Figure_S1.xlsx) | Donor_counts, Pseudobulk_counts, Naive_DE, Tem_DE, CTL_DE | Donor and subset cell counts, with complete donor-paired differential-expression results |
| [Figure_S2.xlsx](../tables/Figure_S2.xlsx) | Donor_counts | Cell counts by donor and median raw spliced/unspliced counts per matched cell |

Figure 2 percentages use a 0–100 scale within each sample and annotation level; group percentages use pooled cell counts. A zero denominator has a blank percentage.

Figure 3 compares each CD8 subset with the other two subsets pooled. `pct.1` and `pct.2` are the fractions of cells expressing the gene in the named subset and the other subsets, respectively (0–1 scale). `p_val_adj` is the Bonferroni-adjusted P value within each comparison. All tested genes are retained.

In Figure S1, each differential-expression sheet contains all 13,249 tested genes. Positive `log2FC` indicates higher expression in the named subset relative to the equally weighted mean of the other two subsets on the log-expression scale. `FDR` is the Benjamini–Hochberg adjustment within that contrast.
