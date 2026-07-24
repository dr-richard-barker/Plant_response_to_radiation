# Legacy reports

The bulky rendered HTML reports that previously lived here — per-contrast DESeq2 reports
(`deg/*.html`), pathway-workflow reports (`pathway/*.html`), and the combined
`deg_stats_report_*.html` — have been **removed to keep the repository lightweight**. They were
rendered views of data that is still preserved in the repository:

- Differential-expression calls: `legacy/data/processed/rnaseq/Up_and_down_genes_DESeq2.csv`
- WGCNA modules / enrichment: `legacy/data/processed/rnaseq/WGCNA_modules.csv`,
  `WGCNA_june_2026_enrichment.csv`
- Analysis scripts: `legacy/scripts/`

The self-contained **Metascape** enrichment report is retained under `metascape/` (including its
GO/PPI enrichment CSVs and figures).

The removed HTML remains recoverable from git history prior to the pruning commit.
