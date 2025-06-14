# RNA-seq Differential Expression and Functional Enrichment Analysis
This project provides a reproducible pipeline for performing differential gene expression analysis (via DESeq2) and functional enrichment analysis (GO and KEGG) using RNA-seq data. The study involves 6 samples divided into two conditions:
- **PAP**: PAP_1, PAP_2, PAP_3  
- **Control**: con_1, con_2, con_3
---
## Differential Expression Analysis (DESeq2)
- R package: `DESeq2`
- Input: gene count matrix (e.g., generated by `HTSeq`), along with sample group annotation
- Comparison: PAP vs Control
- Significance thresholds:
  - |log₂FoldChange| ≥ 1
  - Adjusted p-value (padj) < 0.05 (Benjamini-Hochberg correction)
Script: `scripts/01_deseq2_deg.R`
---
## GO and KEGG Enrichment Analysis
- R package: `clusterProfiler` (v4.6.0)
- Input: list of significantly differentially expressed genes (DEGs)
- Annotation database: `org.Hs.eg.db` (can be changed for other species)
- Threshold: adjusted p-value ≤ 0.05 using BH correction
Script: `scripts/02_go_kegg_enrichment.R`
---
## Quick Start
```bash
# Step 1: Run differential expression analysis
Rscript scripts/01_deseq2_deg.R
# Step 2: Run GO and KEGG enrichment analysis
Rscript scripts/02_go_kegg_enrichment.R
