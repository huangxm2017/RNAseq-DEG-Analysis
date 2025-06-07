# GO and KEGG Enrichment Analysis using clusterProfiler
library(clusterProfiler)
library(org.Hs.eg.db)  # Change for other organisms as needed
library(tidyverse)

# Read significant DEGs
deg <- read.csv("results/deseq2_sig_DEG.csv", row.names = 1)

# Extract gene symbols
gene_symbols <- rownames(deg)

# Convert SYMBOL to ENTREZ ID
gene_df <- bitr(gene_symbols, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

gene_ids <- unique(gene_df$ENTREZID)

# GO enrichment
go_enrich <- enrichGO(gene = gene_ids,
                      OrgDb = org.Hs.eg.db,
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = TRUE)

# KEGG enrichment
kegg_enrich <- enrichKEGG(gene = gene_ids,
                          organism = "hsa",
                          pvalueCutoff = 0.05)

# Output
write.csv(as.data.frame(go_enrich), "results/go_enrichment.csv")
write.csv(as.data.frame(kegg_enrich), "results/kegg_enrichment.csv")