# Differential Gene Expression Analysis using DESeq2
library(DESeq2)

# Load count data (replace with your actual file)
count_data <- read.csv("count_matrix.csv", row.names = 1)
col_data <- read.table("config/sample_info.txt", header = TRUE, row.names = 1)

# Construct DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = ~ group)

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "PAP", "Control"))

# Remove rows with NA adjusted p-values
res <- res[!is.na(res$padj), ]

# Save full results
write.csv(as.data.frame(res), file = "results/deseq2_results.csv")

# Save significant DEGs
sig_res <- res[abs(res$log2FoldChange) >= 1 & res$padj < 0.05, ]
write.csv(as.data.frame(sig_res), file = "results/deseq2_sig_DEG.csv")