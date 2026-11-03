# Differential expression ABA vs mock (2 + 2 replicates) from salmon results
# usage: Rscript deseq2.R <directory with salmon results> <tx2gene table> <output table>
#   e.g. Rscript deseq2.R salmon TAIR10_tx2gene.tsv deseq2_ABA_vs_mock.tsv
args <- commandArgs(trailingOnly = TRUE)
qdir <- args[1]; tx2gene_file <- args[2]; out <- args[3]
suppressMessages({library(tximport); library(DESeq2)})

samples <- data.frame(sample = c("ABA_rep1", "ABA_rep2", "mock_rep1", "mock_rep2"),
                      condition = factor(c("ABA", "ABA", "mock", "mock"),
                                         levels = c("mock", "ABA")))
files <- file.path(qdir, samples$sample, "quant.sf")
names(files) <- samples$sample

# salmon counts reads per transcript; sum them per gene
# (table: transcript ID, gene ID, gene symbol)
tx2gene <- read.delim(tx2gene_file, header = FALSE)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, 1:2])

# DESeq2: normalisation, dispersion estimation, test ABA vs mock
dds <- DESeqDataSetFromTximport(txi, samples, ~ condition)
dds <- DESeq(dds, quiet = TRUE)
res <- as.data.frame(results(dds, contrast = c("condition", "ABA", "mock")))

symbols <- unique(tx2gene[, 2:3])
res <- cbind(gene_id = rownames(res),
             symbol = symbols[match(rownames(res), symbols[, 1]), 2],
             res)
res <- res[order(res$padj), ]
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Genes:", nrow(res), "\n")
cat("padj < 0.05 - higher in ABA:", sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE),
    "- lower in ABA:", sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE), "\n")

pdf(sub("\\.tsv$", "_plots.pdf", out))
print(plotPCA(vst(dds, blind = TRUE, nsub = 1000), intgroup = "condition"))
plotMA(results(dds, contrast = c("condition", "ABA", "mock")), ylim = c(-8, 8))
invisible(dev.off())
cat("Written:", out, "and", sub("\\.tsv$", "_plots.pdf", out), "\n")
