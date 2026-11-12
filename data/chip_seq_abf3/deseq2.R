# Differential expression ABA vs mock from salmon results (salmon/ABA_repN, salmon/mock_repN)
# usage: Rscript deseq2.R <directory with salmon results> <tx2gene table> <output table>
#   e.g. Rscript deseq2.R salmon TAIR10_tx2gene.tsv deseq2_ABA_vs_mock.tsv
args <- commandArgs(trailingOnly = TRUE)
qdir <- args[1]; tx2gene_file <- args[2]; out <- args[3]
suppressMessages({library(tximport); library(DESeq2)})

# one salmon folder per sample (e.g. salmon/ABA_rep1); the condition is the
# part of the folder name before "_rep" (ABA or mock)
sample_names <- sort(basename(dirname(Sys.glob(file.path(qdir, "*", "quant.sf")))))
samples <- data.frame(sample = sample_names,
                      condition = factor(sub("_rep.*$", "", sample_names),
                                         levels = c("mock", "ABA")))
print(samples)
stopifnot(!anyNA(samples$condition), length(unique(samples$condition)) == 2)
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
