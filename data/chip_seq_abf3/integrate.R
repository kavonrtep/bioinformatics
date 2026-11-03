# Are the genes with ABA-gained ABF3 binding induced by ABA?
# usage: Rscript integrate.R <sites with the closest gene> <DESeq2 table>
#   e.g. Rscript integrate.R sites_closest_gene.tsv deseq2_ABA_vs_mock.tsv
# sites_closest_gene.tsv = output of bedtools closest -d (DiffBind sites + TSS):
#   1-6 site (chr, start, end, name, log2 fold change, FDR),
#   7-12 TSS (chr, start, end, gene ID, symbol, strand), 13 distance
args <- commandArgs(trailingOnly = TRUE)
sites <- read.delim(args[1], header = FALSE)
colnames(sites) <- c("chr", "start", "end", "site", "binding_log2fc", "binding_fdr",
                     "tss_chr", "tss_start", "tss_end", "gene_id", "symbol", "strand", "distance")
expr <- read.delim(args[2])

# only sites up to 1 kb from the TSS of a gene
sites <- sites[sites$distance <= 1000, ]
gained <- unique(sites$gene_id[sites$binding_fdr < 0.05 & sites$binding_log2fc > 0])
unchanged <- setdiff(unique(sites$gene_id[sites$binding_fdr >= 0.05]), gained)

expr$class <- "not significant"
expr$class[!is.na(expr$padj) & expr$padj < 0.05 & expr$log2FoldChange > 0] <- "induced by ABA"
expr$class[!is.na(expr$padj) & expr$padj < 0.05 & expr$log2FoldChange < 0] <- "repressed by ABA"

summary_row <- function(name, ids) {
  e <- expr[expr$gene_id %in% ids, ]
  data.frame(genes = name, n = nrow(e),
             induced_pct = round(100 * mean(e$class == "induced by ABA"), 1),
             repressed_pct = round(100 * mean(e$class == "repressed by ABA"), 1))
}
tab <- rbind(summary_row("all genes", expr$gene_id),
             summary_row("ABF3 binding unchanged", unchanged),
             summary_row("ABF3 binding gained in ABA", gained))
print(tab, row.names = FALSE)

# table of genes with gained binding and their expression
g <- sites[sites$binding_fdr < 0.05 & sites$binding_log2fc > 0,
           c("gene_id", "symbol", "distance", "binding_log2fc", "binding_fdr")]
g <- merge(g, expr[, c("gene_id", "baseMean", "log2FoldChange", "padj", "class")], by = "gene_id")
g <- g[order(g$binding_fdr), ]
write.table(g, "genes_gained_binding_expression.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

pdf("integration_plot.pdf")
boxplot(list("all genes" = expr$log2FoldChange,
             "binding unchanged" = expr$log2FoldChange[expr$gene_id %in% unchanged],
             "binding gained" = expr$log2FoldChange[expr$gene_id %in% gained]),
        ylab = "log2 fold change of expression (ABA / mock)", outline = FALSE)
abline(h = 0, lty = 2)
invisible(dev.off())
cat("Written: genes_gained_binding_expression.tsv, integration_plot.pdf\n")
