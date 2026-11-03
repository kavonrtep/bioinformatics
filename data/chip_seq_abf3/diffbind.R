# Differential binding of ABF3: ABA vs mock (3 + 3 replicates)
# usage (in the working directory with bam/ and peaks/):
#   Rscript diffbind.R samples.csv
args <- commandArgs(trailingOnly = TRUE)
sample_sheet <- if (length(args) > 0) args[1] else "samples.csv"
suppressMessages(library(DiffBind))

# 1. read the sample sheet: BAM files, peaks, condition of each sample
db <- dba(sampleSheet = sample_sheet)
print(db)

# 2. consensus peak set (peaks found in at least 2 samples), each peak
#    re-centred on its summit and set to 401 bp; count the reads in each
#    peak in each sample
db <- dba.count(db, summits = 200, bParallel = TRUE)
print(db)

# 3. normalisation (library size) and DESeq2 test ABA vs mock
#    (there is no blacklist of problematic regions for Arabidopsis)
db <- dba.normalize(db)
db <- dba.contrast(db, contrast = c("Condition", "ABA", "mock"))
db <- dba.analyze(db, bBlacklist = FALSE, bGreylist = FALSE)
dba.show(db, bContrasts = TRUE)

# 4. results for all sites (th = 1), sorted by FDR
res <- as.data.frame(dba.report(db, th = 1))
res <- res[order(res$FDR), ]
cat("\nConsensus sites:", nrow(res), "\n")
cat("FDR < 0.05:", sum(res$FDR < 0.05),
    "- stronger in ABA:", sum(res$FDR < 0.05 & res$Fold > 0),
    "- stronger in mock:", sum(res$FDR < 0.05 & res$Fold < 0), "\n")
write.table(res, "diffbind_ABA_vs_mock.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# the same as BED (0-based start), sorted by position, for bedtools:
# chromosome, start, end, name, log2 fold change, FDR
bed <- data.frame(res$seqnames, res$start - 1, res$end,
                  paste0("site", seq_len(nrow(res))),
                  round(res$Fold, 3), signif(res$FDR, 3))
bed <- bed[order(bed[, 1], bed[, 2]), ]
write.table(bed, "diffbind_ABA_vs_mock.bed", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# sites with significantly stronger binding in ABA, 200 bp around the
# centre of the site (= the summit), for the motif search and the profile
# plot; standard BED6: chromosome, start, end, name, log2 fold change,
# strand ("." = no strand)
gained <- bed[bed[, 6] < 0.05 & bed[, 5] > 0, ]
centre <- (gained[, 2] + gained[, 3]) %/% 2
gained[, 2] <- centre - 100
gained[, 3] <- centre + 100
gained[, 6] <- "."
write.table(gained, "sites_gained_200bp.bed", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# 5. plots
pdf("diffbind_plots.pdf")
dba.plotPCA(db, DBA_CONDITION, label = DBA_ID)
dba.plotHeatmap(db)
dba.plotMA(db)
dba.plotVolcano(db)
invisible(dev.off())
cat("Written: diffbind_ABA_vs_mock.tsv, diffbind_ABA_vs_mock.bed,",
    "sites_gained_200bp.bed, diffbind_plots.pdf\n")
