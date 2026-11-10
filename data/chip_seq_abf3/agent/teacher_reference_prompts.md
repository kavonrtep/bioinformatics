# Reference prompts for the teacher

Prompts used in the pilot run of the agent exercise
([chip_seq_abf3_agent.md](../../../docs/chip_seq_abf3_agent.md)). One prompt per step, all in
one opencode session, with `AGENTS.md` and `opencode.json` from this directory in the working
directory. They show the level of detail the free model needs. Do not give them to students
before they have written their own plan.

## Pilot result (2026-11-10, opencode v2.2.1, `opencode/big-pickle`)

- Steps 00-09 finished without intervention in about 15 minutes (mapping 4.5 min).
- Tokens: 77.6k input, 14.1k output, 16.4k reasoning, 1.9M cache read (96 % cached),
  76 model steps, $0.
- Results identical to `../precomputed/`: peaks per replicate 8894 / 8870 / 7549 / 7750 / 4313 /
  6751, 8803 DiffBind consensus sites, 452 gained / 9 lost (FDR < 0.05); top de novo motif
  `CACGTGTC` (80 % of targets), top known motifs bZIP / ABF1.
- The agent followed `AGENTS.md` (scripts in `scripts/`, logs in `logs/`, no installs, used
  `diffbind.R`), corrected two wrong column indices itself, but added biological
  interpretation although told not to.
- Step 09 of the pilot still used `integrate.R`; the integration step below (09b) is the
  current design, where the agent writes the R script from the student's specification.
- Expect students to need more tokens (vaguer prompts, wrong turns, questions).

## Prompts

```text
00  Warm-up: prepare the genome for mapping (decompress, index for samtools and minimap2)
    as described in AGENTS.md. Then tell me how many chromosomes the genome has and how
    long they are.

01  Step 1 - read QC. Run fastp on the 8 ChIP-seq samples (ABF3_* and control_*, not
    RNA-seq), trimmed reads to trimmed/, reports to qc/. Report a table: sample, reads
    before, reads after filtering.

02  Step 2 - mapping. Map the trimmed ChIP-seq reads to TAIR10 with minimap2 (sorted,
    indexed BAM in bam/), then merge the two controls into bam/control.bam. Test on one
    sample first, run the rest in the background. Report a table: sample, total reads,
    % mapped (samtools flagstat), and reads on Mt and Pt (samtools idxstats).

03  Step 3 - coverage tracks. Make CPM-normalised bigWig files for the 6 ABF3 samples and
    the merged control into bigwig/, settings from AGENTS.md. List the files when done.

04  Step 4 - peak calling. Call peaks with MACS3 in each of the 6 ABF3 replicates
    separately against the merged control. Report a table: sample, number of peaks,
    fragment length d.

05  Step 5 - profile around genes. Make a TSS profile and heatmap (+-2 kb) of
    ABF3_ABA_rep1, ABF3_mock_rep1 and control bigWigs over all genes. Tell me which files
    to open.

06  Step 6 - differential binding. Create samples.csv for the 6 ABF3 replicates (format in
    AGENTS.md) and run the DiffBind script. Report: number of consensus sites, number of
    sites with FDR < 0.05 gained and lost in ABA.

07  Step 7 - closest genes. For all DiffBind sites find the closest TSS (bedtools closest,
    HELP/TAIR10_tss.bed). Save to sites_closest_gene.tsv and show the 20 most significant
    sites with gene ID, symbol and distance.

08  Step 8 - motifs. Find DNA motifs in the gained sites (sites_gained_200bp.bed) against
    a shuffled genomic background with HOMER as described in AGENTS.md. Report the best 3
    de novo motifs (consensus, % target, % background) and the top 5 known motifs.

09a Step 9 (optional) - RNA-seq. Quantify the 4 RNA-seq samples with salmon and run the
    DESeq2 script. Report the number of genes with padj < 0.05 up and down in ABA.

09b Integration step: use sites_closest_gene.tsv and deseq2_ABA_vs_mock.tsv. Only sites up
    to 1 kb from a TSS. Groups of genes: gained binding (binding FDR < 0.05 and
    log2FC > 0), unchanged binding (FDR >= 0.05, not in gained), and all genes in the
    DESeq2 table. Induced = padj < 0.05 and log2FoldChange > 0, repressed = padj < 0.05
    and log2FoldChange < 0. Output: table with number of genes and % induced and %
    repressed per group; Fisher test of induced vs not induced, gained vs unchanged;
    boxplot of expression log2FoldChange per group to integration_mine.pdf; table of genes
    with gained binding and their expression. Show me the R script first, then run it.
```

The percentages from 09b can be checked against
`Rscript ../integrate.R sites_closest_gene.tsv deseq2_ABA_vs_mock.tsv`, which uses the same
definitions.
