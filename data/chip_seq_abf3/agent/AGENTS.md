# ChIP-seq of ABF3 - instructions for the agent

You help a student analyse ChIP-seq data. The student decides WHAT to do; you
write and run the commands. Work in small steps and never skip ahead.

## Project

- Working directory: `~/chipseq_abf3` (run everything from here).
- Raw data (read-only, never modify or delete): `data/`
- Helper files (read-only): `~/Desktop/Bioinformatics/data/chip_seq_abf3/`
  (referred to as `HELP/` below)
- Organism: *Arabidopsis thaliana*, genome TAIR10 (chromosomes `1`-`5`, `Mt`,
  `Pt`; names are plain numbers, no `Chr` prefix).

Sequencing data in `data/` (one file per row). All reads are single-end
50 bp - never use paired-end options.

| File | Assay | Plants | Antibody | Treatment | Replicate | Reads |
|---|---|---|---|---|---|---|
| `ABF3_ABA_rep1.fastq.gz` | ChIP | ABF3-YPet | anti-GFP | ABA | 1 | 3 M |
| `ABF3_ABA_rep2.fastq.gz` | ChIP | ABF3-YPet | anti-GFP | ABA | 2 | 3 M |
| `ABF3_ABA_rep3.fastq.gz` | ChIP | ABF3-YPet | anti-GFP | ABA | 3 | 3 M |
| `ABF3_mock_rep1.fastq.gz` | ChIP | ABF3-YPet | anti-GFP | mock | 1 | 3 M |
| `ABF3_mock_rep2.fastq.gz` | ChIP | ABF3-YPet | anti-GFP | mock | 2 | 3 M |
| `ABF3_mock_rep3.fastq.gz` | ChIP | ABF3-YPet | anti-GFP | mock | 3 | 3 M |
| `control_rep1.fastq.gz` | ChIP control | wild type | anti-GFP | not stated | 1 | 2.6 M |
| `control_rep2.fastq.gz` | ChIP control | wild type | anti-GFP | not stated | 2 | 3 M |
| `rnaseq_ABA_rep1.fastq.gz` | RNA-seq | wild type | - | ABA | 1 | 4 M |
| `rnaseq_ABA_rep2.fastq.gz` | RNA-seq | wild type | - | ABA | 2 | 4 M |
| `rnaseq_mock_rep1.fastq.gz` | RNA-seq | wild type | - | mock | 1 | 4 M |
| `rnaseq_mock_rep2.fastq.gz` | RNA-seq | wild type | - | mock | 2 | 4 M |

The control is the same anti-GFP ChIP in plants without the tagged protein;
it is the control (`-c`) for every ChIP replicate. ABA and mock = 4 h
treatment.

| Reference file in `data/` | Content |
|---|---|
| `TAIR10_genome.fa.gz` | genome |
| `TAIR10_annotation.gff3.gz` | annotation (GFF3) |
| `TAIR10_transcripts.fa.gz` | transcripts (for salmon) |

| File in `HELP/` | Content |
|---|---|
| `TAIR10_genes.bed` | protein-coding genes, BED6: chrom, start, end, gene ID, symbol, strand |
| `TAIR10_tss.bed` | 1-bp TSS of the same genes, same columns |
| `TAIR10_tx2gene.tsv` | transcript -> gene table |
| `diffbind.R`, `deseq2.R` | ready R scripts (read the header for usage) |
| `integrate.R` | reference solution for the teacher - do NOT read or use it |
| `samples.csv` | example DiffBind sample sheet - do NOT copy it, write your own |
| `precomputed/` | reference results - use ONLY when the student asks |

## Environment rules

- All tools are already installed and on PATH (conda env `chipseq`): fastp,
  minimap2, samtools, macs3, deeptools, bedtools, seqkit, homer, salmon, R with
  DiffBind, DESeq2, tximport.
- NEVER install anything (no conda/pip/apt/R install.packages, no downloads).
  If a tool is missing, stop and tell the student.
- The computer has 4 CPU threads and 8 GB RAM: use `-t 4` / `-p 4` / `-@ 2`.
- No internet access is needed. Do not download HOMER genomes or databases.

## How to work

1. Before each step, say in 2-3 sentences what you will do and which files you
   will create. Wait for the student if anything is unclear.
2. Write each step as a script `scripts/NN_name.sh` (NN = 01, 02, ...), start
   it with `set -euo pipefail`, and run it with
   `bash scripts/NN_name.sh > logs/NN_name.log 2>&1`.
3. Steps that take longer than 1 minute (mapping, peak calling of all samples):
   first run the command on ONE sample to check it works, then run the whole
   script in the background with
   `nohup bash scripts/NN_name.sh > logs/NN_name.log 2>&1 &`
   and check progress with `tail -n 5 logs/NN_name.log` and `ls -l` of the
   output directory every 30-60 s (`sleep 30`).
4. Each step writes into its own directory (`qc/`, `trimmed/`, `bam/`,
   `bigwig/`, `peaks/`, `diffbind/`, `motifs/`, ...). Never overwrite `data/`.
5. After each step report: the script name, the output files, and the key
   numbers. Take every number from a file or command output (`wc -l`,
   `samtools flagstat`, `grep`, ...) and show the command - never estimate or
   guess a number.
6. If a command fails twice, STOP, show the error and ask the student. Do not
   try random alternatives and do not change the analysis design on your own.
7. Do not interpret the biology unless asked - the student does that.
8. Do not continue with the next step until the student asks.

## Tool notes (use these settings - they are tested)

- **Genome preparation**: `zcat data/TAIR10_genome.fa.gz > TAIR10.fa`,
  `samtools faidx TAIR10.fa`, `minimap2 -x sr -d TAIR10.mmi TAIR10.fa`.
- **fastp** (single-end): `fastp -i data/S.fastq.gz -o trimmed/S.fastq.gz -w 4 -h qc/S.html -j qc/S.json`.
  Summary numbers are in the JSON (`summary.before_filtering`, `after_filtering`).
- **minimap2**: `minimap2 -ax sr -t 4 TAIR10.mmi trimmed/S.fastq.gz | samtools sort -@ 2 -o bam/S.bam`,
  then `samtools index bam/S.bam`. Statistics: `samtools flagstat`,
  `samtools idxstats`. The two controls are merged into one deeper control:
  `samtools merge -@ 4 bam/control.bam bam/control_rep1.bam bam/control_rep2.bam` (+ index).
- **bamCoverage** (deeptools): `--normalizeUsing CPM --binSize 10 --extendReads 200 --minMappingQuality 10 -p 4`,
  output `bigwig/S.bw`.
- **MACS3**: `macs3 callpeak -t bam/S.bam -c bam/control.bam -f BAM -g 1.2e8 -q 0.01 -n S --outdir peaks`.
  Reads are single-end: do NOT use `-f BAMPE`. Run each ChIP replicate
  separately against the merged control. The fragment length is in
  `peaks/S_peaks.xls` (line `# d = `).
- **Profile around genes**: `computeMatrix reference-point -S <bigwigs> -R HELP/TAIR10_genes.bed --referencePoint TSS -b 2000 -a 2000 --skipZeros -p 4 -o M.gz`,
  then `plotProfile -m M.gz -o P.png --perGroup` and `plotHeatmap -m M.gz -o H.png`.
  For regions without strand use `--referencePoint center`.
- **DiffBind**: do NOT write your own DiffBind code. Use
  `Rscript HELP/diffbind.R samples.csv` run in the working directory. You
  must create `samples.csv` (comma-separated) with exactly these columns:
  `SampleID,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller`
  (Condition = `ABA` or `mock`, ControlID = `control`, bamControl =
  `bam/control.bam`, Peaks = the `_peaks.narrowPeak` file, PeakCaller =
  `narrow`; paths relative to the working directory). Outputs are written to
  the working directory: `diffbind_ABA_vs_mock.tsv`, `diffbind_ABA_vs_mock.bed`
  (chrom, start, end, name, log2FC, FDR), `sites_gained_200bp.bed`,
  `diffbind_plots.pdf`. Takes 1-3 minutes.
- **Closest gene**: `bedtools closest -a diffbind_ABA_vs_mock.bed -b HELP/TAIR10_tss.bed -d -t first`
  (both inputs sorted by position; sort with `sort -k1,1 -k2,2n` if needed).
- **Motifs** (HOMER has no Arabidopsis genome installed - do NOT use
  `findMotifsGenome.pl` and do not try to install a genome). Instead:
  `bedtools getfasta -fi TAIR10.fa -bed sites.bed -fo motifs/target.fa`;
  background = 5 rounds of
  `bedtools shuffle -i sites.bed -g TAIR10.fa.fai -chrom -noOverlapping -seed N`
  (N = 1..5, concatenated) -> getfasta -> `motifs/background.fa`; then
  `findMotifs.pl motifs/target.fa fasta motifs/homer -fasta motifs/background.fa -len 8,10 -S 5 -p 4`.
  Results: `motifs/homer/homerResults.html`, `motifs/homer/knownResults.txt`.
- **RNA-seq (optional)**: `salmon index -t data/TAIR10_transcripts.fa.gz -i salmon_index -p 4`;
  `salmon quant -i salmon_index -l A -r data/rnaseq_S.fastq.gz -p 4 --validateMappings -o salmon/S`
  (S without the `rnaseq_` prefix, e.g. `salmon/ABA_rep1`); then
  `Rscript HELP/deseq2.R salmon HELP/TAIR10_tx2gene.tsv deseq2_ABA_vs_mock.tsv`.
- **Binding vs expression (you write this script)**: do NOT use
  `HELP/integrate.R`. Write `scripts/NN_integrate.R` yourself, exactly as the
  student specifies (which sites, which genes, thresholds, plot, test). If
  the specification is incomplete, ask the student instead of choosing
  yourself. Use base R only (`read.delim`, `merge`, `table`, `boxplot`,
  `fisher.test`, `wilcox.test`, `pdf()`); do not load other packages. Inputs:
  - `sites_closest_gene.tsv` - NO header, tab-separated, 13 columns:
    1-6 site (chr, start, end, site name, binding log2 fold change ABA/mock,
    binding FDR), 7-12 TSS (chr, start, end, gene ID, gene symbol, strand),
    13 distance site-TSS in bp. A gene can appear in several rows.
  - `deseq2_ABA_vs_mock.tsv` - with header, tab-separated, one row per gene:
    `gene_id`, `symbol`, `baseMean`, `log2FoldChange` (ABA/mock), `lfcSE`,
    `stat`, `pvalue`, `padj` (`padj` can be `NA` - treat as not significant).
  - Join the two tables on the gene ID (column 10 of the sites = `gene_id`).
  Print every result table to the screen and write it to a `.tsv` file.
- Images (`.png`, `.pdf`, `.html`) and IGV are viewed by the student, not by
  you. Tell the student which file to open.
