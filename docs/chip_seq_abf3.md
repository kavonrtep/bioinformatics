# ChIP-seq - how does the transcription factor ABF3 respond to abscisic acid?

The plant hormone **abscisic acid (ABA)** accumulates in plants exposed to
drought, salt or cold and switches on thousands of stress-response genes. A key
part of this response are the **ABFs** (ABRE-binding factors, also called
AREBs), bZIP transcription factors that bind the **ABA-responsive element
(ABRE)**, a short DNA motif with the core `ACGTG`.

In this exercise you will analyse ChIP-seq of the transcription factor **ABF3**
(AT4G34000) from *Arabidopsis thaliana* seedlings treated with ABA or with a
mock solution, and ask:

1. Where in the genome does ABF3 bind? Close to genes?
2. Does the binding change after ABA treatment? (differential binding)
3. Which DNA motif is in the binding sites?
4. Are the genes with stronger ABF3 binding also more expressed after ABA?
   (RNA-seq from the same experiment)

Everything runs on your computer (or on a VPS with 4 CPUs and 8 GB RAM); all
steps together take about 10 minutes of computing time.

The data come from [Song et al. 2016, *Science* 354: aag1550](https://doi.org/10.1126/science.aag1550)
(GEO [GSE80568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80568)),
who did ChIP-seq of 21 ABA-related transcription factors. Three-day-old
seedlings grown in the dark were treated with 10 µM ABA or with a mock
solution (ethanol) for 4 hours. The plants carry an extra copy of the ABF3
gene (with its own promoter) in which ABF3 is fused with the fluorescent
protein YPet; the fusion protein was precipitated with an anti-GFP antibody,
which also recognises YPet. The control is the same anti-GFP ChIP in
wild-type plants (Col-0) without any tagged protein. To
make the exercise fast, the reads were reduced to 3 million per sample and
shortened to 50 bp.

| Sample | Treatment | Reads |
|---|---|---|
| `ABF3_ABA_rep1`, `rep2`, `rep3` | ChIP of ABF3-YPet, ABA 4 h | 3 M × 50 bp |
| `ABF3_mock_rep1`, `rep2`, `rep3` | ChIP of ABF3-YPet, mock 4 h | 3 M × 50 bp |
| `control_rep1`, `control_rep2` | anti-GFP ChIP in wild type (Col-0) | 2.6 M, 3 M × 50 bp |
| `rnaseq_ABA_rep1`, `rep2` | RNA-seq, ABA 4 h | 4 M × 50 bp |
| `rnaseq_mock_rep1`, `rep2` | RNA-seq, mock 4 h | 4 M × 50 bp |

Contents:

- [Installation and data](#installation-and-data)
- [Part 1 - From reads to binding sites](#part-1---from-reads-to-binding-sites):
  [1.1 QC](#11-quality-control-with-fastp) ·
  [1.2 mapping](#12-mapping-with-minimap2) ·
  [1.3 mapping statistics](#13-mapping-statistics) ·
  [1.4 coverage tracks](#14-coverage-tracks) ·
  [1.5 peaks](#15-peak-calling-with-macs3) ·
  [1.6 IGV](#16-binding-sites-in-igv) ·
  [1.7 binding around genes](#17-where-does-abf3-bind-relative-to-genes)
- [Part 2 - Differential binding](#part-2---differential-binding):
  [2.1 DiffBind](#21-differential-binding-with-diffbind) ·
  [2.2 closest genes](#22-genes-closest-to-the-binding-sites) ·
  [2.3 profile](#23-profile-of-the-gained-sites) ·
  [2.4 motifs](#24-dna-motifs-in-the-gained-sites)
- [Part 3 - Binding and expression](#part-3---binding-and-expression):
  [3.1 salmon](#31-gene-expression-with-salmon) ·
  [3.2 DESeq2](#32-differential-expression-with-deseq2) ·
  [3.3 binding vs expression](#33-are-the-genes-with-gained-binding-induced-by-aba)

## Installation and data

Create the conda environment (a few minutes):

```bash
conda create -n chipseq -y --override-channels -c conda-forge -c bioconda \
    fastp minimap2 samtools macs3 deeptools bedtools seqkit homer salmon \
    bioconductor-diffbind bioconductor-deseq2 bioconductor-tximport bioconductor-genomeinfodb
conda activate chipseq
```

<details>
<summary>What is in the environment</summary>

- `fastp` - quality control and trimming of reads
- `minimap2`, `samtools` - mapping, BAM files
- `macs3` - peak calling
- `deeptools` - coverage tracks (bigWig), profiles and heatmaps
- `bedtools`, `seqkit` - work with intervals and sequences
- `homer` - motif search
- `salmon` - quantification of gene expression from RNA-seq
- R packages `DiffBind` (differential binding), `DESeq2` (differential
  expression), `tximport` (reads salmon results into R)

</details>

The reads, the genome and the annotation are on Zenodo (1.9 GB). The other
files (gene tables, R scripts, results for the case something fails) are in
the course repository, in `~/Desktop/Bioinformatics/data/chip_seq_abf3/`
(if you cloned the repository elsewhere, change the path in `D=...`):

```bash
D=~/Desktop/Bioinformatics/data/chip_seq_abf3
mkdir -p ~/chipseq_abf3/data
cd ~/chipseq_abf3/data
Z=https://zenodo.org/records/XXXXXXX/files
for f in $(cut -d ' ' -f 3 $D/md5sums.txt); do
    wget -c $Z/$f
done
md5sum -c $D/md5sums.txt
cd ~/chipseq_abf3
```

<details>
<summary>Options used</summary>

- `cut -d ' ' -f 3 md5sums.txt` - the third column (separated by spaces) of
  the checksum file = the file names; the `for` loop downloads them one by one
- `wget -c` - continue an interrupted download instead of starting again
- `md5sum -c` - computes the MD5 checksum of each file and compares it with
  the value in `md5sums.txt`. A large download can be damaged without any
  error message (it happened when this exercise was prepared!) - the checksum
  is the only way to be sure. A file reported as `FAILED` must be deleted and
  downloaded again.

</details>

All further commands are run in `~/chipseq_abf3`. Prepare the genome:
decompress it, index it for samtools and for minimap2, and decompress the
annotation for IGV:

```bash
zcat data/TAIR10_genome.fa.gz > TAIR10.fa
samtools faidx TAIR10.fa
cut -f 1,2 TAIR10.fa.fai
minimap2 -x sr -d TAIR10.mmi TAIR10.fa
zcat data/TAIR10_annotation.gff3.gz > TAIR10_annotation.gff3
```

- How many chromosomes does the genome have, and how long are they? What are
  `Mt` and `Pt`?

## Part 1 - From reads to binding sites

### 1.1 Quality control with fastp

```bash
mkdir -p qc trimmed
for s in ABF3_ABA_rep1 ABF3_ABA_rep2 ABF3_ABA_rep3 ABF3_mock_rep1 ABF3_mock_rep2 ABF3_mock_rep3 control_rep1 control_rep2; do
    fastp -i data/$s.fastq.gz -o trimmed/$s.fastq.gz -w 4 -h qc/$s.html -j qc/$s.json
done
```

<details>
<summary>Options used</summary>

- `for s in ...; do ...; done` - repeats the command for each sample; `$s` is
  replaced by the sample name
- `-i`, `-o` - input and output (trimmed) reads
- `-w 4` - number of threads
- `-h`, `-j` - reports in HTML and JSON

</details>

Open some of the reports (`firefox qc/ABF3_ABA_rep1.html`).

- What is the quality of the reads? Are there adapters? How many reads were
  removed?
- Compare the **duplication rate** of the ABF3 ChIP and of the control. Why
  may a control ChIP (without any target protein) contain more duplicates?

### 1.2 Mapping with minimap2

minimap2 with the preset `sr` (short reads) is fast enough to map all eight
samples in about 2 minutes:

```bash
mkdir -p bam
for s in ABF3_ABA_rep1 ABF3_ABA_rep2 ABF3_ABA_rep3 ABF3_mock_rep1 ABF3_mock_rep2 ABF3_mock_rep3 control_rep1 control_rep2; do
    minimap2 -ax sr -t 4 TAIR10.mmi trimmed/$s.fastq.gz | samtools sort -@ 2 -o bam/$s.bam
    samtools index bam/$s.bam
done
# the two controls together give one deeper control
samtools merge -@ 4 bam/control.bam bam/control_rep1.bam bam/control_rep2.bam
samtools index bam/control.bam
```

<details>
<summary>Options used</summary>

- `minimap2 -a` - output in the SAM format; `-x sr` - preset for short
  reads; `-t 4` - threads. `TAIR10.mmi` is the index made above (`-d`).
- `samtools sort -@ 2 -o file.bam` - sorts the alignments by position and
  saves them as BAM
- `samtools index` - index (`.bai`) needed by IGV, MACS3, deeptools, ...
- `samtools merge out.bam in1.bam in2.bam` - joins BAM files

</details>

### 1.3 Mapping statistics

```bash
samtools flagstat bam/ABF3_ABA_rep1.bam
samtools idxstats bam/ABF3_ABA_rep1.bam
samtools idxstats bam/control.bam
```

`samtools idxstats` prints for each chromosome its name, length, number of
mapped reads and number of unmapped reads.

- What fraction of the reads was mapped?
- What fraction of the reads comes from the chloroplast genome (`Pt`) in the
  ABF3 ChIP and in the control? Why can the control contain so many
  chloroplast reads? (Hint: how much DNA does a ChIP without its target
  protein precipitate, and how many copies of the chloroplast genome are in
  one cell?)

### 1.4 Coverage tracks

BAM files are large and slow to view for the whole genome. `bamCoverage`
(deeptools) computes the coverage in 10 bp bins and saves it in the compact
bigWig format. The coverage is normalised to **CPM** (reads per million mapped
reads), so that samples with different numbers of reads can be compared:

```bash
mkdir -p bigwig
for s in ABF3_ABA_rep1 ABF3_ABA_rep2 ABF3_ABA_rep3 ABF3_mock_rep1 ABF3_mock_rep2 ABF3_mock_rep3 control; do
    bamCoverage -b bam/$s.bam -o bigwig/$s.bw --normalizeUsing CPM --binSize 10 \
        --extendReads 200 --minMappingQuality 10 -p 4
done
```

<details>
<summary>Options used</summary>

- `--normalizeUsing CPM` - coverage per million mapped reads
- `--binSize 10` - resolution 10 bp
- `--extendReads 200` - each 50 bp read is extended to 200 bp, the typical
  length of the sequenced DNA fragments (the read is only the beginning of the
  fragment; the protein sits in the middle)
- `--minMappingQuality 10` - ignore reads that map equally well to several
  places (repeats)
- `-p 4` - threads

All options: [bamCoverage](https://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html).

</details>

### 1.5 Peak calling with MACS3

MACS3 looks for regions where the ChIP has significantly more reads than the
control. It first estimates the fragment length from the shift between reads
on the + and - strand around the binding sites, then tests each region
against the local background estimated from the control:

```bash
mkdir -p peaks
for s in ABF3_ABA_rep1 ABF3_ABA_rep2 ABF3_ABA_rep3 ABF3_mock_rep1 ABF3_mock_rep2 ABF3_mock_rep3; do
    macs3 callpeak -t bam/$s.bam -c bam/control.bam -f BAM -g 1.2e8 -q 0.01 -n $s --outdir peaks
done
wc -l peaks/*.narrowPeak
grep "# d = " peaks/*.xls
```

<details>
<summary>Options used</summary>

- `-t` - ChIP (treatment); `-c` - control
- `-f BAM` - input format
- `-g 1.2e8` - effective genome size: the part of the genome where reads can
  be mapped uniquely (for *Arabidopsis* about 120 Mb of the 135 Mb)
- `-q 0.01` - report peaks with FDR (q-value) < 0.01
- `-n`, `--outdir` - name prefix and output directory

Output: `NAME_peaks.narrowPeak` (BED with extra columns: 7 fold enrichment, 8
-log10 p-value, 9 -log10 q-value, 10 position of the summit relative to the
start), `NAME_summits.bed` (summits only), `NAME_peaks.xls` (table with the
parameters of the run). All options: [MACS3 callpeak](https://macs3-project.github.io/MACS/docs/callpeak.html).

</details>

- How many peaks were found in each sample? Are there more peaks in ABA or in
  mock? Which sample is unusual?
- What fragment length (`d`) did MACS3 estimate? Does it agree with the 200 bp
  used in `bamCoverage`?
- Look at the first lines of a `narrowPeak` file
  (`sort -k9,9gr peaks/ABF3_ABA_rep1_peaks.narrowPeak | head`). How wide are
  the peaks? How strong is the enrichment of the best ones?

<details>
<summary>💡 MACS3 failed?</summary>

```bash
mkdir -p peaks
cp $D/precomputed/peaks/*.narrowPeak peaks/
```
</details>

### 1.6 Binding sites in IGV

Start IGV, load `TAIR10.fa` as the genome (`Genomes → Load Genome from
File...`), then `TAIR10_annotation.gff3`, the bigWig files of one ABA and one
mock replicate and of the control, and the peaks of the same replicates
(`File → Load from File...`). Set the same data range for all bigWig tracks
(select the tracks, right-click → "Group Autoscale").

> Working on a VPS without a graphical desktop? Copy the files to your
> computer first, e.g. `scp -r user@server:chipseq_abf3/bigwig .` (the same
> for `peaks/`, `TAIR10.fa`, `TAIR10.fa.fai`, `TAIR10_annotation.gff3`).

Go to these genes (type the gene ID into the search box):

| Gene | ID | Function |
|---|---|---|
| RD29B (`LTI65` in the annotation) | AT5G52300 | ABA-induced protein, classic marker gene of the ABA response |
| RAB18 | AT5G66400 | dehydrin, ABA marker gene (there is a different gene with the same name on chromosome 1!) |
| LEA7 | AT1G52690 | late embryogenesis abundant protein |
| ABI1 | AT4G26080 | protein phosphatase, negative regulator of ABA signalling |
| ACT2 | AT3G18780 | actin, a "housekeeping" gene |
| ABF3 | AT4G34000 | the transcription factor itself |

- Where are the peaks relative to the genes (promoter, gene body, 3' end)?
- At which genes is the binding stronger after ABA? At which is ABF3 bound
  already without ABA?
- Look at the ABF3 gene itself. The signal in all ABF3 ChIP samples is
  extremely high and covers the whole gene, but there is nothing in the
  control. Is it a binding site? (Hint: how were the plants for the ChIP made?
  Where do reads from the extra ABF3-YPet gene copy map?)

### 1.7 Where does ABF3 bind relative to genes?

`computeMatrix` collects the coverage around the transcription start sites
(TSS) of all genes; `plotProfile` and `plotHeatmap` show the average and each
gene separately:

```bash
computeMatrix reference-point -S bigwig/ABF3_ABA_rep1.bw bigwig/ABF3_mock_rep1.bw bigwig/control.bw \
    -R $D/TAIR10_genes.bed --referencePoint TSS -b 2000 -a 2000 --skipZeros -p 4 -o tss_matrix.gz
plotProfile -m tss_matrix.gz -o tss_profile.png --perGroup
plotHeatmap -m tss_matrix.gz -o tss_heatmap.png
```

<details>
<summary>Options used</summary>

- `reference-point --referencePoint TSS` - align all genes at their TSS (the
  start of a `+` gene, the end of a `-` gene, from column 6 of the BED file)
- `-S` - bigWig files; `-R` - regions (all protein-coding genes)
- `-b 2000 -a 2000` - 2 kb before and after the TSS
- `--skipZeros` - skip genes without any signal
- `--perGroup` - all samples in one plot

</details>

- Where is the maximum of the ABF3 signal? Upstream or downstream of the TSS?
- The ABA and mock curves are almost the same. Does it mean that ABA does not
  change the binding? (Hint: how many of the 27,000 genes do you expect to
  change?)
- Look at the heatmap. Is ABF3 bound at all genes or only at some?

## Part 2 - Differential binding

### 2.1 Differential binding with DiffBind

Comparing the numbers of peaks is not enough - a peak may be called in one
sample and missed in another just by chance. **DiffBind** counts the reads in
each binding site in each sample and tests, like RNA-seq, whether the
difference between the conditions is larger than the variability between
replicates (it uses DESeq2 for the test). The sample sheet `samples.csv` tells
DiffBind which files belong to which sample and condition:

```bash
cp $D/samples.csv .
cat samples.csv
Rscript $D/diffbind.R samples.csv
```

Read the commented script (`less $D/diffbind.R`): it makes a consensus set of
sites (peaks found in at least two samples), counts the reads, normalises,
tests ABA vs mock and writes:

- `diffbind_ABA_vs_mock.tsv` - all sites with the log2 fold change (`Fold`),
  p-value and FDR
- `diffbind_ABA_vs_mock.bed` - the same as BED: chromosome, start, end, name,
  log2 fold change, FDR
- `sites_gained_200bp.bed` - sites with significantly stronger binding in ABA
  (FDR < 0.05), 200 bp around the summit
- `diffbind_plots.pdf` - PCA, correlation heatmap, MA plot, volcano plot

<details>
<summary>💡 DiffBind failed?</summary>

```bash
cp $D/precomputed/diffbind_ABA_vs_mock.* $D/precomputed/sites_gained_200bp.bed .
```
</details>

- How many consensus sites are there? How many change significantly after
  ABA, and in which direction?
- Look at the PCA plot. Do the samples separate by the treatment? Which
  sample is different from the others? Compare it with the numbers of peaks
  from step 1.5.
- In the MA plot, where are the significant sites? Why are sites with few
  reads rarely significant even if their fold change is large?

### 2.2 Genes closest to the binding sites

Which genes are regulated from these sites? The simplest guess is the gene
with the closest TSS. `bedtools closest` finds it for each site:

```bash
bedtools closest -a diffbind_ABA_vs_mock.bed -b $D/TAIR10_tss.bed -d -t first > sites_closest_gene.tsv
head sites_closest_gene.tsv
# the most significant sites
sort -k6,6g sites_closest_gene.tsv | head -20
```

<details>
<summary>Options used</summary>

- `bedtools closest -a sites -b TSS` - for each site of A the closest TSS of B
  (both must be sorted by position)
- `-d` - add the distance in bp as the last column (0 = overlapping)
- `-t first` - if two TSS are equally close, report only the first
- columns of the output: 1-6 the site (chromosome, start, end, name, log2 fold
  change, FDR), 7-12 the TSS (chromosome, start, end, gene ID, gene symbol,
  strand), 13 distance
- `sort -k6,6g` - sort by column 6 (FDR) as numbers in the general format
  (`g` understands `1e-10`)

</details>

- Which genes are close to the most significant sites? Look up some of them
  in [TAIR](https://www.arabidopsis.org/) (search by the gene ID). Have you
  seen them in IGV?
- How far from the TSS are the sites? (`cut -f 13 sites_closest_gene.tsv |
  sort -n | uniq -c | head`)
- Find the site at the RAB18 promoter (AT5G66400) - which gene is reported as
  the closest? Look at the region in IGV. What is the problem with the
  "closest gene" assignment for genes in head-to-head orientation?

### 2.3 Profile of the gained sites

Now the average profile only over the sites gained after ABA:

```bash
computeMatrix reference-point -S bigwig/ABF3_ABA_rep*.bw bigwig/ABF3_mock_rep*.bw bigwig/control.bw \
    -R sites_gained_200bp.bed --referencePoint center -b 1000 -a 1000 -p 4 -o gained_matrix.gz
plotProfile -m gained_matrix.gz -o gained_profile.png --perGroup
plotHeatmap -m gained_matrix.gz -o gained_heatmap.png
```

- How large is the difference between ABA and mock at these sites? Are the
  replicates consistent?
- Is ABF3 bound at these sites already without ABA (compare mock with the
  control)? What does it say about how ABA changes the binding - new sites, or
  stronger binding at existing sites?

### 2.4 DNA motifs in the gained sites

Which DNA sequence is common to the gained sites? HOMER compares the frequency
of short sequences (k-mers) in the sites with a background - here random
regions of the genome of the same size (5 times more of them):

```bash
mkdir -p motifs
bedtools getfasta -fi TAIR10.fa -bed sites_gained_200bp.bed -fo motifs/gained.fa
for i in 1 2 3 4 5; do
    bedtools shuffle -i sites_gained_200bp.bed -g TAIR10.fa.fai -chrom -noOverlapping -seed $i
done > motifs/background.bed
bedtools getfasta -fi TAIR10.fa -bed motifs/background.bed -fo motifs/background.fa
findMotifs.pl motifs/gained.fa fasta motifs/homer -fasta motifs/background.fa -len 8,10 -S 5 -p 4
```

<details>
<summary>Options used</summary>

- `bedtools getfasta -fi genome -bed regions -fo out.fa` - sequences of the
  regions
- `bedtools shuffle -g TAIR10.fa.fai -chrom -noOverlapping -seed $i` - the
  same number of regions of the same size at random positions of the same
  chromosome (`-chrom`), not overlapping each other; `-seed` makes the result
  reproducible (each round a different seed)
- `findMotifs.pl target.fa fasta outdir -fasta background.fa` - motif search
  in FASTA sequences against a background
- `-len 8,10` - lengths of the de novo motifs; `-S 5` - report 5 motifs per
  length; `-p 4` - threads (takes about 30 s)

</details>

Open `motifs/homer/homerResults.html` (de novo motifs, with the most similar
known motif) and `motifs/homer/knownResults.html` (enrichment of known motifs
from the HOMER database).

- What is the best de novo motif? In what fraction of the target and of the
  background sequences is it present? Does it contain the ABRE core `ACGTG`?
- Which known motifs are on top? What family of transcription factors do they
  belong to? (Note that the HOMER database contains a motif learned from the
  ChIP-seq of ABF1 from the same study.)
- The ABRE `ACGTG` is in about 80 % of the gained sites. Why does it also
  appear in the random background?

## Part 3 - Binding and expression

A transcription factor binding at a promoter does not necessarily change the
expression of the gene. The same study measured gene expression by RNA-seq
4 hours after ABA or mock treatment (2 replicates each).

### 3.1 Gene expression with salmon

salmon does not map the reads to the genome; it assigns them directly to
transcripts and estimates how many reads come from each transcript. It needs
an index of the transcript sequences:

```bash
salmon index -t data/TAIR10_transcripts.fa.gz -i salmon_index -p 4
mkdir -p salmon
for s in ABA_rep1 ABA_rep2 mock_rep1 mock_rep2; do
    salmon quant -i salmon_index -l A -r data/rnaseq_$s.fastq.gz -p 4 --validateMappings -o salmon/$s
done
grep -i "mapping rate" salmon/*/logs/salmon_quant.log
head salmon/ABA_rep1/quant.sf
```

<details>
<summary>Options used</summary>

- `salmon index -t transcripts.fa -i index` - index of all transcripts
- `salmon quant -l A` - detect the library type (strandedness) automatically;
  `-r` - single-end reads; `--validateMappings` - check each assignment by an
  alignment (more accurate); `-o` - output directory
- `quant.sf` - for each transcript: length, effective length, TPM (transcripts
  per million) and the estimated number of reads

All options: [salmon documentation](https://salmon.readthedocs.io/).

</details>

### 3.2 Differential expression with DESeq2

The script sums the transcripts of each gene (`tximport`), normalises the
counts and tests ABA vs mock with DESeq2:

```bash
Rscript $D/deseq2.R salmon $D/TAIR10_tx2gene.tsv deseq2_ABA_vs_mock.tsv
head deseq2_ABA_vs_mock.tsv
grep -wE "AT5G52300|AT5G66400|AT1G52690|AT3G18780|AT4G34000" deseq2_ABA_vs_mock.tsv
```

The table contains for each gene: `baseMean` (mean normalised count),
`log2FoldChange` (ABA / mock), `pvalue` and `padj` (FDR). The plots are in
`deseq2_ABA_vs_mock_plots.pdf`.

The file `$D/precomputed/deseq2_ABA_vs_mock_full_data.tsv` is the same
analysis with all reads (25-29 million per sample instead of 4 million).

- How many genes are significantly up- and down-regulated by ABA (padj <
  0.05)? How many in the analysis of the full data?
- Why do you find fewer significant genes with fewer reads? Which genes are
  missed - strongly or weakly expressed ones?
- Are the marker genes you looked at in IGV induced? And ABF3 itself?

### 3.3 Are the genes with gained binding induced by ABA?

The script takes the binding sites up to 1 kb from a TSS, splits the genes
into those with ABA-gained binding and those with unchanged binding, and
compares their expression with all genes:

```bash
Rscript $D/integrate.R sites_closest_gene.tsv deseq2_ABA_vs_mock.tsv
Rscript $D/integrate.R sites_closest_gene.tsv $D/precomputed/deseq2_ABA_vs_mock_full_data.tsv
```

It prints the percentage of induced and repressed genes in each group, and
writes `genes_gained_binding_expression.tsv` (genes with gained binding and
their expression) and `integration_plot.pdf` (log2 fold change of expression
in each group).

- What fraction of all genes is induced by ABA? What fraction of the genes
  with gained ABF3 binding? With unchanged binding?
- Is gained binding a good predictor of induction? Name some reasons why a gene
  with gained binding may not be induced, and why an induced gene may have no
  ABF3 binding.
- Compare the results with your RNA-seq (4 M reads) and with the full data.
  Does the conclusion change?
- ABF3 is one of four ABFs (ABF1-4) with very similar DNA-binding domains, and
  the ABFs are themselves induced by ABA. How would you test whether ABF3 is
  *needed* for the induction of a gene, not just bound to it?
