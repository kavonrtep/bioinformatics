# ChIP-seq - where are the centromeres of pea chromosome 6?

The **centromere** is the part of a chromosome where the kinetochore assembles
and the spindle microtubules attach during cell division. Centromeres are not
defined by a particular DNA sequence but by a special histone: in centromeric
nucleosomes, histone H3 is replaced by its variant **CENH3** (called CENP-A in
humans). The position of the centromere can therefore be found by ChIP-seq
with an antibody against CENH3.

Most chromosomes have one small CENH3 domain and a single primary
constriction. Pea (*Pisum sativum*) is different: some of its chromosomes have
a very long primary constriction with **several separate CENH3 domains**
(so-called *metapolycentromeres*). In this exercise you will analyse CENH3
ChIP-seq mapped to the assembled chromosome 6 of pea and ask:

1. How many CENH3 domains are there, where are they and how large are they?
2. Can a standard peak caller (MACS3) find them? If not, why not and what
   works better?
3. What DNA sequences are in the CENH3 domains? Are there genes?
4. How are the satellite arrays in the domains organised, and how could a new
   CENH3 domain have arisen?

Everything runs on your computer (or on a VPS with 4 CPUs and 8 GB RAM); all
computing steps together take about 5 minutes.

The data come from [Macas et al. 2023, *PLoS Genetics* 19: e1010633](https://doi.org/10.1371/journal.pgen.1010633)
(ENA [PRJEB54858](https://www.ebi.ac.uk/ena/browser/view/PRJEB54858)), who
assembled the 177.6 Mb chromosome 6 including its whole centromere (GenBank
[OX352004](https://www.ncbi.nlm.nih.gov/nuccore/OX352004)). The ChIP was done on
chromatin digested with micrococcal nuclease (MNase) to single nucleosomes,
using an antibody against pea CENH3. The **input** is the same MNase-digested
chromatin without immunoprecipitation. To make the exercise fast, only the
first 1.5 million reads of each sample were used, shortened to 110 bp.

| File | Content |
|---|---|
| `CENH3_ChIP.fastq.gz` | ChIP with anti-CENH3 antibody (ERR9981080), 1.5 M × 110 bp |
| `input.fastq.gz` | input chromatin (ERR9981081), 1.5 M × 110 bp |
| `CEN6_OX352004.fa.gz` | assembly of pea chromosome 6 (OX352004) |
| `CEN6_genes.gff3.gz` | gene annotation of chromosome 6 |
| `CEN6_satellites.bed` | satellite DNA arrays (tandem repeats) on chromosome 6 |

Introduction slides - centromeres and CENH3, the experiment, the command line
used in the exercise (commands, paths, variables, loops, pipes, positions),
log2 ChIP/input, domain calling, bedtools, dotplots and PAF:
[chipseq_cenh3_intro.pdf](./slides/chipseq_cenh3_intro.pdf)

Note the differences from the [ABF3 exercise](./chip_seq_abf3.md): a
transcription factor binds short sites (a peak of ~200 bp), but CENH3 covers
regions of hundreds of kilobases; and these regions consist of repetitive DNA.

Contents:

- [Installation and data](#installation-and-data)
- [Part 1 - From reads to CENH3 domains](#part-1---from-reads-to-cenh3-domains):
  [1.1 QC](#11-quality-control-with-fastp) ·
  [1.2 mapping](#12-mapping-with-minimap2) ·
  [1.3 mapping statistics](#13-mapping-statistics---reads-from-repeats) ·
  [1.4 coverage tracks](#14-coverage-tracks) ·
  [1.5 MACS3](#15-peak-calling-with-macs3) ·
  [1.6 domains](#16-cenh3-domains-from-the-chipinput-ratio) ·
  [1.7 IGV](#17-cenh3-domains-in-igv)
- [Part 2 - What is inside the CENH3 domains?](#part-2---what-is-inside-the-cenh3-domains):
  [2.1 satellites](#21-satellite-dna-in-the-domains) ·
  [2.2 genes](#22-genes-in-the-domains) ·
  [2.3 repeat units](#23-structure-of-the-satellite-arrays) ·
  [2.4 similarity of arrays](#24-are-the-arrays-similar-to-each-other) ·
  [2.5 duplication](#25-how-did-a-new-cenh3-domain-arise)

## Installation and data

All tools are in the conda environment `chipseq` from the
[ABF3 exercise](./chip_seq_abf3.md#installation-and-data). If you do not have
it, create a smaller one with only the tools needed here (a few minutes):

```bash
conda create -n chipseq -y --override-channels -c conda-forge -c bioconda \
    fastp minimap2 samtools macs3 deeptools bedtools seqkit
conda activate chipseq
```

<details>
<summary>What is in the environment</summary>

- `fastp` - quality control and trimming of reads
- `minimap2`, `samtools` - mapping, BAM files
- `deeptools` - coverage tracks (bigWig, bedGraph), ChIP/input ratio
- `macs3` - peak calling
- `bedtools`, `seqkit` - work with intervals and sequences

The dotplot programs `dotter` and `Gepard` (part 2) are installed on the lab
computers (Gepard has a shortcut on the Desktop).

</details>

The reads and the genome are on [Zenodo](https://zenodo.org/records/23988194) (210 MB).
The annotation and the checksums are in the course repository, in
`~/Desktop/Bioinformatics/data/chip_seq_cenh3/` (if you cloned the repository
elsewhere, change the path in `D=...`):

```bash
D=~/Desktop/Bioinformatics/data/chip_seq_cenh3
mkdir -p ~/chipseq_cenh3/data
cd ~/chipseq_cenh3/data
Z=https://zenodo.org/records/23988194/files
for f in $(cut -d ' ' -f 3 $D/md5sums.txt); do
    wget -c $Z/$f
done
md5sum -c $D/md5sums.txt
cd ~/chipseq_cenh3
```

<details>
<summary>Options used</summary>

- `cut -d ' ' -f 3 md5sums.txt` - the third column (separated by spaces) of
  the checksum file = the file names; the `for` loop downloads them one by one
- `wget -c` - continue an interrupted download instead of starting again
- `md5sum -c` - checks that the downloaded files are not damaged; a file
  reported as `FAILED` must be deleted and downloaded again

</details>

All further commands are run in `~/chipseq_cenh3`. Prepare the genome:
decompress it, index it for samtools and for minimap2, and decompress the
gene annotation for IGV:

```bash
zcat data/CEN6_OX352004.fa.gz > CEN6.fa
samtools faidx CEN6.fa
cut -f 1,2 CEN6.fa.fai
minimap2 -x sr -d CEN6.mmi CEN6.fa
zcat $D/CEN6_genes.gff3.gz > CEN6_genes.gff3
```

<details>
<summary>Options used</summary>

- `zcat file.gz > file` - decompress a gzipped file; `>` writes the output
  to a new file
- `samtools faidx CEN6.fa` - creates the index `CEN6.fa.fai`: for each
  sequence its name, length and position in the file, so that any region can
  be extracted quickly (used by IGV, samtools, bedtools)
- `cut -f 1,2 CEN6.fa.fai` - prints the first two columns of the index =
  sequence name and length
- `minimap2 -x sr -d CEN6.mmi CEN6.fa` - builds the minimap2 index of the
  genome and saves it to `CEN6.mmi`. A mapper does not search the genome
  sequence directly; it first breaks the genome into short words (k-mers,
  here 21 bases) and records where each occurs, so that reads can be placed
  quickly. `-x sr` - the preset for short reads (sets the k-mer size and
  other parameters; the same preset must be used when mapping); `-d` - the
  file where the index is saved. Indexing takes about 10 seconds and
  1.5 GB of memory; the index is built once and reused for every sample.

</details>

- How long is chromosome 6? The whole pea genome has about 4.5 Gb in 7
  chromosomes. What fraction of the genome is in our reference?

## Part 1 - From reads to CENH3 domains

### 1.1 Quality control with fastp

```bash
mkdir -p qc trimmed
for s in CENH3_ChIP input; do
    fastp -i data/$s.fastq.gz -o trimmed/$s.fastq.gz -f 9 -l 80 -w 4 -h qc/$s.html -j qc/$s.json
done
```

<details>
<summary>Options used</summary>

- `for s in ...; do ...; done` - repeats the command for each sample; `$s` is
  replaced by the sample name
- `-i`, `-o` - input and output (trimmed) reads
- `-f 9` - remove the first 9 bases of each read
- `-l 80` - discard reads shorter than 80 bases after trimming
- `-w 4` - threads; `-h`, `-j` - reports in HTML and JSON

</details>

Open the reports (`firefox qc/CENH3_ChIP.html`) and look at the plot
**Base contents** *before filtering*.

- Is the base composition the same along the whole read? What is unusual
  about the first few positions?
- The DNA fragments were produced by MNase, which cuts DNA between the
  nucleosomes. What does the composition of the first bases tell you about
  the enzyme? Why is it a good idea to remove these bases (`-f 9`)?

### 1.2 Mapping with minimap2

```bash
mkdir -p bam
for s in CENH3_ChIP input; do
    minimap2 -ax sr -t 4 CEN6.mmi trimmed/$s.fastq.gz | samtools sort -@ 2 -o bam/$s.bam
    samtools index bam/$s.bam
done
```

Mapping takes about 1.5 minutes per sample.

<details>
<summary>Options used</summary>

- `minimap2 -a` - output in the SAM format; `-x sr` - preset for short
  reads; `-t 4` - threads. `CEN6.mmi` is the index made above (`-d`).
- `|` - the output of minimap2 goes directly to samtools (no SAM file is
  saved)
- `samtools sort -@ 2 -o file.bam` - sorts the alignments by position and
  saves them as BAM
- `samtools index` - index (`.bai`) needed by IGV, MACS3, deeptools, ...

</details>

### 1.3 Mapping statistics - reads from repeats

```bash
samtools flagstat bam/CENH3_ChIP.bam
samtools flagstat bam/input.bam
```

<details>
<summary>How to read the flagstat output</summary>

The output has this form (`N` and `xx.xx` stand for the numbers you get; the
last lines about mates are left out here):

```
N + 0 in total (QC-passed reads + QC-failed reads)
N + 0 primary
0 + 0 secondary
N + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
N + 0 mapped (xx.xx% : N/A)
N + 0 primary mapped (xx.xx% : N/A)
0 + 0 paired in sequencing
...
```

- `N + 0` - reads that passed QC + reads that failed QC. Nothing was flagged
  as failing QC, so the second number is always 0.
- `in total` - the number of alignment records (lines) in the BAM file, not
  reads: one read can have more than one line.
- `primary` - one line per read, so this is the number of reads that went
  into the mapping (compare with the fastp report).
- `secondary` - alternative placements of a read; minimap2 does not report
  them by default, hence 0. `supplementary` - a read whose parts align to
  different places (split alignment) gets one extra line for each extra
  piece.
- `duplicates` - 0 because duplicates were not marked (no markdup step).
- `mapped` - mapped lines and their fraction of all lines. `primary mapped` -
  mapped reads and their fraction of all reads; this is the fraction to use
  when answering the question below.
- `paired in sequencing`, `read1`, `read2`, `properly paired`, `singletons`,
  `with mate mapped ...` - all 0, the reads are single-end. For paired-end
  data these lines show how many pairs map consistently (both mates mapped,
  at the expected distance and orientation).

</details>

- What fraction of the reads was mapped? The reads come from the **whole**
  genome, but the reference is only chromosome 6. Compare the mapped fraction
  with the fraction of the genome in the reference (step "Installation").
  How is this possible? (Hint: more than 80 % of the pea genome consists of
  repetitive sequences.)

The **mapping quality** (MAPQ, column 5 of the BAM file) says how sure the
mapper is that the read comes from this position. A read that fits equally
well to several places gets MAPQ 0. Count all mapped reads and the reads with
MAPQ ≥ 10, in the whole chromosome and in one region with CENH3
(84.7-85.1 Mb):

```bash
samtools view -c -F 4 bam/CENH3_ChIP.bam
samtools view -c -F 4 -q 10 bam/CENH3_ChIP.bam
samtools view -c -F 4 bam/CENH3_ChIP.bam OX352004:84700000-85100000
samtools view -c -F 4 -q 10 bam/CENH3_ChIP.bam OX352004:84700000-85100000
samtools view -c -F 4 bam/input.bam OX352004:84700000-85100000
```

<details>
<summary>Options used</summary>

- `samtools view -c` - only count the reads
- `-F 4` - skip unmapped reads (flag 4 = unmapped)
- `-q 10` - only reads with MAPQ ≥ 10
- `OX352004:84700000-85100000` - only reads in this region (needs the BAM
  index)

</details>

- What fraction of the mapped reads has MAPQ ≥ 10 in the whole chromosome? And
  in the region 84.7-85.1 Mb?
- In the ABF3 exercise, the reads with low MAPQ were removed
  (`--minMappingQuality 10`). What would happen here if we did the same?
- Compare the number of ChIP and input reads in the region 84.7-85.1 Mb (the
  samples have a similar total number of reads). How many times is the
  region enriched in the ChIP?

**Open IGV now** and keep it open: you will add each result to it as soon as
you produce it. Start IGV (Desktop shortcut, or `igv &` in a second
terminal), load `CEN6.fa` as the genome (`Genomes → Load Genome from
File...`), then load (`File → Load from File...`):

- `CEN6_genes.gff3` (the gene annotation)
- `bam/CENH3_ChIP.bam` and `bam/input.bam` (the `.bai` index from 1.2 must be
  next to each BAM)

Type `OX352004:84800000-84820000` into the search box. Reads are drawn only
when you zoom in to a region shorter than about 30 kb; the grey histogram
above the reads is the coverage, computed by IGV on the fly. IGV draws reads
with MAPQ 0 as transparent (white) boxes.

- How many of the reads in this region are transparent? Why? (Compare with
  the counts above.)
- Compare the ChIP and the input track. Is the difference as large as the
  counts suggested?

### 1.4 Coverage tracks

`bamCoverage` computes the coverage of each sample in 1 kb bins (normalised to
CPM, reads per million mapped reads). `bamCompare` divides the ChIP by the
input in 10 kb bins and saves **log2(ChIP/input)**: 0 = the same as input,
1 = two times more, 3 = eight times more, negative = less than input.

```bash
mkdir -p tracks
for s in CENH3_ChIP input; do
    bamCoverage -b bam/$s.bam -o tracks/$s.bw --normalizeUsing CPM --binSize 1000 -p 4
done
bamCompare -b1 bam/CENH3_ChIP.bam -b2 bam/input.bam --operation log2 --binSize 10000 -p 4 \
    -of bedgraph -o tracks/ChIP_vs_input_10kb.bedGraph
head tracks/ChIP_vs_input_10kb.bedGraph
```

<details>
<summary>Options used</summary>

- `-b1`, `-b2` - ChIP and control BAM; the two samples are first scaled to the
  same number of reads
- `--operation log2` - log2 of the ratio ChIP/input in each bin
- `--binSize 10000` - bins of 10 kb; large bins average out the noise, which
  is fine for domains of hundreds of kb
- `-of bedgraph` - output in the text format bedGraph (chromosome, start, end,
  value), which can be read by MACS3, IGV and also by people
- no `--minMappingQuality`: all mapped reads are used (see 1.3)

All options: [bamCompare](https://deeptools.readthedocs.io/en/latest/content/tools/bamCompare.html).

</details>

Load the three new tracks into IGV (`File → Load from File...`):
`tracks/CENH3_ChIP.bw` and `tracks/input.bw` (select both track names,
right-click → "Group Autoscale", so that they share one scale), and
`tracks/ChIP_vs_input_10kb.bedGraph`. Look at the whole chromosome (type
`OX352004` into the search box), then at `OX352004:132500000-138500000`.

- Where on the chromosome is the ChIP coverage high? Does the
  log2(ChIP/input) track show the same regions?
- In the region 132.5-138.5 Mb both the ChIP **and** the input coverage
  increase strongly. Is the region enriched in CENH3? What would you conclude
  if you looked only at the ChIP coverage, without the input?

### 1.5 Peak calling with MACS3

First try the standard tool - MACS3 in the mode for broad peaks:

```bash
mkdir -p peaks
macs3 callpeak -t bam/CENH3_ChIP.bam -c bam/input.bam -f BAM -g 1.78e8 --broad -n CENH3_macs3 --outdir peaks
wc -l peaks/CENH3_macs3_peaks.broadPeak
head peaks/CENH3_macs3_peaks.broadPeak
```

<details>
<summary>Options used</summary>

- `-t` - ChIP (treatment); `-c` - control; `-f BAM` - input format
- `-g 1.78e8` - genome size (our reference: 178 Mb)
- `--broad` - join neighbouring enriched regions into broad peaks (designed
  for histone modifications that cover a few kb)
- output `NAME_peaks.broadPeak`: BED with extra columns (7 fold enrichment,
  8 -log10 p-value, 9 -log10 q-value). All options:
  [MACS3 callpeak](https://macs3-project.github.io/MACS/docs/callpeak.html).

</details>

- How many peaks did MACS3 find? How long are they (columns 2 and 3)?
- How many CENH3 domains do you expect on one chromosome? Is the MACS3 result
  useful?

<details>
<summary>💡 MACS3 failed?</summary>

```bash
mkdir -p peaks
cp $D/precomputed/CENH3_macs3_peaks.broadPeak peaks/
```
</details>

Load `peaks/CENH3_macs3_peaks.broadPeak` into IGV and look at
`OX352004:84000000-89000000`, next to the log2(ChIP/input) track.

- Where are the MACS3 peaks relative to the enriched regions of the
  log2(ChIP/input) track? Do they cover them?

### 1.6 CENH3 domains from the ChIP/input ratio

MACS3 looks at the reads at the resolution of single DNA fragments. With
1.5 million reads on a 178 Mb chromosome, a CENH3 domain of 500 kb is not a
continuous block of reads but many small piles separated by gaps; MACS3
reports each pile as a separate peak. A simpler approach works much better
for large domains: take the log2(ChIP/input) in 10 kb bins from step 1.4, and
report the regions where it is above a cutoff. `macs3 bdgpeakcall` does
exactly this with a bedGraph file:

```bash
macs3 bdgpeakcall -i tracks/ChIP_vs_input_10kb.bedGraph -c 3 -l 20000 -g 50000 \
    --outdir peaks -o CENH3_domains.narrowPeak
tail -n +2 peaks/CENH3_domains.narrowPeak | cut -f 1-3 > CENH3_domains.bed
cat CENH3_domains.bed
```

<details>
<summary>Options used</summary>

- `-i` - input bedGraph
- `-c 3` - cutoff: bins with log2(ChIP/input) ≥ 3 (8 times more reads than in
  the input)
- `-l 20000` - minimal length of a domain (20 kb = two bins)
- `-g 50000` - join domains separated by a gap of up to 50 kb (e.g. a
  single bin with a low value inside a domain)
- `tail -n +2` - print the file from line 2 (skip the first line, a header
  for genome browsers); `cut -f 1-3` - keep chromosome, start and end (BED)

</details>

Get the lengths of the domains by extracting their sequences:

```bash
bedtools getfasta -fi CEN6.fa -bed CENH3_domains.bed -fo CENH3_domains.fa
seqkit fx2tab -n -l CENH3_domains.fa
seqkit stats CENH3_domains.fa
```

- How many CENH3 domains did you find? How long are the shortest and the
  longest one? How much of the chromosome do they cover together?
- The region from the first to the last domain is the primary constriction
  (centromere) of the chromosome. How long is it? What fraction of the
  chromosome is it?
- Repeat `bdgpeakcall` with the cutoff `-c 2` (and a different output name,
  e.g. `-o domains_c2.narrowPeak`). Which regions are added? Look at column 5
  of the narrowPeak file (10 × the highest log2 ratio in the region). Are the
  added regions as strong as the others? Check them in IGV (1.7).

<details>
<summary>💡 bdgpeakcall failed?</summary>

```bash
cp $D/precomputed/CENH3_domains.bed .
cp $D/precomputed/ChIP_vs_input_10kb.bedGraph tracks/
```
</details>

Load `CENH3_domains.bed` into IGV (and `peaks/domains_c2.narrowPeak`, if you
made it).

### 1.7 CENH3 domains in IGV

All results are now in IGV: the genome `CEN6.fa`, `CEN6_genes.gff3`, the two
BAM files, `tracks/CENH3_ChIP.bw` and `tracks/input.bw` (group-autoscaled),
`tracks/ChIP_vs_input_10kb.bedGraph`, `peaks/CENH3_macs3_peaks.broadPeak` and
`CENH3_domains.bed`. If a track is missing, load it (`File → Load from
File...`). You can drag the track names to reorder them; put the domains and
the MACS3 peaks under the log2(ChIP/input) track.

Add one more track, the **satellite DNA arrays** of chromosome 6
(`File → Load from File...`):
`~/Desktop/Bioinformatics/data/chip_seq_cenh3/CEN6_satellites.bed`. Each
array is one BED feature; its name (shown when you click on it, or with
"Expanded" display mode) is the satellite family. Satellite DNA is explained
in Part 2.

Look at the whole chromosome first, then go to these regions (type them into
the search box):

| Region | What to look at |
|---|---|
| `OX352004` | whole chromosome - all CENH3 domains |
| `OX352004:84000000-89000000` | two domains: MACS3 peaks vs. domains |
| `OX352004:84800000-84820000` | reads inside a domain |
| `OX352004:110500000-113500000` | a large satellite array |
| `OX352004:132500000-138500000` | ChIP and input both high (1.4) |

- Where on the chromosome are the CENH3 domains? Do they coincide with the
  satellite arrays? Is every satellite array bound by CENH3?
- Compare the MACS3 peaks and the domains in the region 84-89 Mb. How many
  MACS3 peaks are inside one domain?
- Look again at the reads in the region 84,800,000-84,820,000 and at the
  domain they are in. Could the domain be found from the reads with
  MAPQ ≥ 10 alone?
- Are the regions added with the cutoff `-c 2` real CENH3 domains?

## Part 2 - What is inside the CENH3 domains?

The following analyses are inspired by
[Macas et al. 2023](https://doi.org/10.1371/journal.pgen.1010633), who used
the same data to describe the structure of the pea chromosome 6 centromere.

### 2.1 Satellite DNA in the domains

**Satellite DNA** consists of short sequences (monomers, typically 150-2,000
bp) repeated head-to-tail many times, forming arrays of up to several Mb.
The file `CEN6_satellites.bed` (the track you loaded into IGV in 1.7)
contains the satellite arrays found on chromosome 6; the 4th column is the
name of the satellite family (FabTR-10, FabTR-48, ...).

```bash
# which satellite arrays overlap each domain
bedtools intersect -a CENH3_domains.bed -b $D/CEN6_satellites.bed -wa -wb
# what fraction of each domain is covered by satellites (last column)
bedtools coverage -a CENH3_domains.bed -b $D/CEN6_satellites.bed
# satellite families in the domains and outside of them
bedtools intersect -a $D/CEN6_satellites.bed -b CENH3_domains.bed -u | cut -f 4 | sort | uniq -c
bedtools intersect -a $D/CEN6_satellites.bed -b CENH3_domains.bed -v | cut -f 4 | sort | uniq -c
```

<details>
<summary>Options used</summary>

- `bedtools intersect -a A -b B` - regions of A that overlap B;
  `-wa -wb` - print both the region of A and the overlapping region of B;
  `-u` - print each region of A once if it overlaps anything in B;
  `-v` - print the regions of A that do **not** overlap B
- `bedtools coverage -a A -b B` - for each region of A: number of overlapping
  regions of B, number of bases covered, length of A and the covered fraction
- `cut -f 4 | sort | uniq -c` - take column 4 (family name), and count how
  many times each name occurs

</details>

- What fraction of the CENH3 domains consists of satellite DNA?
- Which satellite families are in the CENH3 domains? Which families are in
  the centromere (between the first and the last domain), but never bound by
  CENH3?
- One array of FabTR-10 is not in any domain - where is it? (Hint: look at
  its position and length, and in IGV.)
- Is the sequence of satellite DNA enough to make a centromere?

### 2.2 Genes in the domains

```bash
wc -l $D/CEN6_genes.bed
bedtools intersect -a $D/CEN6_genes.bed -b CENH3_domains.bed -u | wc -l
# the whole centromere, from the first to the last domain
echo -e "OX352004\t56520000\t138150000" > centromere.bed
bedtools intersect -a $D/CEN6_genes.bed -b centromere.bed -u | wc -l
```

- Compute the gene density (genes per Mb) in the CENH3 domains, in the whole
  centromere and in the whole chromosome. Are there genes in the centromere?
- Look at the genes inside the domains in IGV. Where in the domains are they?

### 2.3 Structure of the satellite arrays

Extract 20 kb from the middle of several CENH3 domains and make a dotplot of
each against itself with `dotter`:

```bash
mkdir -p arrays
samtools faidx CEN6.fa OX352004:79420001-79440000 > arrays/array_79.4Mb.fa
samtools faidx CEN6.fa OX352004:84800001-84820000 > arrays/array_84.8Mb.fa
samtools faidx CEN6.fa OX352004:88000001-88020000 > arrays/array_88.0Mb.fa
samtools faidx CEN6.fa OX352004:108800001-108820000 > arrays/array_108.8Mb.fa
samtools faidx CEN6.fa OX352004:137400001-137420000 > arrays/array_137.4Mb.fa
samtools faidx CEN6.fa OX352004:138000001-138020000 > arrays/array_138.0Mb.fa
dotter arrays/array_84.8Mb.fa arrays/array_84.8Mb.fa &
```

<details>
<summary>Options used</summary>

- `samtools faidx genome.fa chromosome:start-end` - extract a region of the
  genome (positions are counted from 1, both ends included)
- `dotter seq1 seq2 &` - dotplot of two sequences; `&` runs the program in
  the background so you can use the terminal at the same time

</details>

All these regions are FabTR-10 arrays.

- What do the parallel diagonal lines in the dotplot mean? Measure the
  distance between neighbouring lines (move the cursor and read the positions)
  - this is the length of the repeat unit.
- Is the repeat unit the same in all arrays? Do you see lines at several
  distances in some arrays (e.g. faint lines at a short distance and strong
  lines at a longer one)? What could it mean? (Hint: *higher-order repeats* -
  a block of several monomers repeated as a unit.)

### 2.4 Are the arrays similar to each other?

Compare pairs of arrays with each other:

```bash
dotter arrays/array_137.4Mb.fa arrays/array_138.0Mb.fa &
dotter arrays/array_84.8Mb.fa arrays/array_108.8Mb.fa &
dotter arrays/array_79.4Mb.fa arrays/array_84.8Mb.fa &
```

- Which pairs of arrays are very similar and which are only distantly related?
  (Compare how continuous the lines are.) Where on the chromosome are the
  similar arrays?
- All arrays belong to the same family, FabTR-10. Do the monomers in one array
  look more similar to each other than to the monomers of another array? What
  does it tell you about how satellite arrays evolve?

### 2.5 How did a new CENH3 domain arise?

The domain at 79.4 Mb is small, and its sequence is very similar to the domain
at 84.7 Mb (2.4). Make a dotplot of the 13 Mb region around both domains
with Gepard (Desktop shortcut). Use the same file as both sequences:

```bash
samtools faidx CEN6.fa OX352004:73000001-86000000 > region_73-86Mb.fa
```

In Gepard select `region_73-86Mb.fa` as sequence 1 and sequence 2 and
compute the dotplot (it may take a minute or two).

- Besides the main diagonal, there is a long line parallel to it. What does it
  mean? How long is the duplicated segment and how far apart are the two
  copies?
- Where in the duplicated segment are the CENH3 domains at 79.4 Mb and
  84.7 Mb? Why do you think the region at 79.4 Mb has CENH3?

Check the duplication with minimap2 - align the first part of the region to
the second one:

```bash
samtools faidx CEN6.fa OX352004:73000001-79500000 > part1.fa
samtools faidx CEN6.fa OX352004:79400001-86000000 > part2.fa
minimap2 -c -x asm10 part2.fa part1.fa | sort -k11,11nr | head -3 | cut -f 1-12
```

<details>
<summary>Options used</summary>

- `minimap2 -x asm10 target.fa query.fa` - align two similar long sequences
  (preset for assemblies differing by about 1 %); output in the PAF format
- `-c` - compute the exact base-by-base alignment; without it, columns 10
  and 11 are only estimates (takes about 15 s here)
- PAF columns: 1-4 query name, length, start, end; 5 strand; 6-9 target name,
  length, start, end; 10 number of matching bases; 11 alignment length;
  12 mapping quality. Positions start from 0 and are relative to the
  extracted sequences.
- `sort -k11,11nr` - sort by the alignment length, the longest first

</details>

- Where exactly are the two copies (positions on the chromosome)? How long are
  they? What is their identity (column 10 / column 11)?
- What does the identity say about the age of the duplication?
- Macas et al. 2023 propose that the elongated pea centromeres grow by
  duplications like this one and by the expansion of satellite arrays. How
  could a duplication create a **new** CENH3 domain?
