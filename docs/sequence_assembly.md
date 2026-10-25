# Exercise 3 - Genome assembly of MRSA from Illumina and Nanopore reads

In this exercise you will assemble the genome of the methicillin-resistant
*Staphylococcus aureus* strain USA300 twice - from short Illumina reads and
from long Oxford Nanopore (ONT) reads - and compare the results. All
computations take seconds to a few minutes, so focus on the **interpretation of
the outputs**: what the numbers and plots tell you about the reads, the genome
and the quality of the assembly.

You have Illumina NovaSeq (2×150 bp) and Oxford Nanopore (R10.4.1) reads from a
USA300 isolate. The data contain several problems typical of real sequencing
runs. A closely related
genome - the published genome of USA300 strain FPR3757 - is available, so we
can check what the assemblers got right and what they missed. But be careful:
our isolate is not necessarily identical to the published strain.

Introduction slides - the pipeline, read QC, k-mer spectrum, assembly graphs, QUAST, Flye
output, read pairs in IGV, BUSCO, MRSA, BED files and bedtools:
[assembly_intro.pdf](./slides/assembly_intro.pdf)

Contents:

- [Installation of the conda environments](#installation-of-the-conda-environments)
- [Part 1 - Assembly from short Illumina reads](#part-1---assembly-from-short-illumina-reads):
  [1.1 fastp](#11-quality-control-and-trimming-with-fastp) ·
  [1.2 k-mer spectrum](#12-k-mer-spectrum-and-estimation-of-the-genome-size) ·
  [1.3 MEGAHIT](#13-assembly-with-megahit) ·
  [1.4 k-mer size](#14-the-effect-of-k-mer-size) ·
  [1.5 QUAST](#15-assembly-statistics-with-quast) ·
  [1.6 Bandage](#16-assembly-graph-in-bandage)
- [Part 2 - Assembly from long Nanopore reads](#part-2---assembly-from-long-nanopore-reads):
  [2.1 read statistics](#21-read-statistics) ·
  [2.2 Flye](#22-assembly-with-flye) ·
  [2.3 same assembly?](#23-did-you-get-the-same-assembly-as-your-neighbour) ·
  [2.4 small plasmids](#24-where-are-the-small-plasmids) ·
  [2.5 QUAST vs published genome](#25-comparison-of-the-two-assemblies-with-the-published-genome) ·
  [2.6 reads in IGV](#26-assembly-errors-or-real-differences-reads-mapped-to-the-published-genome) ·
  [2.7 BUSCO](#27-completeness-of-the-assemblies-with-busco)
- [Part 3 - What makes USA300 a dangerous pathogen?](#part-3---what-makes-usa300-a-dangerous-pathogen):
  [3.1 dotplot](#31-whole-genome-dotplot) ·
  [3.2 annotation](#32-genome-annotation) ·
  [3.3 USA300-specific regions](#33-regions-specific-for-usa300) ·
  [3.4 genes in them](#34-genes-in-the-usa300-specific-regions)

Older assembly exercises (Staden, Velvet, hifiasm, MetaCentrum) are in
[sequence_assembly_old.md](./sequence_assembly_old.md).

Data are in `~/Desktop/Bioinformatics/data/genome_assembly2/`:

| File | Content |
|---|---|
| `usa300_R1.fastq.gz`, `usa300_R2.fastq.gz` | Illumina paired-end reads 2×150 bp |
| `usa300_ont.fastq.gz` | Oxford Nanopore reads |
| `usa300_FPR3757_reference.fna.gz`, `usa300_FPR3757_reference.gff.gz` | published genome and annotation of USA300 FPR3757, used **only** for evaluation of the assemblies (Part 2) |
| `precomputed/` | results of the long-running steps, use them only if something fails |

## Installation of the conda environments

Two environments are needed (BUSCO conflicts with some of the other programs).
The installation takes a few minutes:

```bash
conda create -n assembly2 -y --override-channels -c conda-forge -c bioconda \
    python=3.11 fastp kmer-jellyfish megahit flye quast minimap2 samtools seqkit bandage blast bedtools
conda create -n busco -y --override-channels -c conda-forge -c bioconda busco
```

<details>
<summary>Options used</summary>

- `-n assembly2` - name of the new environment
- `-y` - do not ask for confirmation
- `--override-channels -c conda-forge -c bioconda` - take packages only from
  the channels conda-forge and bioconda (bioinformatics software), and ignore
  other channels set in your conda configuration
- the rest - packages to install; `python=3.11` fixes the Python version

</details>

> Note: the k-mer counter Jellyfish is in the package `kmer-jellyfish`. The
> package `jellyfish` is an unrelated Python library.

Prepare a working directory:

```bash
mkdir -p ~/assembly2
cd ~/assembly2
cp ~/Desktop/Bioinformatics/data/genome_assembly2/*.gz .
ls -lh
conda activate assembly2
```

## Part 1 - Assembly from short Illumina reads

### 1.1 Quality control and trimming with fastp

`fastp` checks the quality of the reads and removes adapters and low-quality
ends in one step. The result is an HTML report with statistics before and after
filtering.

```bash
fastp -i usa300_R1.fastq.gz -I usa300_R2.fastq.gz \
      -o trimmed_R1.fastq.gz -O trimmed_R2.fastq.gz \
      --detect_adapter_for_pe -w 4 -h fastp.html -j fastp.json
firefox fastp.html
```

<details>
<summary>Options used</summary>

- `-i`, `-I` - input reads: read 1 and read 2 of each pair
- `-o`, `-O` - output files for the trimmed read 1 and read 2
- `--detect_adapter_for_pe` - find the adapter sequence in the reads
  themselves. Without this option, fastp trims paired-end adapters only
  where the two reads of a pair overlap.
- `-w 4` - number of worker threads (CPU cores) to use (default 3)
- `-h fastp.html` - name of the HTML report (graphs, open it in a browser)
- `-j fastp.json` - name of the JSON report: the same statistics in a
  machine-readable form, e.g. for MultiQC or for your own scripts

All options are described in the [fastp repository](https://github.com/OpenGene/fastp)
or by `fastp --help`.

</details>

Questions:
- How many read pairs are in the dataset? What is the total number of bases?
  What is the approximate sequencing coverage, if the genome size of *S. aureus*
  is about 2.9 Mb?
- What fraction of bases has quality ≥ Q30? What does Q30 mean in terms of the
  probability of an error?
- Which adapter sequence did fastp find? How many reads contained it? Adapters
  are present when the sequenced DNA fragment (insert) is shorter than the
  read - look at the "Insert size estimation" plot. (fastp can estimate the
  insert size only for pairs whose two reads overlap, i.e. for inserts shorter
  than ~270 bp (fastp needs at least 30 bp of overlap); all other pairs are reported as "not overlapped".)
- Look at the "base content" plot after filtering. What is the GC content of
  the reads? (Compare with the human genome, ~41 % GC.)
- Some reads end with a long run of `GGGGG...` (polyG). Count them in the raw
  and in the trimmed reads:
  ```bash
  zcat usa300_R2.fastq.gz | grep -c 'GGGGGGGGGGGGGGGGGGGG$'
  zcat trimmed_R2.fastq.gz | grep -c 'GGGGGGGGGGGGGGGGGGGG$'
  ```
  (`zcat` prints the content of a gzipped file; `grep -c` counts the lines
  matching the pattern; `$` means the end of the line, so the pattern is 20 G
  at the end of a read.)
  The NovaSeq instrument uses two-colour chemistry, in which "no signal" is
  read as G - polyG tails appear when the DNA fragment has ended. fastp
  recognises NovaSeq data from the read names and trims polyG tails
  automatically.

### 1.2 k-mer spectrum and estimation of the genome size

A k-mer spectrum (histogram of k-mer counts) shows how many distinct k-mers
occur once, twice, three times, ... in the reads. K-mers from the genome form a
peak around the average k-mer coverage, k-mers containing sequencing errors
occur only once or twice.

```bash
jellyfish count -C -m 21 -s 50M -t 4 -o reads.jf <(zcat trimmed_R1.fastq.gz trimmed_R2.fastq.gz)
jellyfish histo reads.jf > reads.histo
head -50 reads.histo
```

<details>
<summary>Options used</summary>

- `jellyfish count` - counts all k-mers in the reads and stores the counts in a
  binary file
  - `-C` - canonical k-mers: a k-mer and its reverse complement are counted
    as one, because we do not know from which DNA strand a read comes
  - `-m 21` - k-mer length (k = 21)
  - `-s 50M` - initial size of the hash table (room for 50 million distinct
    k-mers). A good value is about the genome size plus the number of error
    k-mers; if it is too small, Jellyfish still works, only more slowly.
  - `-t 4` - number of threads (CPU cores) to use
  - `-o reads.jf` - output file with the counts
  - `<(zcat ...)` - Jellyfish cannot read gzipped files. `zcat` decompresses
    both read files, and `<( )` passes its output to Jellyfish as if it were
    a file.
- `jellyfish histo` - makes the histogram from the counts: two columns, the
  count (how many times a k-mer was seen) and the number of distinct k-mers
  with this count. By default it stops at count 10,000; the last line
  collects all k-mers seen more often (option `-h` changes this limit).

All options are described in the [Jellyfish repository](https://github.com/gmarcais/Jellyfish)
or by `jellyfish count --help`.

</details>

Open [GenomeScope 2.0](http://genomescope.org/genomescope2.0/), upload
`reads.histo`, set k-mer length 21 and **ploidy 1** (bacteria are haploid).
If the GenomeScope website is not available, use the same model in the
[k-mer Spectrum Lab](https://kavonrtep.github.io/games/kmer-spectrum.html#genomescope)
(tab "Your data: GenomeScope"): open `reads.histo`, set k-mer length 21 and
ploidy 1.

Questions:
- Where is the main peak (k-mer coverage)? What is the estimated genome size?
- What is the meaning of the large number of k-mers with count 1?
- Between the error peak and the main peak there are more k-mers than expected.
  What could they come from? Keep this question in mind until you map the
  reads to the published genome in step 2.6.

### 1.3 Assembly with MEGAHIT

MEGAHIT is a fast and memory-efficient de Bruijn graph assembler. It builds the
graph several times with increasing k-mer size (k = 21, 29, ... 141) and uses
the result of each round in the next one.

```bash
megahit -1 trimmed_R1.fastq.gz -2 trimmed_R2.fastq.gz -t 4 -o megahit
# assembly takes ~1 minute
grep -c ">" megahit/final.contigs.fa
grep ">" megahit/final.contigs.fa | head
```

<details>
<summary>Options used</summary>

- `-1`, `-2` - read 1 and read 2 of the read pairs
- `-t 4` - number of threads (CPU cores)
- `-o megahit` - output directory. MEGAHIT stops with an error if the
  directory already exists - delete it (`rm -r megahit`) before running the
  assembly again.
- Without `--k-list`, MEGAHIT uses k = 21, 29, 39, 59, 79, 99, 119, 141.
- `grep -c ">"` - counts the lines with `>`, i.e. the number of sequences in a
  FASTA file

All options: [MEGAHIT repository](https://github.com/voutcn/megahit) or `megahit --help`.

</details>

The header of each contig contains its length (`len=`) and its average k-mer
coverage (`multi=`).

> Note: the contig names (`k141_27`, ...) are different in every run - your
> neighbour (and your second run) will have the same contigs under different
> names, and a circular contig may start at a different position. Compare
> contigs by their length and coverage, not by their names.

- How many contigs were assembled? What is the length of the longest contig?
- Sort the contigs by coverage:
  `grep ">" megahit/final.contigs.fa | sort -t= -k3,3nr | head`
  (`-t=` splits the header line at `=` signs, `-k3,3nr` sorts by the third
  part, which starts with the coverage, numerically (`n`) from the highest
  (`r`)). Most long
  contigs have similar coverage. Some contigs have a much higher coverage - how
  long are they? What could they be? (Hint: how many copies of the sequence are
  in one bacterial cell?)

### 1.4 The effect of k-mer size

To see why MEGAHIT uses several k-mer sizes, make two assemblies with a single
small and a single large k (each takes a few seconds):

```bash
megahit -1 trimmed_R1.fastq.gz -2 trimmed_R2.fastq.gz -t 4 --k-list 21 -o megahit_k21
megahit -1 trimmed_R1.fastq.gz -2 trimmed_R2.fastq.gz -t 4 --k-list 99 -o megahit_k99
```

<details>
<summary>Options used</summary>

- `--k-list 21` - build the graph with this k only (a comma-separated list,
  e.g. `--k-list 21,41,61`, gives several rounds)

</details>

### 1.5 Assembly statistics with QUAST

```bash
quast megahit_k21/final.contigs.fa megahit_k99/final.contigs.fa megahit/final.contigs.fa \
      -l k21,k99,multi_k -o quast_megahit
firefox quast_megahit/report.html
```

<details>
<summary>Options used</summary>

- the FASTA files - assemblies to compare (any number)
- `-l k21,k99,multi_k` - names of the assemblies in the report, in the same
  order as the files
- `-o quast_megahit` - output directory; `report.html` is the main report

All options: [QUAST manual](https://quast.sourceforge.net/docs/manual.html) or `quast --help`.

</details>

- Compare the number of contigs, the total length, the largest contig and the
  N50 of the three assemblies. Which k-mer size gives the best assembly and
  why? (Think about repeats shorter/longer than k and about k-mers containing
  errors.)
- Is the total length of the assembly in agreement with the genome size
  estimated from the k-mer spectrum?

### 1.6 Assembly graph in Bandage

The contigs are paths in the assembly graph. Convert the final graph of MEGAHIT
(k=141) to a format that Bandage can read and open it:

```bash
megahit_toolkit contig2fastg 141 megahit/intermediate_contigs/k141.contigs.fa > megahit_k141.fastg
Bandage load megahit_k141.fastg &
```

<details>
<summary>Options used</summary>

- `megahit_toolkit contig2fastg 141 ...k141.contigs.fa` - converts the
  contigs of the k = 141 round, together with their connections, to the FASTG
  graph format; the number must be the k of the file. The result is printed,
  `>` saves it to a file.
- `Bandage load file &` - opens the graph in Bandage; `&` runs the program in
  the background, so you can still use the terminal

</details>

Click "Draw graph".
- The graph consists of many separate pieces. Find pieces where several
  contigs meet in one node - these are **repeats** (the assembler cannot
  decide which way to continue). Click on such a node and check its length and
  depth (right panel). Is the depth higher than for the surrounding contigs?
- Other contigs simply end (dead ends). Repeats cannot explain them - what can?
- Can you find a contig that forms a circle? Compare its depth with the rest of
  the graph.

## Part 2 - Assembly from long Nanopore reads

### 2.1 Read statistics

```bash
seqkit stats -a usa300_ont.fastq.gz
```

<details>
<summary>Options used</summary>

- `seqkit stats` - number of sequences, total length, minimal, average and
  maximal length
- `-a` - all statistics: also the quartiles of the read length, N50 and the
  fraction of bases with quality ≥ Q20 and Q30

All seqkit commands: [SeqKit documentation](https://bioinf.shenwei.me/seqkit/) or `seqkit stats --help`.

</details>

- How many reads are there, what is the total number of bases and the
  coverage of the genome?
- What are the mean length, the longest read and the N50 of the reads? Compare
  the read N50 with the length of the longest repeat that you saw in Bandage.

### 2.2 Assembly with Flye

Flye is a long-read assembler based on a repeat graph. `--nano-hq` is the
setting for current ONT reads (R10.4.1, accuracy > 95 %).

```bash
flye --nano-hq usa300_ont.fastq.gz -t 4 -o flye
# assembly takes 3-5 minutes
cat flye/assembly_info.txt
```

<details>
<summary>Options used</summary>

- `--nano-hq` - type of the reads (ONT R10.4 or Guppy5+/Dorado SUP basecalling,
  error rate below 5 %); other types are e.g. `--nano-raw` (older ONT reads),
  `--pacbio-hifi`
- `-t 4` - number of threads
- `-o flye` - output directory; the main results are `assembly.fasta`,
  `assembly_info.txt` and `assembly_graph.gfa`
- Flye does not need the genome size (option `--genome-size` is optional).

All options: [Flye manual](https://github.com/mikolmogorov/Flye/blob/flye/docs/USAGE.md) or `flye --help`.

</details>

<details>
<summary>💡 Flye failed or takes too long?</summary>

Use the precomputed assembly (the same contigs that Flye produces from these
reads; see below why your contig_2 may start at a different position):

```bash
mkdir -p flye
cp ~/Desktop/Bioinformatics/data/genome_assembly2/precomputed/usa300_flye_assembly.fasta flye/assembly.fasta
cp ~/Desktop/Bioinformatics/data/genome_assembly2/precomputed/usa300_flye_assembly_info.txt flye/assembly_info.txt
cp ~/Desktop/Bioinformatics/data/genome_assembly2/precomputed/usa300_flye_assembly_graph.gfa flye/assembly_graph.gfa
```
</details>

The table `assembly_info.txt` contains for each contig: length, mean coverage
(`cov.`), whether it is circular (`circ.`), whether it is a repeat, and its
multiplicity.

- How many contigs were assembled? What is their length and coverage?
- Which contig is the bacterial chromosome? How do you know that it is
  complete?
- What is the second contig? Use its coverage to estimate how many copies of it
  are in one cell.
- Open `flye/assembly_graph.gfa` in Bandage and compare it with the graph of
  the MEGAHIT assembly.

### 2.3 Did you get the same assembly as your neighbour?

Many programs that run on several threads (`-t 4`) are **not deterministic**:
the threads finish their parts of the work in a different order in each run,
and the result can differ slightly. Compare your Flye assembly with your
neighbour's and with the precomputed one:

```bash
cat flye/assembly_info.txt
P=~/Desktop/Bioinformatics/data/genome_assembly2/precomputed
cat $P/usa300_flye_assembly_info.txt
# first 60 bases of each contig
seqkit seq -s -w 0 flye/assembly.fasta | cut -c 1-60
seqkit seq -s -w 0 $P/usa300_flye_assembly.fasta | cut -c 1-60
# align your assembly (query) to the precomputed one (target)
minimap2 -x asm5 $P/usa300_flye_assembly.fasta flye/assembly.fasta | cut -f 1-9
```

The output of minimap2 is a table in the **PAF** format, one line per
alignment (one contig can have several lines). `cut -f 1-9` keeps the first 9
columns:

| column | meaning | here |
|---|---|---|
| 1 | query name | your contig |
| 2 | query length | length of your contig |
| 3, 4 | start and end of the alignment in the query | |
| 5 | strand: `+` same orientation, `-` reverse complement | |
| 6 | target name | contig of the precomputed assembly |
| 7 | target length | |
| 8, 9 | start and end of the alignment in the target | |

Coordinates start at 0 and the end is not included (as in BED), so an
alignment of a whole 1,000 bp contig is `0 1000`. A contig that is the same
circle opened at another position gives two lines: one piece aligns to the
end of the other contig, the rest to its beginning.

<details>
<summary>Options used</summary>

- `P=...` - saves the path to a shell variable, `$P` then stands for it
- `seqkit seq -s -w 0` - prints only the sequences (`-s`), each on one line
  (`-w 0` = no line wrapping); `cut -c 1-60` keeps the first 60 characters
- `minimap2 -x asm5 target query` - aligns one assembly (query) to another
  (target); the preset `asm5` is for very similar genomes (~0.1 % divergence;
  `asm10` and `asm20` allow more). The output is a PAF table, one line per
  alignment; `cut -f 1-9` keeps the first 9 columns.

All options: [minimap2 manual](https://lh3.github.io/minimap2/minimap2.html) or `minimap2` without arguments.

</details>

- Are the chromosomes identical? Do both contig_2 sequences start with the
  same bases? If not, is contig_2 a different sequence, or the same circle
  opened at a different position? (Look at the minimap2 alignment.)
- Why can a circular contig start at any position? Why does it happen to
  contig_2 and not to the chromosome? (Hint: look at the `repeat` and `mult.`
  columns.)
- What may differ between two runs of a non-deterministic program, and what
  must not (sequence, length, completeness of the genome)? What does it mean
  for comparing your results with the teacher's, or for reproducibility of a
  published analysis?
- Optional (4 minutes): Flye has an option `--deterministic`, which runs the
  first step of the assembly on a single thread:
  `flye --nano-hq usa300_ont.fastq.gz -t 4 --deterministic -o flye_det`.
  Run it twice, or compare with your neighbour. Is the result always the same?
  Is it the same as your first assembly?

### 2.4 Where are the small plasmids?

In step 1.3 you found two short contigs (~3.3 kb and ~4.6 kb) with very high
coverage in the MEGAHIT assembly. Are they present in the Flye assembly?

```bash
# extract the two MEGAHIT contigs (replace the names by the names of your contigs)
seqkit grep -p k141_XXX -p k141_YYY megahit/final.contigs.fa > small_contigs.fa
# compare them with the Flye assembly
minimap2 flye/assembly.fasta small_contigs.fa
# how many ONT reads come from them?
minimap2 -x map-ont small_contigs.fa usa300_ont.fastq.gz | cut -f 1,6 | sort -u | cut -f 2 | sort | uniq -c
```

<details>
<summary>Options used</summary>

- `seqkit grep -p NAME` - selects sequences by their name (ID); `-p` can be
  repeated
- `minimap2 target query` - without `-x` minimap2 uses the settings for ONT
  reads; good enough for a quick check
- `minimap2 -x map-ont` - preset for ONT reads
- `cut -f 1,6` - columns 1 (read name) and 6 (contig name) of the PAF table;
  `sort -u` keeps each read only once per contig (one read can have several
  alignments); `cut -f 2 | sort | uniq -c` counts the reads for each contig

</details>

The first minimap2 command prints the full PAF table (no output = no
alignment), one line per alignment:

| column | meaning | here |
|---|---|---|
| 1 | query name | MEGAHIT contig |
| 2 | query length | length of the MEGAHIT contig |
| 3, 4 | start and end of the alignment in the query | |
| 5 | strand (`+` / `-`) | |
| 6 | target name | Flye contig |
| 7 | target length | length of the Flye contig |
| 8, 9 | start and end of the alignment in the target | |
| 10 | number of matching bases (estimate) | |
| 11 | length of the alignment (including gaps) | |
| 12 | mapping quality (0-60; 0 = the query fits equally well elsewhere) | |
| 13+ | extra tags, e.g. `tp:A:P` primary alignment, `dv:f:` estimated divergence | |

Compare columns 3-4 with column 2: how much of the contig is aligned? Without
`-c`, minimap2 does not compute the exact base-by-base alignment, so columns
10-11 are estimates.

The second command (reads vs the small contigs) prints only the summary of
the pipeline: the number of reads and the contig name, e.g. `360 k141_219`.

- Are the small contigs in the Flye assembly? (A short partial alignment of a
  few hundred bp does not mean that the contig is present.) Are they in the
  ONT reads?
- Use [NCBI BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome)
  to find out what they are.
- Long-read assemblers often lose small plasmids (shorter than the typical read
  length). How would you get a complete genome, including all plasmids?

### 2.5 Comparison of the two assemblies with the published genome

Now we will use the published genome of USA300 FPR3757 as a reference. QUAST will
align the contigs to it and report misassemblies, missing parts of the genome
and missing genes:

```bash
quast flye/assembly.fasta megahit/final.contigs.fa -l Flye,MEGAHIT \
      -r usa300_FPR3757_reference.fna.gz -g usa300_FPR3757_reference.gff.gz \
      -o quast_vs_reference
firefox quast_vs_reference/report.html
```

<details>
<summary>Options used</summary>

- `-r` - reference genome; QUAST then reports genome fraction, misassemblies,
  mismatches and indels
- `-g` - annotation of the reference (GFF); QUAST counts the genes that are
  completely or partly present in each assembly
- `-l Flye,MEGAHIT` - names of the assemblies in the report

</details>

- Compare "Genome fraction (%)", "# misassemblies", "# mismatches per 100 kbp",
  "# indels per 100 kbp" and "# genomic features" (genes completely and
  partially present) of both assemblies.
- Open "Icarus" (link at the top of the report) - contig browser. Where are the
  gaps of the MEGAHIT assembly?

Now look at the comparison of the Flye assembly with the reference as a
**dotplot** in Gepard. Gepard reads only uncompressed FASTA files, so
decompress the reference first:

```bash
zcat usa300_FPR3757_reference.fna.gz > usa300_FPR3757_reference.fna
```

<details>
<summary>Options used</summary>

- `zcat file.gz > file` - writes the decompressed content to a new file; the
  gzipped original stays (the same as `gunzip -k`)

</details>

Start Gepard with the shortcut on the Desktop. Load
`usa300_FPR3757_reference.fna` (in `~/assembly2`) as the first sequence and
`flye/assembly.fasta` as the second one, set the word length to 20 and create
the dotplot. Both files contain several sequences (the chromosome and the
plasmids); Gepard puts them one after another, so the plasmids are at the end
of each axis.

- Describe the dotplot: how many pieces does the main diagonal have, and why
  does it not start in the corner of the plot? (Remember step 2.3: where does a
  circular contig start?)
- Which parts of the reference have no line? Compare with the "Genome fraction"
  of Flye in the QUAST report.
- QUAST reported misassemblies in the Flye assembly. Find where they are: the
  file `quast_vs_reference/contigs_reports/all_alignments_Flye.tsv` lists the
  aligned blocks (columns: start and end in the reference, start and end in the
  contig), and the lines between the blocks describe the misassembly. Zoom into
  these regions in Gepard (set the start and end of both sequences to a
  window of ~50-70 kb around the positions and use word length 10). What does
  the dotplot look like at each of them? What kind of rearrangement is it? You
  will check with the reads in step 2.6 whether these are assembly errors or
  real differences.

### 2.6 Assembly errors or real differences? Reads mapped to the published genome

QUAST compares the assemblies with the published genome of USA300 FPR3757 and
reports every difference as an error. But the strain we sequenced is not
necessarily identical to the published one. To decide whether a difference is
an assembly error or a real difference between the genomes, we look at the
reads themselves. Map both Illumina and ONT reads to the published genome:

```bash
# usa300_FPR3757_reference.fna was made in step 2.5 (zcat)
gunzip -k usa300_FPR3757_reference.gff.gz
minimap2 -ax sr -t 4 usa300_FPR3757_reference.fna trimmed_R1.fastq.gz trimmed_R2.fastq.gz \
    | samtools sort -o illumina_to_reference.bam
minimap2 -ax map-ont -t 4 usa300_FPR3757_reference.fna usa300_ont.fastq.gz \
    | samtools sort -o ont_to_reference.bam
samtools index illumina_to_reference.bam
samtools index ont_to_reference.bam
```

<details>
<summary>Options used</summary>

- `gunzip -k` - decompresses the annotation and keeps (`-k`) the gzipped
  original; IGV needs uncompressed files (the reference FASTA was decompressed
  in step 2.5)
- `minimap2 -a` - output in the SAM format (alignments with the read
  sequences) instead of PAF
- `-x sr` - preset for short reads; with two FASTQ files minimap2 aligns them
  as pairs. `-x map-ont` - preset for ONT reads.
- `-t 4` - number of threads
- `| samtools sort -o file.bam` - sorts the alignments by position in the
  reference and saves them in the binary BAM format
- `samtools index` - creates the index (`.bai`) which IGV needs to show a
  region quickly

samtools commands: [samtools documentation](https://www.htslib.org/doc/samtools.html).

</details>

Here minimap2 writes the **SAM** format (`-a`) instead of PAF, and samtools
saves it as BAM. You will look at the alignments in IGV, but you can also see
them as text, e.g. `samtools view illumina_to_reference.bam | head -3`. The
main columns of SAM (one line per read):

| column | meaning |
|---|---|
| 1 | read name (both mates of a pair have the same name) |
| 2 | FLAG - a number that encodes properties of the alignment (paired, mapped, reverse strand, first/second mate, ...); decode it at [Explain SAM flags](https://broadinstitute.github.io/picard/explain-flags.html) |
| 3 | reference sequence name |
| 4 | position of the first aligned base (1-based) |
| 5 | mapping quality (MAPQ, 0-60) |
| 6 | CIGAR - how the read aligns: e.g. `151M` = 151 bases aligned, `100M51S` = 100 aligned and 51 soft-clipped, `D`/`I` = deletion/insertion |
| 7, 8 | reference and position of the mate (`=` = the same reference) |
| 9 | TLEN - insert size: distance between the outer ends of the two mates (negative for the second one) |
| 10, 11 | read sequence and base qualities |

`samtools flagstat illumina_to_reference.bam` summarises the flags (e.g. the
percentage of properly paired reads).

Start IGV, load `usa300_FPR3757_reference.fna` as the genome (`Genomes → Load
Genome from File...`), then the annotation `usa300_FPR3757_reference.gff` and
both BAM files (`File → Load from File...`). In the Illumina track, right-click
and select "View as pairs". IGV colours read pairs with an unexpected distance
(insert size) or orientation of the mates - see
[paired-end alignments](https://igv.org/doc/desktop/#UserGuide/tracks/alignments/paired_end_alignments/)
in the IGV manual. Normal pairs are grey and point towards each other (→ ←).

**A. Regions with low Illumina coverage**
- Zoom out to the whole chromosome. Is the Illumina coverage even? Is the ONT
  coverage even?
- Zoom in to a region where the Illumina coverage drops to (almost) zero and
  look at the sequence. Library preparation of Illumina libraries (PCR
  amplification) is less efficient for AT-rich DNA, and *S. aureus* has a very
  AT-rich genome. Explain now the shoulder in the k-mer spectrum (step 1.2) and
  the gaps in the MEGAHIT assembly.

**B. Structural differences.** QUAST reported misassemblies in the Flye
assembly. Find their positions in the reference: open the Icarus contig
browser from the QUAST report, or look at the alignments of the Flye contigs
to the reference in `quast_vs_reference/contigs_reports/all_alignments_Flye.tsv`
(columns: start and end in the reference, start and end in the contig; lines
between the alignments describe the type of the misassembly). Go to these
positions in IGV:
- One region has no Illumina and no ONT coverage at all. How is it different
  from the AT-rich regions with low coverage from task A? Look at the
  Illumina read pairs around it: what colour are they and what is their insert
  size (click on a read)? Estimate the size of the missing region from the
  insert size of these pairs. Some ONT reads span the whole region - how does
  IGV show it? Which genes of the reference are missing in our isolate?
- At two other positions about 20 kb apart you will find read pairs in which
  both mates point in the same direction (→ → or ← ←), and reads that are
  aligned only partly (soft-clipped, shown in colour after enabling "Show
  soft-clipped bases" in `View → Preferences → Alignments`). The ONT reads are
  split at the same positions (one part of the read is aligned on one strand,
  the other part on the opposite strand). What kind of rearrangement produces
  such pairs? Draw it. How long is the rearranged region?
- Is the Flye assembly wrong, or is our strain different from the published
  one? What would you write in the "Misassemblies" row of QUAST if you were
  sure?

**C. Background noise.** Zoom to any other region of the chromosome (e.g.
10 kb).
- Do you see occasional coloured pairs there as well? How many compared with
  the grey ones?
- Such pairs come from **chimeric fragments** created during library
  preparation (two unrelated DNA pieces ligated together). How can you tell a
  real structural variant from chimeric fragments? (Hint: how many pairs
  support each event?)

### 2.7 Completeness of the assemblies with BUSCO

BUSCO searches the assembly for genes which are present in a single copy in
nearly all genomes of a taxonomic group.

```bash
conda activate busco
busco -i flye/assembly.fasta -o busco_flye -m genome -l staphylococcus_odb12 -c 4
busco -i megahit/final.contigs.fa -o busco_megahit -m genome -l staphylococcus_odb12 -c 4
cat busco_flye/short_summary*.txt busco_megahit/short_summary*.txt
conda activate assembly2
```

<details>
<summary>Options used</summary>

- `-i` - input assembly; `-o` - name of the output directory
- `-m genome` - the input is a genome assembly (other modes: `transcriptome`,
  `proteins`)
- `-l staphylococcus_odb12` - the set of BUSCO genes (lineage dataset) for the
  genus *Staphylococcus*; BUSCO downloads it on the first run. A narrower
  lineage has more genes (`busco --list-datasets` shows all of them).
- `-c 4` - number of CPUs

All options: [BUSCO user guide](https://busco.ezlab.org/busco_userguide.html) or `busco --help`.

</details>

- What is the fraction of complete, fragmented and missing BUSCO genes in each
  assembly?
- According to QUAST, the MEGAHIT assembly misses many more genes than the Flye
  assembly. Why does BUSCO show only a small difference? (Hint: what kind of
  genes does BUSCO use, and what kind of genes are in the missing AT-rich
  regions?)

## Part 3 - What makes USA300 a dangerous pathogen?

USA300 is the most common community-acquired MRSA in North America. We will
compare its genome with the methicillin-sensitive laboratory strain NCTC 8325
and find the genes that USA300 has in addition. Use the Flye assembly - it is
the most complete one.

### 3.1 Whole-genome dotplot

Compare the Flye assembly with NCTC 8325
(`~/Desktop/Bioinformatics/data/genome_assembly/GCF_000013425.1_ASM1342v1_genomic.fna`)
using Gepard or [D-Genies](https://dgenies.toulouse.inra.fr/).

- Is the dotplot line going up or down? Does it start at the corner of the
  plot? Why can the assembly of a circular chromosome start at any position
  and be in any orientation?
- Find places where the line is interrupted - sequences present in only one of
  the two genomes.

### 3.2 Genome annotation

Annotate the Flye assembly with **Bakta** using one of the options:

- **Galaxy**: upload `flye/assembly.fasta` to https://usegalaxy.eu/, run Bakta
  (set "Keep original contig header" to "yes") and download the GFF3 file and
  the annotation table (TSV). Name them `usa300_flye.bakta.gff3` and
  `usa300_flye.bakta.tsv`.
- **Proksee** (https://proksee.ca/): upload `flye/assembly.fasta` and run the
  annotation. You can also add the NCTC 8325 genome as a BLAST track - regions
  without BLAST hits are specific for USA300.
- **Precomputed results** (if the web services are busy): the Bakta annotation
  of the precomputed Flye assembly. Your contig_2 may start at a different
  position (see step 2.3), so use the precomputed assembly together with it:

  ```bash
  cp ~/Desktop/Bioinformatics/data/genome_assembly2/precomputed/usa300_flye.bakta.* .
  cp ~/Desktop/Bioinformatics/data/genome_assembly2/precomputed/usa300_flye_assembly.fasta .
  ```

Look at `usa300_flye.bakta.txt` (summary) and `usa300_flye.bakta.png` (map
of the genome). How many protein-coding genes (CDS) were found?

### 3.3 Regions specific for USA300

Now find the parts of the Flye assembly that have no similar sequence in NCTC
8325. BLASTN finds the similar (shared) parts, bedtools turns them into
intervals, joins them and reports the rest of the genome.

First search the Flye contigs against the NCTC 8325 genome. The output format
`6 qseqid qstart qend` gives one line per alignment: contig name, start and end
of the aligned part of the contig.

```bash
cp ~/Desktop/Bioinformatics/data/genome_assembly/GCF_000013425.1_ASM1342v1_genomic.fna nctc8325.fna
makeblastdb -in nctc8325.fna -dbtype nucl
blastn -task megablast -query flye/assembly.fasta -db nctc8325.fna \
       -evalue 1e-10 -perc_identity 90 -outfmt "6 qseqid qstart qend" > usa300_vs_nctc8325.tsv
head usa300_vs_nctc8325.tsv
wc -l usa300_vs_nctc8325.tsv
```

<details>
<summary>Options used</summary>

- `makeblastdb -in nctc8325.fna -dbtype nucl` - makes a BLAST database from a
  FASTA file (`nucl` = nucleotide sequences); it writes index files next to it
- `blastn -task megablast` - fast search for very similar sequences (word size
  28); good for comparing two strains
- `-query` - sequences to search with; `-db` - database to search in
- `-evalue 1e-10` - report only alignments with an E-value below 10⁻¹⁰
- `-perc_identity 90` - report only alignments with at least 90 % identity
- `-outfmt "6 qseqid qstart qend"` - tabular output (format 6) with the
  selected columns only: query name, start and end of the alignment in the
  query
- `wc -l` - counts lines (here: alignments)

All options: `blastn -help`, or the [BLAST+ manual](https://www.ncbi.nlm.nih.gov/books/NBK279690/).

</details>

bedtools needs the lengths of the contigs (a "genome file"). `samtools faidx`
indexes the assembly; the first two columns of the index are the contig names
and lengths:

```bash
samtools faidx flye/assembly.fasta
cut -f 1,2 flye/assembly.fasta.fai > genome.txt
cat genome.txt
```

<details>
<summary>Options used</summary>

- `samtools faidx` - creates the index `assembly.fasta.fai`: one line per
  sequence with its name, length and position in the file
- `cut -f 1,2` - keeps columns 1 and 2 (name and length)

</details>

Join the alignments into shared regions and take the rest of the genome:

```bash
sort -k1,1 -k2,2n usa300_vs_nctc8325.tsv \
    | bedtools slop -i - -g genome.txt -l 1 -r 0 \
    | bedtools merge -i - -d 1000 > shared.bed
bedtools complement -i shared.bed -g genome.txt > usa300_unique_vs_nctc8325.bed
cat usa300_unique_vs_nctc8325.bed
```

- `sort -k1,1 -k2,2n` - bedtools needs the intervals sorted by contig and start.
- `bedtools slop -l 1 -r 0` - BLAST numbers the bases from 1, the BED format
  from 0 (start 0-based, end the same). Moving the start 1 bp to the left
  converts the BLAST coordinates to BED.
- `bedtools merge -d 1000` - joins overlapping alignments, and also alignments
  less than 1000 bp apart, so short gaps between them disappear.
- `bedtools complement` - everything in `genome.txt` that is not in
  `shared.bed`: the USA300-specific regions longer than 1 kb.

Get the sequences of the regions and their statistics:

```bash
seqkit subseq --bed usa300_unique_vs_nctc8325.bed flye/assembly.fasta > usa300_unique.fasta
seqkit stats usa300_unique.fasta
grep ">" usa300_unique.fasta
```

<details>
<summary>Options used</summary>

- `seqkit subseq --bed file.bed` - extracts the regions listed in a BED file
  from the FASTA file; the names of the new sequences contain the contig name
  and the coordinates of the region

</details>

- How many BLAST alignments are there? How many shared regions
  (`wc -l shared.bed`)?
- How many USA300-specific regions are there, what is their total length and
  which one is the largest?
- Compare the coordinates in `usa300_unique_vs_nctc8325.bed` with the names of
  the sequences in `usa300_unique.fasta`. Why do the start coordinates differ by
  1?

<details>
<summary>💡 BLAST or bedtools failed?</summary>

```bash
cp ~/Desktop/Bioinformatics/data/genome_assembly2/precomputed/usa300_unique_vs_nctc8325.bed .
```
</details>

### 3.4 Genes in the USA300-specific regions

In IGV load `flye/assembly.fasta` as the genome (or
`usa300_flye_assembly.fasta` if you use the precomputed Bakta annotation),
then the Bakta annotation (`usa300_flye.bakta.gff3`),
`usa300_unique_vs_nctc8325.bed` and `shared.bed`. The table
`usa300_flye.bakta.tsv` (columns: contig, type, start, stop, strand, locus
tag, gene, product, ...) is useful for searching, e.g.
`grep -i "leukocidin" usa300_flye.bakta.tsv`.

- Go through the largest regions and look at the annotated genes. Find:
  - **mecA** encoding the penicillin-binding protein **PBP2a**, the cause of
    methicillin resistance. It is part of the mobile element **SCCmec**.
    (Careful: Bakta also reports another, unrelated gene called *mecA* - an
    adapter protein present in all *S. aureus*. Which one is in a
    USA300-specific region?) Which genes are around? Do you see transposases
    or recombinases?
  - Next to SCCmec, genes of the arginine deiminase pathway (arginine
    deiminase, ornithine carbamoyltransferase, carbamate kinase,
    arginine-ornithine antiporter). This is the **arginine catabolic mobile
    element (ACME)**, typical for USA300; it helps the bacteria survive on the
    acidic human skin. (*S. aureus* also has its own chromosomal copy of this
    pathway - find it and check that it is not in a USA300-specific region.)
  - **lukS-PV, lukF-PV** - Panton-Valentine leukocidin, a toxin killing white
    blood cells. What kind of genes surround it outside of the
    USA300-specific region? (Hint: search for "phage" around it.)
- The plasmid (contig_2) is entirely specific for USA300. Which antibiotic
  resistance genes does it carry? Which antibiotics do they affect?
- Methicillin resistance, ACME and PVL are all located on mobile genetic
  elements (SCCmec, ACME, prophage), resistance to other antibiotics on a
  plasmid. What does it tell you about how *S. aureus* strains become
  dangerous?
