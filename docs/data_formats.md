## Sequences & reads

### FASTA (.fa, .fna, .faa, .fasta)
FASTA is a simple text-based format for storing **nucleotide** or **protein** sequences. Each sequence consists of a **header line** (starting with `>`), followed by one or more lines of sequence data.
-   **Header line:** Begins with `>`, followed by an identifier and (optionally) a description.
-   **Sequence lines:** Raw sequence characters (DNA: `A, C, G, T`; protein: 20 amino acid codes).
-   Line wrapping is optional, but many tools expect ≤80 characters per line.

#### Example: single-sequence FASTA

```fasta
>seq1 human beta-globin gene
ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGG
CAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG
```
-   `seq1` is the sequence ID.
-   The description *human beta-globin gene* is optional.
-   The sequence follows on the next lines.

#### Example: multi-FASTA (several sequences in one file)

```fasta
>seq1 Homo sapiens beta-globin
ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTT
>seq2 Mus musculus beta-globin
ATGGTGCACCTGACTGCTGAGGAGAAGTCTGCCGTT
>seq3 Drosophila melanogaster hypothetical protein
MKTAYIAKQRQISFVKSHFSRQDILD
```

-   Each sequence starts with its own `>` header.
-   Nucleotide and protein sequences can both be stored in the same FASTA file (though usually they’re separate).
    

#### Notes
-   Common extensions: `.fa`, `.fasta`, `.fna` (nucleotide), `.faa` (amino acid).
-   Used as input for most bioinformatics tools: BLAST, genome assemblers, aligners, gene annotation pipelines.
-   FASTA has no strict rules for metadata—parsers only guarantee they’ll recognize the sequence ID (up to the first whitespace).
   

 
### FAI (FASTA index, .fai)

A **FASTA index** (`.fai`) is a tab-delimited text file that allows fast random access to sequences in a FASTA file. It's created by tools like `samtools faidx` and contains one line per sequence with essential metadata.

#### Format specification
Each line has 5 tab-separated columns:
1. **NAME** — sequence/chromosome name (exactly as in FASTA header, up to first whitespace)
2. **LENGTH** — total length of the sequence in bases
3. **OFFSET** — byte offset in the file where the sequence starts (after the header line)
4. **LINEBASES** — number of bases per line (not counting newline characters)
5. **LINEWIDTH** — total bytes per line (including newline characters)

#### Example
For a FASTA file `genome.fa`:
```fasta
>chr1
ATGCATGCATGCATGCATGC
GCTAGCTAGCTAGCTAGCTA
>chr2
CGTAGCTAGCTA
```

The corresponding `genome.fa.fai` would look like:
```
chr1	40	6	20	21
chr2	12	53	12	13
```

**Explanation:**
- `chr1`: 40 bases total, starts at byte 6, 20 bases per line, 21 bytes per line (20 + newline)
- `chr2`: 12 bases total, starts at byte 53, all on one line

#### Usage
- **Creating an index:** `samtools faidx genome.fa` (creates `genome.fa.fai`)
- **Extracting regions:** `samtools faidx genome.fa chr1:10-30` (extracts bases 10-30 from chr1)
- **Required by:** IGV, GATK, many variant callers and genome browsers

#### Notes
- The FASTA file must have consistent line lengths within each sequence for the index to work correctly.
- Both the `.fa` and `.fai` files must be present in the same directory.
- Changing the FASTA file invalidates the index—regenerate it if you modify the sequence file.


### Sequence dictionary (.dict)
    
### 2bit (UCSC compact FASTA)
    
### FASTQ (.fastq, .fq, usually .gz)

**FASTQ** is the standard format for storing nucleotide sequences along with their **quality scores**. It's widely used for raw sequencing data from Illumina, PacBio, Nanopore, and other sequencing platforms.

#### Format specification
Each sequence record consists of exactly **4 lines**:
1. **Header line:** Starts with `@`, followed by sequence identifier and optional description
2. **Sequence line:** Raw nucleotide sequence (A, C, G, T, N)
3. **Separator line:** Starts with `+`, optionally followed by the same identifier (usually just `+`)
4. **Quality line:** ASCII-encoded quality scores (same length as sequence)

#### Example: single FASTQ record
```fastq
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9IC
```

**Header components:**
- `SRR001666.1` — sequence identifier
- `071112_SLXA-EAS1_s_7:5:1:817:345` — instrument/flowcell information
- `length=36` — optional metadata

#### Quality scores (Phred scores)
- Each character in line 4 encodes a quality score using ASCII encoding
- **Phred quality score (Q)** = -10 × log₁₀(P), where P = probability of base call error
- **Common encoding:** Phred+33 (Sanger/Illumina 1.8+)
  - ASCII 33-126 represents Q scores 0-93
  - `!` = Q0 (50% error), `I` = Q40 (0.01% error)
- **Example conversions:**
  - `!` (ASCII 33) = Q0 = 100% error probability
  - `+` (ASCII 43) = Q10 = 10% error probability
  - `5` (ASCII 53) = Q20 = 1% error probability
  - `?` (ASCII 63) = Q30 = 0.1% error probability
  - `I` (ASCII 73) = Q40 = 0.01% error probability

| ASCII | Symbol | Phred+33 Score | Error Probability |
|-------|--------|----------------|-------------------|
| 33    | !      | 0              | ~100%            |
| 40    | (      | 7              | ~20%             |
| 53    | 5      | 20             | 1%               |
| 63    | ?      | 30             | 0.1%             |
| 73    | I      | 40             | 0.01%            |

#### Multiple records
```fastq
@SEQ_ID_1
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
@SEQ_ID_2
CCCTTCTTGTCTTCAGCGTTTCTCC
+
;;3;;;;;;;;;;;;7;;;;;;;88
```

#### Paired-end FASTQ files
Paired-end sequencing produces two FASTQ files:
- `sample_R1.fastq.gz` — forward reads
- `sample_R2.fastq.gz` — reverse reads
- Read pairs have matching identifiers with `/1` and `/2` suffixes (or `1:` and `2:` in header)

#### Common operations
- **View compressed FASTQ:** `zcat file.fastq.gz | head -4`
- **Count reads:** `zcat file.fastq.gz | wc -l` then divide by 4
- **Quality control:** `fastqc file.fastq.gz`
- **Convert to FASTA:** `seqkit fq2fa file.fastq.gz > file.fasta`

#### Notes
- Almost always **compressed with gzip** (`.fastq.gz`) to save space
- File sizes: typically 1-50 GB per sample (compressed)
- **Quality filtering** removes low-quality reads before downstream analysis
- **Phred+33** is the current standard; older formats (Phred+64, Illumina 1.3-1.7) are obsolete
- Some tools require uncompressed FASTQ; others can read `.gz` directly


### SRA (.sra; raw runs from ENA/NCBI)
    

## Alignments & mappings

### SAM / BAM / CRAM (+ .bai for BAM index, .crai for CRAM index)

**SAM** (Sequence Alignment/Map), **BAM** (Binary Alignment/Map), and **CRAM** are formats for storing how sequencing reads align to a reference genome. They're essential for most genomics workflows including variant calling, RNA-seq analysis, and genome visualization.

#### The three formats

| Format | Type | Size | Use case |
|--------|------|------|----------|
| **SAM** | Text | Very large | Human-readable, rarely used for storage |
| **BAM** | Binary (compressed SAM) | ~3-5× smaller than SAM | Standard format for most analyses |
| **CRAM** | Binary (reference-based) | ~2× smaller than BAM | Long-term storage, requires reference genome |

**In practice:** Almost everyone uses BAM. SAM is mainly for viewing/debugging, CRAM for archiving.

---

## SAM Format (Sequence Alignment/Map)

SAM is a **tab-delimited text format** with two sections:
1. **Header lines** (optional) — start with `@`, contain metadata
2. **Alignment lines** — one line per aligned read

#### Example SAM file
```sam
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:248956422
@PG	ID:bwa	PN:bwa	VN:0.7.17
READ1	99	chr1	10000	60	76M	=	10200	276	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:0	AS:i:76
READ1	147	chr1	10200	60	76M	=	10000	-276	TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:0	AS:i:76
```

#### Header lines (start with @)
- `@HD` — file format version and sort order
- `@SQ` — reference sequence (chromosome) names and lengths
- `@PG` — programs used to create the file

#### Alignment lines: 11 mandatory fields

Each alignment has **11 required tab-separated fields**:

| # | Field | Name | Example | Meaning |
|---|-------|------|---------|---------|
| 1 | QNAME | Query name | READ1 | Read identifier |
| 2 | FLAG | Bitwise flag | 99 | Alignment properties (see below) |
| 3 | RNAME | Reference name | chr1 | Chromosome/contig name |
| 4 | POS | Position | 10000 | Leftmost mapping position (1-based) |
| 5 | MAPQ | Mapping quality | 60 | Confidence score (0-60) |
| 6 | CIGAR | CIGAR string | 76M | Alignment structure |
| 7 | RNEXT | Mate reference | = | Mate's chromosome (= means same) |
| 8 | PNEXT | Mate position | 10200 | Mate's mapping position |
| 9 | TLEN | Template length | 276 | Insert size (for paired-end) |
| 10 | SEQ | Sequence | ACGT... | Read sequence |
| 11 | QUAL | Quality | IIII... | Base quality scores (Phred+33) |

**Plus optional fields** (tags like `NM:i:0`, `AS:i:76`) with additional information.

---

#### Understanding FLAG (column 2)

The **FLAG** is a number that encodes multiple properties as **binary bits**. Common flags:

| Flag | Meaning |
|------|---------|
| 1 | Read is paired |
| 2 | Read mapped in proper pair |
| 4 | **Read unmapped** |
| 8 | Mate unmapped |
| 16 | **Read reverse strand** |
| 32 | Mate reverse strand |
| 64 | **First in pair** |
| 128 | **Second in pair** |
| 256 | Secondary alignment |
| 512 | Failed quality checks |
| 1024 | PCR duplicate |
| 2048 | Supplementary alignment |

**Example:** FLAG = 99
- 99 = 1 + 2 + 32 + 64
- Meaning: paired, proper pair, mate reverse strand, first in pair

**Tip:** Use online FLAG explainers (search "SAM flag decoder") or `samtools flags 99`

---

#### Understanding CIGAR (column 6)

**CIGAR** (Compact Idiosyncratic Gapped Alignment Report) describes how the read aligns to the reference.

**Format:** number + operation letter

| Operation | Letter | Meaning |
|-----------|--------|---------|
| Match/mismatch | M | Aligned bases (match or mismatch) |
| Insertion | I | Bases in read, not in reference |
| Deletion | D | Bases in reference, not in read |
| Skipped | N | Skipped region (for RNA-seq splicing) |
| Soft clip | S | Bases in read, not aligned |
| Hard clip | H | Bases not in read sequence |
| Padding | P | Silent deletion in padded alignment |

**Examples:**
- `76M` — 76 bases match/mismatch (simple perfect alignment)
- `50M2I24M` — 50 matches, 2 inserted bases, 24 matches
- `30M1000N46M` — 30 matches, 1000 bp skipped (intron), 46 matches (RNA-seq)
- `5S71M` — 5 soft-clipped bases, 71 aligned bases

---

#### Mapping quality (MAPQ, column 5)

**MAPQ** = -10 × log₁₀(P), where P = probability the mapping is wrong

| MAPQ | Meaning |
|------|---------|
| 0 | Unreliable (read maps to multiple locations) |
| 1-10 | Low confidence |
| 20 | 1% error probability |
| 30 | 0.1% error probability |
| 40+ | High confidence (99.99%+ correct) |
| 60 | Maximum quality (unique mapping) |

**Filtering:** `samtools view -q 20` keeps only reads with MAPQ ≥ 20

---

## BAM Format (Binary Alignment/Map)

BAM is a **binary compressed version of SAM**. It's much smaller and faster to process.

#### Key features
- **Binary format** — not human-readable, requires tools to view
- **Compressed** with BGZF (blocked gzip) — typically 3-5× smaller than SAM
- **Indexed** — allows fast random access to genomic regions
- **Sorted** — usually sorted by coordinate for efficient processing

#### Common operations

**Convert SAM to BAM:**
```bash
samtools view -b file.sam > file.bam
```

**View BAM as SAM (first 10 reads):**
```bash
samtools view file.bam | head -10
```

**View BAM header:**
```bash
samtools view -H file.bam
```

**Sort BAM by coordinate:**
```bash
samtools sort file.bam -o sorted.bam
```

**Index BAM (creates .bai file):**
```bash
samtools index sorted.bam
# Creates sorted.bam.bai
```

**Extract reads from specific region:**
```bash
samtools view sorted.bam chr1:1000000-2000000
```

**Count mapped reads:**
```bash
samtools view -c -F 4 file.bam
```

**Get alignment statistics:**
```bash
samtools flagstat file.bam
samtools stats file.bam
```

**Filter by mapping quality:**
```bash
samtools view -q 20 -b file.bam > filtered.bam
```

---

## BAM Index (.bai)

A **BAM index** (`.bai` file) enables fast random access to specific genomic regions without reading the entire file.

**Creating an index:**
```bash
samtools index sorted.bam
# Creates sorted.bam.bai (or sorted.bai)
```

**Requirements:**
- BAM must be **sorted by coordinate** (not by read name)
- Both `.bam` and `.bai` files must be in the same directory
- Most tools (IGV, GATK, variant callers) require indexed BAMs

---

## CRAM Format

**CRAM** is a more compressed format that stores differences from a reference genome rather than full sequences.

#### Advantages
- **2× smaller** than BAM (sometimes more)
- Lossless compression (can perfectly reconstruct original data)
- Growing adoption for long-term storage

#### Disadvantages
- **Requires reference genome** to read/write
- Slightly slower to process than BAM
- Less universal tool support (though improving)

#### Usage
```bash
# Convert BAM to CRAM
samtools view -C -T reference.fa file.bam -o file.cram

# Convert CRAM to BAM
samtools view -b -T reference.fa file.cram -o file.bam

# Index CRAM
samtools index file.cram
# Creates file.cram.crai
```

**When to use CRAM:**
- Long-term archiving of alignment data
- When storage space is critical
- For submission to public databases (ENA/SRA prefer CRAM)

---

## Paired-end reads in SAM/BAM

For paired-end sequencing, each read pair appears as **two separate lines** with:
- Same **QNAME** (read identifier)
- Matching **FLAG** bits (64 for first, 128 for second)
- **RNEXT/PNEXT** pointing to mate's location
- **TLEN** showing insert size

**Example paired-end reads:**
```
READ1  99  chr1  10000  60  76M  =  10200  276  [sequence1]  [quality1]
READ1 147  chr1  10200  60  76M  =  10000 -276  [sequence2]  [quality2]
```
- FLAG 99 (64+32+2+1): paired, proper pair, mate reverse, first in pair
- FLAG 147 (128+16+2+1): paired, proper pair, reverse strand, second in pair

---

## Common workflows using SAM/BAM

### 1. Read mapping workflow
```bash
# Map reads to reference genome
bwa mem reference.fa reads_R1.fq.gz reads_R2.fq.gz > aligned.sam

# Convert to BAM, sort, and index
samtools view -b aligned.sam | samtools sort -o sorted.bam
samtools index sorted.bam

# Get alignment statistics
samtools flagstat sorted.bam
```

### 2. Viewing in genome browser
```bash
# BAM must be sorted and indexed
samtools sort file.bam -o sorted.bam
samtools index sorted.bam

# Open sorted.bam and sorted.bam.bai in IGV
```

### 3. Variant calling
```bash
# Most variant callers need sorted, indexed BAM
bcftools mpileup -f reference.fa sorted.bam | bcftools call -mv -o variants.vcf
```

### 4. RNA-seq quantification
```bash
# Count reads per gene
featureCounts -a genes.gtf -o counts.txt sorted.bam
```

---

## Tools that use SAM/BAM/CRAM

**Aligners (create SAM/BAM):**
- **BWA, Bowtie2** — DNA alignment
- **STAR, HISAT2** — RNA-seq alignment
- **minimap2** — long-read alignment

**Manipulation tools:**
- **samtools** — Swiss Army knife (view, sort, index, filter)
- **Picard** — Java tools for BAM manipulation
- **sambamba** — faster alternative to samtools

**Downstream analysis:**
- **GATK, FreeBayes, bcftools** — variant calling
- **featureCounts, HTSeq** — read counting for RNA-seq
- **deepTools, bedtools** — coverage and enrichment analysis
- **IGV, UCSC Genome Browser** — visualization

---

## File sizes

**Example for 30× human genome sequencing:**
- **FASTQ** (raw reads): ~100 GB
- **SAM** (aligned): ~300 GB
- **BAM** (aligned, compressed): ~50-70 GB
- **CRAM** (aligned, reference-compressed): ~25-35 GB

**Storage recommendation:** Keep BAM for active analysis, convert to CRAM for long-term storage.

---

## Important notes

1. **Always sort BAM files by coordinate** before indexing or using with most tools
2. **Always create an index** (.bai or .crai) for sorted alignment files
3. **BAM is the standard** — most tools expect BAM, not SAM
4. **Check alignment statistics** with `samtools flagstat` to verify mapping quality
5. **Coordinate system:** SAM/BAM uses **1-based coordinates** (like GFF, unlike BED)
6. **File names matter:** Keep `.bam` and `.bam.bai` together with matching names
7. **CIGAR N vs D:** N is for long skips (introns in RNA-seq), D is for small deletions

---

## Quick reference: samtools commands

```bash
# View BAM as text
samtools view file.bam

# View header only
samtools view -H file.bam

# Convert SAM to BAM
samtools view -b file.sam > file.bam

# Sort BAM
samtools sort file.bam -o sorted.bam

# Index BAM
samtools index sorted.bam

# Extract region
samtools view sorted.bam chr1:1000-2000

# Count reads
samtools view -c file.bam

# Alignment statistics
samtools flagstat file.bam
samtools stats file.bam

# Filter by quality
samtools view -q 30 -b file.bam > filtered.bam

# Extract mapped reads only
samtools view -F 4 -b file.bam > mapped.bam

# Extract unmapped reads
samtools view -f 4 -b file.bam > unmapped.bam

# Merge BAM files
samtools merge output.bam input1.bam input2.bam

# Convert to CRAM
samtools view -C -T ref.fa file.bam -o file.cram
```


### PAF (minimap2 pairwise mapping)
    
### PSL (BLAT alignments)
    
### BEDPE (paired‑end features)
    
### MAF (Multiple Alignment Format) — *whole‑genome multi‑alignments*
    
### AXT (pairwise genome alignment blocks)
    
### Chain / Net (UCSC liftover alignments)
    

## Variants & genotypes

### VCF / BCF (+ .tbi/.csi index; gVCF for per‑base blocks)
    
### PLINK 1: .bed/.bim/.fam; PLINK 2: .pgen/.pvar/.psam
    
### MAF (Mutation Annotation Format) — *cancer variant summaries* *(name clash with the MAF above)*
    

## Annotations & genomic features

### BED (BED3/6/12), bed12

**BED** (Browser Extensible Data) is a flexible tab-delimited text format for defining genomic regions/features. It's one of the most common formats in genomics, used for genes, peaks, regulatory elements, and any coordinate-based annotations.

#### Key features
- **0-based, half-open coordinates:** start is inclusive, end is exclusive
  - Example: `chr1  100  200` represents bases 101-200 (100 bases total)
- **Tab-delimited** plain text
- **Flexible:** supports 3 to 12+ columns (BED3, BED6, BED12 are most common)

#### BED3 (minimal format)
The minimum required fields:
```
chrom  chromStart  chromEnd
```

**Example:**
```
chr1	1000	2000
chr2	5000	6000
chrX	10000	10500
```

#### BED6 (standard format)
Adds strand, name, and score:
```
chrom  chromStart  chromEnd  name  score  strand
```

**Example:**
```
chr1	1000	2000	feature1	100	+
chr2	5000	6000	feature2	200	-
chrX	10000	10500	feature3	50	.
```

**Field descriptions:**
1. **chrom** — chromosome name (chr1, chr2, etc.)
2. **chromStart** — start position (0-based)
3. **chromEnd** — end position (exclusive)
4. **name** — feature name/label
5. **score** — score (0-1000), often used for visualization intensity
6. **strand** — `+`, `-`, or `.` (unknown/not applicable)

#### BED12 (gene/transcript structure)
Full format including exon structure (used for genes, transcripts):
```
chrom  chromStart  chromEnd  name  score  strand  thickStart  thickEnd  itemRgb  blockCount  blockSizes  blockStarts
```

**Example (gene with 3 exons):**
```
chr1	1000	5000	GENE1	0	+	1200	4800	255,0,0	3	400,300,500	0,1500,3500
```

**Additional fields (7-12):**
7. **thickStart** — start of coding region (CDS start, for display)
8. **thickEnd** — end of coding region (CDS end)
9. **itemRgb** — RGB color for display (e.g., `255,0,0` for red)
10. **blockCount** — number of exons/blocks
11. **blockSizes** — comma-separated list of exon sizes
12. **blockStarts** — comma-separated list of exon start positions (relative to chromStart)

**Interpreting the example:**
- Gene spans chr1:1000-5000 (4000 bp total)
- 3 exons:
  - Exon 1: positions 1000-1400 (400 bp)
  - Exon 2: positions 2500-2800 (300 bp)
  - Exon 3: positions 4500-5000 (500 bp)
- CDS: 1200-4800
- Display in red (255,0,0)

#### Common use cases
- **Peak calling:** ChIP-seq, ATAC-seq peaks (BED6 or narrowPeak)
- **Gene annotations:** genes, exons, transcripts (BED12)
- **Genomic regions:** enhancers, promoters, repeats (BED3-6)
- **Coverage tracks:** depth/coverage information
- **Variant regions:** regions to include/exclude in variant calling

#### Coordinate system gotcha: 0-based vs 1-based
BED uses **0-based, half-open** coordinates:
- The first base of chromosome 1 is position 0
- A region from 0-100 includes bases 1-100 (100 bases)

Compare with **1-based, closed** (GFF/GTF, VCF, SAM):
- The first base is position 1
- A region from 1-100 includes bases 1-100 (100 bases)

**Example conversion:**
| Format | Coordinates | Meaning |
|--------|------------|---------|
| BED | chr1:100-200 | Bases 101-200 |
| GFF | chr1:101-200 | Bases 101-200 |

#### Common operations
```bash
# Sort BED file by coordinate
sort -k1,1 -k2,2n file.bed > sorted.bed

# Merge overlapping regions
bedtools merge -i sorted.bed

# Find overlaps between two BED files
bedtools intersect -a regions.bed -b genes.bed

# Extract regions from genome
bedtools getfasta -fi genome.fa -bed regions.bed

# Convert to other formats
bedtools bed12tobed6 -i genes.bed12
```

#### Tools that use BED
- **bedtools** — Swiss Army knife for BED manipulation
- **IGV, UCSC Genome Browser** — visualization
- **MACS2, HOMER** — peak calling (output BED/narrowPeak)
- **deepTools, bedGraphToBigWig** — coverage analysis

#### Notes
- **No header line** by default (though some tools accept `track` or `browser` lines)
- Fields beyond column 12 are allowed (custom BED extensions)
- **narrowPeak** (BED6+4) and **broadPeak** (BED6+3) are BED variants used for ChIP-seq
- Always specify if your coordinates are 0-based (BED) or 1-based (GFF) to avoid off-by-one errors


### bedGraph
    
### GFF3 / GTF

**GFF3** (General Feature Format version 3) and **GTF** (Gene Transfer Format, also called GFF2) are tab-delimited formats for storing gene annotations and genomic features. They use **1-based, closed** coordinates (unlike BED's 0-based).

---

## GFF3 (.gff3, .gff)

GFF3 is the modern standard for gene annotations, supporting hierarchical relationships between features (gene → transcript → exon → CDS).

#### Format specification
**9 tab-separated columns:**
```
seqid  source  type  start  end  score  strand  phase  attributes
```

#### Column descriptions

1. **seqid** — chromosome/contig name (e.g., `chr1`, `scaffold_1`)
2. **source** — annotation source (e.g., `Ensembl`, `RefSeq`, `GenBank`, `.` for unknown)
3. **type** — feature type (e.g., `gene`, `mRNA`, `exon`, `CDS`)
4. **start** — start position (**1-based, inclusive**)
5. **end** — end position (**1-based, inclusive**)
6. **score** — numeric score (`.` if not applicable)
7. **strand** — `+`, `-`, or `.` (unknown/not applicable)
8. **phase** — reading frame for CDS (0, 1, 2, or `.`)
   - **0** = first base is the start of a codon
   - **1** = second base is the start of a codon
   - **2** = third base is the start of a codon
9. **attributes** — semicolon-separated key=value pairs (e.g., `ID=gene1;Name=ABC1`)

#### Example: complete gene structure
```gff3
##gff-version 3
chr1	Ensembl	gene	1000	5000	.	+	.	ID=gene001;Name=ABC1;biotype=protein_coding
chr1	Ensembl	mRNA	1000	5000	.	+	.	ID=transcript001;Parent=gene001;Name=ABC1-201
chr1	Ensembl	exon	1000	1500	.	+	.	ID=exon001;Parent=transcript001
chr1	Ensembl	exon	3000	3500	.	+	.	ID=exon002;Parent=transcript001
chr1	Ensembl	exon	4500	5000	.	+	.	ID=exon003;Parent=transcript001
chr1	Ensembl	CDS	1200	1500	.	+	0	ID=cds001;Parent=transcript001
chr1	Ensembl	CDS	3000	3500	.	+	0	ID=cds002;Parent=transcript001
chr1	Ensembl	CDS	4500	4800	.	+	0	ID=cds003;Parent=transcript001
```

**Hierarchy:**
- **gene001** has child **transcript001**
- **transcript001** has children: **exon001-003** and **cds001-003**
- Relationships defined by `ID` and `Parent` attributes

#### Key GFF3 attributes
- **ID** — unique identifier for this feature
- **Name** — human-readable name
- **Parent** — ID of parent feature (defines hierarchy)
- **Alias** — alternative names
- **Note** — free-text description
- **Dbxref** — cross-references to databases (e.g., `UniProtKB:P12345`)
- **Ontology_term** — links to ontologies (e.g., `GO:0006915`)

#### Common feature types (column 3)
- **gene** — gene locus
- **mRNA** — messenger RNA transcript
- **exon** — transcribed region
- **CDS** — coding sequence (translated region)
- **five_prime_UTR / three_prime_UTR** — untranslated regions
- **ncRNA, tRNA, rRNA** — non-coding RNAs
- **pseudogene** — non-functional gene copy

#### GFF3 directives (lines starting with ##)
```gff3
##gff-version 3
##sequence-region chr1 1 248956422
##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9606
```
- **##gff-version 3** — required header
- **##sequence-region** — defines sequence boundaries
- **##FASTA** — can embed reference sequences at end of file

---

## GTF (Gene Transfer Format / GFF2)

GTF is an older format still widely used (especially by RNA-seq tools like STAR, StringTie, Cufflinks). It's similar to GFF3 but has different attribute syntax and no formal hierarchies.

#### Format specification
**Same 9 columns as GFF3**, but:
- **Attributes (column 9):** space-separated `key "value";` pairs (not `key=value`)
- No `ID/Parent` system—relationships inferred from shared `gene_id` and `transcript_id`

#### Example: GTF gene structure
```gtf
chr1	HAVANA	gene	1000	5000	.	+	.	gene_id "ENSG00000000001"; gene_name "ABC1"; gene_biotype "protein_coding";
chr1	HAVANA	transcript	1000	5000	.	+	.	gene_id "ENSG00000000001"; transcript_id "ENST00000000001"; gene_name "ABC1"; transcript_name "ABC1-201";
chr1	HAVANA	exon	1000	1500	.	+	.	gene_id "ENSG00000000001"; transcript_id "ENST00000000001"; exon_number "1";
chr1	HAVANA	exon	3000	3500	.	+	.	gene_id "ENSG00000000001"; transcript_id "ENST00000000001"; exon_number "2";
chr1	HAVANA	exon	4500	5000	.	+	.	gene_id "ENSG00000000001"; transcript_id "ENST00000000001"; exon_number "3";
chr1	HAVANA	CDS	1200	1500	.	+	0	gene_id "ENSG00000000001"; transcript_id "ENST00000000001"; exon_number "1";
chr1	HAVANA	CDS	3000	3500	.	+	0	gene_id "ENSG00000000001"; transcript_id "ENST00000000001"; exon_number "2";
chr1	HAVANA	CDS	4500	4800	.	+	0	gene_id "ENSG00000000001"; transcript_id "ENST00000000001"; exon_number "3";
```

#### Key GTF attributes
- **gene_id** — unique gene identifier (required)
- **transcript_id** — unique transcript identifier (required)
- **gene_name** — human-readable gene name
- **transcript_name** — human-readable transcript name
- **exon_number** — exon order in transcript
- **gene_biotype / transcript_biotype** — feature type (protein_coding, lincRNA, etc.)

---

## GFF3 vs GTF: Key differences

| Feature | GFF3 | GTF (GFF2) |
|---------|------|------------|
| Version | 3 | 2 |
| Attributes | `key=value` | `key "value";` |
| Hierarchy | `ID/Parent` system | Inferred via `gene_id`/`transcript_id` |
| Flexibility | More feature types, formal ontology | Limited to gene-centric features |
| Modern support | ✓ Preferred for new annotations | Still common in RNA-seq tools |

---

## Common operations

```bash
# Extract all genes
grep -P "\tgene\t" annotation.gff3

# Extract genes on chr1
awk '$1=="chr1" && $3=="gene"' annotation.gff3

# Convert GFF3 to GTF
gffread annotation.gff3 -T -o annotation.gtf

# Convert GTF to GFF3
gffread annotation.gtf -o annotation.gff3

# Extract gene names
grep -P "\tgene\t" annotation.gff3 | grep -oP 'Name=\K[^;]+'

# Sort by coordinate
sort -k1,1 -k4,4n annotation.gff3 > sorted.gff3

# Extract features in a region
awk '$1=="chr1" && $4>=1000000 && $5<=2000000' annotation.gff3

# Get CDS sequences from genome
gffread -x cds.fa -g genome.fa annotation.gff3

# Get protein sequences
gffread -y proteins.fa -g genome.fa annotation.gff3
```

---

## Tools that use GFF3/GTF

**Annotation sources:**
- **Ensembl, RefSeq, NCBI, GENCODE** — provide GFF3/GTF gene annotations
- **Prokka, Bakta** — bacterial genome annotation (GFF3 output)
- **Augustus, MAKER** — eukaryotic gene prediction

**Analysis tools:**
- **featureCounts, HTSeq** — count RNA-seq reads per gene (need GTF)
- **STAR, HISAT2** — splice-aware RNA-seq aligners (use GTF for known junctions)
- **StringTie, Cufflinks** — transcript assembly (GTF input/output)
- **bedtools** — intersect GFF with BED/other GFF files
- **gffread** — convert, extract, validate GFF/GTF files
- **AGAT** — GFF/GTF manipulation toolkit

**Genome browsers:**
- **IGV, UCSC, Ensembl, JBrowse** — visualize annotations

---

## Important notes

1. **1-based coordinates** (unlike BED's 0-based)
   - chr1:100-200 in GFF = 100 bases (positions 100-200 inclusive)
   - chr1:100-200 in BED = 100 bases (positions 101-200)

2. **Phase (column 8) only matters for CDS**
   - Phase adjusts for split codons across exons
   - Always `.` for non-CDS features

3. **Multiple transcripts per gene**
   - Alternative splicing → same gene, multiple mRNA features
   - Each transcript has its own exon/CDS children

4. **GTF vs GFF3 compatibility**
   - Many tools accept both but may prefer one format
   - Use `gffread` to convert between formats
   - Check tool documentation for format requirements

5. **Comments**
   - Lines starting with `#` are comments (except `##` directives in GFF3)

---

## Which format should I use?

- **GFF3** — modern standard, more flexible, better for complex annotations
- **GTF** — required by many RNA-seq tools (STAR, StringTie, Cufflinks, featureCounts)
- **Both** — keep both versions or convert as needed with `gffread`


### AGP (assembly scaffolding)
    
### genePred (UCSC)
    
### RefFlat
    
### narrowPeak / broadPeak (ChIP/ATAC peaks)
    

## Big/track‑optimized browser formats

### bigWig, bigBed, bigInteract
    
### WIG (legacy)
    
### 2bit (again; reference sequences for browsers)
    
### .sizes (chrom sizes)
    

## Expression & single‑cell

### Counts tables (.tsv/.csv from featureCounts/HTSeq)
    
### Matrix Market: .mtx + barcodes.tsv + features.tsv (10x)
    
### HDF5‑based: .h5ad (AnnData), .loom, .h5 (kallisto/alevin)
    
### quant.sf (Salmon), genes.results / isoforms.results (RSEM)
    

## Comparative genomics / pangenomes

### GFA1/GFA2 (assembly/variation graphs)
    
### .mmi (minimap2 index), .fmd/.bwt (BWA/SA indexes)
    
### msh (Mash MinHash sketches), .jf (Jellyfish k‑mer DB), .kmc\_pre/.kmc\_suf (KMC)
    

## Metagenomics

### Kraken/Kraken2 DBs (.k2d, etc.)
    
### Centrifuge DBs
    
### CAMI profile formats (taxonomic profiles; TSV variants)
    
### BIOM (feature tables for microbiome)
    

## Motifs, profiles & homology search

### BLAST outputs: tabular (outfmt 6), XML (outfmt 5), ASN.1 (outfmt 11); BLAST DBs
    
### DIAMOND: .daa (binary alignments), .dmnd (DB)
    
### HMMER: .hmm (profile HMM), Stockholm (.sto) alignments
    
### MEME/JASPAR/TRANSFAC motif formats (PFM/PWM)
    

## Multiple sequence alignment & phylogenetics

### Aligned FASTA (.afa/.fasta)
    
### Clustal (.aln), Stockholm (.sto), PHYLIP (.phy), MSF
    
### NEXUS (.nex/.nxs)
    
### Newick (.nwk/.tree), PhyloXML (.xml)
    

## 3D genome / conformation capture

### .hic (Juicer)
    
### .cool / .mcool (Cooler)
    
### .pairs / .pairsam (pair lists)
    

## Epigenomics & methylation

### bedMethyl
    
### Bismark cytosine reports, CGmap
    
### bigWig (methylation/coverage tracks)
    

## Structural biology & docking

### PDB (.pdb), mmCIF (.cif), MMTF
    
### MRC/CCP4 (.mrc/.map; EM maps)
    
### SDF / MOL2 (ligands), PDBQT (AutoDock)
    

## Pathways, ontologies & gene sets

### SBML, BioPAX, KGML (KEGG)
    
### GO: OBO (ontology), GAF/GPAD (annotations)
    
### GSEA: GCT, CLS, GMT/GMX (gene sets)
    

## Indices, compression & tabular standards (practical must‑knows)

### bgzip vs gzip; .gzi (bgzip index)
    
### tabix indexes: .tbi/.csi for coordinate‑sorted TAB/VCF/BED‑like files
    
### .fai/.bai/.crai as above
    
### Parquet/Arrow (increasingly used for large genomics tables)
    

---