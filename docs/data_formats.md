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