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
    
### Sequence dictionary (.dict)
    
### 2bit (UCSC compact FASTA)
    
### FASTQ (.fastq, .fq, usually .gz)
    
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
    
### bedGraph
    
### GFF3 / GTF
    
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