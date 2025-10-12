# Sequence Alignment

## Dotplots

Dotter is a graphical dotplot program for detailed comparison of two sequences.
It generates a dotplot where the first sequence to be compared is set along the
x-axis, and the second sequence along the y-axis. If two sequences are similar
in two regions, a dot is plotted with an intensity proportional to the
similarity.

### Dotplot programs:
- Desktop applications:
  - dotter (https://www.sanger.ac.uk/science/tools/seqtools), linux only
  - gepard (http://cube.univie.ac.at/gepard)  - suitable for larger sequences, platform independent
  - re-DOT-able https://www.bioinformatics.babraham.ac.uk/projects/redotable/
- Web applications
  - YASS server (https://bioinfo.lifl.fr/yass/index.php), suitable for large sequences or whole chromosomes
  - D-genies (http://dgenies.toulouse.inra.fr/) - whole genome comparison

#### Dotter
To run *dotter*, open terminal window and  type "dotter" followed by two file names:
```bash
dotter file1.fasta file2.fasta
```
- Sequences from dotter program must be provided in FASTA format
- Dotter can be used for protein-proten, DNA-DNA or protein-DNA comparison

#### Dotter keyboard shortcuts
to navigate in `dotter` program use mouse click to select region
- To move cursor in dotter use *arrows* for vertical or diagonal direction
- to move diagonally use `<`, `>`, `[` or `]`. 
- To zoom select area with mouse while holding middle mouse button.

#### Gepard 
Gepard is Java based program, also available on Windows or MacOS. See https://cube.univie.ac.at/gepard
Gepard can be run either of Desktop or using command line:
To run gepard:
```
conda biotools gepard
gepard
```
In Graphical user interface select pair of sequences you want to compare. Test different word size settings

### Exercise 1.1 - Simple self-comparison using dotplot
Make self comparison of the following sequences and identify repetitive sequences.
- How long is repetitive element?
- Is it in forward or reverse-complement orientation.
- How many copies have repetitive element?
- Draw schematic of sequence structure
- What is the position (coordinates) of repetitive region(s)?

Sequences for comparison are located in  `~/Desktop/Bioinformatics/data/dotter_sequences2`
To run self comparison run `dotter` as:
```bash
cd ~/Desktop/Bioinformatics/data/dotter_sequences2
# to get information about sequence run:
seqkit stat seq1.fasta
# to show dotter
dotter seq1.fasta seq1.fasta
...
```

```txt
>seq1
TACCTTGACTTCAGATCTTCCTGAAATAGTGCCGTTTTTATTTTTTGTGATATTTGATAT
ACCTTTAGCTCTGGGCACTATAGCAATTCTATGCATCGATTGCTATAATGGCAGAACATG
GATTCCCCCCCTCCAGACTTAAAGGAATTCGAGAAGACTGGGACTCAAAAAATGTAGAAG
ACCTTGAAGATGGCTACGGACATGCTACTGTTTGTATAACTCCAATTTGCAAATAAGCCA
TTAGAATTAACTTTTTATTTACCAAACGGTCTTCAAACGGATCTCTTGGCATACGCGCCA
TGAATTCAAGAACCAAAGCAAAGTTAAGTACATGATTGCCCATACCTTGTTGCAATATTG
AGTTTCGACGAGTTTTGCATATTAAAAGGTCAAAAACCTGCATGCTACTGTTTGTATAAC
TCCAATTTGCAAAGAGCCGATATTGGCGTTGCAATGGGTATTTCTGGATCTGA
```

```txt
>seq2
TACCTTGACTTCAGATCTTCCTGAAATAGTGCCATGCTACTGTTTGTATAACTCCAATTT
GCAACGTTTTTATTTTTTGTGATATTTGATATACCTTTAGCTCTGGGCACTATAGCAATT
CTATGCATCGATTGCTATAATGGCAGAACATGGATTCCCCCCCTCCAGACTTAAAGGAAT
TCGAGAAGACTGGGACTCAAAAAATGTAGAAGACCTTGAAGATGGCTACGGACATGCTAC
TGTTTGTATAACTCCAATTTGCAAATAAGCCATTAGAATTAACTTTTTATTTACCAAACG
GTCTTCAAACGGATCTCTTGGCATACGCGCCATGAATTCAAGAACCAAAGCAAAGTTAAG
TACATGATTGCCCATACCTTGTTGCAATATTGAGTTTCGACGAGTTTTGCATATTAAAAG
GTCAAAAACCTGCATGCTACTGTTTGTATAACTCCAATTTGCAAAGAGCCGATATTGGCG
TTGCAATGGGTATTTCTGGATCTGACGTTTCTAAGCAGCCGGCAGATATGATTCTATTTG
ATGACAACTTTGCATCAAGTGTTG
```

```txt
>seq3 
ATACCTTTAGCTCTGGGCACTATAGCAATTCTATGCATCGATATCGGCGATGACAACTTT
GCATCAAGTGTTGTACCTTGACTTCAGATCTTCCTGAAATAGTGCCATGCTACTGTTTGT
ATAACTCCAATTTGCAACGTTTTTATTTTTTGTGATATTTGATATACCTTTAGCTCTGGG
CACTATAGCAATTCTATGCATCGATTGCTATAATGGCAGAACATGGATTCCCCCCCTCCA
GACTTAAAGGAATTCGAGAAGACTGGGACTCAAAAAATGTAGAAGACCTTGAAGATGGCT
ACGGACATGCTACTGTTTGTATAACTCCAATTTGCAAATAAGCCATTAGAATTAACTTTT
TATTTACCAAACGGTCTTCAAACGGATCTCTTGGCATACGCGCCATGAATTCAAGAACCA
AAGCAAAGTTAAGTACATGATTGCCCATACCTTGTTGCAATATTGAGTTTCGACGAGTTT
TGCATATTAAAAGGTCAAAAACCTGCATGCTACTGTTTGTATAACTCCAATTTGCAAAGA
GCCGATATTGGCGTTGCAATGGGTATTTCTGGATCTGACGTTTCTAAGCAGCCGGCAGAT
ATGATTCTATTTGATGACAACTTTGCATCAAGTGTTGTCCCTGTAGGTCCATTCTTGTCC
GTAGCCATCTTC
```

```txt
>seq4
ACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAA
GTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCACTGCCCCAACAAACTAATGCC
ATGCAGGACATGTTTTATTTGGGCAAATTCCTGATCGACGAAAGTTTTCAATTGCGCCAG
CGGGAACCCCGGCCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGT
GGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAA
AACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGA
ACTTTTGACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTT
CGTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGGGCA
GTGCCCGGATAGCATCAACGCTGCGCTGATTTGCCGTGGCGAGAAAATGTCGATCGCCAT
TATGGCCGGCGTATTAGAAGCGCGCGGTCACAACGTTACTGTTATCGATCCGGTCGAAAA
ACTGCTGGCAGTGGGGCATTACCTCGAATCTACCCACTGCCCCAACAAACTAATGCCATG
CAGGACATGTTTTATTTGGGCAAATTCCTGATCGACGAAAGTTTTCAATTGCGCCAGCGG
GAACCCCGGCCTCGGACGCTTTGCCGATAAGCTGCCGTCAGAACCACGGGAAAATATCGT
TTATCAGTGCTGGGAGCGTTTTTGCCAGGAACTGGGTAAGCAAATTCCAGTGGCGATGAC
CCTGGAAAAGAATATGCCGATCGGTTCGGGCTTAGGCTCCAGTGCCTGTTCGGTGGTCGC
GGCGCTGATGGCGA
```

### Exercise 1.2 - Identification of repetitive motifs using dotplot

#### example repeats in DNA sequences 
##### Inverted repeat in MITE element
 Triticum aestivum DNA, mobile element MITE contains inverted repeat. Compare the
 sequences of MITE element against itself. 
- What is the position of inverted repeat?
- Can you identify palindromes in the sequence? 
- Is repeat perfect or imperfect?

```bash
cd ~/Desktop/Bioinformatics/data/dotter_sequences/
# make dotplot
dotter inverted_repeat.fasta inverted_repeat.fasta
```

##### Direct repeat in LTR retrotransposon
Boundary of transposable element is defined by long terminal repeat (LTR).
Make dotplot of sequence which contain LTR retrotransposon against itself. 
- What is the length of the whole LTR retrotransposon?
- What is the length of LTR?
- Are LTR sequences identical?
```bash
cd ~/Desktop/Bioinformatics/data/dotter_sequences/
dotter ltr_rt2.fasta ltr_rt2.fasta
```
##### Tandem repeat
Make dotplot of sequence containing tandem repeat.
- how do you interpret dotplot?
- what is a monomer length of tandem repeat(s)
- what is unusual about tandem repeat in `tandem_repeat2.fasta`?

Hint: distances between diagonal parallel lines can be used to estimate length of monomer in tandem repeat.
```bash
cd ~/Desktop/Bioinformatics/data/dotter_sequences
dotter tandem_repeat.fasta tandem_repeat.fasta
dotter tandem_repeat2.fasta tandem_repeat2.fasta
```

### Exercise 1.3 - Comparison of sequences with insertions, deletions, inversions
Use dotter to visualize alignments and identify insertions or deletions in sequences AX02 and AX03 relative to AX01. Sequence are located in `/Desktop/Bioinformatics/data/dotter_sequences/` directory. 
To create dotplot use:
```bash
dotter AX01.fasta AX02.fasta
dotter AX01.fasta AX03.fasta
```
- Examine the dot plots to identify any insertions or deletions (indels) or iversions in AX02 and AX03 with respect to the AX01 sequence. Describe the location and nature of these indels.
- Based on what you see in the dot plots, answer the following:
  - When comparing AX01 and AX02, would you use local or global alignment? Justify your choice.
  - What type of alignment (local or global) would you use to compare AX01 with AX03? Explain your reasoning.

- Performing Global and Local Alignments:
  - For global alignment use: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&BLAST_SPEC=GlobalAln
    (change gap cost in algorithm parameter to 'Existence:5, extension: 2)
  - For local alignemnt use: https://blast.ncbi.nlm.nih.gov/Blast.cgi?BLAST_SPEC=blast2seq&LINK_LOC=align2seq&PAGE_TYPE=BlastSearch 
-  Perform both types of alignments for:
  - AX01 vs AX02
  - AX01 vs AX03

 Compare the results of the dot plots with the global and local alignments. Which method (global or local) gave you a clearer understanding of the differences (e.g., indels) between the sequences? How the choice of alignment (local vs. global) impacts the interpretation of the sequence similarities and differences?

<details>
<summary>💡 FASTA sequences</summary>

```fasta
>AX01
GCGAGCGCCTCGTTCAGCTTGTTGGTGATGATATCTCCCCAGAATTGATACAGATCTTTC
CCTCGGGCATTCTCAAGACGGATCCCCATTTCCAGACGACTTTTGCCAAACTGGCGGATG
TAGCGAAACTGCGATAAGGCTGCATTAAATCGAGCGGGCGGAGTACGCCATACAAGCCGG
AAAGCATTCGCAAATGCTGTTGGGCAAAATCGAAATCGTCTTCGCTGAAGGTTTCGGCCT
GCAAGCCGGTGTAGACATCACCTTTAAACGCCAGAATCGCCTGGCGGGCATTCGCCGGCG
TGAAATCTGGCTGCCAGTCATGAAAGCGAGCGGCGTTGATACCCGCCAGTTTGTCGCTGA
TGCGCATCAGCGTGCTAATCTGCGGAGGCGTCAGTTTCCGCGCCTCATGGATCAACTGCT
GGGAATTGTCTAACAGCTCCGGCAGCGTATAGCGCGTGGTGGTCAACGGGCTTTGGTAAT
CAAGCGTTTTCGCAGGTGAAATAAGAATCAGCATATCCAGTCCTTGCAGGAAATTTATGC
CGACTTTAGCAAAAAATGAGAATGAGTTGATCGATAGTTGTGATTACTCCTGCGAAACAT
CATCCCACGCGTCCGGAGAAAGCTGGCGACCGATATCCGGATAACGCAATGGATCAAACA
CCGGGCGCACGCCGAGTTTACGCTGGCGTAGATAATCACTGGCAATGGTATGAACCACAG
GCGAGAGCAGTAAAATGGCGGTCAAATTGGTAATAGCCATGCAGGCCATTATGATATCTG
CCAGTTGCCACATCAGCGGAAGGCTTAGCAAGGTGCCGCCGATGACCGTTGCGAAGGTGC
AGATCCGCAAACACCAGATCGCTTTAGGGTTGTTCAGGCGTAAAAAGAAGAGATTGTTTT
CGGCATAAATGTAGTTGGCAACGATGGAGCTGAAGGCAAACAGAATAACCACAAG
>AX02
GCGAGCGCCTCGTTCAGCTTGTTGGTGATGATATCTCCCCAGAATTGATACAGATCTTTC
CCTCGGGCATTCTCAAGACGGATCCCCATTTCCAGACGATAAGGCTGCATTAAATCGAGC
GGGCGGAGTACGCCATACAAGCCGGAAAGCATTCGCAAATGCTGTTGGGCAAAATCGAAA
TCGTCTTCGCTGAAGGTTTCGGCCTGCAAGCCGGTGTAGACATCACCTTTAAACGCCAGA
ATCGCCTGGCGGGCATTCGCCGGCGTGAAATCTGGCTGCCAGTCATGAAAGCGAGCGGCG
TTGATACCCGCCAGTTTGTCGCTGATGCGCATCAGCGTGCTAATCTGCGGAGGCGTCAGT
TTCCGCGCCTCATGGATCAACTGCTGGGAATTGTCTAACAGCTCCGGCAGCGTATAGCGC
GTGGTGGTCAACGGGCTTTGGTAATCAAGCGTTTTCGCAGGTGAAATAAGAATCAGCATA
TCCAGTCCTTGCAGGAAATTTATGCCGACTTTAGCAAAAAATGAGAATGAGTTGATCGAT
AGTTGTGATTACTCCTGCGAAACATCATCCCACGCGTCCGGAGAAAGCTGGCGACCGATA
TCCGGATAACGCAATGGATCAAACACCGGGCGCACGCCGAGTTTACGCTGGCGTAGATAA
TCACTGGCAATGGTATGAACCACAGGCGAGAGCAGTAAAATGGCGGTCAAATTGGTAATA
GCCATGCAGGCCATTATGATATCTGCCAGTTGCCACATCAGCGGAAGGCTTAGCAAGGTG
CCGCCGATGACCGTTGCGAAGGTGCAGATCCGCAAACACCAGATCGCTTTAGGGTTGTTC
AGGCGTAAAAAGAAGAGATTGTTTTCGGCATAAATGTAGTTGGCAACGATGGAGCTGAAG
GCAAACAGAATAACCACAAG
>AX03
GCGAGCGCCTCGTTCAGCTTGTTGGTGATGATATCTCCCCAGAATTGATACAGATCTTTC
CCTCGGGCATTCTCAAGACGGATCCCCATTTCCAGACGACTTTTGCCAAACTGGCGGATG
TAGCGAAACTGCGATAAGGCTGCATTAAATCGAGCGGGCGGAGTACGCCATACAAGCCGG
AAAGCATTCGCAAATGCTGTTGGGCAAAATCGAAATCGTCTTCGCTGAAGGTTTCGGCCT
GCATAAATTTCCTGCAAGGACTGGATATGCTGATTCTTATTTCACCTGCGAAAACGCTTG
ATTACCAAAGCCCGTTGACCACCACGCGCTATACGCTGCCGGAGCTGTTAGACAATTCCC
AGCAGTTGATCCATGAGGCGCGGAAACTGACGCCTCCGCAGATTAGCACGCTGATGCGCA
TCAGCGACAAACTGGCGGGTATCAACGCCGCTCGCTTTCATGACTGGCAGCCAGATTTCA
CGCCGGCGAATGCCCGCCAGGCGATTCTGGCGTTTAAAGGTGATGTCTACACCGGCTTGC
CGACTTTAGCAAAAAATGAGAATGAGTTGATCGATAGTTGTGATTACTCCTGCGAAACAT
CATCCCACGCGTCCGGAGAAAGCTGGCGACCGATATCCGGATAACGCAATGGATCAAACA
CCGGGCGCACGCCGAGTTTACGCTGGCGTAGATAATCACTGGCAATGGTATGAACCACAG
GCGAGAGCAGTAAAATGGCGGTCAAATTGGTAATAGCCATGCAGGCCATTATGATATCTG
CCAGTTGCCACATCAGCGGAAGGCTTAGCAAGGTGCCGCCGATGACCGTTGCGAAGGTGC
AGATCCGCAAACACCAGATCGCTTTAGGGTTGTTCAGGCGTAAAAAGAAGAGATTGTTTT
CGGCATAAATGTAGTTGGCAACGATGGAGCTGAAGGCAAACAGAATAACCACAAG
```
</details>

### Exercise 1.4 - Comparison of HER proteins using dotplot
compare sequence of HER proteins - Human epidermal growth receptors using dotter program. 
- download protein sequences of receptor protein-tyrosine kinase from Uniprot:

| accessions | name        |
|------------|-------------|
| P00533     | ERBB2       |
| P21860     | ERBB3       |
| Q15303     | ERBB4       |
| O18735     | ERBB2 (Dog) |
|------------|-------------|

- data for this exercise could  be downloaded from directly from uniprot using the wget command:
```sh 
mkdir -p ~/data/dotter_sequences
cd ~/data/dotter_sequences
wget https://www.uniprot.org/uniprot/P00533.fasta -O ERBB2.fasta
wget https://www.uniprot.org/uniprot/P21860.fasta -O ERBB3.fasta
wget https://www.uniprot.org/uniprot/Q15303.fasta -O ERBB4.fasta
wget https://www.uniprot.org/uniprot/O18735.fasta -O ERBB2_dog.fasta
# check sequence statistics:
seqkit stat *.fasta
```
- compare ErbB2 against itself
- compare *ErbB2* against *ErB3*. Notice the patterns
  in the dotplot and try to find functional domains, for example cysteine rich
  regions B
- try different dotplot stringency using slider on gray-scale strip
- sequences below dotplot correspond to position of blue cross, you can change
  the position of cross either using mouse or  by keyboard shortcuts.
- click on the line and then use arrows to find a good alignment.
- when you identify match, move along diagonal.
- compare *ErbB2* with *ErbB2-dog*. Do you see the same pattern? 
- compare all proteins to all proteins (concatenate all four FASTA file into one sequence using `cat` command)

```sh 
# selfcomparison:
dotter ERBB2.fasta ERBB2.fasta
# compare two sequences against each other:
dotter ERBB2.fasta ERBB3.fasta # human ERBB2 vs human ERB3  (paralogs)
dotter ERBB2.fasta ERBB2_dog.fasta  # human ERBB2 vs dog ERBB2  (orthologs)
# all to all comparison:
cat E*.fasta > all_erb.fasta # first we need concatenated sequences in single fasta file
dotter all_erb.fasta all_erb.fasta
```

Domain structure of *ERBB2* protein:
![ERBB2](../fig/ERBB2.png)

Domain structure of *ERBB3* protein:
![ERBB3](../fig/ERBB3.png)

### Exercise 1.5 - comparison of HOX proteins using dotplot
Download sequence for protein from uniprote, concatenate all HOX protein
sequence into single FASTA file and make all-to-all comparison using dotplot.

```
wget https://rest.uniprot.org/uniprotkb/P49639.fasta -O HOXA1.fasta
wget https://rest.uniprot.org/uniprotkb/P20719.fasta -O HOXA5.fasta
wget https://rest.uniprot.org/uniprotkb/P09067.fasta -O HOXB5.fasta
wget https://rest.uniprot.org/uniprotkb/P14653.fasta -O HOXB1.fasta

cat H*.fasta > all_hox.fasta
dotter all_hox.fasta all_hox.fasta
```
Which pairs of protein are more similar to each other? What part of the proteins is conserved (N or C end)?

### Exercise 1.6 - Locate exon/intron boundaries using dotter (splice sites).
- download sequence AC108130.3 from genbank, save only region from 60000 to 119999 
- download cDNA sequence of GABA A receptor: https://www.ncbi.nlm.nih.gov/nuccore/21265167?report=fasta
- run dotter on these two sequences, identify exon/intron structure
- Are the splice sites consensus splice sites? In vertebrates, the intron starts
  with GT and ends with AG, which are called consensus splice sites.
- What's going on at the 3' end of the cDNA?
- download GABAA1 protein sequence - https://www.ncbi.nlm.nih.gov/protein/27808653?report=fasta
- make dotter of genomic dna vs protein sequence
- what is different, why is the protein alignment shorter than cDNA

```sh
dotter AC108130.3.fna BC030696.1.fna  # genome vs cDNA
dotter AC108130.3.fna GBRA1_HUMAN.fna   # genome vs protein 
```
- Data are also available in ~/Desktop/bioinformatics/data/dotter_sequences
- When using `dotter` to compare  DNA to protein, DNA sequence must be in forward orientation!

### Exercise 1.7 - Identification of insertions, deletions, duplications  - more complex example
Compare two genomic regions `a_region` and `b_region`
- first do self comparison for each sequence
- then compare `a_region` against `b_region`
- What you can say about these genomic regions? Are there any insertions,
  duplications or deletions?

```bash
cd ~/Desktop/Bioinformatics/data/dotter_sequences/
dotter a_region.fasta a_region.fasta
dotter b_region.fasta b_region.fasta
dotter a_region.fasta b_region.fasta
```

### Exercise 1.8 - Identifying overlaps and creating a "sequence assembly" using dot plots
Make dotplot from following sequences stored in file:
```
~/Desktop/Bioinformatics/data/dotter_sequences/dna_examples/overlaping_sequences.fasta
```

- Sequence file contain SeqA, SeqB and SeqC. 
- Get information about sequences in file using `seqkit` program:
```bash
seqkit stat overlaping_sequences.fasta
seqkit fx2tab --length --name overlaping_sequences.fasta
```

- Generate an all-to-all dot plot using dotter. This will compare each sequence (seqA, seqB, and seqC) against the others in a single analysis.
- Examine the resulting dot plot to identify overlaps among seqA, seqB, and seqC.
- Using the overlaps you identified, Create a simple schematic to show how the sequences fit together into a longer continuous sequence (contig).

<details>
<summary>💡 Hint</summary>

![scheme_dotter_overlap](../fig/scheme_dotter_overlap.png)

</details>

### Exercise 1.9 - Compare two genomic regions with dotter
Use a dot plot to compare two large genomic regions and identify structural variation.

```
cd ~/Desktop/Bioinformatics/data/dotter_sequences/
dotter genomeA_part.fasta genomeB_part.fasta
```
- How would characterize difference between genomes?
- For larger sequences, the dot plot may appear "noisy" due to the presence of short repetitive sequences scattered throughout the genome. Identify regions in the dot plot that appear messy or contain many small dots. Zoom into this areas and identify what kind of sequences are causing these 'noise'

### Exercise 1.10 - Whole genome comparison with `Gepard` program
 *Mycoplasma hyopneumoniae* is a bacterial pathogen that primarily affects pigs, causing enzootic pneumonia, a chronic respiratory disease.The genome of Mycoplasma hyopneumoniae is relatively small can be analyzed using `Gepard` dotplot program. Download genomic sequences of two strains, characterize them with `seqkit` program  and compare them using `Gepard`. You will find shortcut to `Gepard` program on Desktop. 

```txt
https://zenodo.org/record/4485547/files/mycoplasma-232.fasta
https://zenodo.org/record/4485547/files/mycoplasma-7422.fasta
```

<details>
<summary>💡 Details</summary>

To download sequences use `wget` command:
```bash
cd ~/Desktop/
mkdir mycoplasma
cd mycoplasma
wget https://zenodo.org/record/4485547/files/mycoplasma-232.fasta
wget https://zenodo.org/record/4485547/files/mycoplasma-7422.fasta
# check sequence statistics:
# use seqkit program to get information about sequences
seqkit stat *.fasta
```
</details>


- You can  use G and H (slow movement) or J and K (fast navigation) to slide along the current diagonal
- Do you observe any diagonal lines in the dot plot other than the main diagonal? What do these lines indicate about the genomic organization (e.g., inversions, duplications)?
 - Recalculate dotplot with word length 100
 - Are there any regions of divergence (areas with no or few dots) between the two genomes? What could these regions represent (e.g., strain-specific genes, deletions, insertions)?
 - How many deletions in you can observe in strain 232?

### Exercise 1.11 - Identification of problems in sequences from SRA database

#### Illumina data : SRR2911427  (Migratory locust WGS)

- download sequences from SRA database using `fastq-dump` command line program 
- for documentation see https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump )
- inspect data and create dotplot

```bash
mkdir ~/tmp
cd ~/tmp
fastq-dump -X 20 --split-files --fasta SRR2911427
# -X 20 = download 20 sequences only
# --fasta = convert sequences to fasta format
# --split-file = create two files one for each pair
# SRR291142 - accession ID
ls -l
cat SRR2911427_1.fasta
dotter SRR2911427_1.fasta SRR2911427_2.fasta
```

what does it mean? 


#### Illumina data : SRR453021 (Nicotian repanda - WGS)

```bash
fastq-dump -X 50 --split-files  --fasta SRR453021
dotter SRR453021_2.fasta SRR453021_2.fasta 
dotter SRR453021_1.fasta SRR453021_1.fasta 
dotter SRR453021_1.fasta SRR453021_2.fasta 
```
select repeated sequences using dotter and search with NCBI blast:

https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome

## Pairwise alignment
### Sequences for alignment:
-  `ERB2_HUMAN` : http://www.uniprot.org/uniprot/P04626.fasta   
-  `EGFR_DROME` : http://www.uniprot.org/uniprot/P04412.fasta   
-  `Unknown protein` : http://www.uniprot.org/uniprot/Q8SZW1.fasta
(Receptor tyrosine-protein kinase erbB-2, Epidermal growth factor receptor)

To download sequences use either web browser of try to use `wget` command in terminal:
```bash
mkdir ~/data/erb
cd ~/data/erb
wget http://www.uniprot.org/uniprot/P04626.fasta
wget http://www.uniprot.org/uniprot/P04412.fasta
wget http://www.uniprot.org/uniprot/Q8SZW1.fasta
```

### Exercise 2.1 - Compare global and local alignments
- global alignment is performed by program `needle`
  - http://www.bioinformatics.nl/cgi-bin/emboss/help/needle
- for local alignment use program `water`, 
  - http://www.bioinformatics.nl/cgi-bin/emboss/help/water

- Programs `needle` and `water` are available from command line or from EBI web interface: http://www.ebi.ac.uk/Tools/emboss/
- Sequences for alignments are located in directory `~/Desktop/bioinformatics/data/alignment_sequences`
- compare ERB2 (P04626.fasta) vs EGFR (P04412.fasta) using `needle` and then using `water` using command lne programs:
```bash
# command example:
needle P04626.fasta P04412.fasta
water P04626.fasta P04412.fasta
```
same programs are also available from web interface:
- https://www.ebi.ac.uk/Tools/psa/emboss_water/
- https://www.ebi.ac.uk/Tools/psa/emboss_needle/

- compare ERB2 (P04626.fasta) vs Unknown protein (Q8SZW1.fasta) using `needle` and then using `water`
- what is difference between local and global alignments?
- what happened what gap penalty is increased to 20 and extend_penalty to 5 when using local alignment
- what happened with global alignment if you change `end gap panalty` setting.
- by default BLOSUM62 scoring matrix is used, what happend when you use PAM10?
- compare these protein sequence using `dotter`

```bash
# command line example using PAM10
water P04626.fasta P04412.fasta -datafile EPAM10
```
#### differences between PAM10 and BLOSUM62 matrices
PAM10 : ftp://ftp.ncbi.nih.gov/blast/matrices/PAM10
BLOSUM62 : ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62

#### using blast (blast2seq) to create local aligment for two sequences:
https://blast.ncbi.nlm.nih.gov/Blast.cgi?BLAST_SPEC=blast2seq&LINK_LOC=align2seq&PAGE_TYPE=BlastSearch
blast2seq can be used instead of `needle`. It also provide graphical view of alignment and non-interactive dotplot. Use blast2 seq on  `P04626.fasta` and  `P04412.fasta` sequences and explore results. Compare alignments and dotplot.

You can paste either AA sequences to the blast form or you can use just accession ID (P046256, P04412).
### Exercise 2.2 - Pairwise alignment using NCBI blast
Compare two sequences of human Hexokinase and yeast Hexokinase using NCBI blast (You will have to use option Align two or more sequences). 

Compare sequences using global and local alignment:
- For global alignment use: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&BLAST_SPEC=GlobalAln
- For local alignemnt use: https://blast.ncbi.nlm.nih.gov/Blast.cgi?BLAST_SPEC=blast2seq&LINK_LOC=align2seq&PAGE_TYPE=BlastSearch 

Note: Switch to Protein tab!

Describe results of BLAST and answer following questions:
+ Sequence Identity: What is the percent identity between human Hexokinase and yeast Hexokinase?
+ Alignment Coverage: What percentage of the human Hexokinase and yeast Hexokinase sequences is covered by the alignment
+ Sequence of human Hexokinase is longer that yeast Hexokinase,  explain why is it so. Hint: explore Dot Plot tab in BLAST output.
```text
>HXK1_HUMAN Hexokinase-1
MIAAQLLAYYFTELKDDQVKKIDKYLYAMRLSDETLIDIMTRFRKEMKNGLSRDFNPTAT
VKMLPTFVRSIPDGSEKGDFIALDLGGSSFRILRVQVNHEKNQNVHMESEVYDTPENIVH
GSGSQLFDHVAECLGDFMEKRKIKDKKLPVGFTFSFPCQQSKIDEAILITWTKRFKASGV
EGADVVKLLNKAIKKRGDYDANIVAVVNDTVGTMMTCGYDDQHCEVGLIIGTGTNACYME
ELRHIDLVEGDEGRMCINTEWGAFGDDGSLEDIRTEFDREIDRGSLNPGKQLFEKMVSGM
YLGELVRLILVKMAKEGLLFEGRITPELLTRGKFNTSDVSAIEKNKEGLHNAKEILTRLG
VEPSDDDCVSVQHVCTIVSFRSANLVAATLGAILNRLRDNKGTPRLRTTVGVDGSLYKTH
PQYSRRFHKTLRRLVPDSDVRFLLSESGSGKGAAMVTAVAYRLAEQHRQIEETLAHFHLT
KDMLLEVKKRMRAEMELGLRKQTHNNAVVKMLPSFVRRTPDGTENGDFLALDLGGTNFRV
LLVKIRSGKKRTVEMHNKIYAIPIEIMQGTGEELFDHIVSCISDFLDYMGIKGPRMPLGF
TFSFPCQQTSLDAGILITWTKGFKATDCVGHDVVTLLRDAIKRREEFDLDVVAVVNDTVG
TMMTCAYEEPTCEVGLIVGTGSNACYMEEMKNVEMVEGDQGQMCINMEWGAFGDNGCLDD
IRTHYDRLVDEYSLNAGKQRYEKMISGMYLGEIVRNILIDFTKKGFLFRGQISETLKTRG
IFETKFLSQIESDRLALLQVRAILQQLGLNSTCDDSILVKTVCGVVSRRAAQLCGAGMAA
VVDKIRENRGLDRLNVTVGVDGTLYKLHPHFSRIMHQTVKELSPKCNVSFLLSEDGSGKG
AALITAVGVRLRTEASS
>HXKA_YEAST Hexokinase-1
MVHLGPKKPQARKGSMADVPKELMDEIHQLEDMFTVDSETLRKVVKHFIDELNKGLTKKG
GNIPMIPGWVMEFPTGKESGNYLAIDLGGTNLRVVLVKLSGNHTFDTTQSKYKLPHDMRT
TKHQEELWSFIADSLKDFMVEQELLNTKDTLPLGFTFSYPASQNKINEGILQRWTKGFDI
PNVEGHDVVPLLQNEISKRELPIEIVALINDTVGTLIASYYTDPETKMGVIFGTGVNGAF
YDVVSDIEKLEGKLADDIPSNSPMAINCEYGSFDNEHLVLPRTKYDVAVDEQSPRPGQQA
FEKMTSGYYLGELLRLVLLELNEKGLMLKDQDLSKLKQPYIMDTSYPARIEDDPFENLED
TDDIFQKDFGVKTTLPERKLIRRLCELIGTRAARLAVCGIAAICQKRGYKTGHIAADGSV
YNKYPGFKEAAAKGLRDIYGWTGDASKDPITIVPAEDGSGAGAAVIAALSEKRIAEGKSL
GIIGA
```

<details>
<summary>💡 Hint</summary>

- This size difference reflects a fundamental evolutionary relationship: mammalian hexokinases, including human HK1, evolved through **gene duplication** and **fusion** of an ancestral 50 kDa hexokinase gene similar to modern yeast hexokinase. The mammalian enzyme essentially represents a "doubled" version of the yeast enzyme, with both halves sharing extensive sequence homology to yeast hexokinase.
- Human HK1: Only the C-terminal domain retains catalytic activity, while the N-terminal domain has lost its catalytic function despite maintaining the structural framework. The N-terminal domain contains a critical amino acid substitution - lysine to glutamate in the ATP-binding site
- Allosteric Regulation Mechanisms, The regulatory mechanisms differ significantly between the two enzymes:
  - Human HK1 exhibits dual inhibition mechanisms by glucose-6-phosphate (G6P):
  - Active site inhibition: G6P competes directly with ATP at the C-terminal catalytic domain
  - Allosteric inhibition: G6P binds to the regulatory N-terminal domain, causing conformational changes that indirectly displace ATP from the active site
</details>

### Exercise 2.3 - Comparing P53 and P63 Protein Sequences
Use NCBI BLAST to compare the protein sequences of human P53 and P63. These proteins share the DNA-binding domain and oligomerization domains but have distinct transactivation domains, reflecting their different roles in cellular processes.

Use the NCBI BLAST's "Align two or more sequences" (https://blast.ncbi.nlm.nih.gov/Blast.cgi?BLAST_SPEC=blast2seq&LINK_LOC=align2seq&PAGE_TYPE=BlastSearch ) option to align these two sequences. From BLAST results:
- Identify the shared  domain in the alignment. What is the percent identity/similarity within this domain?
- What part of the P53 and P63  are likely transactivation and oligomerization domains
- Run search with p53 and p63 as queries against Conserved Domain Database https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi. Was you conclusion about position of shared domains correct?

P53 and P63 sequence in fasta format:
```text
>sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens OX=9606 GN=TP53 PE=1 SV=4
MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK
SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE
RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS
SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP
PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG
GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD
>sp|Q9H3D4|P63_HUMAN Tumor protein 63 OS=Homo sapiens OX=9606 GN=TP63 PE=1 SV=1
MNFETSRCATLQYCPDPYIQRFVETPAHFSWKESYYRSTMSQSTQTNEFLSPEVFQHIWD
FLEQPICSVQPIDLNFVDEPSEDGATNKIEISMDCIRMQDSDLSDPMWPQYTNLGLLNSM
DQQIQNGSSSTSPYNTDHAQNSVTAPSPYAQPSSTFDALSPSPAIPSNTDYPGPHSFDVS
FQQSSTAKSATWTYSTELKKLYCQIAKTCPIQIKVMTPPPQGAVIRAMPVYKKAEHVTEV
VKRCPNHELSREFNEGQIAPPSHLIRVEGNSHAQYVEDPITGRQSVLVPYEPPQVGTEFT
TVLYNFMCNSSCVGGMNRRPILIIVTLETRDGQVLGRRCFEARICACPGRDRKADEDSIR
KQQVSDSTKNGDGTKRPFRQNTHGIQMTSIKKRRSPDDELLYLPVRGRETYEMLLKIKES
LELMQYLPQHTIETYRQQQQQQHQHLLQKQTSIQSPSSYGNSSPPLNKMNSMNKLPSVSQ
LINPQQRNALTPTTIPDGMGANIPMMGTHMPMAGDMNGLSPTQALPPPLSMPSTSHCTPP
PPYPTDCSIVSFLARLGCSSCLDYFTTQGLTTIYQIEHYSMDDLASLKIPEQFRHAIWKG
ILDHRQLHEFSSPSHLLRTPSSASTVSVGSSETRGERVIDAVRFTLRQTISFPPRDEWND
FNFDMDARRNKQQRIKEEGE
```
## Multiple sequence alignment
### Exercise 3.1 - Multiple sequence alignment - Cyclin-dependent kinase
Cyclin-dependent kinases (CDKs) are a group of enzymes that regulate the
progression of the cell cycle by adding phosphate groups to other proteins, a
process called phosphorylation. They are activated by binding to regulatory
proteins called cyclins, which undergo cyclic changes in abundance and activity
throughout the cell cycle.

Create multiple sequence alignment for group of CDKs from human and mouse. Use
program `mafft`.  use default setting. Before running
`mafft` check help documentation using `mafft --help`

```bash
mkdir -p ~/data/cdk
cd ~/data/cdk
cp ~/Desktop/Bioinformatics/data/alignment_sequences/CDK/cdk.fasta .
# inspect sequence using less command
less cdk.fasta
# alignment using mafft program
mafft --help
mafft cdk.fasta > cdk_mafft_aligned.fasta
# alterantive alignment using muscle program
muscle --help
muscle -align cdk.fasta -output cdk_muscle_aligned.fasta
# some muscle version have different options:
# muscle -in cdk.fasta -out cdk_muscle_aligned.fast
seqkit stat *.fasta

# evaluate alignment using transitive consistency score
t_coffee -infile cdk_mafft_aligned.fasta  -evaluate
t_coffee -infile cdk_muscle_aligned.fasta  -evaluate
# the above steps will generate html file with TCS results, inspect them in firefox
```

Inspect alignment from mafft using *Jalview* program. 
- Try different coloring schemes - clustal, percentage identity, hydrophobicity
  and by conservation
  - *Clustal*:  It colors the sequences according to the amino acid type and the level of conservation within the alignment. Highly conserved residues are highlighted in specific colors, making it easier to spot conserved motifs or domains across different sequences. This scheme is derived from the Clustal series of programs, which are widely used for sequence alignment.
  - *Zappo*: The Zappo coloring scheme colors amino acids based on their physicochemical properties, such as charge, polarity, and hydrophobicity.
  - *Taylor*: , like Zappo, colors amino acids based on their physicochemical properties but uses a different set of colors. It is designed to make it easier to distinguish between different types of amino acids by using a broader and more intuitive color palette.
  - The *Hydrophobicity* coloring scheme highlights amino acids based on their hydrophobic or hydrophilic properties.
  - Percentage Identity: This scheme colors residues based on the percentage of sequences in the alignment that have the same amino acid at a given position. Highly conserved positions (where most sequences share the same amino acid) are colored differently from variable positions.
- Jalview is using several metrics calculated from score:
  - *Conservation*: How similar or different are the amino acids at this position?" If all the amino acids are the same, it's considered highly conserved (very important and likely critical for the protein's function). If the amino acids are different but still similar in their chemical properties, it's also considered conserved, but to a lesser extent
  - Quality: measure of base on blosum62 score
  - Consensus - by default it shows histogram, can be switched to sequence logo
  - Occupancy
- *Overview window* (View->Overview window) shows complete alignment
- How are sequences in the alignment sorted? Reorder sequences using Calculate->Sort->By Pairwise Identity.
- Try to identify the most conserved regions.
- What are the coordinates of most conserved region related to *CDK14* sequence.
- Compare this conserved regions with conserved regions which can be identified
  using *conserved domain database*.  Use this link for search
  https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi and CDK14 sequence below.
```text
>CDK14
MCDLIEPQPAEKIGKMKKLRRTLSESFSRIALKKDDTTFDEICVTKMSTRNCQGMDSVIK
PLDTIPEDKKVRVQRTQSTFDPFEKPANQVKRVHSENNACINFKTSSTGKESPKVRRHSS
PSSPTSPKFGKADSYEKLEKLGEGSYATVYKGKSKVNGKLVALKVIRLQEEEGTPFTAIR
EASLLKGLKHANIVLLHDIIHTKETLTLVFEYVHTDLCQYMDKHPGGLHPDNVKLFLFQL
LRGLSYIHQRYILHRDLKPQNLLISDTGELKLADFGLARAKSVPSHTYSNEVVTLWYRPP
DVLLGSTEYSTCLDMWGVGCIFVEMIQGVAAFPGMKDIQDQLERIFLVLGTPNEDTWPGV
HSLPHFKPERFTLYSSKNLRQAWNKLSYVNHAEDLASKLLQCSPKNRLSAQAALSHEYFS
DLPPRLWELTDMSSIFTVPNVRLQPEAGESMRAFGKNNSYGKSLSNSKH
```
*Transitive consistency score* : TCS is an alignment evaluation score that makes it possible to identify in an MSA the most correct positions. It has been shown that these positions are the most likely to be structurally correct and also the most informative when estimating phylogenetic trees. The TCS evaluation and filtering procedure is implemented in the T-Coffee package and can be used to evaluate and filter any third party multiple sequence alignment  
### Exercise 3.2 - Multiple alignment from HSPB8 proteins
Create MSA for set of orthologs of HSPB8 protein (Heat shock protein beta-8) and identify conserved regions.

Make copy of fasta file and then rename fasta headers:
```bash
cd 
mkdir -p data/hspb8
cd data/hsbb8
cp ~/Desktop/Bioinformatics/data/alignment_sequences/HSP8.fasta .

```
Create alignment using `mafft` program.

```bash
mafft HSP8.fasta > HSP8_aln.fasta
```

Open resulting alignment in `Jalview` program.
- Inspect alignment, Try different coloring schemes. (see https://www.jalview.org/help/html/colourSchemes/index.html)
- What part of proteins is conserved? Use `Colour -> By Conservation`  to change coloring threshold.
- Compare conserved part with domains annotation
  - Go to https://www.ncbi.nlm.nih.gov/protein/NP_055180.1
  - Select `analyze this sequence/identify conserved domains`
  - Will you be able to identify conserved domain if you use only mouse, cow, pig and human sequences?
  - In Jalview. select subset of sequences(mammals) and create alignment again from 'Web Service -> Alignment -> mafft with defaults'. Is conserved domain still visible in the new alignment?

### Exercise 3.3 - Alignment of protein isoforms, alignment editing
Investigate the alignment of 11 alternatively-spliced gene products from the human erythrocyte membrane protein band 4.1 (EPB41) gene, focusing on how different multiple sequence alignment (MSA) programs handle the dataset. The aim is to compare the performance of three popular MSA programs—MAFFT, MUSCLE, and ClustalW—when aligning sequences that differ only by deletions due to alternative splicing.

Correct alignment of isoforms will contain only matches and gaps, no mismatches!

- Sequences can be obtained from `../data/alignment_sequences/epb41.fasta`
- Open the JalView desktop application and load the unaligned sequences to visually inspect their similarities and differences.
- Use the JalView web services menu to access the MAFFT, MUSCLE, and ClustalW alignment services. Perform an MSA with each program using the EPB41 isoform dataset.
- Keep the results accessible for comparison, either by keeping the tabs/windows open or by saving the output files.
- Compare the three alignments using JalView to assess how each program performed. Remember, an ideal alignment for isoforms should contain only matches and gaps, with no mismatches, since the sequences are identical except for the presence of deletions.
- Set coloring and sort sequences based on pairwise comparison.
- Evaluate the alignments to determine which program produced the best results based on the criteria of correctly handling deletions without introducing mismatches. In Consensus histogram set option *Ignore Gaps in Consensus*. This will help you identify problematic regions.
- Choose one of the alignments you believe has the potential to be corrected easily manually. Using JalView, edit this alignment to resolve any issues, aiming to create an ideal alignment that reflects the expected pattern of matches and gaps for alternatively spliced isoforms.
- Inspect structural annotation of this gene in UCSC genome browser - https://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=2356525077_wGXvGEoZ9M0QWUR4riuRK94xqF93&db=hg38&position=lastDbPos
  (enter EPB41 to the search window)

    
<details>
<summary>💡 Using Jalview</summary>

## Jalview: quick orientation

-   Two edit modes:
    -   **Normal mode** (default): mostly mouse-based editing; use Shift/Ctrl + drag to add/remove gaps.
    -   **Cursor mode**: keyboard-centric; toggle with **F2**. 

### Mode & help
-   **F2** – Toggle Cursor mode on/off. 
-   **F1** – Open built-in help. 
-   **Esc** – Clear selection (and cancel partial commands in Cursor mode). 

### Gaps (Cursor mode)
-   **Space** – Insert a gap at cursor; prefix with a number to insert *n* gaps (e.g., `12` then **Space**)
-   **Delete** or **Backspace** – Remove a gap at cursor; prefix number removes *n* gaps
-   **Ctrl+Space** or **Shift+Space** – Insert gaps across a **group** (all sequences in the group); same for group delete with **Ctrl/Shift+Delete/Backspace**. (Win uses Shift+Space; macOS supports Ctrl+Space too.)
    
### Navigation (Cursor mode)
-   **Arrow keys** – Move cursor; hold **Shift** to jump between aligned region ↔ gap runs. 
-   **p S** – Jump to sequence *p* (e.g., `7 S`).
-   **p C** – Jump to column *p*.
-   **p P** – Jump to residue *p* in current sequence.
-   **p1, p2 ↩︎** – Jump to column *p1*, sequence *p2* (e.g., `13,5` then **Enter**). 

### Making a rectangular selection (Cursor mode)

-   Move to top-left → **Q** (marks start)
-   Move to bottom-right → **M** (marks end)  
    Now edits apply within the box. 

### Column/sequence operations

-   **Ctrl+L** – Remove columns left of leftmost column marker.
-   **Ctrl+R** – Remove columns right of rightmost column marker.
-   **Ctrl+E** – Remove gapped columns (columns containing only gaps).
-   **Ctrl+Shift+E** – Remove *all* gaps.

### Clipboard & undo/redo

-   **Ctrl+X / Ctrl+C / Ctrl+V** – Cut / Copy / Paste to current alignment.
-   **Ctrl+Shift+V** – Paste to **new** alignment window.
-   **Ctrl+Z / Ctrl+Y** – Undo / Redo sequence edits (gap insertions/deletions, trims, etc.).

### Hiding & selecting columns (handy while editing)

-   **B / Alt+B / Ctrl+B** – Add highlighted columns / invert add / toggle column selection marks.
-   **H / Shift+H / Ctrl+H / Ctrl+Shift+H** – Hide/reveal columns or sequences; “hide all but selection”.

</details>

### Exercise 3.4 - Multiple sequence alignment of the mitochondrial 16S gene
The mitochondrial 16S gene is a widely studied genetic marker in molecular
biology, which is used for species identification, phylogenetic analysis, and
evolutionary studies. 16S gene codes for a RNA subunit of the mitochondrial
ribosome and contains many regions with high substitution rates.
We will use *MAFFT* to align the two sets of sequences, and visualize the
resulting alignments with program called *AliView*. Alignments can be edited
manually or automatically with the software *BMGE*, which determines the most
reliable alignment positions based on the proportion of missing data and their
entropy score.

make new directory and copy sequences. Each sequence is iin one file. We will
concatenate to single multi FASTA file using `cat` command

```bash
mkdir MSA_16s
cd MSA_16s
cat ~/Desktop/bioinformatics/data/alignment_sequences/16s/*.fasta > 16s.fasta
# get information about sequences
seqkit stat 16s.fasta
```
Inspect resulting file with `less` command. 

Align sequence using `mafft` program, at first use default setting. Before running
`mafft` check help documentation using `mafft --help`

```bash
mafft --help
# align 16s with defaults
mafft 16s.fasta > 16s_aln.fasta
# explore output with less command
less 16s_aln.fasta
seqkit stats 16s_aln.fasta
```

Inspect alignment using `Aliview` program.
```bash
~/Desktop/bioinformatics/bin/aliview 16s_aln.fasta
```
Inspect alignment. By default, `mafft` keep order of sequences in the alignment
same as in input file. Close `aliview` and rerun `mafft` with `reorder` option. 
```
mafft --reorder 16s.fasta > 16s_aln.fasta
```
inspect alignment in `Aliview`

#### Manual editing of alignment : 

In the AliView window of 16s_aln.fasta, place the cursor on the sequence that
bridges the first of the two large gaps near the end of the alignment (around
bp 2000) and zoom in (ctrl + mouse wheel) until you can see the labels of the
individual bases. You'll see that the taxon responsible for these gaps is called
*Balistecaprisc*. It appears that the sequence alignment for this taxon is correct
up to this gap , but that the sequence is not homologous to other taxa between
bp ~ 2000 and the end of the alignment.

Use the cursor to select position begining of gap around 2000 bp of the
'Balistecaprisc' sequence. Use 'Expand selection Right' in the 'Selection' menu.

Remove this part of the 'Balistecaprisc' sequence using 'Clear selected bases'
in the 'Edit' menu or just press delete

After removing this part of the 'Balistecaprisc' sequence, the two large gaps
near the end are not bridged by any sequence anymore. Remove these gaps entirely
using 'Delete gap-only columns' in the 'Edit' menu.

Have a look at the regions which appears to be poorly aligned. Use the cursor to
click in the ruler area (above the alignment) and select the regions delimited
by boundary sites which appear to be reliably aligned, in contrast to the
alignment block between these boundaries.

In the 'Align' menu, click 'Change default Aligner program > for realigning all
(or selected blocks)'.

Click the third radio button to select 'Mafft-globalpair' as the default
algorithm for realignment. Make sure that the specified MAFFT installation path
is correct, and confirm with 'OK'.

Click 'Realign selected block' in the 'Align' menu.

Does the alignment look more reliable now? Once more, remove gap-only columns,
and save the alignment file.

check the lengths of resulting alignments.

#### Automatic evaluation of alignment
##### using `t_coffee` program 
(this can take several hours to finish)

TCS is an alignment evaluation score that makes it possible to identify the most
correct positions in an MSA. It has been shown that these positions are the most
likely to be structuraly correct and also the most informative when estimating
phylogenetic trees. The TCS evaluation and filtering procedure is implemented in
the T-Coffee package and can be used to evaluate and filter any third party MSA

The TCS is most informative when used to identify low-scoring portions within an
MSA. It is also worth noting that the TCS is not informative when aligning less
than five sequences.

check is program t_coffee is installed, if not run:
```bash
sudo apt-get install t_coffee
```
evaluate alignment:

```bash
# this can take several minutes to finish
t_coffee -infile  16s_aln.fasta -evaluate
```

##### using `bmge` program
BMGE is able to perform biologically relevant trimming on a multiple alignment
of DNA, codon or amino acid sequences. BMGE is designed to select regions in a
multiple sequence alignment that are suited for phylogenetic inference. For each
character, BMGE computes a score closely related to an entropy value.
Calculation of these entropy-like scores is weighted with BLOSUM or PAM
similarity matrices in order to distinguish among biologically expected and
unexpected variability for each aligned character

```
conda create -n bmge -c conda-forge -c biconda bmge
conda activate bmge
```

run `bmgi`
```
bmge  -i 16s_aln.fasta -t DNA -of 16s_filtered.fasta -oh 16s_filtered.html
```
`bmge` generates two output - html is suitable for viewing in web browser.
fitered.fasta can be used for phylogenetic analysis.

### Exercise 3.5 - MSA of globins proteins
 Create and analyze a multiple sequence alignment (MSA) of proteins from the globin family. Globins are oxygen-binding proteins found in many organisms, including humans. They play a crucial role in oxygen transport and storage. Some well-known globin family members are hemoglobin, myoglobin, and neuroglobin.

Protein Sequence Set: Select protein sequences from different organisms representing hemoglobin, myoglobin, and neuroglobin. For example:

- Human Hemoglobin Subunit Alpha (P69905)
- Human Hemoglobin Subunit Beta (P68871)
- Human Myoglobin (P02144)
- Human Neuroglobin (Q9NPG2)
- Mouse Hemoglobin Subunit Alpha (P01942)
- Mouse Hemoglobin Subunit Beta (P02088)
- Mouse Myoglobin (P04247)
- Mouse Neuroglobin (Q9ER97)

P69905; P68871; P02144; Q9NPG2; P01942; P02088; P04247; Q9ER97

1. Which regions of the aligned sequences are conserved across all the globin
   family members? What might be the functional significance of these conserved
   regions?
2. Can you identify any organism-specific or protein-specific sequence
   variations? What might be the evolutionary or functional implications of
   these differences?

Inspect 3D structure of neuroglobin protein  - for example
https://www.rcsb.org/3d-sequence/1OJ6?assemblyId=1

### Exercise 3.6 - Identification of catalytic triad residues in serine proteases
Create and analyze an MSA of serine proteases, a family of enzymes that cleave
peptide bonds in proteins. They play essential roles in digestion, blood
clotting, and immune responses. Some well-known serine proteases include
trypsin, chymotrypsin, and elastase.

- Human Trypsin-1 (P07477)
- Human Chymotrypsinogen B (P17538)
- Human Neutrophil Elastase (P08246)
- Mouse Trypsin-2 (P07146)
- Mouse Chymotrypsinogen B (Q9CR35)
- Mouse Neutrophil Elastase (Q3UP87)
- Drosophila melanogaster Serine proteinase stubble (Q05319)
-  Drosophila melanogaster Chymotrypsin (Q9VVA6)
- Xenopus laevis Complement C3 (Q91701)
- Manduca sexta Chymotrypsinogen (Q25503)

P07477; P17538; P08246; P07146; Q9CR35;Q3UP87;Q05319; P42280;Q91701;Q25503

Sequences can be obtained either 
  - from UniProt database (https://www.uniprot.org/) using above accessions with 'list' option
  - or directly from Jalview program using 'Fetch sequences' command from File menu

Identify the conserved catalytic triad residues in the aligned sequences. The
conserved catalytic triad residues in serine proteases are typically histidine (H),
aspartate (D) , and serine (S). These residues form a charge relay system that enables
the nucleophilic attack by the serine residue on the peptide bond of the
substrate. The conservation of this triad is essential for the
catalytic mechanism of serine proteases.

In serine proteases, the catalytic triad and other essential structural elements
are typically under negative selection, while surface loops and other flexible
regions might be under positive selection, reflecting the adaptation to
different substrates and physiological conditions.

In the alignment, there are two conserved Serines, Two identify which one is part of the catalytic triad, you will need to check the function of Human Trypsin in Uniprot database - https://www.uniprot.org/uniprotkb/P07477

### Exercise 3.7 - Identification of Bacterial Homologs of Human Neuroglobin and Analysis of Heme-Binding Pocket Conservation
*Objective*:  In this assignment, you will identify bacterial homologs of human
neuroglobin, create a multiple sequence alignment (MSA) with human and mouse
globin homologs (hemoglobin, myoglobin, and neuroglobin), and analyze the
conservation of the heme-binding pocket residues, particularly the two histidine
residues, in the bacterial sequences.

*Background*: Globins are a family of heme-containing proteins found in various
organisms, including bacteria, plants, and animals. They are involved in various
functions, such as oxygen transport and storage. The heme-binding pocket in
globins typically contains two conserved histidine residues that coordinate with
the iron atom in the heme group, allowing the protein to reversibly bind oxygen
or other small ligands.

Some of the key functions of bacterial globins include:
- Oxygen transport and storage: Bacterial globins maintain oxygen supply for
  cellular respiration in microaerophilic or facultative anaerobic bacteria.
- Oxygen sensing and regulation: They function as oxygen sensors, helping
  bacteria adapt to changing oxygen levels and modulating gene expression.
- Nitric oxide detoxification: Flavohemoglobins detoxify nitric oxide,
  converting it to a less toxic form.
- Oxidative stress protection: Bacterial globins scavenge reactive oxygen
  species to protect cells from oxidative stress.
- Terminal oxidases: Cytochrome bd-type oxidase globins act as terminal oxidases
  in the respiratory chain, transferring electrons to oxygen.
- Sulfide oxidation: Sulfide:quinone oxidoreductases (SQR) globins oxidize
  hydrogen sulfide for energy in sulfur bacteria.
- Sensing and signaling: Bacterial globins act as sensors and signal
  transducers, detecting environmental changes and triggering cellular
  responses.

#### Tasks:

1. Obtain the amino acid sequence of human neuroglobin from a database, such as
   UniProt or NCBI Protein.
2. Use a sequence similarity search tool  BLASTP, to identify bacterial homologs
   of human neuroglobin.
   - limit search to 50 sequences and e-value threshold to 1e-5
3. Retrieve FASTA sequence for the 10 best bacterial hits based on E-value,
   sequence identity, and query coverage.
4. Import FASTA with bacterial globins to Jalview program. 
5. Collect the amino acid sequences of human and mouse globin homologs
   (hemoglobin, myoglobin, and neuroglobin) from UniProt directly from Jalview
   program. Use following accessions: P69905; P68871; P02144; Q9NPG2; P01942;
   P02088; P04247; Q9ER97 (use Fetch sequences command from FIle menu)
6. Create a multiple sequence alignment (MSA) using the Jalview program with
   bacterial, human and mouse sequences.
7. Analyze the MSA to assess the conservation of the two histidine residues in
   the heme-binding pocket among the bacterial globin sequences.
8. Clearly describe your observations and conclusions regarding the conservation
   of the two histidine residues in the heme-binding pocket.

#### Example of bacterial globin with known 3D structure:
- https://www.rcsb.org/3d-sequence/1TU9
#### Example of neuroglobin
https://www.rcsb.org/3d-sequence/1OJ6?assemblyId=1

In 3d viewer, highlight *metal coordination hem* for both structures, in bacterial globin highlight binding site PPI (Q)

<details>
<summary>💡 Hint</summary>

The proximal His is conserved to ligate the heme iron.
The distal position is what varies—bacteria often use Gln (plus Tyr/Trp neighbors) to re-shape the H-bonding landscape,
slow O₂ off-rates, limit autoxidation, and favor NO-dioxygenase/sensing functions instead of the vertebrate-style His
gate used for reversible gas transport

</details>

### Exercise 3.8 - Alignment-Based Primer Design

#### Motivation
Cuscuta([[https://en.wikipedia.org/wiki/Cuscuta]]), commonly known as dodder, is a parasitic plant that affects a wide range of host species. Proper identification of different Cuscuta species is important for understanding their ecological impacts and for managing affected plants. The 5.8 rDNA gene is part of the ribosomal DNA (rDNA) gene cluster, which includes the 18S, 5.8S, and 28S rDNA genes separated by internal transcribed spacers (ITS1 and ITS2). The 5.8 rDNA gene is commonly used for this type of task because it contains both conserved and variable regions, making it suitable for distinguishing between species while still being conserved enough for primer design. Your task is to design primers that could be used for species identification in further experiments. These primers will help amplify specific conserved regions, aiding in the confirmation of species identity.

#### Organization of rDNA
 The ITS1 and ITS2 regions are highly variable, which makes them useful for distinguishing between closely related species, whereas the 18S, 5.8S, and 28S regions are more conserved, providing stable targets for primer design.
![rDNA](./rDNA.png)

The set of rDNA sequences from Cuscuta species can be found in file `../data/alignment_sequences/5.8S_Cuscuta.fasta`. Your task is to identify conserved regions and design primers accordingly. The exercise will guide you through the following steps:

#### Sequence Alignment with Jalview

Open the provided set of sequences using Jalview.

Perform a multiple sequence alignment to visualize conserved and variable regions among the different Cuscuta species. Use `mafft` program with L-INS-i setting.

Inspect the alignment to identify suitable regions for primers, and save the alignment for use in the next step.

#### Generating a Consensus Sequence with EMBOSS Cons

Use the EMBOSS Cons program available on the web to create a consensus sequence from your alignment. [Link to EMBOSS Cons](http://www.ebi.ac.uk/jdispatcher/msa/emboss_cons?stype=dna&matrix=EDNAFULL)

The consensus sequence should highlight the regions that are conserved across all species, which will be key for primer design. Set a suitable identity threshold ( you have to unfold parameters setting). Here the identity parameter mean "The required number of identities at a site for it to give a consensus at that position".

#### Primer Design with Primer3

Use the consensus sequence as input for Primer3 (web version: [Primer3](https://primer3.ut.ee/)) to design primers.

Select a suitable product size range, aiming for as large a product as possible, and click on "Pick Primers."

Do the suggested primers correspond to the conserved regions observed in the alignment from step 1?

#### Deliverable
Once your primers are designed, they could potentially serve to confirm species identity in further experiments. 

#### Tools Needed

Jalview: For sequence alignment.

EMBOSS Cons: To create a consensus sequence from the alignment.

Primer3: To design primers based on the consensus sequence.

This exercise will give you hands-on experience in using bioinformatics tools to identify conserved regions and design primers, skills that are crucial for species identification and molecular biology research.

## Amino Acid codes

![Amino Acid Codes](../fig/aa_codes.png)

## Jalview

- Jalview has *two navigation and editing modes*: _normal mode_, where editing and navigation is
  performed using the mouse, and _cursor mode_ where editing and navigation are performed using
  the keyboard. The *F2* key is used to switch between these two modes.
- Navigation in cursor mode:
  - Jump to Sequence n: Type a number n then press [S] to move to sequence (row) n.
  - Jump to Column n: Type a number n then press [C] to move to column n in the alignment.
  - Jump to Residue n: Type a number n then press [P] to move to residue number
    n in the current sequence.
- Overview of the whole alignment, especially when it is large. Select *View*
⇒ *Overview Window*
- Find sequence - *ctrl-F*

- *Clustal X* color scheme
https://www.jalview.org/help/html/colourSchemes/clustal.html 
- *Blosum62* Gaps are coloured white. If a residue matches the consensus sequence residue at that
position it is coloured dark blue. If it does not match the consensus residue
but the 2 residues have a positive Blosum62 score, it is coloured light blue.