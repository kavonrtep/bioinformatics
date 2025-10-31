# BLAST

## Exercise 1.1- BLASTN - NCBI web interface
Make nucleotide blast  search against following sequence:
```
>gi|202471|gb|M57671.1|OCOINS Octodon degus insulin mRNA, complete cds
GCATTCTGAGGCATTCTCTAACAGGTTCTCGACCCTCCGCCATGGCCCCGTGGATGCATCTCCTCACCGT
GCTGGCCCTGCTGGCCCTCTGGGGACCCAACTCTGTTCAGGCCTATTCCAGCCAGCACCTGTGCGGCTCC
AACCTAGTGGAGGCACTGTACATGACATGTGGACGGAGTGGCTTCTATAGACCCCACGACCGCCGAGAGC
TGGAGGACCTCCAGGTGGAGCAGGCAGAACTGGGTCTGGAGGCAGGCGGCCTGCAGCCTTCGGCCCTGGA
GATGATTCTGCAGAAGCGCGGCATTGTGGATCAGTGCTGTAATAACATTTGCACATTTAACCAGCTGCAG
AACTACTGCAATGTCCCTTAGACACCTGCCTTGGGCCTGGCCTGCTGCTCTGCCCTGGCAACCAATAAAC
CCCTTGAATGAG
```

### search I
setting of the search:
- In the section "Program Selection" select the option "Somewhat similar
  sequences (blastn)"
- Choose "Nucleotide collection (nr/nt)" as the search database. NR is the "Non
  Redundant" database, which contains all non-redundant (non-identical)
  sequences from GenBank and the full genome databases.
- Check "Show results in a new window" and  Click the BLAST button to launch the search.

After the search has completed, make yourself familiar with the BLAST output page. After a header with some information about the search, there are four main parts:

#### Graphic Summary
each hit is represented by a line showing which part of the query sequence the alignment covers. The lines are coloured according to alignment score.

#### Descriptions
a table with a one-line description of each hit with some alignment statistics, The columns in the Descriptions table are:
- Description — the description line from the database
- Max score — the alignment score of the best match (local alignment) between the query and the database hit
- Total score — the sum of alignment scores for all matches (alignments) between the query and the database hit (if there is only one match per hit, these two scores are identical)
- Query cover — the percentage of the query sequence that is covered by the alignment(s)
- E value — the Expect value calculated from the Max score (i.e. the number of unrelated hits with that score or better you would expect to find for random reasons)
- Ident — the percent identity in the alignment(s)
- Accession — the accession number of the database hit.

#### Alignments
the actual alignments between the query and the database hits.

#### Taxonomy
Taxonomic summary of all positive hits

First, take a look at the best hit. Since our search sequence (the query) was taken from GenBank which is part of NR, we should find an identical sequence in the search. Make sure this is the case!
Find the best hit from human (Homo sapiens) that is not a synthetic construct
and answer following questions:
- what is the identifier (Accession)?
- what is the alignment score ("max score")?
- what is the percent identity and query coverage?
- what is the E-value?
- are there any gaps in the alignment?

### Search II
Make a new BLASTN search with the same query sequence, this time with select database
 to Human genomic + transcript (Human G+T). For the best hits answer following:
- what is the identifier (Accession)?
- what is the alignment score ("max score")?
- what is the percent identity and query coverage?
- what is the E-value?
- are there any gaps in the alignment?

Compare hit of this search with hit from previous search.  You may have noticed
that the E-value changed, while the alignment score did not. Why?

**Hint**:  Expand the "Search summary" section near the top by clicking it. Compare
E values and  the sizes (in basepairs) of the databases we used for the two BLAST searches

## Exercise 1.2 - BLASTP - NCBI web interface
The following is the product of a gene in C. elegans that seems to be important in development.
```
>C.elegans protein
MFHPGMTSQPSTSNQMYYDPLYGAEQIVQCNPMDYHQANILCGMQYFNNSHNRYPLLPQMPPQFTNDHPY
DFPNVPTISTLDEASSFNGFLIPSQPSSYNNNNISCVFTPTPCTSSQASSQPPPTPTVNPTPIPPNAGAV
LTTAMDSCQQISHVLQCYQQGGEDSDFVRKAIESLVKKLKDKRIELDALITAVTSNGKQPTGCVTIQRSL
DGRLQVAGRKGVPHVVYARIWRWPKVSKNELVKLVQCQTSSDHPDNICINPYHYERVVSNRITSADQSLH
VENSPMKSEYLGDAGVIDSCSDWPNTPPDNNFNGGFAPDQPQLVTPIISDIPIDLNQIYVPTPPQLLDNW
CSIIYYELDTPIGETFKVSARDHGKVIVDGGMDPHGENEGRLCLGALSNVHRTEASEKARIHIGRGVELT
AHADGNISITSNCKIFVRSGYLDYTHGSEYSSKAHRFTPNESSFTVFDIRWAYMQMLRRSRSSNEAVRAQ
AAAVAGYAPMSVMPAIMPDSGVDRMRRDFCTIAISFVKAWGDVYQRKTIKETPCWIEVTLHRPLQILDQL
LKNSSQFGS
```
Use the protein-protein blast page to perform search against the non redundant
protein database (nr) to identify protein and find its homologs.

blast link : https://blast.ncbi.nlm.nih.gov/Blast.cgi

- Based on the identical hit to C. elegans what is the identity of this protein?
- Aside from the C. elegans proteins, what is the most significant hit? What is
  its identity and E-value?
- switch for "graphic summary" - for blastp there is additional information
  about identified conserved domains.
- Explore link to other databases  - find coresponding entry in uniprot database.
- Check also taxonomy report. How many hits was found total. How many hits was
  found in "chordata" taxonomy group? If no hits are found in Chordata, try to modify search
  parameters to identify similar sequence in Chordata group.

## Exercies 1.3 - BLASTX - NCBI web interface
The part of genomic sequence of Twort bacteriophage is bellow. Use BLAST to
search protein database  (Non-redundant protein sequences  - nr) to find
homologous protein from other bacteriophages and viruses.
- For search try to use Blosum62 and Blosum45 scoring matrix. What is the
  difference between searches with different matrices
- What proteins (Accessions) has alignment with E value < 1e-10

```
>Staphylococcus phage Twort, partial sequence
AAAGATGCAGAGTTAGCTGTAATGGAAATCAACAAAAAACAATTGGAGGACTAATCT
TAATGAGTAAATTTTCAAATATTCTAGAAGAATATAATAAATTACAATCACAAGATG
TTGAAAAATCATTAGAAGAAAATAAAGATGAAGAACCTAAAGAAGAGGCTACTGTAG
AATCAGTTACGGAAGAACAAGTAGTTGAAACAGATGCACCGCAAAAAGAAGAACCAC
AACAAGTATCTGAAGAAGACGCTAAAAAAGCACAAGAAGAATCTAAGAAATTAGAAT
CAGAAAAGCAAGAAGAAGATAAAGAAGTAGAGAAGTCTGTTAAAGATTCTAAAGACC
CAGTAGACCATAAAGATACTAAAACTGAAGACAAAGACAATGAGAAACGTAAAAACA
AAAAAGAAGATAAAGAAGACGAGTCTAAAGAAGAAGATGAAAAAGAATCTAAGAAAG
ATAAAGACAAAGAAGATAAAAAGTCTGAAAACAAAAAAGATTCTGAAAAAGTTAAAA
AGTCAGCTTTATCTGATGAAGATATTGTAGAAGGATTTAGTACAGTATTAAAATCTT
TACAAGACTTACCTAAACAATTTGCTACTAAAGATGATGTCAAGGAAATTAAAAAAT
CTTTAGAAGAATTACAAGATGCTTTTGCTGAAAAAGAAAAGAAACAAGAAGAAAAAG
TAGAAACTATTAAAGAAGAAGTTAACAAAGAACAAGAAGATAAAGAAGAAGAAAATA
CTGATGAATCAGTAGAAAAATCAGTAACAACTTCTAACACTGCACAGCAAGATGATG
TTAATTATGTTTCTAAATCAGCAGTTGTAGAAGAAGAAGTACAAGAGGAACAACCTG
AGGAAGATAAGCAAGAGGTTAATACAATTACACAAGAAGACCGTGAGGCTTTCATGA
ATAAATTCAAATCAGAATCTCAACGTCGTGACAAACCTACACGTCAATTAAATGATG
CATATTTAGCATATATGGACGTTCGTAATAATGGTGAAAATGCAAGTCCAAGTTCTT
TAAAAACTGTTAAAGATTTTATTAAGTAATACAAAGTAGTTGTGTTATATTATACAT
GAAATTAAATTAATAAAA
```


## Exercise 1.4 - Identification of coding sequence using BLASTX vs BLASTN
characterization of unknown DNA PCR amplified fragment from unknown
non-cultivatable microorganism
```
>clone12
AACGGGCACGGGACGCATGTAGCTGGAACAGTGGCAGCCGTAAATAATAATGGTATCGGA
GTTGCCGGGGTTGCAGGAGGAAACGGCTCTACCAATAGTGGAGCAAGGTTAATGTCCACA
CAAATTTTTAATAGTGATGGGGATTATACAAATAGCGAAACTCTTGTGTACAGAGCCATT
GTTTATGGTGCAGATAACGGAGCTGTGATCTCGCAAAATAGCTGGGGTAGTCAGTCTCTG
ACTATTAAGGAGTTGCAGAAAGCTGCGATCGACTATTTCATTGATTATGCAGGAATGGAC
GAAACAGGAGAAATACAGACAGGCCCTATGAGGGGAGGTATATTTATAGCTGCCGCCGGA
AACGATAACGTTTCCACTCCAAATATGCCTTCAGCTTATGAACGGGTTTTAGCTGTGGCC
TCAATGGGACCAGATTTTACTAAGGCAAGCTATAGCACTTTTGGAACATGGACTGATATT
ACTGCTCCTGGCGGAGATATTGACAAATTTGATTTGTCAGAATACGGAGTTCTCAGCACT
TATGCCGATAATTATTATGCTTATGGAGAGGGAACATCCATGGCTTGTCCACATGTCGCC
GGCGCCGCC
```


- use `blastn` and `blastx` to characterize clone12. What tool is relevant to use?
- check DNA sequence using ORFfinder https://www.ncbi.nlm.nih.gov/orffinder
- use `blastp` against longest ORF, run blastb directly from orffinder, compare results with blastx search.
- What kind of enzyme is coded by `clone12`?
- does the clone12 represent complete CDS of putative protein?

<details>
<summary>Details</summary>
In blastx AA alignment is truncated probably because of low complexity
filtering, try to turn off filtering and run in again.
When using BLASTP from orf finder - select non-redundant protein database (nr).
</details>

- compare best hit from blastx directly using tblast against DNA sequence of
  clone12
  - click on accession - this will take you to protein entry.
  - on protein entry page, click on left panel on *run BLAST*
  - select * tblastn* - query is a protein accession (this represent query)
  - select `Align two or more sequences`
  - second sequence is clone12 nucleotide sequence
  - run blast and explore results
  - inspect also dotplot page

## Exercise 1.5 - Using BLAST (blast2seq) to create local aligment for two sequences:
Two sequences can be also compared using BLAST web interface. For comparison use:
(two align two sequences to each other - click on checkbox *Align two or more sequences*  in NCBI BLAST form)
- `ERB2_HUMAN` : http://www.uniprot.org/uniprot/P04626.fasta
- `EGFR_DROME` : http://www.uniprot.org/uniprot/P04412.fasta
paste these sequences to BLAST form. Use blastp.

https://blast.ncbi.nlm.nih.gov/Blast.cgi?BLAST_SPEC=blast2seq&LINK_LOC=align2seq&PAGE_TYPE=BlastSearch
blast2seq can be used instead of `needle`. It also provide graphical view of alignment and non-interactive dotplot. Use blast2 seq on  `P04626.fasta` and  `P04412.fasta` sequences and explore results. Compare alignments and dotplot.

Alternativelly, it is possible to enter can use just accession ID to search window (P04626, P04412).

## Exercise 1.6 - Identification of species using NCBI BLAST
In your experiment, you are working with plant species Cuscuta campestris,
Cuscuta californica, and Cuscuta japonica. These plants are difficult to
distinguish based solely on morphology. To confirm the correct species
identification, you performed PCR amplification of the Internal Transcribed
Spacer (ITS) genes for rRNA, which are regions between rRNA genes that are often
used for species identification. You obtained the following sequences:
```
>speciesA
GGAGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTCGAAACCTG
CCCAGCAGAACGACTCGAGAACCTGTTTCACATACAACACATTCATTGGGGGCTGTTCTC
TCGGGCACGCGCCTCCAATGATCAACGAACCCCCGGCGCGGAACGCGCCAAGGACTACTC
AAACGAGATCGTCGGGCCATCGTGCCCCGTCCGCGGGTGCATGGGTGGCGTTGGCGTCTT
TAATAACATAAACGACTCTCGGCAACGGATATCTCGGCTCTCGCATCGATGAAGAACGTA
GCGAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGC
AAGTTGCGCCCAAAGCCGTCAGGCCGAGGGCACGTCTGCCTGGGCGTCACGCATCGCGTC
GCTCCCCTCCCGTTGCGGAGCGGGGAGCGGATGATGGCCTCCCGTGCCCGACCTTGGATG
CGGCTGGCTGAAATGTTGGTCCTTGACGACTGACGTCACGGCGAGTGGTGGTCGTACCTA
GTGTGCTTATCGTCGCGTCGTGCCCAGTCATCTTGGGATTTTGACCCTTTTGAGCTGGTG
TGAGCTGGCTCTCTGACCGCGACCCCAGGTCAGGCGGGACTACCCGCTGAGTTTAAGCAT
ATCAATAAGCGGAGGAAAAGAAACT
>speciesB
GGAGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTCGAAACCTG
CCCAGCAGAACGACTCGAGAACATGTTTCACATACAACACATTCATTGGGGGTTGTGCTC
TCGGGCACACGCCCCCAATGATCAACGAACCCCCGGCGCGGAACGCGCCAAGGATTACTC
AAATGAGATCGTCGGGCCATCGTGCCCCGTCCGCGGGTGCATTGGTGGCATTGGCGTCTT
TAATAACATAAACGACTCTCGGCAACGGATATCTCGGCTCTCGCATCGATGAAGAACGTA
GCGAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTCTTTGAACGC
AAGTTGCGCCCAAAGCCGTCAGGCCGAGGGCACGTCTGCCTGGGCGTCACGCATCGCGTC
GCCCCCCTCCCATTGCGGAGCGGGGAGCGGATGATGGCCTCCCGTGCCCGACCTTGGATG
CGGTTGGCCGAAATGTTGGTCCTTGACGATAGACGTCACGGCGAGTGGTGGTCGTACCTA
GTGTGCTTATCGTCGCGTCGTGCCCTATCGTCTTGCGATTTTGACCCTTTTGAGTTGGTG
TGAGCCGGCTCTCTGACCGCGACCCCAGGTCAGGCGGGACTACCCGCTGAGTTTAAGCAT
ATCAATAAGCGGAGGAAAAGAAACT
>speciesC
GGAGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTCGAAGCCTC
GCTGAGAAAGACTTGTTAACCTGTACCAATTCATGATTCGAGTGTCGTGGTCATTTTCTG
ATTTGCCCATGACGAACACAAAAACACCGGCGCGGCAGCGCCAAGGAATTTCGTGATGAG
TATGCTGCCTCATATAGCTCGTACTCTACTGCTTGTGAGGTTGGCTTCCTTTAAGAAAAA
TGACTCTCGGCAATGGATATCTCGGCTCTTGCAACGATGAAGAACGTAGCGAAATGCGAT
ACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAAACTTTGAACGCAAGTTGCGCCTC
ATGCCATTAGGTTGAGGGCACGTTTGCTTGGGTGTCATGCGTTATGTCTTCCCTCTCGTG
CGTGGAGTGGGAATAGATTGTGGCCTCCTGGGCCCTTCCTTGGGCGTGGTTGGCCGAAAA
AGTTGTCCTTGACTCTGTCGATGCCTTGGTGTGTGGTGGACGTACCAAGTGTGCATGATT
GCCAGCCTTGCTCGGCTTCATTGTGGCGTTCGGATCCTATGAGGCTGTCGGTTTTGGCTC
TTTGATTGCGGCCCCAAGTCAGGCGAGACCACCCGCTGAGTTTAAGCATATCAATAAGCG
GAGGAGAAGAAACT
```
You will use the NCBI BLAST (Basic Local Alignment Search Tool) web service to
compare these sequences with a database of known sequences. This step will allow
you to identify the closest matches between your sequences and those stored in
the database, which is crucial for species identification. For each sequence,
check the list of top hits and identify the source species. Evaluate the
percentage identity and sequence coverage, which are important for assessing the
quality of the match. Records with high identity and good sequence coverage are
likely to correspond to the correct species matching your sequences.

Questions:

- What are the BLAST results for each sequence, including percentage identity and
   coverage? To which species of Cuscuta were these sequences assigned?
- How would you evaluate the reliability of the results? Did you encounter any
ambiguities or problems during the analysis? What were these problems?

## Exercise 1.7 - Identification of mutations in gyrA gene
Neisseria meningitidis is a Gram-negative encapsulated bacterium isolated only
from humans, where it can cause serious invasive infections (primarily
septicemia and meningitis). Treatment of meningococcal disease requires
treatment of patients and chemoprophylaxis for contacts. Currently, antibiotics
such as ciprofloxacin are recommended for chemoprophylaxis. Resistance of
meningococci to ciprofloxacin is associated with mutations in the quinolone
resistance-determining region (QRDR) of the gyrA gene.

Your task will be to determine how many mutations are present in the gyrA gene
of the sequenced N. meningitidis strains compared to the reference genome.

Reference genome: NZ_CP021520.1

The sequenced strains code for the following gyrA proteins:
```
>strain_1
MTDATIRNDHKFALETLPVSLEDEMRKSYLDYAMSVIVGRALPDVRDGLKPVHRRVLYAMHELKNNWNAAYK
KSARIVGDVIGKYHPHGDIAVYDTIVRMAQDFAMRYVLVDGQGNFGSIDGLAAAAMRYTEIRMAKISHEMLA
DIEEETVNFGPNYDGSEHEPLVLPTRFPTLLVNGSSGIAVGMATNIPPHNLSDTINACLRLLDAPDTEIDEL
IDIIQAPDFPTGATIYGLSGVREGYKTGRGRVVIRAKTHTEPIGKNGEREAIVIDEIPYQVNKAKLVEKIGE
LVREKTLEGISELRDESDKSGMRVVIELKRNENAEVVLNQLYKLTPLQDSFGINMVVLVDGQPRLLNLKQIL
SEFLRHRREVVTRRTLFRLKKARHEGHIAEGKAVALSNIDEIIRLIKESPNAAEAKEKLLARPWRSSLVEEM
LTRSGLDLEMMRPEGLTANIGLKEQGYYLSEIQADAILRMSLRKLTGLDQEEIVESYKNLMGKIIDFVDILS
KPERITQIIRDELEEIKTNYGDERRSEINPFGGDIADEDLIPQREMVVTLTHGGYIKTQPTTDYQAQRRGGR
GKQAAATKDEDFIETLFVANTHDYLMCFTNLGKCHWIKVYKLPEGGRNSRGRPINNVIQLEEGEKVSAILAV
REFPEDQYVFFATAQGIVKKVQLSAFKNVRSQGIKAIALKEGDYLVGAAQTGGSDDIMLFSNLGKAIRFNEY
WEKSGNDEAEDADIETEISDDLEDETADNENALPSGKHGVRPSGRGSGGLRGMRLPADGKIVSLITFAPEAA
QSDLQVLTATANGYGKRTPIADYSRKNKGGQGNIAINTGERNGDLVAATLVSETDDLMLITSGGVLIRTKVE
QIRETGRAAAGVKLINLDEGETLVSLERVAEDESELSDASVISNVTEPEVEN
>strain_2
MTDATIRHDHKFALETLPVSLEDEMRKSYLDYAMSVIVGRALPDVRDGLKPVHRRVLYAMHELKNNWNAAYK
KSARIVGDVIGKYHPHGDTAVYDTIVRMAQNFAMRYVLIDGQGNFGSVDGLAAAAMRYTEIRMAKISHEMLA
DIEEETVNFGPNYDGSEHEPLVLPTRFPTLLVNGSSGIAVGMATNIPPHNLSDTVNACLRLLDAPDTEIDEL
IDIIQAPDFPTGATIYGLSGVREGYKTGRGRVVMRGKTHIEPIGRNGEREAIVIDEIPYQVNKAKLVEKIGD
LVREKTLEGISELRDESDKSGMRVVIELKRNENAEVVLNQLYKLTPLQDSFGINMVVLVDGQPRLLNLKQIL
SEFLRHRREVVTRRTLFRLKKARHEGHIAEGKAVALSNIDEIIKLIKESPNAAEAKDKLLAHPWRSSLVEEM
LTRSGLDLEMMRPEGLAANIGLKEQGYYLSEIQADAILRMSLRNLTGLDQEEIVESYKNLMGKIIDFVDILS
KPERITQIIRDELEEIKTNYGDERRSEINPFGGDIADEDLIPQREMVVTLTHGGYIKTQPTTDYQAQRRGGR
GKQAAATKDEDFIETLFVANTHDYLMCFTNLGKCHWIKVYKLPEGGRNSRGRPINNVIQLEEGEKVSAILAV
REFPEDQYVFFATAQGLVKKVQLSAFKNVRAQGIKAIALKEGDYLVGAAQTGGADDIMLFSNLGKAIRFNEY
WEKSGNDEAEDADIETEILDGIEDETADSENALPSGKHGVRPSGRGSGGLRGMRLPADGKIVSLITFAPETE
ESGLQVLTATANGYGKRTPIADYSRKNKGGQGNIAINTGERNGDLVAATLVGETDDLMLITSGGVLIRTKVE
QIRETGRAAAGVKLINLDEGETLVSLERVAEDESELSDASVISNVTEPEVEN
>strain_3
MTDPTIRHDHKFALETLPVSLEDEMRKSYLDYAMSVIVGRALPDVRDGLKPVHRRVLYAMHELKNNWNAAYK
KSARIVGDVIGKYHPHGDTAVYDTIVRMAQNFAMRYVLIDGQGNFGSVDGLAAAAMRYTEIRMAKISHEMLA
DIEEETVNFGPNYDGSEHEPLVLPTRFPTLLVNGSSGIAVGMATNIPPHNLSDTVNACLRLLDAPDTEIDEL
IDIIQAPDFPTGATIYGLSGVREGYKTGRGRVIMRGKTHIEPIGKNGEREAIVIDEIPYQVNKAKLVEKIGD
LVREKTLEGISELRDESDKSGMRVVIELKRNENAEVVLNQLYKLTPLQDSFGINMVVLVDGQPRLLNLKQIL
SEFLRHRREVVTRRTLFRLKKARHEGHIAEGKAVALSNIDEIIKLIKESPNAAEAKDKLLAHPWRSSLVEEM
LTRSGLDLEMMRPEGLAANIGLKEQGYYLSEIQADAILRMSLRNLTGLDQEEIVESYKNLMGKIIDFVDILS
KPERITQIIRDELEEIKTNYGDERRSEINPFGGDIADEDLIPQREMVVTLTHGGYIKTQPTTDYQAQRRGGR
GKQAAATKDEDFIETLFVANTHDYLMCFTNLGKCHWIKVYKLPEGGRNSRGRPINNVIQLEEGEKVSAILAV
REFPEDQYVFFATAQGLVKKVQLSAFKNVRAQGIKAIALKEGDYLVGAAQTGGADDIMLFSNLGKAIRFNEY
WEKSGNDEAEDADIETEISDGIEDETADSENALPSGKHGVRPSGRGSGGLRGMRLPADGKIVSLITFAPETE
ESGLQVLTATANGYGKRTPIADYSRKNKGGQGNIAINTGERNGDLVAATLVGETDDLMLITSGGVLIRTKVE
QIRETGRAAAGVKLINLDEGETLVSLERVAEDESELSDASVISNVTEPEAEN
>strain_4
MTDATIRHDHKFALETLPVSLEDEMRKSYLDYAMSVIVGRALPDVRDGLKPVHRRVLYAMHELKNNWNAAYK
KSARIVGDVIGKYHPHGDIAVYDTIVRMAQDFAMRYVLVDGQGNFGSVDGLAAAAMRYTEIRMAKISHEMLA
DIEEETVNFGPNYDGSEHEPLVLPTRFPTLLVNGSSGIAVGMATNIPPHNLSDTVNACLRLLDAPDTEIDEL
IDIIQAPDFPTGATIYGLSGVREGYKTGRGRVVMRGKTHIEPIGRNGEREAIVIDEIPYQVNKAKLVEKIGD
LVREKTLEGISELRDESDKSGMRVVIELKRNENAEVVLNQLYKLTPLQDSFGINMVVLVDGQPRLLNLKQIL
SEFLRHRREVVTRRTLFRLKKARHEGHIAEGKAVALSNIDEIIRLIKESPNAVEAKDKLLARPWRSSLVEEM
LTRSGLDLEMMRPEGLAANIGLKEQGYYLSEIQADALRMSLRNLTGLDREEIVESYKNLMGKIIDFVDILS
KPERITRIIRDELEEIKTNYGDERRSEINPFGGDIADEDLIPQREMVVTLTHGGYIKTQPTTDYQAQRRGGR
GKQAAATKDEDFIETLFVANTHDYLMCFTNLGKCHWIKVYKLPEGGRNSRGRPINNVIQLEEGEKVSAILAV
REFPEDQYVFFATAQGIVKKVQLSAFKNVRSQGIKAIALKEGDYLVGAAQTGGSDDIMLFSNLGKAIRFNEY
WEKSGNDEAEDADIETEISDDLEDETADNENALPSGKHGVRPSGRGSGGLRGMRLPADGKIVSLITFAPEAA
QSDLQVLTATANGYGKRTPIADYSRKNKGGQGNIAINTGERNGDLVAATLVSETDDLMLITSGGVLIRTKVE
QIRETGRAAAGVKLINLDEGETLVSLERVAEDESELSDASVISNVTEPEVEN
```
Use the TBLASTN tool on NCBI to compare each of the protein sequences with the
reference genome NZ_CP021520.1. Set the mode to "Align two sequences" for a
direct comparison of your sequence (query) with the reference sequence
(subject). Enter the reference genome ID in the subject field and run the
search for each strain's protein sequence.

Inspect different alignments views (pairwise, query-anchored, ...).

Questions:

- How many mutations did you find in each strain's protein compared to the
reference genome?
- Did any of the strains contain an insertion or deletion?
- Can you use BLASTX instead of TBLASTN for this analysis? Why or why not?


# Use of BLAST from command line

When to use blast from CLI:
- If you have large number of queries (it is possible do download whore `nr` database and run blast localy)
- When your sequence database is not part of public databases (NCBI, EBI,...)
- If you need to automate your similarity search
- More detailed manual can be found at https://www.ncbi.nlm.nih.gov/books/NBK279684/


## Basic commands:

`makeblastdb`, `blastn`, `blastp`, `blastx`, `tblastx`

<details>
<summary>Details</summary>
explan differences in commands
</details>

Basic use of blast commands
```bash
  makeblastdb -help
  blastn -help
  # create database:
  makeblastdb -dbtype nucl -in  dna_sequences.fasta
  # or for proteins:
  makeblastdb -dbtype prot -in  prot_sequences.fasta
  # nucleotide - nucleotide search
  blastn -db database_file -query query_sequences.fasta -out output_file
```


The most used blast options:
```
-db <String>
   BLAST database name

-out <File_Out>
   Output file name
   Default = `-'

-evalue <Real>
   Expectation value (E) threshold for saving hits
   Default = `10'

-word_size <Integer, >=4>
   Word size for wordfinder algorithm (length of best perfect match)

-outfmt <String>
   output format
```
for complete options type `blastn -help`



## Exercise 2.1

files in exercise:
- query : `~/Desktop/bioinformatics/data/blast_data/proteins.fasta`
- database : `~/Desktop/bioinformatics/data/blast_data/db/pdbaa`

<details>
<summary>Details</summary>
input sequenc contain two proteins - sequence1 and sequence2
sequence1: is cytochome c oxidase subunit )
sequence2: - HIV1 envelope protein
db is fraction of PDB protein database
</details>

run protein blast with default parameters in terminal
```bash
cd
mkdir blast_search
cd blast_search
# copy query and database to directory with data:
cp ~/Desktop/Bioinformatics/data/blast_data/proteins.fasta .
cp ~/Desktop/Bioinformatics/data/blast_data/db/pdbaa .

# inspect the query file  protein.fasta
cat proteins.fasta
seqkit stats proteins.fasta
# inspect fasta file we will use as database
seqkit stats pdbaa

# fasta file db/pdbaa will be used as database, it must be formated using
# makeblastdb command to make data blast compatible
makeblastdb -in pdbaa -dbtype prot
# after succesfull creation of database, information about size of database is printet to stdout
# Additional files in db directory were created, what are these files?
ls -l
# run blastp with default settings:
blastp -query proteins.fasta -db pdbaa -out proteins_blastp_default.txt
# inspect output with less command or text editor
less proteins_blastp_default.txt
```


try command line blast with different parameters:
```bash
# see all passoble blast options:
blastp -h
# or
blast -help
# blast documentation is long, it can be more convenient pipe it to less
blast -help | less
```

The most important blastp/blastn options:
- `-task` type of algorithm for search (blast/megablast for blastn or blastp/blastp-fast for blastp)
- `-outfmt` output format
- `-word_size`
- `-evalue`
- `-num_alignments` Number of database sequences to show alignments for (Default = 250)


```bash
# default:
  blastp -query proteins.fasta -db pdbaa -out proteins_blastp_align_all.txt
# limit output to 10 alignments
  blastp -query proteins.fasta -db pdbaa -out proteins_blastp_align_top10.txt -num_alignments 10
# return max 10 alingments, hits must be e-value 1e-30 or less
  blastp -query proteins.fasta -db pdbaa -out proteins_blastp_align_1e-30.txt -num_alignments 10 -evalue 1e-30
# tabular ouput with descriptions
  blastp -query proteins.fasta -db pdbaa -out proteins_blastp_1e-30_table.csv -evalue 1e-30 -outfmt 7
# plain tabular output
  blastp -query proteins.fasta -db pdbaa -out proteins_blastp_1e-30_table2.csv -evalue 1e-30 -outfmt 6
# simple html output
  blastp -query proteins.fasta -db pdbaa -out proteins_blastp_1e-30.html -evalue 1e-30 -outfmt 2 -html
# inspect all output using less command
# html output should be viewed in firefox!
```

<details>
<summary>Details</summary>
Tabular output can be easily viewed using spreadsheet software (LibreOffice Calc, MS Excel).
</details>


## Exercise 2.2 - extract hits from database and create alignment with query:
The aim of this exercise is to perform a BLASTP search against a protein
database using a single query sequence and subsequently analyze the
results, extracting hit information, retrieving and ordering the hit sequences,
and finally, creating a multiple sequence alignment.

```bash
# Run BLASTP search
cd
mkdir blast_search2
cd blast_search2
# copy query and database to directory with data:
cp ~/Desktop/Bioinformatics/data/blast_data/query2.fasta .
cp ~/Desktop/Bioinformatics/data/blast_data/db/pdbaa .
makeblastdb -in pdbaa -dbtype prot


blastp -query query2.fasta -db pdbaa -out query_blastp_1e-10_table.txt -evalue 1e-10 -outfmt 6
# -query: Specifies the input query file (in FASTA format).
# -db: Specifies the protein database to search against.
# -out: Specifies the output file to save results.
# -evalue 1e-10: Sets the e-value threshold for reporting hits (only hits with an e-value less than 1e-10 will be reported).
# -outfmt 6: Specifies the output format (6 is a tabular format).


# Inspect the structure of the output from BLAST using the 'less' command
less query_blastp_1e-10_table.txt
# The columns in the tabular output are:
# 1.qaccver 2.saccver 3.pident 4.length 5.mismatch 6.gapopen
# 7.qstart 8.qend 9.sstart 10.send 11.evalue 12.bitscore

# Count the number of hits for the query sequence
wc query_blastp_1e-10_table.txt
# wc: Counts the number of lines, words, and characters in a file (use -l to only count lines).


# Extract a list of sequence IDs:
# Extract the second column from the BLAST output, which contains the database IDs
cut -f 2 query_blastp_1e-10_table.txt > all_hits_id.txt
# -f 2: Extracts the second field (subject IDs).

seqkit grep -f all_hits_id.txt pdbaa -o all_hits.fasta
# -f all_hits_id.txt: Specifies the file with the list of sequence IDs.
# -o all_hits.fasta: Specifies the output file for the extracted sequences.

seqkit stats  all_hits.fasta


# Explore the extracted sequences with Dotter
dotter query2.fasta all_hits.fasta

# If cdbfasta is not installed, use the following command to install it:
sudo apt install cdbfasta


# Extract sequences in the order of hit significance
cdbfasta pdbaa
# cdbfasta creates an index (.cidx file) for the FASTA file to allow rapid access to specific sequences.

# Create a binary index of the FASTA file for fast sequence retrieval
cat all_hits_id.txt | cdbyank pdbaa.cidx > all_hits_in_order.fasta
# cdbyank retrieves sequences from a .cidx indexed FASTA file.

# Check the ordered sequences using the Dotter program
dotter query2.fasta all_hits_in_order.fasta

# Concatenate the query sequence with the ordered hit sequences to create a combined FASTA file:
cat query2.fasta all_hits_in_order.fasta > query_with_hits.fasta
# cat concatenates the files into a single FASTA file.

# Create a multiple sequence alignment using the MAFFT program
mafft --help
mafft query_with_hits.fasta > query_with_hits_aligned.fasta

# Inspect the aligned sequences with the 'less' command
less query_with_hits_aligned.fasta

# View the alignment with the Jalview program

```

The whole process can be automated with simple bash script. Create text file with following content:
```bash
#!/bin/bash
# Check if correct number of arguments is given
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <query_fasta> <database_fasta> <output_dir>"
    exit 1
fi
# Assign input arguments to variables
QUERY=$1
DATABASE=$2
OUTPUT_DIR=$3
# Create output directory and move to it
mkdir -p $OUTPUT_DIR

# Copy input files to the working directory
cp $QUERY $OUTPUT_DIR/query.fasta
cp $DATABASE $OUTPUT_DIR/database.fasta

cd $OUTPUT_DIR
# Create BLAST database
makeblastdb -in database.fasta -dbtype prot
# Run BLASTP search
blastp -query query.fasta -db database.fasta -out query_blastp_1e-10_table.txt -evalue 1e-10 -outfmt 6
# Extract sequence IDs of hits
cut -f 2 query_blastp_1e-10_table.txt > all_hits_id.txt
# Retrieve sequences of hits
seqkit grep -f all_hits_id.txt database.fasta -o all_hits.fasta
# Create binary index for database
cdbfasta database.fasta
# Extract sequences in order of significance
cat all_hits_id.txt | cdbyank database.fasta.cidx > all_hits_in_order.fasta
# Concatenate query and hits
cat query.fasta all_hits_in_order.fasta > query_with_hits.fasta
# Create multiple sequence alignment with MAFFT
mafft query_with_hits.fasta > query_with_hits_aligned.fasta
# Notify the user of completion
echo "BLASTP search and alignment complete. Results are in $OUTPUT_DIR."
```
 and save it as `blast2alignment.sh` and run following command to make script executable
 ```bash
chmod +x blast2alignment.sh
 ```

 This script automates the process of conducting a BLASTP search followed by a
 multiple sequence alignment. It takes three inputs: a query sequence (FASTA
 format), a protein database (FASTA format), and an output directory for storing
 the results. The automation offers several benefits:

- Reproducibility: The script ensures that the analysis steps are easily
   repeatable, enabling consistent results and reliable workflows.
- Efficiency: By
   integrating all commands, the script saves time and minimizes human error
   compared to running each step manually.
- Scalability: Automating the process
   makes it feasible to handle larger datasets and batch analyses, which are common
   in bioinformatics research.

The script carries out several steps, including creating a BLAST
database, running a BLASTP search, extracting and retrieving significant hits,
and aligning the sequences using MAFFT. The results are then stored in the
specified output directory.

To run script:
```bash
cp ~/Desktop/Bioinformatics/data/blast_data/query3.fasta .
./blast2alignment.sh query3.fasta pdbaa output3
# inspect output3 directory
```

## Exercise 2.3 - BLAST against remote database
Instead of having to download the entirety of NR or other NCBI databases, we can BLAST against the version held on the website. This ensures we have the most up to date version but is also significantly slower. We use the -remote command to do this. Lets BLAST out sequences against NR held on the NCBI website by typing:

<details>
<summary>Details</summary>
take  too long or does not work at all
</details>

```bash
blastp -query proteins.fasta -remote -db nr -out proteins_nr.txt -outfmt 6 -evalue 1e-30
```

```bash
# install  NCBI Entrez utilities on the command line
sudo apt install ncbi-entrez-direct

```

## Exercies 2.4
- visit ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot
- there are several file with protein sequences named `*protein.faa.gz`
- download these sequences, you can use web browser of you can try to use `wget` command.
  Note that  `wget` command can accept wildcards like {1..3}
- unzip downloaded fasta  ( use `gunzip` command)
- concatenate all fasta files into one (use `cat` commad)
- create blast database (`makeblastdb`)
- inspect restuls
- query sequence is protein from Danio rerio, you fill find it in `~/Desktop/bioinformatics/data/blast_data/danio_rerio_proteins.fasta`
- use `blastp` and find best hit to Danio rerio protein in Mus musculus refseq sequences
- what are these sequences, what is identity and e-value of the best hits
- which protein has more hits by blast

# Localization of sequences in genome using BLAT in Ensemble genome

Use BLAT to find location of the sequence in Caenorabditis elegans genome. For
search, use Ensemble database: https://www.ensembl.org/index.html

```
>C.elegans unkown sequence
GAATATTTAGGAGATGCAGGAGTTATTGATAGCTGCAGTGATTGGCCGAACACACCTCCT
GATAACAATTTTAATGGTAAGAGTTGAACTCCAAAACTGTAAGTAGAGGTGGCTGCTCTC
TCTCTCTGACTTTTATGCCTGCCTACGTACCTTCTAATACTTATTTGTTTGATATGGATG
TTTAGTGAAGATAAAGGGTAGATAGAGGCATTTCTCATCTGCCCAAGATGAGCATGAATA
TATTTAATACAAAATCAACACTGAGAATTTTAGAGACCGATTTTAAATGTGACCCAATTT
TTTTCAGGAGGATTTGCACCAGATCAACCTCAGCTAGTCACACCGATTATTTCTGATATT
CCGATAGATCTCAATCAAATATATGTTCCAACACCTCCACAATTACTTGATAATTGGTGT
TCAATCATTTATTATGAACTGGATACACCCATTGGTGAAACCTTTAAGGTATGTTTTTCT
ATGAAATCTGATGACTATTCATTCATGGTGCAAATCGCCTAGAAATTTTTGTGAAAGAGC
```

- What is the genomic location of sequences. Is the sequence part of any gene?
