# BLAST

## Exercise 1.1 - BLASTN - NCBI web interface
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

Introduction slides - BLAST databases (core_nt, ClusteredNR) and what the hit list does not show (slides 2-3):
[week4_intro.pdf](./slides/week4_intro.pdf)

### search I
setting of the search:
- In the section "Program Selection" select the option "Somewhat similar
  sequences (blastn)"
- As the search database keep the default "Core nucleotide database (core_nt)".
  `core_nt` contains GenBank nucleotide sequences, except most of the eukaryotic
  chromosome-level (whole genome) sequences, which makes it less than half the size
  of the full "Nucleotide collection (nr/nt)". Despite its name, `nt` is not
  really non-redundant - identical sequences from different submissions are
  kept as separate records.
- Expand "Algorithm parameters" at the bottom of the form and change "Max target
  sequences" from the default 100 to **250**. The database contains many insulin
  sequences from other mammals, and with only 100 hits shown the human sequence
  may not be on the list at all.
- Check "Show results in a new window" and  Click the BLAST button to launch the search.

> **Tip:** The BLAST form remembers some settings (e.g. the organism filter)
> from your previous search, even in a new browser tab. Before each new search,
> click **"Reset page"** (top right corner of the form) and check all settings.

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

First, take a look at the best hit. Since our search sequence (the query) was taken from GenBank, which is part of core_nt, we should find an identical sequence in the search. Make sure this is the case!
Find the best hit from human (Homo sapiens) that is not a synthetic construct
(tip: use Ctrl+F and search for "Homo sapiens" in the Descriptions table)
and answer following questions:
- what is the identifier (Accession)?
- what is the alignment score ("max score")?
- what is the percent identity and query coverage?
- what is the E-value?
- are there any gaps in the alignment?

### Search II
Make a new BLASTN search with the same query sequence, this time in the
"Database" section switch from "Standard databases" to "Genomic + transcript
databases" and select "Human genomic plus transcript (Human G+T)". For the best
hits answer following:
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

> **Note on protein databases:** the default protein database on the NCBI BLAST
> page is now **ClusteredNR**. In ClusteredNR, proteins that are >90 % identical
> (and of similar length) are grouped into clusters and only one representative
> per cluster is searched - the search is faster and the hit list is more
> diverse, but you will not see the individual (e.g. identical) sequences. For
> this exercise select **"Non-redundant protein sequences (nr)"** in the Database
> menu. Optionally, repeat the search with ClusteredNR and compare the hit lists.

- Based on the identical hit to C. elegans what is the identity of this protein?
- Aside from the C. elegans proteins, what is the most significant hit? What is
  its identity and E-value?
- switch for "graphic summary" - for blastp there is additional information
  about identified conserved domains.
- Explore link to other databases  - find coresponding entry in uniprot database.
- Check also taxonomy report. How many hits was found total. How many hits was
  found in "chordata" taxonomy group? If no hits are found in Chordata, try to modify search
  parameters to identify similar sequence in Chordata group.

<details>
<summary>💡 Hint</summary>

By default, BLAST reports only 100 best hits ("Max target sequences" in
"Algorithm parameters"), and these can all come from nematodes. You can either
increase "Max target sequences", or restrict the search to a taxonomic group by
typing e.g. "Chordata" into the "Organism" field (or use "exclude" to remove
Nematoda from the search).
</details>

## Exercise 1.3 - BLASTX - NCBI web interface
The part of genomic sequence of Twort bacteriophage is below. Use BLASTX to
search protein database (select "Non-redundant protein sequences (nr)", not the
default ClusteredNR) to find homologous protein from other bacteriophages and
viruses.
- For search try to use BLOSUM62 and BLOSUM45 scoring matrix ("Algorithm
  parameters" → "Matrix"). What is the difference between searches with
  different matrices?
- What proteins (Accessions) has alignment with E value < 1e-10?
- The encoded protein is very rich in lysine (K) and glutamate/aspartate (E/D).
  Look at the alignments - are there regions masked by low complexity filter
  (shown as `X` in the query)? Repeat the search with the filter switched off
  ("Filters and Masking" → uncheck "Low complexity regions"). How did the
  results change? Are the new hits trustworthy?

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
  - first run ORFfinder with default settings ("ORF start codon to use": "ATG
    only"), then repeat it with "Any sense codon". Why is the longest ORF
    different? (Hint: our sequence is a PCR fragment - does it contain start
    and stop codon?)
- use `blastp` against longest ORF, run blastp directly from ORFfinder, compare results with blastx search.
- What kind of enzyme is coded by `clone12`?
- does the clone12 represent complete CDS of putative protein?

<details>
<summary>Details</summary>
In blastx AA alignment is truncated probably because of low complexity
filtering, try to turn off filtering and run in again.
When using BLASTP from ORFfinder - select non-redundant protein database (nr).
</details>

- compare best hit from blastx directly using tblastn against DNA sequence of
  clone12
  - click on accession - this will take you to protein entry.
  - on protein entry page, click on left panel on *Run BLAST*
  - select *tblastn* - query is a protein accession (this represent query)
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
Spacer (ITS) region of nuclear ribosomal DNA. ITS1 and ITS2 are non-coding
spacers located between the 18S, 5.8S and 26S rRNA genes; they evolve faster
than the rRNA genes and are therefore often used for species identification
(DNA barcoding). You obtained the following sequences:
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
likely to correspond to the correct species matching your sequences. Use
`blastn` (megablast) with the default database core_nt.

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

Introduction slides - DNA gyrase, quinolones and the QRDR (slide 5):
[week4_intro.pdf](./slides/week4_intro.pdf)

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
search for each strain's protein sequence. (You can paste all four sequences
into the query box at once; on the results page, use the "Results for" menu to
switch between queries.)

Inspect different alignments views (pairwise, query-anchored, ...).

Questions:

- How many mutations did you find in each strain's protein compared to the
reference genome?
- Did any of the strains contain an insertion or deletion?
- Ciprofloxacin resistance in *N. meningitidis* is mostly caused by
  substitutions of **Thr91** (typically T91I) and, less often, **Asp95** of
  GyrA. Find these two positions in your alignments (in TBLASTN output, the
  query coordinates are the positions in the GyrA protein). Which strains are
  probably resistant to ciprofloxacin?
- Some strains differ from the reference at many positions, but only a few of
  them are relevant for resistance. Where do so many differences come from?
  (Hint: *Neisseria* species are naturally competent and readily take up DNA
  from related commensal *Neisseria* species living in the same host.)
- Can you use BLASTX instead of TBLASTN for this analysis? Why or why not?

<details>
<summary>Backup - if NCBI BLAST is not available</summary>

The same analysis can be done with local BLAST+, the reference genome is
available in `~/Desktop/Bioinformatics/data/blast_data/gyrA/`:

```bash
cd ~/Desktop/Bioinformatics/data/blast_data/gyrA/
# -subject is used instead of -db for comparison of two sequences (no database needed)
tblastn -query gyrA_strains.fasta -subject NZ_CP021520.1.fasta -max_hsps 1 | less
# summary table (pident, mismatches, gap openings, coordinates):
tblastn -query gyrA_strains.fasta -subject NZ_CP021520.1.fasta -max_hsps 1 -outfmt 7
```
</details>


## Exercise 1.8 - Spread of the colistin resistance gene mcr-1
Colistin is an old antibiotic that is used as a "last resort" drug against
multi-resistant Gram-negative bacteria. Resistance to colistin used to be rare
and was caused only by chromosomal mutations, which cannot spread between
bacteria. In 2015, Liu et al. (*Lancet Infect Dis* 16:161) described the first
**plasmid-borne** colistin resistance gene, *mcr-1*, in *E. coli* from pigs and
from patients in China. The gene encodes a phosphoethanolamine transferase,
an enzyme that modifies lipid A in the bacterial outer membrane so that the
positively charged colistin binds less well.

Introduction slides - colistin and lipid A, plasmids and mobile elements (slides 6-8):
[week4_intro.pdf](./slides/week4_intro.pdf)

The sequence of the *mcr-1* gene from the original plasmid pHNSHP45
(`KP347127.1`) is in
`~/Desktop/Bioinformatics/data/blast_data/mcr1/mcr1_cds.fasta`, the protein
sequence (MCR-1, `AKF16168.1`) is in
`~/Desktop/Bioinformatics/data/blast_data/mcr1/mcr1_protein.fasta`.

<details>
<summary>📄 Sequence <code>mcr1_cds.fasta</code> (mcr-1 gene, 1626 bp) - click to show, then copy and paste into the BLAST query box</summary>

```
>mcr-1 KP347127.1:22413-24038 Escherichia coli SHP45 plasmid pHNSHP45, mcr-1 gene (phosphoethanolamine transferase)
ATGATGCAGCATACTTCTGTGTGGTACCGACGCTCGGTCAGTCCGTTTGTTCTTGTGGCG
AGTGTTGCCGTTTTCTTGACCGCGACCGCCAATCTTACCTTTTTTGATAAAATCAGCCAA
ACCTATCCCATCGCGGACAATCTCGGCTTTGTGCTGACGATCGCTGTCGTGCTCTTTGGC
GCGATGCTACTGATCACCACGCTGTTATCATCGTATCGCTATGTGCTAAAGCCTGTGTTG
ATTTTGCTATTAATCATGGGCGCGGTGACCAGTTATTTTACTGACACTTATGGCACGGTC
TATGATACGACCATGCTCCAAAATGCCCTACAGACCGACCAAGCCGAGACCAAGGATCTA
TTAAACGCAGCGTTTATCATGCGTATCATTGGTTTGGGTGTGCTACCAAGTTTGCTTGTG
GCTTTTGTTAAGGTGGATTATCCGACTTGGGGCAAGGGTTTGATGCGCCGATTGGGCTTG
ATCGTGGCAAGTCTTGCGCTGATTTTACTGCCTGTGGTGGCGTTCAGCAGTCATTATGCC
AGTTTCTTTCGCGTGCATAAGCCGCTGCGTAGCTATGTCAATCCGATCATGCCAATCTAC
TCGGTGGGTAAGCTTGCCAGTATTGAGTATAAAAAAGCCAGTGCGCCAAAAGATACCATT
TATCACGCCAAAGACGCGGTACAAGCAACCAAGCCTGATATGCGTAAGCCACGCCTAGTG
GTGTTCGTCGTCGGTGAGACGGCACGCGCCGATCATGTCAGCTTCAATGGCTATGAGCGC
GATACTTTCCCACAGCTTGCCAAGATCGATGGCGTGACCAATTTTAGCAATGTCACATCG
TGCGGCACATCGACGGCGTATTCTGTGCCGTGTATGTTCAGCTATCTGGGCGCGGATGAG
TATGATGTCGATACCGCCAAATACCAAGAAAATGTGCTGGATACGCTGGATCGCTTGGGC
GTAAGTATCTTGTGGCGTGATAATAATTCGGACTCAAAAGGCGTGATGGATAAGCTGCCA
AAAGCGCAATTTGCCGATTATAAATCCGCGACCAACAACGCCATCTGCAACACCAATCCT
TATAACGAATGCCGCGATGTCGGTATGCTCGTTGGCTTAGATGACTTTGTCGCTGCCAAT
AACGGCAAAGATATGCTGATCATGCTGCACCAAATGGGCAATCACGGGCCTGCGTATTTT
AAGCGATATGATGAAAAGTTTGCCAAATTCACGCCAGTGTGTGAAGGTAATGAGCTTGCC
AAGTGCGAACATCAGTCCTTGATCAATGCTTATGACAATGCCTTGCTTGCCACCGATGAT
TTCATCGCTCAAAGTATCCAGTGGCTGCAGACGCACAGCAATGCCTATGATGTCTCAATG
CTGTATGTCAGCGATCATGGCGAAAGTCTGGGTGAGAACGGTGTCTATCTACATGGTATG
CCAAATGCCTTTGCACCAAAAGAACAGCGCAGTGTGCCTGCATTTTTCTGGACGGATAAG
CAAACTGGCATCACGCCAATGGCAACCGATACCGTCCTGACCCATGACGCGATCACGCCG
ACATTATTAAAGCTGTTTGATGTCACCGCGGACAAAGTCAAAGACCGCACCGCATTCATC
CGCTGA
```
</details>

<details>
<summary>📄 Sequence <code>mcr1_protein.fasta</code> (MCR-1 protein, 541 aa) - click to show, then copy and paste into the BLAST query box</summary>

```
>MCR-1 AKF16168.1 phosphoethanolamine transferase MCR-1 [Escherichia coli]
MMQHTSVWYRRSVSPFVLVASVAVFLTATANLTFFDKISQTYPIADNLGFVLTIAVVLFGAMLLITTLLS
SYRYVLKPVLILLLIMGAVTSYFTDTYGTVYDTTMLQNALQTDQAETKDLLNAAFIMRIIGLGVLPSLLV
AFVKVDYPTWGKGLMRRLGLIVASLALILLPVVAFSSHYASFFRVHKPLRSYVNPIMPIYSVGKLASIEY
KKASAPKDTIYHAKDAVQATKPDMRKPRLVVFVVGETARADHVSFNGYERDTFPQLAKIDGVTNFSNVTS
CGTSTAYSVPCMFSYLGADEYDVDTAKYQENVLDTLDRLGVSILWRDNNSDSKGVMDKLPKAQFADYKSA
TNNAICNTNPYNECRDVGMLVGLDDFVAANNGKDMLIMLHQMGNHGPAYFKRYDEKFAKFTPVCEGNELA
KCEHQSLINAYDNALLATDDFIAQSIQWLQTHSNAYDVSMLYVSDHGESLGENGVYLHGMPNAFAPKEQR
SVPAFFWTDKQTGITPMATDTVLTHDAITPTLLKLFDVTADKVKDRTAFIR
```
</details>

### Part 1 - Where is mcr-1 found today?
- Run **blastn** (megablast) with `mcr1_cds.fasta` against "Core nucleotide
  database (core_nt)".
- In "Algorithm parameters", increase "Max target sequences" to 500.

Questions:
- How many hits have ~100 % identity and 100 % query coverage? Does this
  number tell you how many sequences in the database carry *mcr-1*? (Hint:
  what limits the number of hits shown?)
- Open the "Taxonomy" tab (and its "Organism" report). In how many different
  bacterial genera and species was the identical gene found? Which species are
  most common?
- Look at the descriptions of the hits. Is the gene located on chromosomes or
  on plasmids?
- The bacterial species you found are separated by millions of years of
  evolution, yet they carry an (almost) identical gene. What does it tell you
  about how the gene spreads? Compare with Exercise 1.7: how would you expect a
  chromosomal gene like *gyrA* to look in different species?

### Part 2 - Where did mcr-1 come from?
A gene that spreads on plasmids usually originates on the chromosome of some
other bacterium. To find the possible source, we need to search for more
distant homologs at the protein level.

- Run **blastp** with `mcr1_protein.fasta` (or just type the accession
  `AKF16168.1`). Keep the default database **ClusteredNR** and set "Max target
  sequences" to 250.

In ClusteredNR, proteins that are >90 % identical are grouped into one
cluster. Thousands of identical MCR-1 proteins from Part 1 are therefore
represented by a single hit, and more distant relatives can appear in the hit
list. (In the full "nr" database, the first several thousand hits would all be
MCR-1 from Enterobacterales.)

Questions:
- The results are shown in the "Clusters" tab. Look at the first hit: how
  many members (proteins) and organisms does it contain (column "Cluster
  Composition")? Click on ⊞ to see the cluster contents.
- The column "Cluster Ancestor" shows the common taxonomic group of all
  proteins in the cluster. Find the best hits from bacteria that do **not**
  belong to Enterobacterales (*Escherichia*, *Salmonella*, *Klebsiella*, ...).
  To which genus and family do they belong? What is their percent identity?
- Some hits are named MCR-2, MCR-6 etc. These are other plasmid-borne *mcr*
  variants discovered after *mcr-1*. In which organism was the MCR-6 hit found?
- Read the descriptions of the hits - what is their function? What does it
  mean that a similar protein is found on the chromosome of other bacteria?
- Look at the "Graphic Summary" - which conserved domains does MCR-1 contain?

<details>
<summary>💡 Hint</summary>

Bacteria of the genus *Moraxella* (family Moraxellaceae, together with
*Psychrobacter* and *Faucicola* which you also find among the hits) carry
chromosomal genes for MCR-like phosphoethanolamine transferases, which are
the most similar known proteins to MCR-1 outside Enterobacterales.
*Moraxella* is therefore considered to be the most likely original source of
*mcr-1* (and of other *mcr* genes), which was later mobilized by the insertion
sequence IS*Apl1* to plasmids in Enterobacterales.
</details>

<details>
<summary>Backup - if NCBI BLAST is not available</summary>

**1. EMBL-EBI BLAST** (https://www.ebi.ac.uk/jdispatcher/sss/ncbiblast) -
Part 1: program BLASTN, database "ENA Sequence Prokaryote". Part 2: use
**UniProt BLAST** (https://www.uniprot.org/blast) with `mcr1_protein.fasta`,
target database "UniProtKB", and in "Restrict by taxonomy" enter
`Moraxella` (UniProt BLAST does not allow excluding taxa, so we restrict the
search to the candidate source group instead).
Compare the identity of these hits with the hits from Enterobacterales.

**2. Local BLAST** - the backup databases are in
`~/Desktop/Bioinformatics/data/blast_data/mcr1/backup_db/`:
`mcr1_hits.fasta` contains regions (*mcr-1* ± 1000 bp) of 36 sequences found
by NCBI BLAST in Part 1, `mcr_like_proteins.fasta` contains the proteins found
by NCBI BLASTP in Part 2.

```bash
mkdir -p ~/blast_backup/mcr1
cd ~/blast_backup/mcr1
cp ~/Desktop/Bioinformatics/data/blast_data/mcr1/*.fasta .
cp ~/Desktop/Bioinformatics/data/blast_data/mcr1/backup_db/*.fasta .
# Part 1
makeblastdb -in mcr1_hits.fasta -dbtype nucl
blastn -query mcr1_cds.fasta -db mcr1_hits.fasta \
  -outfmt "6 sacc pident length qcovs stitle" | less -S
# count the hits by genus (2nd word of the description)
blastn -query mcr1_cds.fasta -db mcr1_hits.fasta -outfmt "6 stitle" \
  | awk '{print $2}' | sort | uniq -c
# how many hits are on plasmids?
blastn -query mcr1_cds.fasta -db mcr1_hits.fasta -outfmt "6 stitle" | grep -c -i plasmid
# Part 2
makeblastdb -in mcr_like_proteins.fasta -dbtype prot
blastp -query mcr1_protein.fasta -db mcr_like_proteins.fasta -max_target_seqs 50 \
  -outfmt "6 sacc pident length qcovs evalue stitle" | less -S
```

**3. Pre-computed NCBI results** - if nothing else works, the results of the
NCBI searches (run in October 2026) are in
`~/Desktop/Bioinformatics/data/blast_data/mcr1/precomputed/`.
</details>

## Exercise 1.9 - Where did SARS-CoV-2 come from?
Soon after the first SARS-CoV-2 genome was published in January 2020
(isolate Wuhan-Hu-1, accession `NC_045512.2`), scientists all over the world
used BLAST to search for its closest relatives. The answer is important for
understanding the origin of the virus: which animal viruses are most similar,
and how similar are they? In this exercise, you will repeat this analysis.

Introduction slides - SARS-CoV-2 genome, spike RBD and recombination (slides 10-12):
[week4_intro.pdf](./slides/week4_intro.pdf)

### Part 1 - whole genome
- Open the nucleotide BLAST page (blastn). Instead of pasting a sequence, you
  can type the accession `NC_045512.2` into the query box.
- Database: "Core nucleotide database (core_nt)"
- Organism: type `Severe acute respiratory syndrome coronavirus 2` and select
  "taxid:2697049", then check **exclude**. GenBank contains millions of
  SARS-CoV-2 genomes from patients - without excluding them, all hits would be
  SARS-CoV-2 itself.
- Click "Add organism", type `synthetic construct`, select "taxid:32630" and
  check **exclude** as well. Laboratories have deposited many synthetic copies
  of the SARS-CoV-2 genome (e.g. clones used for reverse genetics).
- Program: "Highly similar sequences (megablast)"
- In "Algorithm parameters", set "Max target sequences" to 250
- Run BLAST

Questions:
- Even after the exclusions, the hits at the top of the list are >99.5 %
  identical to the query. Read their descriptions - what are they, and why
  were they not removed by the exclusion? (Hint: click on the accession and
  check the "ORGANISM" line of the record.) What does it tell you about the
  quality of annotations in public databases?
- Scroll down to the first hits with identity below 99 %. These are the real
  relatives of SARS-CoV-2. What are the five best of them? From which host
  species were these viruses isolated and where (check the hit descriptions or
  click the accession)?
- What is the percent identity and query coverage of the best hit? How many
  nucleotide differences does it mean over the 30 kb genome?
- Some hits (e.g. RpYN06 or BANAL-20-103) have higher percent identity than the
  best hit, but they are ranked lower. Compare their "Max score" and "Total
  score". Click "Graphics" next to such a hit (in the "Alignments" tab) - the
  alignment is split into two parts. Which region of the SARS-CoV-2 genome is
  missing in between? Which gene is located there? (Use the genome annotation
  of `NC_045512.2` to find out.)
- Look at the "Taxonomy" tab. Which viruses are found most often among the
  hits? Judging from their names, in which animals were they found?
- How similar is SARS-CoV-2 to the SARS coronavirus from the 2002-2003
  epidemic (SARS-CoV-1)? Use blastn "Align two or more sequences" with
  `NC_045512.2` as query and `NC_004718.3` (SARS-CoV-1 Tor2) as subject, and
  compare the identity and query coverage with the bat viruses found above.

### Part 2 - receptor-binding domain (RBD) of the spike protein
The receptor-binding domain of the spike protein binds to the human ACE2
receptor and determines which hosts the virus can infect. In the Wuhan-Hu-1
genome it is located at positions 22517-23185.

- Repeat the search with the same settings, but enter `22517` and `23185` into
  the "Query subrange" From/To fields. Because the region is short and more
  divergent, change the Program to "Somewhat similar sequences (blastn)".
  (If you use "Edit Search" on the results page, check that both excluded
  organisms are still set - the form does not always keep them.)

Questions:
- Which virus is the closest relative in the RBD region? Is it the same virus
  as in the whole-genome search?
- Find RaTG13 (`MN996532.2`) and the pangolin coronavirus MP789
  (`MT121216.1`) in the hit list. Compare their identity to SARS-CoV-2 in the
  whole genome (Part 1) and in the RBD. Which part of the genome is more
  divergent?
- The function of the RBD depends on its amino acid sequence. Compare the RBD
  at the protein level: use **tblastx** with "Align two or more sequences",
  query `NC_045512.2` with the same query subrange (22517-23185), and as
  subjects enter `MN996532.2`, `MT121216.1` and `MZ937000.1` (BANAL-20-52, one
  accession per line). Which virus has the RBD most similar to SARS-CoV-2 at
  the protein level? Is the order the same as at the nucleotide level?

<details>
<summary>💡 Hint</summary>

Different parts of the genome have different closest relatives. Coronaviruses
frequently **recombine** - when two viruses infect the same animal, the progeny
can carry a mosaic genome. The phylogenetic history of a recombinant genome
cannot be described by a single tree, and the "closest relative" depends on
which region you look at. RaTG13 is one of the closest viruses over the whole
genome (~96 %), but its RBD is much more different (~90 % amino acid
identity). The pangolin virus MP789 is more distant over the whole genome, but
its RBD is ~97 % identical to SARS-CoV-2 at the protein level. The BANAL
viruses from bats in Laos (collected in 2020) are close relatives both over
the whole genome and in the RBD, and their spike protein can bind human ACE2.
</details>

<details>
<summary>Backup - if NCBI BLAST is not available</summary>

**1. EMBL-EBI BLAST** (https://www.ebi.ac.uk/jdispatcher/sss/ncbiblast) -
the EBI service does not allow excluding taxa. If you search "ENA Sequence
Viral", all the best hits will be SARS-CoV-2 genomes. Try it, and you will see
why the "exclude" option was so useful. This backup is therefore not suitable
for this exercise.

**2. Local BLAST** - a database of 21 coronavirus genomes (SARS-CoV-2
relatives from bats and pangolins, SARS-CoV-1, and MERS-CoV as an outgroup)
is in `~/Desktop/Bioinformatics/data/blast_data/sars_cov2/`:

```bash
mkdir -p ~/blast_backup/sars_cov2
cd ~/blast_backup/sars_cov2
cp ~/Desktop/Bioinformatics/data/blast_data/sars_cov2/*.fasta .
cp ~/Desktop/Bioinformatics/data/blast_data/sars_cov2/backup_db/coronavirus_genomes.fasta .
# list of genomes in the database
grep ">" coronavirus_genomes.fasta
makeblastdb -in coronavirus_genomes.fasta -dbtype nucl
# Part 1 - whole genome (megablast is the default task of blastn)
# only the best HSP (local alignment) for each virus is reported (-max_hsps 1)
# columns: subject, % identity, alignment length, query coverage, E-value, bit score, description
blastn -query sars_cov2_genome.fasta -db coronavirus_genomes.fasta -max_hsps 1 \
  -outfmt "6 sacc pident length qcovs evalue bitscore stitle" | sort -k6,6nr | less -S
# all HSPs for one virus (RpYN06) - which part of the SARS-CoV-2 genome is not covered?
blastn -query sars_cov2_genome.fasta -db coronavirus_genomes.fasta \
  -outfmt "6 sacc pident length qstart qend" | grep MZ081381.1 | sort -k4,4n
# Part 2 - RBD
blastn -task blastn -query sars_cov2_rbd.fasta -db coronavirus_genomes.fasta -max_hsps 1 \
  -outfmt "6 sacc pident length qcovs evalue bitscore stitle" | sort -k6,6nr | less -S
# Part 2 - RBD at the protein level (tblastx translates both query and database)
tblastx -query sars_cov2_rbd.fasta -db coronavirus_genomes.fasta -max_hsps 1 \
  -outfmt "6 sacc pident length bitscore stitle" | sort -k4,4nr | less -S
# SARS-CoV-1 (Tor2) vs SARS-CoV-2 - direct comparison of two genomes
seqkit grep -p NC_004718.3 coronavirus_genomes.fasta > sars_cov1.fasta
blastn -query sars_cov2_genome.fasta -subject sars_cov1.fasta -outfmt 7
```

**3. Pre-computed NCBI results** - if nothing else works, the results of the
NCBI searches (run in October 2026) are in
`~/Desktop/Bioinformatics/data/blast_data/sars_cov2/precomputed/`.
</details>

## Exercise 1.10 - Syncytin: a viral gene that makes the placenta
About 8 % of the human genome consists of endogenous retroviruses (ERVs) -
remnants of ancient retroviral infections of germ-line cells, inherited by all
descendants. Most of them have accumulated mutations and are no longer
functional. A few viral genes, however, were "domesticated" and now serve the
host. The best known example is **syncytin-1**, the envelope (*env*) gene of
the HERV-W endogenous retrovirus. The retroviral Env protein normally fuses
the viral membrane with the membrane of the host cell. Syncytin-1 uses the same
ability to fuse placental cells (trophoblasts) into the syncytiotrophoblast, a
multinucleate layer separating maternal and fetal blood.

Introduction slides - endogenous retroviruses, syncytiotrophoblast and dating an insertion (slides 14-16):
[week4_intro.pdf](./slides/week4_intro.pdf)

Human syncytin-1 (UniProt `Q9UQF0`) is in
`~/Desktop/Bioinformatics/data/blast_data/syncytin/syncytin1_human.fasta`.

<details>
<summary>📄 Sequence <code>syncytin1_human.fasta</code> (human syncytin-1, 538 aa) - click to show, then copy and paste into the BLAST query box</summary>

```
>sp|Q9UQF0|SYCY1_HUMAN Syncytin-1 OS=Homo sapiens OX=9606 GN=ERVW-1 PE=1 SV=1
MALPYHIFLFTVLLPSFTLTAPPPCRCMTSSSPYQEFLWRMQRPGNIDAPSYRSLSKGTP
TFTAHTHMPRNCYHSATLCMHANTHYWTGKMINPSCPGGLGVTVCWTYFTQTGMSDGGGV
QDQAREKHVKEVISQLTRVHGTSSPYKGLDLSKLHETLRTHTRLVSLFNTTLTGLHEVSA
QNPTNCWICLPLNFRPYVSIPVPEQWNNFSTEINTTSVLVGPLVSNLEITHTSNLTCVKF
SNTTYTTNSQCIRWVTPPTQIVCLPSGIFFVCGTSAYRCLNGSSESMCFLSFLVPPMTIY
TEQDLYSYVISKPRNKRVPILPFVIGAGVLGALGTGIGGITTSTQFYYKLSQELNGDMER
VADSLVTLQDQLNSLAAVVLQNRRALDLLTAERGGTCLFLGEECCYYVNQSGIVTEKVKE
IRDRIQRRAEELRNTGPWGLLSQWMPWILPFLGPLAAIILLLLFGPCIFNLLVNFVSSRI
EAVKLQMEPKMQSKTKIYRRPLDRPASPRSDVNDIKGTPPEEISAAQPLLRPNSAGSS
```
</details>

<details>
<summary>📄 Sequence <code>syncytinA_mouse.fasta</code> (mouse syncytin-A, 617 aa, used in Part 2) - click to show, then copy and paste into the BLAST query box</summary>

```
>sp|Q5G5D5|SYNA_MOUSE Syncytin-A OS=Mus musculus OX=10090 GN=Syna PE=1 SV=1
MVRPWVFCLLLFPCSSAYSDSWMPLVNLTQHLLQEANSSFSSNCWVCLSIQTQRSLAMPA
PLRTWTETPMKLRIMYSARTLSGPYPITDLERRLQNFQPLTPHSSFVNPDQRAIAFLQIT
SVTGILPILSRITSVRYPDDHVYESAQRPIWGSLSTQTILTSQAPLCISRFFKNSNHATF
VGKLPASLCNHTFQLSPSANHQSIDLSSSYAFAPLMAMPGSKWRNPLRFSGPPSLNSGMP
HYSCPIDDIHCHTYPTTPWRSCPSFPASTCYNLTLFEPDNSSHPITLSVDTTYFKIKLQG
HKDPYPLFQYQPLMGAALSGQYSIWEYEPTVKKNGGITPNIFSHLVSLTYSFCLNSSGVF
FLCGNSTYVCLPANWSGVCTLVFQYPDIELLPNNQTISVPLFATVPSSVPASRRKRALPL
LPLLAGLGIASALGLGIAGITTSTVYFQQLSKALSDSLDEIATSIISLQDQIDSLAGVVL
QNRRALDLIVAERGGTCLFLQEECCFYINQSGVVRHAARKLRERASELGTSSSSWIQWLG
LGPWLPSWLTSLMAPILFILVLLVFRPCLLNCLTHSVSRRMSSFIHTTTEGHVDKILLLR
ESQYKRLPQEPPEEDAV
```
</details>

### Part 1 - How many copies of HERV-W env are in the human genome?
- Open the **tblastn** page (protein query vs. translated nucleotide database).
  You can type the accession `Q9UQF0` into the query box.
- Database: "RefSeq Genome Database (refseq_genomes)", Organism: type
  `Homo sapiens` and select "Homo sapiens (taxid:9606)" (not the entries with
  a gene name in brackets, such as "Homo sapiens (PARIS)")
- Run BLAST

Questions:
- Each hit in the list is a whole chromosome (or a scaffold), and each
  chromosome contains many separate alignments (HSPs). On which chromosomes
  did you find hits? Open the alignments for one chromosome - how many HSPs
  does it contain?
- How many HSPs cover (almost) the whole query (~538 aa) with high identity?
  On which chromosomes are they?
- Inspect the alignments of a few weaker hits. Look for stop codons (`*`) and
  frameshifts in the translated subject sequence. Could these copies produce a
  functional protein?
- Where is the copy that encodes the functional syncytin-1 located? (Hint:
  UniProt entry `Q9UQF0` lists the gene name *ERVW-1*.)

### Part 2 - Which organisms have syncytin-1?
- Run **blastp** with the same query. Keep the default database
  "ClusteredNR" - each hit represents a cluster of very similar proteins, so
  the hit list is not flooded by hundreds of almost identical sequences and
  the taxonomic overview is easier to read. Set "Max target sequences" to 250.

Questions:
- Look at the "Taxonomy" tab and at the identity of the hits. In which group
  of mammals do you find proteins with >80 % identity to syncytin-1? What is
  the identity of hits from New World monkeys (e.g. *Aotus*, *Callithrix*) and
  from other mammals (carnivores, ruminants, rodents)?
- Old World monkeys and apes separated from New World monkeys about 40 million
  years ago. What does your result tell you about when the HERV-W retrovirus
  entered the genome of our ancestors?
- Mice also have a placenta with a syncytiotrophoblast. Find the hits from
  *Mus musculus* in the list. What are they called, what is their identity
  and query coverage? Are they orthologs of syncytin-1?
- Mouse placenta uses syncytin-A (UniProt `Q5G5D5`,
  `~/Desktop/Bioinformatics/data/blast_data/syncytin/syncytinA_mouse.fasta`).
  Compare it with human syncytin-1 using blastp "Align two or more sequences".
  How similar are they? Which part of syncytin-1 is aligned? Syncytin-1 is
  cleaved into the surface subunit SU (residues 21-317, binds the receptor on
  the target cell) and the transmembrane subunit TM (residues 318-538, performs
  membrane fusion) - which subunit is more conserved and why?

<details>
<summary>💡 Hint</summary>

HERV-W integrated into the genome of the common ancestor of Old World
monkeys and apes (catarrhine primates), after their separation from New World
monkeys. Syncytin genes were captured independently from different
retroviruses in different mammalian lineages (primates, rodents, carnivores,
ruminants, ...). This is an example of **convergent evolution**: the same
function (cell-cell fusion in the placenta) is performed by proteins of
independent viral origin. The proteins are recognizably related (all are
retroviral Env proteins, the fusion machinery in the TM subunit is the most
conserved part), but they are not orthologs.
</details>

<details>
<summary>Backup - if NCBI BLAST is not available</summary>

**1. Other web services** - Part 1: **Ensembl BLAST**
(https://www.ensembl.org/tools/blast) - paste the protein sequence, select
species "Human", database "Genomic sequence" and program "TBLASTN". You can
open each hit in the Ensembl genome browser and see which genes are nearby.
Part 2: **UniProt BLAST** (https://www.uniprot.org/blast), query `Q9UQF0`,
database "UniProtKB". In the results, use the "Taxonomy" filter on the left
side to see in which organisms the hits are found. The "Align two or more
sequences" task can be done with EMBOSS Needle
(https://www.ebi.ac.uk/jdispatcher/psa/emboss_needle) or with local BLAST
(see below).

**2. Local BLAST** - the backup databases are in
`~/Desktop/Bioinformatics/data/blast_data/syncytin/backup_db/`:
`herv_w_env_loci.fasta` contains the regions of the human genome (GRCh38) found
by NCBI TBLASTN in Part 1 (± 500 bp), `syncytin1_blastp_hits.fasta` contains
the proteins found by NCBI BLASTP in Part 2.

```bash
mkdir -p ~/blast_backup/syncytin
cd ~/blast_backup/syncytin
cp ~/Desktop/Bioinformatics/data/blast_data/syncytin/*.fasta .
cp ~/Desktop/Bioinformatics/data/blast_data/syncytin/backup_db/*.fasta .
# Part 1 - tblastn against HERV-W env loci from the human genome
makeblastdb -in herv_w_env_loci.fasta -dbtype nucl
tblastn -query syncytin1_human.fasta -db herv_w_env_loci.fasta -max_hsps 1 \
  -outfmt "6 sacc pident length qstart qend qcovs evalue bitscore" | sort -k8,8nr | less -S
# alignments - look for stop codons (*) in the translated subject sequences
tblastn -query syncytin1_human.fasta -db herv_w_env_loci.fasta -max_hsps 1 | less
# Part 2 - blastp against proteins similar to syncytin-1
makeblastdb -in syncytin1_blastp_hits.fasta -dbtype prot
blastp -query syncytin1_human.fasta -db syncytin1_blastp_hits.fasta -max_target_seqs 250 \
  -outfmt "6 sacc pident length qcovs evalue stitle" | less -S
# mouse syncytin-A vs human syncytin-1
blastp -query syncytin1_human.fasta -subject syncytinA_mouse.fasta
```

**3. Pre-computed NCBI results** - if nothing else works, the results of the
NCBI searches (run in October 2026) are in
`~/Desktop/Bioinformatics/data/blast_data/syncytin/precomputed/`.
</details>

## Exercise 1.11 - Is my genome assembly contaminated?
You sequenced the genome of a *Staphylococcus aureus* isolate and assembled it
(you will do this yourself in the genome assembly practicals). Before
submitting the assembly to a database, you should always check it for
contamination. DNA from other organisms can get into your sample from many
sources: the person handling the sample, laboratory reagents ("kitome"), other
samples sequenced on the same run, cloning vectors, or control DNA added by the
sequencing facility. Assembly programs can also make errors and join unrelated
sequences together (chimeric contigs).

Introduction slides - sources and signs of contamination (slides 18-19):
[week4_intro.pdf](./slides/week4_intro.pdf)

Your draft assembly contains 11 contigs. The FASTA header of each contig shows
its length and its average sequencing coverage (`cov`):
`~/Desktop/Bioinformatics/data/blast_data/contamination/draft_assembly_contigs.fasta`

```bash
cd ~/Desktop/Bioinformatics/data/blast_data/contamination/
grep ">" draft_assembly_contigs.fasta
seqkit stats draft_assembly_contigs.fasta
```

If you work on a computer without the course data, copy the contigs from here
and save them as `draft_assembly_contigs.fasta` (or paste them directly into
the BLAST query box):

<details>
<summary>📄 Sequence <code>draft_assembly_contigs.fasta</code> (11 contigs, 88 kb) - click to show, then copy and paste into the BLAST query box</summary>

```
>contig_01 len=7000 cov=130x
TTGCAACTGAAATTATCAATTTAGAATCAATTTTAGAATTACCTAAAGGTACGGAACATTTCGTCAGCGA
TTTACATGGTGAATATGAAGCTTTCCAACACGTATTACGCAACGGTTCTGGGAACGTGCGAGCGAAAATC
AATGATATTTTCAAAGAGAGACTTTCAACTAAGGAGCTTAATGACTTAACTGCTCTTGTCTACTATCCAG
AAGACAAATTAAAATTGATTAAAAGTGATTTCCAAAATTGCGGTCAACTTAATGTCTGGTATATCACAAC
AATCGAACATTTAATTGAGTTAATTAAATATTGTTCTTCAAAATATACGCGTTCAAAACTGCGTAAAGCA
TTGCCAAAACAATACGTTTATATTATTGAGGAACTACTGTATAAAAGTAATGAATATCAAAATAAAAAAT
CATATTACGAAACACTTGTTAACCAAGTAATTGAACTTAAACAGGCAGATGATTTAATTATTGGACTCGC
TTATTCCGTACAACGCTTAGTCGTCGATCATTTACACGTTGTCGGTGATATTTATGATCGTGGACCACAA
CCAGATAAAATTATGGATACACTGATTAATTATCATTCCCTAGATATTCAATGGGGTAATCATGATGTGC
TTTGGGTTGGAGCCTATGCTGGGTCAAAAGTATGCTTAGCAAACTTACTTCGAATTTGTGCACGTTATGA
TAATTTAGATATTATCGAAGACGCTTATGGCATTAATTTAAGACCACTGCTTACTTTAGCTGAAAAATAC
TATGACGCAGATAATCCTGCTTTTAAGCCTAAAAAAAGACCTGACAAACACGAACGTTTAACTCAACGTG
AAGAAAGTCAAATTACTAAAATTCATCAAGCTATTGCGATGATTCAATTCAAGTTAGAAATACCAATTAT
TAAACGTCGTCCAAATTTCGAAATGGAAGAACGTCTTGTGCTTGAAAAGGTTAATTATGATACAAATGAA
ATTACAGTTTATGGTAATACATACCCATTGAAAGACACATGTTTCCAAACTGTCAATCGTGATAATCCAG
CAGAATTACTACCTGAAGAAGAAGAAGTCATGAATAAACTATTATTGTCATTCCAACAATCTGAAAAATT
ACGTCGTCATATGTCTTTCTTGATGCGTAAAGGCTCTCTTTACTTACCATATAATGGCAATTTACTCATT
CATGGTTGTATTCCAGTTGATGAAAATGGTGAGATGGAATCATTTGAAATTGATGGTCATACTTACAGCG
GCCAAGAATTATTAGATGTGTTTGAGTATCATGTCCGTAAATCATTTGATGAAAAAGAAAATACTGATGA
CTTATCGACGGATTTAGTTTGGTATTTATGGACTGGGAAATATTCGTCACTATTTGGTAAACGTGCCATG
ACTACGTTTGAGCGATACTTTATTGCAGATAAAGCTTCTCATAAAGAAGAAAAGAATCCGTACTATCATC
TTCGTGAAGATGTGAATATGGTTCGTAAAATGCTCAGTGATTTCGGATTAAATCCAGATGAAGGACGCAT
TATTAATGGTCACACACCAGTGAAAGAAATCAATGGCGAAGATCCTATCAAGGCTGATGGAAAGATGCTT
GTCATTGATGGTGGCTTTTCAAAAGCTTATCAGTCAACAACTGGCATTGCAGGCTATACACTATTGTATA
ATTCATTCGGTATGCAACTCGTTGCACACCAACAATTTAACGCTAAAGAAAAAATACTTTCCGAAGGTAT
CGATGAACTTTCTATAAAGCGCGTCGTAGATAAAGAATTACAACGTAAGAAGATTCGTGATACCAATATT
GGTAAAGATCTTCAAGCTCAAATAGATATTTTGAAAATGTTAATGCATGATCGTTACTTAGACTAAAGAC
GACATTTGATGCATGTATCGTGCACTTATTTAGTGGTACTGTGTTTACGCATTGGCTACGCAAACACGAT
GATAAATAGGCGACCTAGCTACACCGAAAAAAATCCCCTCACTATTGCGAGTAGTGAGGGGATTGGTGTA
TCTATGCTAATTCCTGTTAACTGTATTATAACAACTTTACAGCAAAAAAATGCAAAATCAACTAATAAAT
TAATGAAATGACCGACACCCTCGCTTCTACTCCAATTTAACAAAATACTATATTAGTACTATATTTTTCA
TTAACATGTATGTTACATTTTCTATATTACTTAAAATTACGAAGGGACATTTTTAAATCATGCGTATCTT
AAAAGAGTCCATTATTGTGGCATTTGCCTTTGTTGGTGTTGTCGTTGGTGCCGGCTTTGCTACTGGTCAA
GAAATTTTCCAGTTTTTCACAAGTCATGGCGCATATAGCATTTCAGGCATTATTGTAACAGGACTATTGA
TTACTTTAGGTGGAATGGTTGTCATGCATACAGGTCATCATCTAAAGTCCAGAAATCATTCTGATTCAAT
TAACTATTTCTTATACCCCTCTATTGCAAGAGGTTTTGATATTATTTTAACAATGTTTATGTTGTCTTTA
GCTATTATTATGACTGCAGGTGGTGCGTCAACCATTCATCAAAGTTTCAACTTACCGTATTGGCTGAGCG
CACTCATATTAGTCGCCTTTATTTTAGCAACACTGTTTCTAAAATTCGATCGTTTAATTGCTGTGCTTGG
CGGTGTTACCCCATTTTTAATTGCGATTGTCATTATGATTGCGGTCTACTATTTCACAACAAGTCATCTT
GATTTTACTGCCGCTAATAATGATGCTCAGATTCATAAGCAGAAATCATTATCACCTGGATGGTGGTTTG
ATGCGATTAACTATGCAAGCTTGCAAATTGCTGCTGCCTTCAGCTTCTTATCAGTGATGGGTAGTAAAGT
TAAATATCGTGACTCAACGTTATACGGGGGCTTGATTGGCGGTTTAATCATTACATTTTTACTCATGATG
ATTAATCTAGGTTTAATTTCTCAATTCGATAAAATTAAACACGTAGATCTACCTACATTAAAATTAGCGA
CACAAATGTCTCCGTCAATTGGTATTATTATGTCTGTCATTATGATACTTGTCATCTACAATACTGTTGT
TGGATTAATGTATGCATTTGCGTCACGTTTCAGCGTTCCATTCAGCAGACGTTACTTCATCATTATTATT
ACAATGGCTGTCATCACTTATATTAGTACATTTATCGGTTTCATTTCATTAATTGGAAAAGTATTCCCTA
TTATGGGATTGTTCGGTTTCATCTTACTCATACCTGTACTCTATAAAGGTTTAATTAAGCGTATTACCGG
CAAATCTCATATCGATTAAAATATCTATACTATTGTAAAATTTAGTACAGCACTTCAAAATCAAAAAAGC
GGGCAGATGAAAAGTAAAAATCTTTTCCATCTGTCCGCTTCTTAATAGAACTGTTATTTCCCAAGTATTT
GTTGACCAGCCGTTAATCCAGTTTCTGTAATTTCATGGCCCTTCACCCAAACTTCTTCGACTTGTGCCCC
ACGTGTATTAAACAAGTTAATGACTTGTTCACTTGCAGCTAATGGCACAATTGGATCATGTTTCCCCATA
GAAAGCAACACACTGACATCTGATAAATCCTTTGTTGACGTTACTTCAACTGGGTATAACGGTGCATATA
ACAATGCTTTTTTAAATGGTGCTTCTGAACGCAACATTAAGTTAATCGCTATATTTGATCCATTTGAAAA
TCCAACAAGTACTGCTTTTTCAATATCAAAATCATAACGTTCAGCAGCTTCTTTAATGAACGTCAACAAT
TCTTGTCCACGAAATGCCAAATCTTCTTCATCATAAACACCTTCACCAAGACGTTTGAAATAACGGTTCA
TCCCATTTTCTGAAACTTGTCCTCTAATACTTAACAAGTGATAATTTTCATTCAATGCTTCGCCTAACGG
TAATAAATCGAACTCATCACCACCTGTACCATGCAATAATATTAGTGTTGGCGCACCATTTTGTCCTTCT
CTAAAAATATGTTCCATCACAAATCCTCCTCATTCCAATTAACCATGTTGACGCTTCGTATTAAAAGGTC
TAACTTCCGATTCAATATATTCTCTTTTATTTTCTAAAAATGGTGGTAAGGATAACCCTTCGCCTAATGT
TTCATAAGGTTCATCTTCCATAAATCCTGGTCCATCTGTTGAAATTTCTATTAAAATATGCCCCACACGT
GCATATAATGCTTCAAAATAGAAACGATTAACGATGCCTGAGTTATTAATACCTACCTCTTTATATTTCG
TTGCCCACGCTTCTATTGCATCATGATCTTTCACACGAAATGACACATGATGTACCTCACCATAACCTTG
ACGTGCTGCTGGCCCTTTATCATCTTTTATTAAGATTACCTGGCCACCATTGCCTCCTTCGCCAACTTCA
AGTAATGCGACATTATCTTCATGCGCAATAGTTGTCATACCGTAAACAGTCTCTAAAATATTTTTAAAGT
CGTCAAAATAACTTACTTTAATTTCAATGGGGCCTAATCCATAAATCGCTTTATCTACTGGAACCGGTCC
ATTCTTCCAAGGTACACCAGGTGCTACCCCTTCATTTAACTCATCTGAAATTAATTGATACACTTGGCCA
TCGACTTCTTCAAATGGCAACACTTTTTTACCAAATAATTCTTGAATACCTTCGTGTTTAACACCAAACT
CATCAAAGCGCTGTTCATAATATGTTAATGCGTCATCGTTAGGCACTCTAAAAGACGGTCTTGTAATGGA
ATTTGTTCCTGCCTGCCCTTTTGTAATATTTGGAAAATCAAAGAACGTCATGTCTGTACCTGCCGAACCT
ACATCATCTGCAAAAAAAGTATGATACGTATAAATATCATCTTGATTGACTGTCTTTTTAACTAAACGCA
TGCCTAGTACTTCTGTAAAAAATTTATAATTACGTTCTGCATCATCTGTCATTGCAGTAACATGGTGAAT
ACCTAGTAATTCATGATTATTTATCATTTTATATACCCTCTTTATATTTTTAGTGCATTTATCTCGAAAT
CGAGATACTTTAATGTAAATTTTTAGTGCACGCCTTACACTTCTGTAGATTGTGCACTTAGTTTCTTAAA
CGCTTGTTGTAAGATTGTTAATTCATCCTTTGTTAACACATCAAACGCTTTTGTTAATGTCTCAGCATGC
TTAGGGAAAATATCTGCCATTTGACTTTGACCTTTTTCAGTTAAACAAGCCATATATACACGTTTATCAT
CTTTATCCTTTTCACGTGTAATCCAACCTTTGTCCTCTAATTGACTTACAACATATGAAATGCTGCTACT
TGCAATTAATACGCGGTCTCTAATACGTTGAATTGGTTGCGGACCTTTATTATAAAGCAACTCGAGCACT
GCAAATTCAGTAATATTTAAGCCATATCGCTTTACGTCTTCTTTTGTAATTTGCTCTAATGTGTCTAACG
CCCTATTCATTCCGACAAAAACATTGAGAGATTGTTTCGTTCGATCCATATGATCACCAACTTTTTTATT
TTATTTCGATTTCGAGATAAATTATATAGCTGAATAGATTTAATGTCAAGATACGTGTTTTCATCTTCAT
TTCAAACTGTGCTTAATAAAAACGACTACATAGATGTGGCAATGTCGTCACATGACTGTCATGCACCTTT
GACTTCGTTTAAGTTAAAATAGAATAGGGAAATACTATAATCAATCATGATAGGAGTTGATATGTATGAG
TAACCTTGAAATCAAACAAGGCGAGAACAAATTCTATATTGGTGATGATGAAAATAATGCTTTAGCTGAA
ATCACATACCGTTTTGTGGATAATAATGAAATTAACATTGATCATACAGGCGTATCTGATGAACTTGGTG
GTCAAGGTGTTGGCAAAAAACTACTTAAAGCAGTTGTTGAACACGCTCGAGAAAATAATTTGAAAATTAT
TGCCTCATGTTCATTTGCCAAACATATGTTAGAAAAAGAAGATTCATATCAAGATGTATATCTTGGTTAA
ATAAATTTTGCGATGGTGCAACTGTACCATCGCTTTTTTGTTCGCATCAATATAAAATTGCAAAACGACT
AGTTATTTTGCATTTCAAACTTTATGCTATCAATTCATTTTCTTTTCTGGTATACTTTTTATACTAACTT
ATAGAAAAATAGAAAAAGAGGCGATGCATAATGGCATTTCACGATAAAACAGCAACACAAGTGACAAACA
TTGTACTGAACGTAAGAGACTTAGACTTAATGACAACATTTTACAAAAACATTTTAGGACTATCCGTTAA
AAGTTCTGACGATAATACAACCGTATTATCTGTTGGGACTGGTGGCCATACTCTAACGTTACATTTATTA
GAAGACGGCCGTCAGACTTCCCCACGTGAAGCAGGGCTTTTTCATATAGCATTTTTATTACCAACTACTG
AGGATCTAGCTAACTTCTTATATTTCGTGGCACAAAAAAATATGGGCATCGGCGCTGGTGATCATTTAGT
AAGTGAAGCTTTATATTTCAACGACCCCGAAGGTAATGGTATTGAAGTCTATCGCGATAGACCTTCATCT
TCATGGGAATGGCAAAATGGCAAAGTTAAAATGGATACTTTAGAAGTTGATAGCCAAACCTTATTAACAC
ATCGTACTGATGAAGGTTGGCAAGGAATGCCAGCAAAAGGCATGATTGGACACTTACATTTAAAAACACA
TGATTTAGACGCAGCATATCAATGTTACATTGAACAATTAGGATTCCAACATGTGTCTGATTTTCCACGT
GCGCTATTTATGTCAACGAATCATTATCACCATCATATTGCTGCTAATACTTGGCAATCCAATCAAGCAA
GACAAAACAATTCACAAAGTTATGGATTAACACATGTTGATATTTATCAACCCGACGCCTCAACACATAC
ATTTACGGCACCAGAAGGGTTTGACATCACAGTTCATTCAAATACTGACCTTGTTCCAGAAAAATAATCA
>contig_02 len=2686 cov=25x
TCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCT
GTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCTGG
CTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGCACCATATGCGGTGTGAAATACCGCACAGAT
GCGTAAGGAGAAAATACCGCATCAGGCGCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATC
GGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTA
ACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCGAGCTCGGTACCCGGGGAT
CCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATT
GTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATG
AGTGAGCTAACTCACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAG
CTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGC
TCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACG
GTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAAC
CGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGAC
GCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCT
CGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTG
GCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTG
TGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGT
AAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGT
GCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTC
TGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAG
CGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATC
TTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAA
AAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTA
AACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCA
TCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTG
CTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAG
GGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCT
AGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCAC
GCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCAT
GTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTA
TCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGA
CTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTC
AATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGG
CGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGAT
CTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAA
GGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTAT
CAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGC
GCACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACCATTATTATCATGACATTAACCTATAAAAA
TAGGCGTATCACGAGGCCCTTTCGTC
>contig_03 len=9000 cov=121x
GCATTTTTAATCATTTGCCCAATACTTATAAGTTCATTCAATGGATAAAAGAAAGGGGGCAGTGACTTAA
CACTTGTGCCCTCTTCAAAATTATCTTTAATGACTTGTTCAGTAATCAATGCCTGAATTGCTTCATAAAC
ATCGTCGGTTTTTACATTTAATTCTCTGGCAATTGATAAATTACTTGGAAGTTGGTCACCATATTCCAAT
TTATTTTGCTTAATCAAATTTAAAATGTATTGTTTAACGTTATATACCGCATTCATTTCTGACATTTTTT
TACACCTCGTCATTATTAAGATCTATTAAAACTTGTCTAGGCTTGCTTCCTTTTTGTGGACCGATTACCT
GATTGCGTTCTAAATCATCCATCAACCTAGATGCTCTATTATAACCAATTCTAAATTGGCGTTGTAACAA
TGATGTACTTGCCTTTTGTTGTTCAACAACAAACAAATACGCTTCATCATATAAAGCATCTTCACTTTTC
ATTTCCGATTTATCCACTGGTGCATCTGGTTCCATTTCTTTTACATAATTTGCCTGTTGTTGTTCTACTA
CATAATTTACAACATCTTGCACCTCTTGGTCACTTAAAAACGCCCCTTGAATACGTGTTTGTGATGAGTC
ACCATTTCCAACGTATAACATGTCACCTTTACCAAGTAACTTTTCTGCGCCGCCAGTACCAATAATAGTT
CTTGAATCTGTTTGAGAACTCACAGCAAAAGCTATTCTAGATGGAATATTATTTTTAATGATACCTGTAA
TTACATCCACAGAAGGTCTTTGTGTCGCTACAATTAAATGTATACCTGCTGCACGTGCCATTTGTGTAAT
ACGTTGAATCGCATTTTCAACTTCTTTACCAGCTACCATCATTAAATCTGCAAGCTCATCTACAATAACA
ACGATATAAGGTAATTCTGGTTGCTTCTCATCTAATTCTTGATTTTGCTTACGGATTAATTCGTTATAAC
CTTTAATGTTTCTAGTTGATGAATGTTGGAATAAATCATAACGTCTTTCCATCTCAGCTACAATTTTTTC
TAAAGCTTGAGCAGCTTTATGAGGATTTGTAACAACCGGAATTAATAAGTGTGGAATTCCGTTATAAACA
TTTAGTTCAACCATTTTCGGATCGATTAACATAAGTTTAACTTCATGCGGCTTAGCATTTAATAAAATAC
TTGTAATAATACCATTTATACAAACAGATTTACCACTACCCGTCGATCCTGCCACCAATAAGTGTGGCAT
TTCATTTAGTGGAACAGTAATTGGATCACCTGATATATCTCTTCCTAATCCAACTTCTAGTTTATTATTA
GACGGGAATTTTTCATCTAAAACTTCTTTTAGTGAAACTAATGAAATTTTCTCATTTGGCACTTCAATAC
CTACTGCAGAACGACCAGGTATTGGCGCTTCGATTCTAACATCTTTTGCTGCTAAAGCTAATGCAATATC
ATTATGCAAGTTTACAATTTTACTCACTTTAACCCCTTGAGCTGGTTGAATTTCATATTGAGTTACTGCA
GGACCAATTTTAATTTGTGTCACTTTTGCATTTACCCCAAAATCTTTTAATGTATTCTCTAGTACTTGTC
CTTTACGTTGCACTTCAGCTTTAGATGTTGCTTTTTGTTTTGCAGGTTGATTAAGTAAAGTTAACGGTGG
AACAACATACGATACATTCGTAACTTCACCAGCTTCTTCAATAGTGTTTTCACTTTCAACACTGTTATGA
GTTTGTTCTGTTAATTGTTCTTGCTGATCTGCTTGATGATGATTTACGATGTTATTCGAACTATTCTCTG
CATCAAACACTCGTTTTTTTCGACTTGGTTGACTCTGGCTTTTACTTTCATTTTCAGTATGCCCATAAAT
TGGAATATCTCTTTCTTGAGGCACTTCCGTAAAATCACTCACATCTTTTATTTGTGGCTGTTCATTTTGA
CGTGCCTTTTGTTCTTCTTTAAGTCTTGCTTTTTCTTCACGCTTCAATTGTTTTTCTTGGTTTCTTTCCG
ACATTTTTTCATTAAATGAACCAAACCAAGCTTTTATATTTTCCAGTGCAACTTTTGCAACTTCACGATG
TTGATGGTTTGTTAATAAAATAACACTTGAGCATAATAATAAAATAGTAATAATACATACACCAAATAAT
GAAATTAAAGGTACGCTTAACTCTAATAAATAAAAGCCTAATACACCGCCACCAAAATTTGGGAAATGAC
TGTGTTGGTATGACTGATACACATAAGAAAGTACAGGTTCTCTTTCAGCTTTGATACCACTATTAAAATG
AAAAACTAACTGTGATACAAATAGCAATGCAATTTGCAATACAATCGAACCAGCCGTTCGTCTAGTTTTA
GGAATACGTTTAGAGTATGTAATAAAACCAGTTGCTAAGAGTACTAAAATATATGTTAAATATCTACTGT
ACCCAAATAAATAATTAAAGAAGCTGTCAATTAGACGACCTATTATTCCTAATTGGAAAACACCCAACAC
CATTAATACAACTACTAAAATAGCTATGACATAACGTATCGGATTATCATTCTTCTTTTTCCTCGAATTT
GTTCTTTTTTTTGATGTTGTTTTTTTCTTAGCTGTCGATTTCTTTTTTGCTTGTGCCAATTGACGCACCT
CCTTATGTTGCCTTAGTTAAAACGACGTTATTATATATAATTTCCTATTTGATGTTTAGAATATGTAGTC
CAAACAGACGATGAAAGAGGTTAGAACAATATCATCTATGTATCAATTTACTTAAGATAACAGATTAATC
GTTACCCCTTCAATTGACTATCATAAAATATAAGTAAATTGTCCATGACAACCATTATATTTCAGTGACT
TGTTCTGACCTCTTTTATGTTAATGACTTTTTGATTAAATTTCAGAAATTACTGGAATAATCATAGGACG
ACGTTTTGTACTTTCGAATAATAGTTTACTAATTTGATCACGCATATTTTGTTTAATTTCAGACCATTCT
ATGCGTTTTTCTTGTAAACCAGCCTCTACTATTTCACGTACTTTCTCTTCTGCTTCACGTAATAAGTCTT
CACTTTCACGTACATATACAAACCCACGAGATTGAATTTCAGGTCCCGCAGCTATACGTCTATTTTTAGG
ATCTAACGTTACAACAGCAATAAAGATACCATCTTCTGCTAAAAGATGACGGTCTCTCAACACGATATTT
CCTACATCCCCAATACCAATGCCATCTATTAAAATATTTCCTGAATTTACCTTTTCATTTAATATCATAT
CTTTACCGTTGTAATTAATGACATCTCCTTTTTCCACAAGGAAAATCTTTTCTGGTGCAACACCTGCTTC
AGCTGCTAGCTTCGCATGTGCTATCTGCATTTTAAATTCACCTTGTACAGGAATAAAGTATTCAGGTTTC
ATAATATTAATCATCATTTTTAATTCTTCCATGCAACCATGACTTGAAGCATGAATCTTTTTGTTATTTG
GAATAATATGTGCGCCAGCACGTACAAGCTCATTTAATGTATTCGCAATGATAACTTCCATATTAGCAGA
AGCCGTAATTGCTAAAAATACAGAATCGCCTTCTTCGATATTCATAATTTTATGCTTATGTTGCGCCATT
TGACTTAAGGCTTCTACAGGTTCTCCTTGCATACCAGTAGCTATAATTATCACTTCATTTTTAGGATAAT
TATCAACTTCTGTTATAGGAATTAGCAAATCTTTAGGAATGTCGAAATACCCCATTTTACGAGCAATATT
AAATGAACTTTCAAGTGATCTTCCTAAAAATGACACTTTACGATTTAGCTTGCTAGCAATATTTAAAACT
TGCTGAATACGTATAAAGTTCGAAGCATAACATGAAACTATCAAGCGACCTCGCACTTTTGCAAAAGCAT
CATACATATGATGTTCAATCACATTTTCCGGAGTATTATATCCAGGTTTCTCTGCCTCAGTAGAATCACT
GATTAAGACAAATACGCCTTCTTCACCAATCTCTGCCATACGTTTAATATCTGGTGCATAATGTCCATGT
AAACTTTGGTCAAACTTAAATTCACCTGTATACACAATGGCACCATATGAAGTGTGAATACAAACACCTA
AACTATCAGGAATACTGTGTGTCGTATTAAAGAAACTAATATTCACGTTTTTGAATCTCATAATTGAATC
ATTATTAACTGTATAGTAGCGAACTTTTTTATCAATATTACGGGCTTTCATATTTTCTTTAATTAACGCT
ATTGTCAATTTAGATCCATATACTGGTGCATCTAATTGTTCTAAAACATAACTCACTGCACCAATCGCGT
GCTCATGTCCGTGTGTAAGGAATATACCCTTCAATTTATCTTTATTTTCAAGTACGTATGAAATGTCTGG
TATAACAATATCAATACCTAGCATTTCGTCTTCTGGAAACATAAGTCCAGCATCTAACATAAACATTTCA
TCGTCTACTTCAACGATATACATATTTTTAGCAATTTCGCCAACACCGCCTAATGGTATAATGCGAATAT
CTTTATTCTTTTTCTTTATTAAACTCAAAATCTTACCTCCTAAATTTGTTACCCGTCCATATATAAACTC
ACTATTCATTATAAGTTAAAAAACACTTGATGTACACTATTAATCACTATACATACAGGATAATTAACAG
TTGCACACAAAATAATAAAGCGTTTTAAAATGTTCAATACAACAAATTCAAAACGGGTCACTTCGTAAAA
AGTGGCATAAAAAAATCACATGACGATACATTGTCGTCATATGACTTTAAATATTATTCTTCTAATGCTC
TATGTGAAGCATTTACTCGACCTTGTTTATCAATTTCAGTAATCTTAACTTCAATTGTGTCACCGATTTT
TAATACATCTTCCACTTTTTCAATTCTATTTTTTGAAATTTGTGAAATGTGAAGCAACGCATCTTTACCT
GGGAATAGGCCTACAAACGCACCGTATTTTTCAATACGTTTAACAGTGGCTTGATAAGTTTGACCTACTT
CCGCTTCACGTGTAATTTCCTCAATGATTTCACGAGCACGATTTATCATAGCTTGATCAACAGCACCAAT
AAAGATAGTACCATCTTGTTCAATATCTAATTTAACACCTGTTTCATCAATAATTTCGTTAATTTTTTTA
CCACCAGGTCCGATAACATCTCTAATCTTATCTGGTTTAATTGTCATAGTTACAACTTTTGGCGCGTAAG
CACTTAATTCAGTACGTGGTTGATCAATTGTTTGTAACATATGATTCATTATTTCTAAACGACCACGTCT
CGCTTGTTCTAGAGCCTCTTCGATAATTTCACGCGTTAAACCGTCAATTTTAATATCCATTTGGATTGCT
GTAATACCTTCTTTAGTACCAGCGACTTTAAAGTCCATATCACCTAATGCATCTTCCATACCTTGGATAT
CAGTTAAAATCGTATAGCTATCTTCACGTGTAACAAGGCCCATAGCAATACCAGCAACTGGTGCTTTAAT
CGGTACGCCCGCATCCATTAATGCTAATGTTGATCCACAAATTGACGCTTGAGATGATGAACCATTTGAT
TCAAGTACCTCACTTACAATACGAATTGTATATGGGAAATCAGCAGTATCAGGAATAATATATTTTAATG
CTCTTTCACCTAACGCACCATGTCCAATTTCACGACGACCTGGCGCACGTACTGGACCAGTTTCACCTAC
TGAAAAATTCGGGAAGTTGTAATGATGCATGAATCTTTTTTCTTCTTCAGGTCCTAAACCATCAATTAAT
TGATAATCGCCTAAAGCACCTAATGTTAAAACTGAAAGTGCTTGAGTCTGACCACGTGTAAATAGACCTG
AACCATGCGTTCTAGGTAAAATACCAACTTCAGAATCTAATGGACGGATTTCATCAGGTTTACGGCCGTC
TGGTCTAATTTTTTCATCTGCAATTAAACGTCGAACTTCTTCTTTCACTAATTCATTTAAAATTGCATAA
ACTTCTTTAATAAGTAATTCATTCTCTGGATCTTCTTCATCGATAAATTCATTGACGATTTCTTCTTTTA
AGTTATCAAGATTTTCATCTCGTTGTTGTTTATCAAATGTTAAAACTGTTTCTTTAAGTCCTTTTTCTTC
GGTTAAAGATTTTACACGTTCAACTAGCGCTTCATCACGCTCTGCTGGAATAAATTCTTGTTTAACAGGT
TGAATGTGGTCGACGATTTGTTGTTGGAAATCAACTAAACGTTGAATCTCTTCATGACCAAAGAAAATCG
CCTCTAACATTTCTTGTTCAGTAATCTCACTAGCGCCTGCCTCTACCATGTTTACCGCATCTTTATGACC
AGCTACTTCAAGGTCTAAACGAGAAACTTCTTTTTCTTCTACTGTTGGGTTAATGATATATTTACCGTCA
ATATAACCCACATTTACACCGGCGATTGGCCCTTGGAATGGAATATCCGACACACTAAGCGCCATAGATG
AACCAATCATTGCAGCCATTTGTGGTGAACAATCAGGATCTGCACTTAATACCATGTTCATAATTTGAAC
ATCATGCTTATATCCTTTAGGGAATAAAGGTCTAATTGGTCTATCAATTAATCGCGCAGTTAATGTTGCA
TCGTCACCAGGACGTCCTTCTCTCTTTTTAAATCCACCAGGAATTTTACCCGCAGCGTACATTTTTTCTT
CATAGTTCACTGTTAATGGGAAGAAATCTCCATCACGAGGTTCTTTTGATGCAGTTGCCGTCGATAACAC
GACTGTATCTCCATAACGAACCAATACAGCGCCATTTGCTTGTTTAGCTAATTGCCCTGTTTCAATCGTT
AAAGATCTTCCTGCCCATTCAGTTTTAAAAACTTTCTTTTCTTGAGACATTATAAATCTCCTCTCTTGAA
TCTAATAGTTATATCATATCATTTATTATCGTATTTTAATATTAGTTATAAAATATTGAATTGGAAAAAA
AGAGAGCGGGATATCCATCTATTATTTACAAAATAGACTTATCGTCCCACCCTTATCAGTATAAATATAA
GTTATATTTTTTAAAAAAACATACTTTATTTTTAATATCTTGATTATTACTCTAGGTGTACTTATATCAT
CTTTATCTCATTGTAACCTTATTATTCACTCAATATTAGCTTAACATATAATAAATGTTAATATTTGTAA
AAAGCTTTTTTAATATAAATTAAGGTTGGAACATAAATATGCCCCAACCTCAAAGACGTTATATTAAGAT
TAACGACGGATACCAAGTGATTTAATTAATTCACGGTAACGTTGAATATCTTTACTACGTAAGTAGTTTA
ATAAATGTCTACGACGACCTACCATTTTTAATAATCCACGACGTGAATGGTGGTCTTTTTTGTGTGTACG
TAAGTGTTCGTTTACTGCGTTGATTTCTGCAGTAAGTACAGCGATTTGTACTTCTGGTGAACCAGTATCA
GTTTCGTGTACACGGTATTCTTTAATGATTTCGTTTTTACGTTCTTGTGAAATTGCCATAATCAATTTCC
TCCTTAAAATGTTTGAATACTCCTTAATCCGAGTTAGGCGTTGGAGATTCGACCTACTAAGAAAAAGGTG
TGGTATAGATTGTTTATATAACTATACGCAACTCTAGATATTATAAGCTACTTCATCACCAAAATCAACT
GCTAATAAATATTTAGCACGCGATTTATCATCGTTCATTTGTTTAACTAATGGGTCGATACCATCAAATT
TAATCTCAGGACGTAAGAAATGATGCCAATTCACGGTCACTCGTTCACCATAAATATTATCCTCAAAGTC
AAAGATATTCACTTCGATGACAACTTCTGCTTTGTTAGGATCATGAAATGTTGGCTTTACACCTATGTTA
GCTACCCCTCGATATAATTTATTTTCAGTGCCGATTTCAATACTAACAGCATAAACACCTTTACGAGGTA
ACAAATAATCATCACTAGGTTGAATGTTAGCTGTTGGGAAGCCAATAGTTCTTCCCCTTTTTTCACCTTG
CACTACAGTGCCTTTAATAGAATATATATAGCCTAAAGCATCATTCGCTTTTTGCAACTCACCATTGATT
AAATCTTGACGAATAGAAGTTGTAGAAATTTTTTCATTTTCAATTTCTTGTTTACTCACAATTGTCGTAT
TAAACGCATCATATTCTTGAAGTACAGTCATATTACCTTTTCCAAATTTACCAAAAGTAAAATCAAAACC
AGCAATGACTTCTTTTACATTATTTTTAATTATATAATTTTCAACAAAATCTTCTACGCTCACATTAGCA
AACCTAGATGAAAAATTAACCACTATACAATAATCAATATCATGTTGGCTAATTTTTTCGATTTTATCTG
AAAGTGGCGTTAAATACGTTGTTCGTTTTCTTTTAGGATTCAACACGACAGACGGATGCGGATCAAATGT
CATCACCGCTTTTTTTAAACTGCGTGCCTCAGCTATTTCGTTTAATATATCAAAGACTTTGTCATGACCT
TTATGCATGCCATCGAAAAATCCGAATGCCATTGCAACATCCTCTGTAATATACTGTTTAGATTGTATAG
GATGTGTCACTTCTATGACTTTCATAAATTCTATCTCCTTTAATTAAAGACTTTTTTAGGTTTAATTTCT
GATTCTTTTGTAGGGTGTACCATATAAATTGCTAATACTTTTTCTGAATCATCATCAATAAATACAATTT
GGTCTTTAATTTTGTTATCAAATTCATTTTTATTAAATTTCTGCCCATTTAAAATACGTTTTTTTATGTG
CGAATCTTTAATTTTAATGCTTGGCAAACCCTTTAATCCATATTCTAAAGGAAACAATTTATTTTGCAAT
GAATCCTGCTCATGAAGTTCTTTTATTTGTTCTAATGTAAGGCTATCTTTCAACACAAATCCACCAGACT
CGATTCGTGTTAATTTCGACATATGTGCCGGAAAGCCTAA
>contig_04 len=12000 cov=124x
AGTATTATAAATTATCTTTAAATAAATATATGTGAAGTAAAAATTACACGTTAGCATATCGATTATGTCA
TTTCTTTTAACATATTAACTAGGGAAACGTTAAAAGTTAACGGTTGATATCTAACTAAAAACAAGGTCAC
AGTAGTATGTTTTAATCTGGCGTCTATTACAAATAAAAATTACATCTATAATTATTCGTTTTCTTTTTTG
AAAGTAATAGCCAATTAATATCATACATACTGGAGTGACTATAAGGAGGACATTATTATGAGAGCAGCAG
TTGTAACGAAAGATCACAAAGTAAGTATTGAGGACAAAAAGTTAAGAGCTTTAAAACCTGGTGAAGCGTT
GGTACAAACGGAATATTGTGGCGTTTGTCATACCGATTTACATGTTAAGAATGCTGATTTTGGTGATGTT
ACAGGCGTTACTTTAGGTCATGAAGGTATTGGTAAAGTCATCGAAGTTGCGGAAGATGTAGAATCATTAA
AAATTGGAGACCGTGTGTCTATCGCTTGGATGTTCGAAAGCTGTGGAAGATGTGAATATTGTACAACAGG
TCGTGAAACACTTTGCCGTAGTGTGAAAAATGCTGGTTATACAGTAGATGGTGCAATGGCTGAACAAGTT
ATTGTTACTGCAGACTATGCTGTGAAAGTACCTGAAAAATTAGATCCAGCAGCAGCGTCTTCTATTACAT
GCGCAGGTGTGACAACTTATAAAGCTGTAAAAGTAAGTAATGTAAAACCTGGACAATGGTTAGGTGTTTT
TGGTATAGGTGGTTTAGGTAACCTAGCTTTACAATATGCTAAAAACGTTATGGGGGCTAAAATTGTTGCC
TTCGACATCAATGATGATAAATTAGCATTCGCGAAAGAATTAGGTGCTGATGCTATTATTAATTCTAAAG
ATGTTGATCCAGTTGCAGAAGTTATGAAATTAACTGATAACAAAGGATTAGATGCAACAGTGGTAACTTC
AGTTGCTAAGACGCCATTTAACCAAGCGGTTGATGTTGTAAAAGCTGGTGCAAGAGTTGTTGCCGTTGGT
TTACCTGTTGATAAAATGAACTTAGATATCCCAAGATTAGTGCTTGATGGTATTGAAGTAGTAGGTTCAC
TTGTTGGTACAAGACAAGACTTACGTGAAGCGTTTGAATTTGCTGCTGAAAATAAAGTAACACCTAAAGT
TCAATTAAGAAAATTAGAAGAAATCAATGATATTTTTGAAGAAATGGAAAATGGTACTATAACTGGTAGA
ATGGTTATTAAATTTTAAAAATATCAACTGACTATATAGATAAAGAAGGTAGTGCTCTGAACACTATCAT
TATAATCAAACACGAGGTTTTCATGAAAGATAGTGAAAATCTCGTGTTTTTTGGTTTTGAGGTGTTGTTT
GTATTTTATAAAATGGCTTACATATATGAAGCGTTGATTAAGTATGGAATTGTTAATTAATTGAACCTAT
TTAGCTTTAAGAAGGCATAACAAGATGACCTTATTTTATGCTATAATATTTCTATTATGCGAAGATTAAG
GTGAGTAGTAAATTGGATAAAAAAGTAAGTATTCAAACAAAGCAAGTGTTGAAACAGCACAACGAAAAAG
AAAAATTTGAATTTACTACTGAAGGAACTTGGCAACAAAGGCAATCTAACTTTATTCGGTATGTAGAACA
AATTGAGGATGCAACAGTTAATGTTACAATAAAAGTGGATGATGATAGCGTTAAGTTGATTCGTAAAGGC
GACATTAATATGAATTTGCATTTTGTTGAAGGACAAACGACAACAACTTTTTACGATATATCGGCTGGAC
GAATTCCACTAGAAGTTAAAACATTACGCATTTTACATTTCGTAAGTGGAGACGGTGGCAAGCTAAAGAT
TCATTATGAATTATATCAAGATAATGAAAAAATGGGTTCTTATCAATATGAAATTAACTATAAGGAGATA
GGCGAATGAATATTATTGATCAAGTGAAACAAACATTAGTAGAAGAAATTGCAGCAAGTATTAACAAAGC
AGGATTAGCAGATGAGATTCCTGATATTAAAATTGAAGTTCCTAAAGATACAAAAAATGGAGATTATGCT
ACTAATATTGCGATGGTACTGACTAAGATTGCAAAGCGTAATCCTCGTGAAATTGCTCAAGCGATTGTTG
ATAACTTAGATACTGAAAAAGCACATGTAAAACAAATTGACATTGCTGGTCCAGGATTCATTAATTTTTA
CTTAGATAATCAGTATTTAACAGCAATTATTCCTGAAGCAATTGAAAAAGGTGATCAATTTGGACATGTA
AATGAATCAAAAGGTCAAAATGTATTGCTTGAGTATGTTTCAGCTAACCCTACAGGAGATTTACATATTG
GTCATGCTAGAAATGCAGCAGTTGGTGATGCTTTAGCTAATATTTTAACTGCAGCTGGCTATAATGTAAC
ACGTGAATATTATATTAATGATGCTGGTAATCAAATTACTAACTTAGCGCGTTCGATTGAAACACGTTTC
TTTGAAGCTTTAGGTGACAATAGTTATTCAATGCCAGAAGATGGCTATAATGGAAAAGATATTATTGAAA
TAGGTAAAGATTTAGCAGAGAAACACCCTGAAATTAAAGATTATTCTGAAGAAGCACGTTTGAAAGAATT
TAGAAAATTAGGCGTAGAATACGAAATGGCTAAATTGAAAAATGATTTAGCAGAGTTCAATACGCATTTT
GATAATTGGTTTAGTGAAACATCTTTATATGAAAAAGGAGAAATTCTTGAAGTTTTAGCAAAAATGAAAG
AATTAGGTTATACGTATGAAGCTGATGGCGCTACATGGTTACGTACAACTGATTTTAAAGACGACAAAGA
CAGAGTATTAATTAAAAATGACGGTACATATACGTATTTCTTACCAGATATTGCGTACCACTTCGATAAA
GTAAAACGTGGTAATGACATTTTAATCGATTTATTTGGTGCTGATCATCATGGTTATATTAATCGTTTGA
AAGCATCTCTTGAAACGTTTGGTGTAGATAGTAATCGTTTAGAAATTCAAATCATGCAAATGGTTCGTTT
AATGGAAAATGGTAAAGAAGTGAAGATGAGTAAACGTACTGGTAATGCGATTACATTAAGAGAAATTATG
GACGAAGTTGGCGTTGACGCTGCACGTTATTTCTTAACTATGCGTAGTCCTGATAGTCACTTTGATTTTG
ATATGGAATTAGCGAAAGAGCAATCTCAAGACAATCCAGTTTACTATGCTCAATATGCACATGCGCGTAT
TTGTTCAATTTTAAAACAAGCGAAAGAGCAAGGTATTGAAGTGACTGCTGCGAATGATTTTACAACGATT
ACTAATGAAAAAGCGATTGAATTGTTGAAAAAAGTAGCTGATTTCGAACCTACAATTGAAAGTGCTGCTG
AGCATAGATCGGCACATAGAATTACTAATTATATTCAAGATTTAGCTTCTCATTTCCATAAATTCTATAA
TGCTGAAAAAGTGTTAACAGATGATATTGAAAAAACAAAAGCACATGTTGCTATGATTGAAGCGGTCAGA
ATTACATTGAAAAATGCATTGGCAATGGTCGGTGTAAGCGCACCTGAATCAATGTAAGAACATTTATATA
CACTCCAACGTAGAGTTTCTCGAAAGATACTTTGTGTTGGAGTGTTTTTTTTAGGTATGTGACATATTGG
GGAATGCTTAGTATGTGAATAAGGTTAAGAGGAACACAGTTGGATGCTCTGCACAACTGCATAAGAGAGC
CTGAGACATAAATCAATGTTCTATGCTCTACAAAGTTATAATGGCAGTAGTTGACTGAACGAAAATTCGC
TTGTAACAAGCTTTTTTCAATTCTAGTCAACCTTGCCGGCGGGGCCCCAACAAAGAGAAATTGGATTCCC
AATTTCTACAGACAATGCAAGTTGGGGTGGGACGACGAAATAAATTTTACGATAATATCATTTCTGTCCC
ACTCCCTCTAAAATGGAGGGTGTAAATGTTAGGAACTGATGAATTATATAAAGTTTTATATGAACATCTC
GGACCACAATTTTGGTGGCCTGCTGATAATGACATTGAAATGATGTTAGGTGCAATTTTAGTTCAAAATA
CTAGATGGCGAAATGCAGAAATTGCATTGAATCAGATTAAAGAACATACGCATTTTAATCCAAATCATAT
ATTAGAACTACCTATTGAAACGTTACAATCATTGATACATTCAAGTGGCTTTTATAAAAGTAAATCACTG
ACGATTAAAACATTATTAACATGGTTAGCACGACATCATTTCAATTATCAAGAGATTAATGAGCGATATA
AAGGTGGATTAAGAAAAGAATTATTATCTTTGAAAGGTATTGGAAGTGAAACAGCAGATGTCTTACTTGT
TTATATATTCGGACGTATTGAATTTATTCCAGATAGCTATACAAGAAAAATATATGATAAATTAGGATAT
GAAAACACTAAAAATTATGATCAATTAAAAAAAGTAGTCACATTACCAAATCATTTTACAAATCAAGATG
CTAATGAATTTCATGCTCTGTTAGATGTATTTGGTAAACATTACTTTAGAGACAAAGATATAAAGAATTA
TGATTTTTTAGAACCTTACTTTAAAAAGTAAACGCTGTGAAGTTAGATAGATGAGTTTATATGAAATATA
AAAAATAATTTACTATTTTCTTTTAGTATGTGGACTTATATAATAAATAGAAGCATATAAAGAAAAAAAC
AGTTGTTTGTTTGTGCAGCAACTGCATAAGAGCCCCTAATCGCTAAAGCTCAAGGGGAGTAAAGGAATAC
AGTTGTTTGTGCAGCAACTGCATAAAAGCCTCTAATCACTAAAGGTGAAGAGGAACGCAGTTGGATGCTA
AGGCACAACTGCATAAAAGCCTCTAATCGCTAAAGATGAAGAGGAACGCAGTTGGATGCTACCGCACAAC
TGCATAAATCCCTCTAATCGCTAAAGCGAAAAGTGGGATTAAAAAGGAGATGTGATAGTGTGAAGAAATC
GTTAATTGCTTTTATTTTGATTTTTATGCTTGTCCTGAGTGGCTGTGGTATGAAAGATAATGATAAACAA
GGTAGCAATGATAATGGCTCGTCTAAATCGCCGTACCATAGAATTGTTTCGTTAATGCCTAGTAATACTG
AAATTTTATATGAATTAGGATTAGGTAAATACATAGTTGGTGTTTCAACGGTTGATGATTATCCAAAAGA
TGTGAAAAAGGGTAAGAAACAATTTGATGCTTTGAATCTAAATAAAGAGGAACTTTTAAAGGCAAAGCCA
GATCTAATTCTTGCGCATGAGTCGCAAAAGGCAACTGCTAATAAAGTATTGTCATCATTAGAGAAACAAG
GCATCAAAGTAGTGTATGTTAAAGATGCACAATCAATTGATGAAACTTACAACACATTTAAGCAAATTGG
GAAATTAACGCATCATGATAAGCAGGCTGAACAACTTGTTGAGGAAACTAAAGATAATATCGATAAAGTC
ATAGATTCAATTCCTGCTCATCATAAAAAATCAAAAGTATTTATTGAGGTTTCATCAAAGCCTGAAATAT
ATACAGCAGGGAAGCATACATTTTTTAATGATATGTTAGAAAAATTAGAAGCCCAAAATGTGTATAGTGA
CATTAATGGTTGGAACCCTGTAACGAAGGAAAGTATTATTAAAAAGAACCCAGATATATTAATTTCGACG
GAAGCTAAGACAAGATCAGATTATATGGATATCATCAAAAAAAGAGGTGGATTCAATAAAATTAATGCTG
TCAAGAATACACGTATTGAAGTTGTAAATGGTGATGAAGTATCAAGACCAGGTCCACGTATTGATGAAGG
ATTAAAAGAATTAAGAGATGCAATTTATAGAAAATAAACCATTCTAATTATGCCCCTTATTGCTACATGT
AAAAAATACATGTTTGAGATAAGGGGTTTTTAAAATATATTTAGTGAATGATAGCAACGCGAGTATGTGA
TTGCTATAATGAATGTAATTATCGATGAAAGAAAAGAGAATGCTATGACATTTAATAAAGTATTATTGAG
CTGGATAGTCATATTGATTATAACAACTAGCATATATCTATTTTGGCAGTTGGGCGATATCAATGATGTA
TTTAACCAGTCTATTTTAATCAATGTTAGATTACCGAGATTATTAGAAGCATTGTTGACAGGTATGATAT
TAACTGTTGCAGGCCTTATATTTCAAACAGTTTTAAATAATGCATTGGCAGATAGCTTTACATTAGGATT
GGCAAGCGGCGCTACATTTGGTTCAGGATTAGCATTATTTTTAGGTTTAACAACGTTATGGATTCCTGTA
TTTTCAATAACATTTAGTTTGATAACATTAATAACTGTATTAGTCATTACGTCGGTATTGAGCCAAGGCT
ATCCAGTTAGAATCTTAATATTAAGTGGTTTAATGATTGGTGCGTTATTCAATTCACTTCTATATTTTTT
GATTTTATTAAAACCTCGCAAATTAAATACAATTGCCAATTATCTGTTTGGTGGTTTTGGTGATGCAGAA
TACTCAAATGTATCTATAATAGCAATCACATTTATCATTGCATTGTTTGGTATATTTATCATTCTTAATC
AACTAAAGTTATTGCAATTAGGAGAACTAAAAAGTCAGTCACTAGGCTTAAATGTTCAATTGATTACATA
TATCGCGTTATGTATAGCTTCTATGATAACGGCGATAAATGTCGCATATGTTGGCATCATTGGATTCATT
GGTATGGTGATACCGCAACTCATTAGAAAATGGCAGTGGAAACAATCATTAGGAAGACAATTGGCTTTAA
ATATTGTAACTGGAGGACAAATAATGGTTATGGCAGATTTTATTGGTAGCCATATATTGTCACCAGTACA
AATACCGGCAAGTATTATCATTGCATTAATTGGTATACCAGTGTTATTTTACATGCTAATATCTCAGTCG
AAACGGTTACACTAGCACACGACATTTGCTAAAATAAAAATAACTATAAACATAAAGAGGGCATAAGCGA
TGGATTTGAATCAAATTAAAGCAGTTGTATTTGATTTAGAAGGTACGTTGTTGGACAGAGTTAAATCTCG
AGAGAAATTTATCGAAGAGCAATATGAACGATTTCATGACTACTTAATTCATGTTCAACTGGCAGATTTT
AAAAAAGCATTTATTGAGCTAGATGACGATGAAGATAATGATAAACCTGATTTATATAAAGAAATCATTA
AACGTTTCCATGTAGATAGGTTAACTTGGAAAGACTTATTTAATGATTTTGAAATGCATTTTTATCGTTA
TGTATTTCCTTATTACGATACTTTGTATACACTAGAAAAGCTATCGCAAAAAGGCTTTCAAATTGGTGTT
ATCGCAAATGGTAAATCTAAGATTAAACAATTTCGATTACATTCACTTGGTTTGATGCATGTTATTAATT
ATTTATCAACATCAGAAACAGTTGGTTTTCGTAAACCACATCCTAAAATTTTTGAAGATATGATTGATCA
ACTAGGGGTATTACCTGAGCAAATTATGTATGTTGGCGATGATGCGTTAAATGATGTAGCTCCAGCACGA
GCTATGGGCATGGTTAGTGTATGGTATAAACAAGAAGATGCTGAAATTGAACCACTCGAAGAAGAAGTTG
ATTTTACAATTACAACAGTGGAAGAATTATTAACCATTTTACCAATAAAAAATGATAATAAAGGAGAAAA
TTATGGATCTATTTACTAGAAAAGATGGAACATCGATACATTACAGTACATTAGGTGAAGGCTATCCTAT
CGTATTGATTCATACTGTACTTGATAATTATTCTGTGTTTAATAAATTAGCAGCAGAATTAGCAAAATCA
TTTCAAGTTGTGTTAATTGATTTACGTGGACATGGCTATTCTGATAAACCTCGTCACATTGAAATAAAAG
ATTTTTCTGATGACATTGTTGAATTACTTAAATATTTATATATTGAAGAAGTTGCATTTGTATGCCATGA
AATGGGTGGAATCATTGGTGCGGATATTTCAGTACGTTATCCTGAATTTACATCATCACTTACGTTGGTA
AATCCAACATCTATTGAAGGTGAATTACCGGAAGAACGTTTATTTAGAAAATATGCCCATATTATTCGAA
ACTGGGATCCTGAAAAACAAGATAAATTTTTAAATAAGCGTAAGTATTATCGTCCGAGAAAAATGAATCG
ATTCCTCAAACATGTCGTAGATACAAATGAAATATCAACTAAAGAAGAAATTCAAGCAGTTAAAGAGGTA
TTCAAAAACGCTGATATTTCTCAAACTTATAGAAATGTCGTAGTACCGACAAAAATTATTGCAGGAGAAT
TCGGTGAAAGAACAACAAGATTAGAAGCTAAAGAAGTAGCTGATTTAATCCAAAATGCGGACTTTGAAGT
ATATCAAGAATCAAGTGCATTCCCATTTGTTGAAGAGCAAGAAAGATTCGTCGAAGATACAGCTGCATTT
ATCAACAAACATCACGATGAAAAGCATGTTTAATTTTTAGTGATAAGTGAAATATAAAAGTAGTAATTAT
TATTTGTAATATAATTGTAATATGACTGTTGTTTTAGAAATGATTGTTGTATAATGGTTGTGAGCACAAG
ACGATATAAATAAAAGGAGAGAGTGCTTATGAAAAAGCTACTAACGGCAAGTATAATTGCATGTTCTGTT
GTAATGGGAGTAGGCTTAGTGAACACTAGTGCCGAAGCAGCAAGTGGCAACTCTATTGATACTGTTAAAC
AATTAATTAAGGGTGATCAGTCATTAGAAAATGTGAAAATTGGCGAATCTATTAAAGATGTTTTAACTAA
GTACAAAAATCCGATGTATTCTTACAATGAAGATGGAACTGAACATTATTACGAATTCCATACTAAAAAA
GGTATGTTATTAGTAACTACTGATGGTAAGAAAAACAATGGTAAAGTAACTCATATTTCAATGATGTATA
ATGATGCTAATGGTCCAACATATCAAGCTGTTAAAAATTATGTTGGCAAAGCAGTAACACATACGGAATA
TAGCAAAGTTGCTGGTAATTTTGGATATATTGAAAAAGGCAAAACGACTTATCAATTTGCCTCAGCACCA
AAAGATAAAAACATAAAATTATATCGTATTGATTTGGAAAAATAATAATTAAATTAATCGTGCTTAACCG
ACAAGTTAAATCGCAATACTTGTCGGTTATTTTATTACATGGCATAATACAGTGAATATTAATTTTTATG
ATTTGCGAAGGGGAAATTAAACCAATACAATAAATAAAGATGTTTCATTTTAAAGGAGTACTACCAATGC
GTATTAACGATAAAATTTTACTAGAAAATATAGAAGATTACTTTAATCATAAAGGTTTATCACCGCATTT
GATTGATGATATTAAAGAGAAAGTAATTACTGATATAAAAAATTCTGAAAAGAAAGATCAAGATTATATT
GAATATAAAAGAAAATCTCCAGCACAAATCATATTAATGATTCAAAGGAATTTGTTTGCTTTACAAATGA
ATCCAGTTATTTTCTTTATTATAAACTTCATTCTCATATCCTATTTATATGATAAACAGTATGTTCAGTT
TCAAGCTATTACTGGAATGAGTCTATTTTATTGTTTAGTGATTTTTCCAATGACTATTGTTGTATACCTA
AGGGTGTCTCAAAAGAATTACTTGCGTAGTAATAAAATAGAAATGATTATGGGTACAATTATCGCTATCA
TATCCTTGTTATTAATTATATTACAAGCATTTAATATTACTTGGGGCGTTATACCAATTACAAATTTTGG
ACATCAATTTTTCTTTTTCATTGGTATTATTTTAGTAATTGCCGGCATATTTTATAAGCGACTTGAGTTT
TCGGGAATCGGGTTATTATTTTGTCAAAAAACCGTCGATGCAATGATTCATAATCCACAATCAGCCCAGA
TTTTTTCATTAATTATATGGATATTATTAGTAGTTCTAGTTATATATTTCACAATTAGATTATCTTCACG
TACAAGATTATAAATATGATAAAACTATTCACTTGATTAATTGTATTAATTGAGATGAATAGTTTTTTTA
TTGTTGGAATAACTTTTGGTAATTTATAAATAATTTAAAAAAATTGTTTATAAAATGGAAGCGTATATAG
AATGAAGGTTGGGTATATAGTTTATTGAGGGAGGTGTCACAATGAATAAAGTCACAATTAATCCTCAAAT
CCAATTAACTTATCAAATTGAAGGTAAAGGGGATCCTATAATATTACTTCATGGATTGGATGGTAATTTA
GCTGGATTTGAAGATTTGCAACATCAACTAGCATCATCATATAAAGTACTTACTTACGATTTAAGAGGTC
ATGGCAAGTCTTCTAAAAGTGAATCATACGATTTAAACGATCACGTTGAGGATTTAAAAATTCTAATGGA
GAAGTTAAATATTCATGAGGCACATATTCTAGGACATGATTTAGGTGGGGTAGTTGCTAAGTTATTTACA
GATAAATATGCTTATCGTGTAAAATCATTAACTACCATTGCATCGAAGAAAGATGACTTAATACACAGCT
TTACTCAATTGTTAATACAATATCAAGATGATATAGCGGGTTTTAATAAGTCTGAAGCGTATATTCTTTT
ATTTTCTAAATTGTTTAGAAATCAAGAGAAGACGATGAAATGGTATCAAAAACAAAGAATATATAGCATT
AAGTCTGAGGATGATAGTGCGGTGGCAATTCGTTCATTAATTTTGCATAAAGATGAACCTATGTATTTAA
AAAAACGTACATGTGTACCTACTTTGTTAATTAATGGGGAACATGATCCTTTGATTAAAGATAAGAATCA
TTTTAAATTGGAAGCGCATTTTTTAAATGTTACGAAAAAAATCTTCGAACATTCAGGACATGCACCGCAT
ATTGAAGAACCAGAAGCATTTATGAATTATTATTTAAAATTTTTAAAAAGCGTATCATAATATGTGATAT
ATAAACCTAGGGCATAAAGTCCTTAGGCAATGTGAAAAAGCTGATTACTATTCATTATTTGATAGAAATC
AGCTTTTTTTGAAATGTATTTGATATATACTGCTCGTTATGCGGCTATCTTCCTTATATTAAGTGCCATT
AGTGCAAAACCTCTTAACAATTAGGTAAAAAGAGCATAAAAAAAGGAAGTTTAATAGAATGTATCATCTA
TCAAACTTCACCAAATTGCGCTAAACAAAATTATAGTTCAATTTCGTTGTTTGCTTCAGTGATTCGTTTA
TTTACTCGACTCAATAATGATTCGATTTTTTTACGTTGTTGTGCATTAACAAGAATTAATACAGTTCTTT
CATCATGCTCATTACGTTTTTTATCGAAGTAATCTTCTTGAGATAAAATTTTAACTGCTTTAACAACTTG
TGGTTGTTTGTAGTTTAAATGATTAATAATATCTTTAAGATAGTATTCTTTCTCTTTGTTTTCGCTGATG
TATGTCAATACAGCGAATTCTTCAAAGCTAATTGAAAATTCCTTTTTAATTAAACTTTTTAATTTGTCAG
CATAAGTGACCATTGATAACAACTCAAAGCAATCATTGATTTTTGTAATTGCCATGTTTAAAACCTCCCT
ATTTGATGCATCTTGCTCGATACATTTGCCCGATAATATATTGATATCTAATCTTTATTTATTATAGATA
TGTTAGTCATAATTTTGCATTAAATAAGTTTTATTAAATATATTTAATGCTCTATTATTTAGTTAATTAT
AACTAATTAAAAATGAGAAGTAAACAAAAAAGTGTTTATAAAACAAATTATTCATTAGACACAGTGATTG
TATTTCTGGGTTAGCATTTGGTTTAGTCAAAAATATCAGCATTTTTTATAATTTACCTTAATTTAATCGA
CAATACAAGGTAATTTATTTAATGAAACGA
>contig_05 len=7501 cov=67x
TGATAATTCACATCATCTTCTCTTTCAGTATCTCCATCACGTAATGGAAGTTTCAGATTTGCAGGAACTT
TGCCTCCAAACTTACTAACTTTTCCTTCTTCTTTAGCAGCTTCTATAGCTTGTTCAATGGCTTTTATCGT
ACTTGTATCTGATTTAGGAATGATTAAACTGATTGAATACTTTGCTTCTTGCCCTTCTTGCATACTGTGA
GGTTCAAAAATATGTGCATATGATGCTCTTACTTTTCCTGTAATCACTTTAGTTTTATTTAATACTTTTG
CTTTCATGTTTATATACCGTCCTTTTTAATTTTTATAGTTTGTCAAAATCATCTTCAGCAGATTGCTTTA
TAGCTGGTCGTTTATCAGACTCGGTAGCAAGTGTTAATTTACCTTGCGGCTTTTCTATAAAGCCTTCTGC
AATTTTAGAAAATGCTTTTTTGCCGATTAATTTTTCTAAATTCGTAATGCTAAGTAACTTGGTTTCTGTA
ATATCTTCAGGTTTATAACCCGCTTCAACTAACCTGTCTCGGACTGCTTTTGTATCAGTTATCACTCTCC
TTGAACGTCCCGTGACTAGCTTCCAACCTGGATACGTTTTATCATTCTCTTTCGCTTGTTCTAACGCATA
TTTCTCTACTTCATCAGCCCATTTTTTGATATCAGGCAGTTTATATAAAAGTTCTGCAATCTCTTCATCA
CTCAACAAATGTGGTGGCTTTTGAGGCACATTTTGCATGTATTCTGCACGTGTTCTACATGAATGCTTTA
TCTTACAGAATCTACAATGACTACCTGCTTTAAACTCACCTTCACCGTTATAAGCAAGTCTGGCTAAGGG
TTTAACAAAATCGGTTCCCCATTGAAGTAATCTTGAGATTGGTAACTCTTCAGTAGAAAAGTTATCTATT
CGAGGTTGTATGATAGTCATGCGAACTGTATGAATGTCATACACTAAACTAAGCAGTTCATATGCGCCCA
AGCCATATAATCTAAGTTGAGGATTATCTATAGCTGAAACTTCAATGCCTTTACCGTATTTAAGGTCAAT
AATTTCAAGTACACCACCTGAAAATATAATGACATCACCAGTACCAAAAGATTCAGGGACGTATTTACCT
AAATCCAATTTTGTTTCAAATAAAGCTATTACATCGTCATCTCTACTCAAAGCCTCGTTATACTTTTCTT
CTACATTAGCTACGTACTCTTCAACATATTCGCGCAACTCTTCACTGTAATATTGATTTCGCTTATAATT
TTGAAAAGCTTTATTAAACTCAAACTGTGTTAGGCCTTCATATTTAAGACTGAAATATAACTCACTTAAC
TCATGAGCGAATGTACCTTCTTCAGCAAAAACTGAACTTTTATCTGCAATACCTTCACTTGCCTTAATAC
TCGGTGGACAGTTTAGCCATTGTTTTGCTCCACTTGCACTAAGCTTTGCATGAGCTCTATTTGAGTGATC
TAGCTTCATGCATTTATTCTCGCATTCATAAAATCAACAATTTTTTCATAATGCTCTTCTTTGATAGTAG
ATAGCTTATCCGCACCAAGTTCGTTAAGTTTATTTCTAAATTCTTTCTTATCAGAAGTGTCTGCTTTTTT
AAGGAACTCTTTTCCTACTGATAAAACATAATCTTTAGTCAAATCAGCAGAAGTTTCCTTAACTTCTTCA
ATTGATTCCAGTTGAGCTGTTTCTTCTTTTGGCATTGGTGCTTCTTTAACTTTCTCTTGTACAATTGATG
AATCTACAGTTGATAGTTCAGTGTTTAACACACGTAAGTTCTTATTTAATAGTTTTAATTCTTCAAAAAT
ACCTTCTAATACTGCCATTGATTAACTCCTCCTTAAAATTGGTTAGCTAGACGAATCATTAACTTGATAC
GTTCTTCTATTTCTCTAGGGTCATCACTTTGTTCATTTAATCTTGCTAACAATTCAAATTGTTCTTCTAA
AATTTCTTTTTTACGTTCGACTACACTTAAATGCAATTGCGATTCAATAACACGCCAGACACCCCAACTT
TCCATTGCAATCTTTCCTTTTTTCTTAAGTTTTGAAAGTGTGGATTTTGCATGTGTTTTAGATATCCCAA
AAGTTTCAACAACATCATCAGAATTGAAATTGTCATATGTTGCAAAATGTGATAGTATTTTTTGTTGTAA
GGTCATATTAATAACTCCTTATATAATTATTTAAGACAATTGCTCATCTTGCACTGTTACTTGCTCCAAC
AAGTAGCAGTTTTTTTATTCTTCATAAAAGTATTCCTTATAAAATATGAATGTCACTATACTTGCGAATC
CTGCAATTGACCATGCTGTAGTGAAGTACAGCAATGGCATGAGTACAATCGCTAAGACCGTGAAGCACAA
CACTGCTAATAGGTAGCTTTTATATGTGTCGCTCATTTGATAATCCTCCTAATACCATTTTTTATGCTTT
CTGATCAAATACTCTTCTAATTTAGAAATATTAATCAATGTTCCCGTTGCTGAATAATCAATGTATAAAT
TTTCTACACCTAAATTATCTTCACGGTAATATTTCAACCAGTTGTATACTGTACTTCTACATACTCCAAA
CAATTGATGGATTTGTGTAGGTGTTGCGTATAACTTTTTCACAAATTTTTCTTCGCCTCGATATGTGTTT
TCTGGTGTTGGTGGTATTATGATTTTTGGCATCTCTATCACTCCTTTAGATAAATGTTAAAGTTTGTTAT
TATTCGCCCTGTATTGAAGTTCTCTATCTAATGCATAGAAAACTTTGTTTATTTCTAAGTAGCTGTAATC
ACTTTTTTTAATAAGTTCTAATATTTCCGCTCCTAAGTTACGTTCCTTTTCCGTTAAATAGGATGAAGAA
GCATCAGCTTTGCTAGAAACTTGTGGGACGCCTATACGCAATCCTTCTGATCTTGTGTTCATTTGTTTAT
GCTCCTTTCGTGTATAATGTTGTTATCAACCTAAGGAGGTGATAACATGCCCTTGACATCTGATGAATTT
GATACACTTACTAAAGACCAACAATATATCTTGTCCGTACTCTACAAAGATTATTTAGAATGTGTAAAGT
TAGGTTCGATTAAATTAACCTGCAATAATTTTGGAAGTGCTAAAGATATACATACAAAGTATTTTCAAAA
ACTACATTTCGAAGATGTAAAATACGATTTAAATAAACTTAAAAACTCTGGGTTCCTAAACGGCGTGTAT
GCTAGTAACACTATTTATCATGTAACAATTTCAGACAAGACTGTTGTTTACTTTGAAAATGAGTTTAAAA
ACAATTTAAAAAGTATCATTGATAGCATTTCTATAATTGCTTCAATAATTCCTGGTCTCTAGTTGGGTTT
ATAACTTCCCAATCATTTGCCATGAGGTCATCGGCTGAAGGTTGCCAATATCTGATAAGGTTTGTCCCAT
CGCTATTTGAAATGATGCATTGTAAAAAACTATCATTTGTTGGTAATATCTTAGTTCGATGACTTTCTTT
CCAATCTTTCCGTGTCATAGAGACAAGATTTTTTGTAGCTATCTTAGTTGCTTCTTGAATGTTCATTTGT
TATTCCTCCTTTTAAGATGTTTGTTTTTCTTTAAATGCTAAAATAATTGATTTCTTTTTATCATTCGTGA
ATACGAAATTTTCGTATTCATTACCTAAAAAAATATCATCATATTTAACATTAAAAGCACTCATATACTT
AGAAAGTAAACTATCTTTAATGTTTGTAGAGTCTTTTTCCATATTTTGAATTGTACGTGATGAGACCTTA
AATAAATCTCCTAACTCTTTTTGAGTCAATCCGTAATCAGTCCTCAACTCTTTTAATGTTTTCATGTTGT
CACCGCCTTTCGTAAACCTAATATAATACGAAATTTTCGTATTGTCAACATTAAATACGTTTTTTTCGTA
AAAAACTTTACTATGATATGAAAATTTCGTATAATAAGAAAAAAGGAGGTAAGTAATATGAACAAAGAAA
GAAATATTATATACCAGAAATGGTGTACGGGGTTGAGAAATCGTATTTTTTCTTGCGCTCATCAGAAATG
GTGACCTGATTAATCACCACATCAATACGTTTAGAGTCCAGCGACGCCAGCATACCGTCCCATTTGGTCG
GTTTTAGTGACGCCTCAACGCCAAGATGTTTTGCCAGCTGTTGGGCAAATTCCACTTCAAAACCGGTTAA
TTTGCCGTCATCTCCCTGAAAACTGAACGGCGGATAAGTTCCTTCCAGCCCTACCAGCAGCGTGCCGCGC
TCTTTAACTTTATTAAGCAGACCTTCATCTGCAAAACTTTTAACGCTCATGCCCGCAACCAGCGCCACGG
CCATCACACCCATCAATGCCTGACGTCCCAGATGTGCTAATTTCATATTCACCCCGAATGTTGTTATGTC
TGTTTGCAGTGTAGAGCCATTGTTTGTAAACACAAAAACAACTCCGCTACATCTTATTCTTATTTAATAT
ATATCAGAAGAAGGCAGGCTGGAGGATTTCTGCACAAGTCTGGTACAGGTTTGGCGCTGATAGTGTTGAT
ACTTCCGTAACGCGATGCGGTAATTGTTGGTGCTGGTGGTGGGCAGCCACGGGGCGAGGATTTCATCAAG
GAAACCGTCCTGCAGCTGATCGAGGGAAATATTTTGCTCGTGTAGATGATTCCCCAAACGGCGCAAACGA
ACGACATATTCACGCACCGTACCTGGGCTCATTTCAGTTTGCTCAAAAAGAAATTGCTTAAAGCCGATAA
TATCGAAAAAATCACTCTGCGTTTTGCAATGTAAATCACCGCAAAATCGGCACAATGCCGCCCATGATTT
TTGTTCCGTTTGCCAGCCATTTTCATCAAGCAGCGTGTCCAGGCGGGAAATCTCGATTTTATTCACTATT
TTGCCGTCGCGAACTAAGGTAATGCGATCGAGTAATTTACGGCAATGCGCGCAATGGGTCTGGCTGTGTT
TAAAGTCTTTAAGATAGCGGCTTAATGGCCGTCTTTTCAGGTGCTGCACCATCATTAAGAACTCCTGGTA
GTCAAAGTTAAAGTGCGGCATTTACTGACGTTATAACTTACCCAGTTTAGTGCGTAACCGTTTAATAGCC
TGGCTGTGTAACTGACTGACCCGCGATTCCCCGACCTCCAGCACCGCGCCAATCTCTTTGAGATTCAGCT
CTTCCTGGTAATAGAGGGTTAATACCAGTTTTTCGCGCTCCGGCAACGTTTCGATGGCTTCCATCACCCG
CTGGCGCAGATTACTGTCCAGTAGTTGTTGTAGCGGGTTTTCTCGCTGATGATCATCAGTAACCAGTTCG
ATGCTATCGCCGTGCTCTTCGCGCCACTCATCGTAGGAGAAGAGCTGGCTGTTATTGGTGTCGAGCAACA
TTTGGCGATAATCGGCAATATCGATCCCTAAACGTTCCGCTACCTCAGTTTCCGTGGCGTTGCGGCCAAG
TTCCTGCTCCAGTTGCCCTATTGCCTGTGCCACTTCACGCGCGTTGCGTCGCACGCTGCGCGGCACCCAG
TCACGGCTGCGAAGTTCATCCAGCATAGCGCCACGGATACGCTGCACTGCGTAAGTTGTAAATGCCGTTC
CTTGTAGGGCGTCATAGCGTTCGACGGCATTAAGTAACCCAATGCCGCCCGCCTGTAGCAGATCGTCAAG
TTCCACGCTCGCGGGCAGTCGAACCTGCAGGCGCAATGCTTCGTGACGCACCAGCGGGACATAACGCTGC
CACAGCGAGTGTTTATCCATTACACCTTCAGCGGTATAGAGTGAATTCACGATAAACAGCCCTGCGTTAT
ATGAGTTATCGGCATGATTATCCGTTTCTGCAGGGTTTTTAATCGGACGATTAGTGGGTGAAATGAGGGG
TTATTTGGGGGTTACAGGTAAATTCCAGGCAGAAAAAAACCCCGCCGGTGGCGGGGAAGCACGTTGCTGA
CAAATTGCGCTTTATGTTGCCGGATGCGGCGTAAACGCCTTATCCGGCCTACAAAAATGTGCAAATTCAA
TAAATTGCAATTCAACTTGTAGGCCTGATAAGCGCAGCGCATCAGGCAATTTGGCGTTGCCGTCAGTCTC
AGTTAATCAGGTTACAACGATTAACCCTGCAGCAGAGACAGAACCTGCTGCGGTACCTGGTTAGCTTTTG
CCAACACGGAGTTACCGGCCTGCTGGATGATCTGCGCTTTCGACATATTGGACACTTCGGTCGCATAGTC
GGCGTCCTGAATACGGGACTGCGCTTCAGACAGGTTGGTAGTGGTGTTGTTCAGGTTGGTAACCGCGGAA
TCCAGACGGTTTTGCACCGCACCGAGGGAAGAACGGAATTTGTCTACAGATGCGATAGCATCGTCCAGCG
CTTTCAGCGGATCCGTGGTTTTACCATTTGCAACAGCAGTCAGAGCCTCACCACCAGCAGTCAAACCTGT
TTGCAGATTACCGCCATTTAAATCGGCAGAATCGTATGTTTTACCATCAATATCGACCACTTCTGTTTTG
CCATCATCTCCGCCCAGTTTGACCGCGGTTGGAGAACTGGCGGCACCGGAAGAGTCAGTATAGGTAATAG
TTTTAACAGAAACAGCACCAGTAGTTTCATTCACATCCGCAGCGTAAAGATTGCCATTTGTATCTTTAAG
CGCATATGTATCGGTATCATTACCCTTGGAATCCTGCAGTTTTACTAAGCTAACAGCACCAAGGTTGGCA
GTTGCGGAACCTGCAGTATTATCAATCTGAACAGGTGTACCGCCTGAAGTGATAGTTGTAGCTTTAGTAG
TATTTGCATCAGTTACAGTTGCATTTGCCGTTGCTCCAGTCGCCATTGTCACTGTACCATCATTAGCAAC
TGTTACTGCGTAATACTTCCCATCGTTATCACCACCGGTGATTTTCGCATAGTAATCATTACCATTATCA
GTATAAACACCCTCAATTGAAGCTGGGTTAGTTCCGCCAGTATCAGTGGCTGCTTCCGTAGAAAGGGTAA
TTCCAGTAAGTTTAATATTGTTTGTGGTGGTAGCACCAAAAGCAGTTACTGGAGCACTAGTGGTAACTGT
ATCGTTATTTTTAACGCTAAAACCATCAAGGCCAAGAGTTTTAGCATCAATCTGCTTCAGATCGATAGTG
ATAGTCTGGTTATCATTTGCGCCAACCTGGATTTTCATGGAGCCATTTTTTGCCAGCACGTTCACGCCGT
TGAACTGGGTCTGACCAGATACGCGGTCAATTTCATCCAGACGGGATTTAATTTCGTCCTGGATAGAAGA
CAGATCAGACTCAGAGTTAGTACCGGTAGTGGCCTGTACCGTCAGTTCACGCACACGCTGTAAGTTGTTG
TTGATTTCGGACAGCGCGCCTTCGGTGGTCTGCGCAACGGAGATACCGTCGTTGGCGTTACGGGCCGCCT
GAGTCAGGCCTTTAATGTTAGAGGTGAAACGGTTAGCAATCGCCTGACCCGCTGCGTCATCCTTCGCGCT
GTTAATACGCA
>contig_06 len=8001 cov=6x
CCGTGGTTTTCTGGCGTGCTGCGTAGTATCAGCGGCAATTACCTGATGGACTGGGCATTTCATGCGCTGT
GCGTCGTTGGGCTGATGATCATAACCCTGCGTTTTGCACCAGTACGTTTTCCGCAGCTGTGGGTCAAAGA
GGCATGATGCGACGCTTGTTCCTGCGCTTTGTTCATGCCGGATGCGGCTAATGTAGATCGCTGAACTTGT
AGGCCTGATAAGCGCAGCGTATCAGGCAATTTTTATAATTTAAACTGACGATTCAACTTTATAATCTTTG
AAATAATAGTGCTTATCCCGGTCGTTTATTTCGCGAATAACCCGACAAGGAACGCCAGCCGCCACGACGT
TTGGTGGAATGTCTTTTGTGACGATACTACCCGCGCCAATAACAGAATTATCCCCGATGGTGACGCCTGG
ATTAATAACCACATGACTTCCGATCCAGACGTTATTGCCAATCGTTATCGGAAAAGAGTACATCTCGCCG
TTTTTTCTCAATTCATGGTGTACAGGGTGTCCCGTAACGGAAAGAGTAACGTTGGGTGCAATCAGTACGT
TATCACCGATTGTTACCGTGTAGTCATCGACAATGGTTAAATTGAAATTTGCATAAAAATTGCGGCCTAT
ATGGATGTTGGAACCGTAAGAGAAATAGACAGGCGGTTCTACCCAGGCGTTTTCCCCTACCGTGGCAAAC
ATTTCTTTAATCAGGCTTTCTCTTTTTTCAACTTCTGATGGATGCGAGTGATTAAACTCATACATTAACG
TTTTCCCACGAAGTCTTTTTTCCGGTAAGCCTTCGCACATATCGGTAAATAGCTTGCCTGCTCTTATTCT
TTCGGTCATTGGCATGTTCAATGCGATCACTCCGTTATGATATGTTGGTCGGATAAGGCGCTCGCGCCGC
ATCCGACATTGATTGCTTAAGCGACTTCATTCACCTGACGACGCAGCAGGGAAAGCGGGCCGGGGCCGCT
AAGCGTGAACACGGAAATTAAGGTGAAGCCCAGCGCCACCAGACCCAGCACCAGATAAGCGCCCTGGAAA
CCGATGCTTTCATACATATTGCCCGCCAGTACAGACATAAAAATCATCGCCAGTTGCTTAAAGAAGCAGA
AACAGACCAGATAAATCGTCGCTGAAAAACGCACTTCAAACTGGCTGGTAATATATTTAAAGCAGCCCAC
CAGCAGGAACGGTACTTCAAACATATGCAGCGTTTTCAGAATAACCACTTCCAGCGCTGAGGTGGCGAAC
GATGAGCCAATAATACGTACAGACATAATAGTGCCAGCCAGCAGCAGGGCGTTTTTCCCACCGATGCGAT
TAATGATCAGTGGCGCAAAGAACATAATCGAGGCGTTAAGTAATTCGCCCATTGTCGTTACGTAGCCAAA
TACCCGCGTACCCTGTTCACCGGTAGCAAAGAACGAAGTAAAGAAATTAGCAAACTGTTGGTCAAAAACA
TCGTAGGTGCAGGAAACGCCAATAACATACAGTGACAAAAACCACAGTTTTGGCTGTCTGAACAGTTCCA
GTGCCAGCTTAAGGCTAAATGCCGAATGGTTGGCACCTACCGCATTGGCAACCGTGGCAGAAGAGGGCGC
ATCCGTTTTGGCGAAAAAGAGTAAAACGGCGAGGATGAGTGCACAGCCAGAGCCCAGCCAGAAAACAAAC
TGATTATTGATGGTGAACATGATGCCGACAATCGAGGCACACAGCGCCCAGCCAACACAGCCAAACATCC
GCGCGCGACCAAATTCGAAATTACTGCGACGGCTGACTTTCTCAATAAATGCCTCTACTGCTGGCGCACC
GGCGTTAAAACAAAAGCCTAGATAAATACCACCAACAATCGATCCTACTAAAATGTTGTATTGTAACAGT
GGCCCGAAGATAAAAATAAAGAACGGCGCAAACATCACTAACATGCCGGTAATAATCCACAGCAGGTATT
TGCGCAGCCCGAGTTTGTCAGAAAGCAGACCAAACAGCGGTTGGAATAATAGCGAGAACAGAGAAATAGC
GGCAAAAATAATACCCGTATCACTTTTGCTGATATGGTTGATGTCATGTAGCCAAATCGGGAAAAACGGG
AAGTAGGCTCCCATGATAAAAAAGTAAAAGAAAAAGAATAAACCGAACATCCAAAAGTTTGTGTTTTTTA
AATAGTACATAATGGATTTCCTTACGCGAAATACGGGCAGACATGGCCTGCCCGGTTATTATTATTTTTG
ACACCAGACCAACTGGTAATGGTAGCGACCGGCGCTCAGCTGGAATTCCGCCGATACTGACGGGCTCCAG
GAGTCGTCGCCACCAATCCCCATATGGAAACCGTCGATATTCAGCCATGTGCCTTCTTCCGCGTGCAGCA
GATGGCGATGGCTGGTTTCCATCAGTTGCTGTTGACTGTAGCGGCTGATGTTGAACTGGAAGTCGCCGCG
CCACTGGTGTGGGCCATAATTCAATTCGCGCGTCCCGCAGCGCAGACCGTTTTCGCTCGGGAAGACGTAC
GGGGTATACATGTCTGACAATGGCAGATCCCAGCGGTCAAAACAGGCGGCAGTAAGGCGGTCGGGATAGT
TTTCTTGCGGCCCTAATCCGAGCCAGTTTACCCGCTCTGCTACCTGCGCCAGCTGGCAGTTCAGGCCAAT
CCGCGCCGGATGCGGTGTATCGCTCGCCACTTCAACATCAACGGTAATCGCCATTTGACCACTACCATCA
ATCCGGTAGGTTTTCCGGCTGATAAATAAGGTTTTCCCCTGATGCTGCCACGCGTGAGCGGTCGTAATCA
GCACCGCATCAGCAAGTGTATCTGCCGTGCACTGCAACAACGCTGCTTCGGCCTGGTAATGGCCCGCCGC
CTTCCAGCGTTCGACCCAGGCGTTAGGGTCAATGCGGGTCGCTTCACTTACGCCAATGTCGTTATCCAGC
GGTGCACGGGTGAACTGATCGCGCAGCGGCGTCAGCAGTTGTTTTTTATCGCCAATCCACATCTGTGAAA
GAAAGCCTGACTGGCGGTTAAATTGCCAACGCTTATTACCCAGCTCGATGCAAAAATCCATTTCGCTGGT
GGTCAGATGCGGGATGGCGTGGGACGCGGCGGGGAGCGTCACACTGAGGTTTTCCGCCAGACGCCACTGC
TGCCAGGCGCTGATGTGCCCGGCTTCTGACCATGCGGTCGCGTTCGGTTGCACTACGCGTACTGTGAGCC
AGAGTTGCCCGGCGCTCTCCGGCTGCGGTAGTTCAGGCAGTTCAATCAACTGTTTACCTTGTGGAGCGAC
ATCCAGAGGCACTTCACCGCTTGCCAGCGGCTTACCATCCAGCGCCACCATCCAGTGCAGGAGCTCGTTA
TCGCTATGACGGAACAGGTATTCGCTGGTCACTTCGATGGTTTGCCCGGATAAACGGAACTGGAAAAACT
GCTGCTGGTGTTTTGCTTCCGTCAGCGCTGGATGCGGCGTGCGGTCGGCAAAGACCAGACCGTTCATACA
GAACTGGCGATCGTTCGGCGTATCGCCAAAATCACCGCCGTAAGCCGACCACGGGTTGCCGTTTTCATCA
TATTTAATCAGCGACTGATCCACCCAGTCCCAGACGAAGCCGCCCTGTAAACGGGGATACTGACGAAACG
CCTGCCAGTATTTAGCGAAACCGCCAAGACTGTTACCCATCGCGTGGGCGTATTCGCAAAGGATCAGCGG
GCGCGTCTCTCCAGGTAGCGAAAGCCATTTTTTGATGGACCATTTCGGCACAGCCGGGAAGGGCTGGTCT
TCATCCACGCGCGCGTACATCGGGCAAATAATATCGGTGGCCGTGGTGTCGGCTCCGCCGCCTTCATACT
GCACCGGGCGGGAAGGATCGACAGATTTGATCCAGCGATACAGCGCGTCGTGATTAGCGCCGTGGCCTGA
TTCATTCCCCAGCGACCAGATGATCACACTCGGGTGATTACGATCGCGCTGCACCATTCGCGTTACGCGT
TCGCTCATCGCCGGTAGCCAGCGCGGATCATCGGTCAGACGATTCATTGGCACCATGCCGTGGGTTTCAA
TATTGGCTTCATCCACCACATACAGGCCGTAGCGGTCGCACAGCGTGTACCACAGCGGATGGTTCGGATA
ATGCGAACAGCGCACGGCGTTAAAGTTGTTCTGCTTCATCAGCAGGATATCCTGCACCATCGTCTGCTCA
TCCATGACCTGACCATGCAGAGGATGATGCTCGTGACGGTTAACGCCTCGAATCAGCAACGGCTTGCCGT
TCAGCAGCAGCAGACCATTTTCAATCCGCACCTCGCGGAAACCGACATCGCAGGCTTCTGCTTCAATCAG
CGTGCCGTCGGCGGTGTGCAGTTCAACCACCGCACGATAGAGATTCGGGATTTCGGCGCTCCACAGTTTC
GGGTTTTCGACGTTCAGACGTAGTGTGACGCGATCGGCATAACCACCACGCTCATCGATAATTTCACCGC
CGAAAGGCGCGGTGCCGCTGGCGACCTGCGTTTCACCCTGCCATAAAGAAACTGTTACCCGTAGGTAGTC
ACGCAACTCGCCGCACATCTGAACTTCAGCCTCCAGTACAGCGCGGCTGAAATCATCATTAAAGCGAGTG
GCAACATGGAAATCGCTGATTTGTGTAGTCGGTTTATGCAGCAACGAGACGTCACGGAAAATGCCGCTCA
TCCGCCACATATCCTGATCTTCCAGATAACTGCCGTCACTCCAGCGCAGCACCATCACCGCGAGGCGGTT
TTCTCCGGCGCGTAAAAATGCGCTCAGGTCAAATTCAGACGGCAAACGACTGTCCTGGCCGTAACCGACC
CAGCGCCCGTTGCACCACAGATGAAACGCCGAGTTAACGCCATCAAAAATAATTCGCGTCTGGCCTTCCT
GTAGCCAGCTTTCATCAACATTAAATGTGAGCGAGTAACAACCCGTCGGATTCTCCGTGGGAACAAACGG
CGGATTGACCGTAATGGGATAGGTCACGTTGGTGTAGATGGGCGCATCGTAACCGTGCATCTGCCAGTTT
GAGGGGACGACGACAGTATCGGCCTCAGGAAGATCGCACTCCAGCCAGCTTTCCGGCACCGCTTCTGGTG
CCGGAAACCAGGCAAAGCGCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGC
CTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGG
TTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATCCGTAATCATGGTCATAGCTGTTTCCTGTG
TGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTG
CCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTC
GTGCCAGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCCAGGGTGGT
TTTTCTTTTCACCAGTGAGACGGGCAACAGCTGATTGCCCTTCACCGCCTGGCCCTGAGAGAGTTGCAGC
AAGCGGTCCACGCTGGTTTGCCCCAGCAGGCGAAAATCCTGTTTGATGGTGGTTAACGGCGGGATATAAC
ATGAGCTGTCTTCGGTATCGTCGTATCCCACTACCGAGATATCCGCACCAACGCGCAGCCCGGACTCGGT
AATGGCGCGCATTGCGCCCAGCGCCATCTGATCGTTGGCAACCAGCATCGCAGTGGGAACGATGCCCTCA
TTCAGCATTTGCATGGTTTGTTGAAAACCGGACATGGCACTCCAGTCGCCTTCCCGTTCCGCTATCGGCT
GAATTTGATTGCGAGTGAGATATTTATGCCAGCCAGCCAGACGCAGACGCGCCGAGACAGAACTTAATGG
GCCCGCTAACAGCGCGATTTGCTGGTGACCCAATGCGACCAGATGCTCCACGCCCAGTCGCGTACCGTCT
TCATGGGAGAAAATAATACTGTTGATGGGTGTCTGGTCAGAGACATCAAGAAATAACGCCGGAACATTAG
TGCAGGCAGCTTCCACAGCAATGGCATCCTGGTCATCCAGCGGATAGTTAATGATCAGCCCACTGACGCG
TTGCGCGAGAAGATTGTGCACCGCCGCTTTACAGGCTTCGACGCCGCTTCGTTCTACCATCGACACCACC
ACGCTGGCACCCAGTTGATCGGCGCGAGATTTAATCGCCGCGACAATTTGCGACGGCGCGTGCAGGGCCA
GACTGGAGGTGGCAACGCCAATCAGCAACGACTGTTTGCCCGCCAGTTGTTGTGCCACGCGGTTGGGAAT
GTAATTCAGCTCCGCCATCGCCGCTTCCACTTTTTCCCGCGTTTTCGCAGAAACGTGGCTGGCCTGGTTC
ACCACGCGGGAAACGGTCTGATAAGAGACACCGGCATACTCTGCGACATCGTATAACGTTACTGGTTTCA
CATTCACCACCCTGAATTGACTCTCTTCCGGGCGCTATCATGCCATACCGCGAAAGGTTTTGCGCCATTC
GATGGTGTCAACGTAAATGCATGCCGCTTCGCCTTCCGGCCACCAGAATAGCCTGCGATTCAACCCCTTC
TTCGATCTGTTTTGCTACCCGTTGTAGCGCCGGAAGATGCTTTTCCGCTGCCTGTTCAATGGTCATTGCG
CTCGCCATATACACCAGATTCAGACAGCCAATCACCCGTTGTTCACTGCGCAGCGGTACGGCGATAGAGG
CGATCTTCTCCTCCTGATCCCAGCCGCGGTAGTTCTGTCCGTAACCCTCTTTGCGCGCGCGCGCCAGAAT
GGCTTCCAGCTTTAACGGTTCCCGTGCCAGTTGATAGTCATCACCGGGGCGGGAGGCTAACATTTCGATT
AATTCCTTGCGGTCTTGTTCCGGGCAAAAGGCCAGCCAGGTCAGGCCCGAGGCGGTTTTCAGAAGCGGCA
AACGTCGCCCGACCATTGCCCGGTGAAAGGATAAGCGGCTGAAACGGTGAGTGGTTTCGCGTACCACCAT
TGCATCAACATCCAGCGTGGACACATCTGTCGGCCATACCACTTCGCGCAACAGATCGCCCAGCAGTGGG
GCCGCCAGTGCAGAAATCCACTGTTCGTCACGAAATCCTTCGCTTAATTGCCGCACTTTGATGGTCAGTC
GAAAACTATCATCGGAGGGGCTACGGCGGACATATCCCTCTTCCTGCAGCGTCTCCAGCAGTCGCCGCAC
AGTGGTGCGATGCAGGCCGCTGAGTTCCGCCAGCAGCCCGACGCTGGCACCGCCATCAAGTTTATTTAAC
ATATTTAATAACATTAGACCGCGGGTTAAGCCGCGCACGGTTTTGTATTCCGTCTGCTCATTGTTCTGCA
TATTAATTGACATTTCTATAGTTAAAACAACGTGGTGCACCTGGTGCACATTCGGGCATGTTTTGATTGT
AGCCGAAAACACCCTTCCTATACTGAGCGCACAATAAAAAATCATTTACATGTTTTTAACAAAATAAGTT
GCGCTGTACTGTGCGCGCAACGACATTTTGTCCGAGTCGTGAGGTACTGAAATGGCAATACAACACCCTG
ACATCCAGCCTGCTGTTAACCATAGCGTTCAGGTGGCGATCGCTGGTGCCGGCCCGGTTGGGCTGATGAT
GGCGAACTATCTCGGCCAGATGGGCATTGACGTGCTGGTGGTGGAGAAACTCGATAAGTTGATCGACTAC
CCGCGTGCGATTGGTATTGATGACGAGGCGCTGCGCACCATGCAGTCGGTCGGCCTGGTCGATGATGTTC
TGCCGCACACTACGCCGTGGCACGCGATGCGTTTTCTCACCCCGAAAGGCCGCTGTTTTGCTGATATTCA
GCCAATGACCGATGAATTTGGCTGGCCGCGCCGTAACGCCTTTATTCAGCCGCAGGTCGATGCGGTGATG
CTGGAAGGGGTGTCGCGTTTT
>contig_07 len=10000 cov=118x
CAAAATAAACGTCCAGACTATATGACTGCATTTTGGAACATTGTTAACTGGAAAAAAGTTGATGAATTAT
ACCAAGCAGCAAAATAATATAACTTAATATAATTGAGGTGGAGCATCTACAAGGTGTTCCACTTTTTTGT
ACTTAATATCTTTACTTATTTGTTTTTGTGTAATGGGATAGCACGTAAAAGTGGAAAAATAATAAGATAA
ATTTCAACTCGAAATACCTTATTATTGTTGACTATAACGATGACCTAAATATTATATACATTAGACAACA
CTATTAATTAATCAACATATTGACATTAAATAAATTGACAAAATAAGTAATTATTGTGAATCTAAAGTGA
AATTTTTATAAAAAAATGTAATGATTCAAAAATTTTGTTGCATTTCTTTTGTAATCGTATGATAATGTAA
ATGTAATCAAATTGTAATATAAGGGGACAAGACAATGAAAAAATTAGCAACAGTAGGTTCTTTAATTGTA
ACAAGCACTTTAGTATTCTCAAGTATGCCTTTTCAAAATGCGCATGCCGACACAACTTCAATGAATGTGT
CGAATAAACAAAGCCAAAATGTACAAAATCATCGTCCTTATGGCGGAGTAGTACCACAAGGAATGACGCA
AGCACAATATACTGAATTAGAGAAAGCTTTACCCCAATTAAGCGCTGGCAGTAATATGCAAGACTATAAT
ATGAAATTGTATGATGCGACGCAAAATATTGCTGATAAATACAATGTGATAATTACAACTAATGTAGGGG
TATTTAAACCACATGCTGTTAGAGATATGAATGGCCATGCGTTACCTTTAACAAAAGATGGCAATTTTTA
TCAAACGAATGTAGATGCAAATGGTGTTAATCATGGTGGTAGTGAAATGGTGCAAAATAAAACAGGTCAT
ATGAGTCAACAAGGCCATATGAATCAGAACACACATGAACCAACAGCCACACATGCAACAAGGTCATATG
CAATCATCAAACCATCAAATGATGAGTCCAAAAGCAAATATGCATTCATCAAATCATCAAATGAACCAAA
GTAACAAAAAAGTTTTACCAGCTGCTGGTGAAAGTATGACATCAAGTATTCTTACTGCAAGTATTGCCGC
ACTACTATTAGTATCTGGGTTATTCTTAGCATTTAGACGACGTTCAACAAATAAATAAACATAATACGAT
TAATAATAGAAAAATCGTGTGATTATCTGAGCGAGCCTAGGACATAAATCAATGTCCTAGGCTCGCTAAT
GTTATATTGGCAGTAGTTGACTGAATGAAATTGCGCTTGTAACAAGCTTTTCCATTTCTTACCAACTCCC
TAAATAGTCATCAAAAAATTTCTTATATTTTAATAATTTTTAATAATCCGATTGTCTTATACGTGTCAGT
GTTAATTCAGATATTTCCTGTGGAATATACCACTTATTAATCATAATTGGATAAGGTGTTTGTGCGTACA
GTGTTTCAATAATCAGCCAACAATGTGTATCACCATCAAACACGTGACTATGATTTTTGAAGTGGGGCGC
TTTGGTAATAGACATTTTTAAATCTGATTGATATGCATTGCTATAAATCGTTTGCTCAACGAATGTCTTC
ATGTCGTCTTCGTTTTGTGTATTCACTTTAAATGTGTCAATGACATTTAACGGTATAAAGGTAAAGCAAA
ATGCATCAGCTTGCTTAGAATGATTGTCCTTTTTTTGATAATAGCGTTCCATTGCAATGACGGCAGAAGG
ATGGTTTGCAAACAAATGATTTGTATATTCACTTTCTAAATCAACACGATAATTAATTGATGACATAGAT
ACGCGAGCTAGCAATATTTGATCAAGTGGATGCTTAAATTGATCCATACTTGAAGCGTGTTGGGCATTTG
TTTGTGGAATAACAAAGTGTCCCTTCCCTCTTGTACTCTCTACGATGCCATCTTCGGCTAACAATTTTAT
AGCTTGGCGCAAAGTCATACGACTGACATCAAAGCGCGCACAAAGTTCCTTTTCAGTAGGTAATGCATGG
CCACTCGGATATTTTCCTATTTGAATTTCTTTATATAACGTATTATAAATCGTTAAAAATTTTGGTTGTG
TTTGCGTCACGTAGACAACCTCCATAAAGTTACTTAATCACTCTCATCATACAATAATTTTTACTCAAAT
TGGAAAAATTATAAAAATTAAATATAGATAGGCTTTGAAAATTAGTTTTATACAAGGTTAGTAGCTGTAA
CTGTAAAATGTTCTTAATATTGTCAAAATGTAATGCTTGAAAGCGCTTTTAAAAAATATTATTATATACA
TGGTTAGACAAATAGACAAATCACTATACAAATATTGGGAGGAATATTTTATGAAATCAACACCACACAT
TAAACCAATGAATGACGTCGAAATTGCAGAAACGGTTCTATTGCCAGGAGATCCGTTAAGAGCTAAGTTC
ATTGCAGAAACTTATTTGGATGATGTGGAACAGTTCAATACAGTGCGAAACATGTTTGGTTTTACCGGAA
CATATAAAGGTAAAAAAGTTTCTGTCATGGGTTCAGGTATGGGTATGCCATCTATTGGCATTTACTCTTA
TGAATTAATTCATACATTTGGTTGTAAAAAATTAATTCGCGTTGGCTCTTGTGGCGCGATGCAAGAAAAC
ATTGATTTATATGATGTGATTATTGCACAAGGTGCCTCTACTGATTCAAATTACGTTCAACAATATCAAT
TACCAGGTCATTTTGCGCCAATTGCTTCTTATCAATTATTAGAAAAAGCAGTTGAAACAGCACGTGACAA
AGGTGTACGTCATCATGTAGGTAATGTGTTATCAAGTGATATTTTCTATAACGCGGATACAACAGCGAGT
GAACGTTGGATGCGTATGGGTATTTTAGGTGTAGAAATGGAATCAGCTGCATTATACATGAATGCAATTT
ACGCTGGTGTCGAAGCATTAGGTGTGTTCACAGTGAGCGATCATTTAATTCATGAAACGTCAACAACACC
TGAGGAAAGGGAACGTGCATTTACAGATATGATTGAAATTGCACTGTCATTGGTGTAGATGATTATGAAT
GTTGAATATTCTAAAATAAAGAAAGCAGTACCTATTTTATTATTCTTATTTGTATTCAGTTTGGTTATAG
ACAACTCATTTAAATTGATTTCTGTAGCCATTGCTGATGACTTAAACATATCTGTAACGACAGTAAGTTG
GCAAGCGACATTAGCCGGTTTAGTAATTGGTATTGGCGCTGTAGTATACGCTTCATTATCTGATGCCATT
AGTATACGCACACTATTTATTTATGGCGTGATATTAATCATTATCGGATCAATTATTGGTTACATTTTCC
AACATCAATTCCCATTACTTTTAGTTGGACGTATTATTCAAACTGCCGGTTTAGCTGCTGCAGAGACATT
ATATGTGATATATGTTGCAAAGTATCTTTCTAAAGAGGACCAGAAGACTTACCTTGGCTTAAGTACGAGC
AGTTATTCCTTGTCATTAGTTATCGGTACATTATCAGGTGGATTTATTTCTACGTATTTACACTGGACAA
ATATGTTTTTAATTGCATTAATCGTAGTATTTACGTTGCCATTCCTATTTAAATTATTACCAAAAGAAAA
TAATACGAATAAAGCTCATTTAGATTTTGTTGGCTTAATTCTAGTGGCAACTATTGCTACAACAGTCATG
CTGTTTATTACGAACTTTAATTGGTTATATATGATTGGTGCCTTAATTGCGATTATCGTTTTTGCGCTAT
ATATTAAAAATGCGCAACGTCCATTAGTAAATAAATCATTTTTCCAAAATAAACGTTATGCTTCATTTTT
ATTTATAGTATTTGTAATGTATGCTATCCAATTGGGTTATATTTTTACGTTCCCATTCATAATGGAGCAA
ATTTATCATCTGCAACTAGACACAACATCACTGTTATTAGTACCGGGTTATATAGTAGCAGTCATTGTTG
GTGCACTAAGTGGTAAAATCGGCGAATATCTGAATTCAAAACAAGCGATTATCACAGCAATTATTTTAAT
AGCACTGAGCTTGATTTTACCTGCATTTGCAGTAGGTAATCACATTTCAATCTTCGTCATTTCTATGATA
TTCTTTGCAGGTAGCTTTGCTTTAATGTATGCACCTTTACTTAACGAAGCCATTAAAACAATAGATCTTA
ATATGACAGGTGTGGCTATTGGTTTTTATAATTTAATTATTAATGTGGCGGTATCTGTAGGTATTGCGAT
TGCTGCGGCTCTAATCGATTTTAAAGCATTAAATTTCCCAGGCAATGATGCATTAAGTTCACATTTCGGT
ATTATTTTAATTATTTTAGGTTTAATGAGTATTGTCGGATTAGTTTTATTCGTCAGCTTAAATCGTTGGA
CACAATCTGAAAAATAAATAGATATAAATTCGCGAGATATATTCGTATTTATAGTAAAATTAAATAAAGA
GATTATATAACACGAGGAGTAGTAAGTATGAAATTTGAGAAATATATAGATCACACTTTATTGAAGCCTG
AGTCAACACGTACGCAAATCGATCAAATCATCGATGAAGCGAAAGCATACAATTTTAAATCTGTATGTGT
GAATCCAACACATGTTAAATATGCAGCAGAGCGACTAGCTGATTCAGAGGTGCTCGTTTGTACGGTAATA
GGATTCCCATTAGGTGCGTCGACAACTGCAACGAAAGCATTTGAAACAGAAGATGCAATTCAAAATGGTG
CAGATGAAATTGACATGGTCATCAACATCGGCGCATTAAAAGATGGACGTTTTGATGATGTACAACAAGA
CATTGAAGCAGTGGTTAAAGCTGCGAAAGGTCACACAGTAAAAGTGATTATTGAGACGGTATTGTTGGAC
CATGACGAAATTGTAAAAGCGAGTGAATTAACAAAAGCGGCTGGTGCGGACTTCGTTAAAACTTCAACAG
GTTTTGCAGGTGGCGGTGCGACTGCAGAAGACGTTAAATTAATGAAAGATACAGTAGGTGCTGATGTAGA
AGTAAAAGCATCAGGTGGCGTACGTAATTTAGAAGATTTCAATAAAATGGTTGAAGCAGGTGCGACACGT
ATTGGTGCGAGCGCAGGTGTTCAAATTATGCAAGGTTTAGAAGCAGATTCAGATTACTAATATATATAAA
TTTTGGGAGTGATAGCTATGACAAGACCATTTAATCGTGTACATTTAATCGTAATGGATTCAGTAGGTAT
TGGTGAAGCGCCAGACGCAGCTGATTTTAAAGATGAAGGTTCACATACTTTAAGACATACCTTAGAAGGT
TTCGATCAAACTTTACCAAACCTTGAAAAGTTAGGTCTAGGGAACATCGATAAATTACCAGTAGTAAATG
CAGTTGAACAACCAGAAGCATACTATACTAAATTGAGTGAAGCTTCAGTTGGTAAAGATACAATGACTGG
TCACTGGGAAATTATGGGATTAAATATTATGCAACCTTTTAAAGTATACCCTAATGGATTCCCTGAAGAG
TTAATTCAACAAATTGAAGAAATGACAGGTCGTAAAGTTGTTGCTAACAAACCGGCATCGGGTACGCAAA
TTATCGATGAGTGGGGCGAGCACCAAATGAAAACTGGTGACTTAATTGTTTATACAAGTGCAGACCCAGT
ATTGCAAATTGCTGCACATGAAGACATTATCCCATTAGAAGAGTTATATGATATTTGTGAAAAGGTTCGT
GAGTTGACAAAAGACCCTAAATATTTAATTGGTCGTATTATCGCACGTCCATATGTTGGTGAACCAGGAA
ACTTTACACGTACATCTAATCGACATGACTATGCGTTAAAACCTTTTGGTAAAACTGTCTTAGATCATTT
GAAAGACGGTGGTTATGATGTTATTGCCATCGGTAAAATTAATGACATTTATGATGGTGAAGGTGTAACA
GAAGCGGTTCGTACGAAGAGTAACATGGACGGTATGGATCAATTGATGAAAATTGTTAAGAAAGATTTCA
CAGGTATTAGCTTCTTAAACTTAGTAGACTTTGATGCATTATACGGTCATCGTCGTGATAAACCAGGTTA
TGCACAAGCAATTAAAGATTTCGATGATCGCTTGCCAGAACTGTTTAGCAACTTAAAAGAAGACGATTTA
GTAATTATTACAGCAGACCATGGTAATGACCCGACAGCGCCAGGTACGGACCATACGAGAGAATATATCC
CAGTAATTATGTACAGTCCGAAATTTAAAGGTGGTCATGCACTAGAAAGTGATACTACATTCAGTTCTAT
CGGTGCAACTATAGCAGATAATTTCAACGTAACATTACCAGAGTTCGGTAAAAGTTATTTAAAGGAATTG
AAATAGAATAAATTTAGATATTATAAAAACAGCAGTGAAGTTAACTATAACAATAGTTTTCTTCACTGCT
GTTTTTATTATAATAGAGAAACGTAAGACGGTAGGACTTCTTATTTAGGAGTATCCTGATTTAATGTTAA
ACAATACGTTTTCGGATTGAACCGGAAATTAAATCGACAATTGCGACCATTAGTACTAAACCGATTAATA
TAATACCTACACGGTCCCAAGAACGTGTTTGAATGGCAAATATGAGTGGTGTCCCGATACCACCAGCCCC
AATTAGCCCCAGTATAGAAGCTGAACGTAAGTTTAGTTCAAAGCGATAAAGTATGAGTGATAGAAAGGCA
GGCATAATTTGTGGTATGACTGCAAATACGAGTGTTTTAATCTTATTCGCACCACTGGCCTTTAATGATT
CTACAGCACTGAAATCTAGACCTTCAATATCTTCAGCTAAAAGTTTCCCAAGCATACCTACGGAATGGAT
ACCTAAAGCTAATACACCTGAAAATGAACCTGGGCCAACAGCTTTGATAAATATAAGTGCCATTACAATT
TCTGGGAAGACACGTATAACACTTAAAATAAATTTGCTAACACCTGAAACTGGGCGTAGCTTTACCATAT
TATTTGCACCTAGAAATGCTAATGGAATACAGATAATTGCGGCGATGAAAGTACCTACAACGGCTATCGC
AAAGGTTTCAAGTAAACCACGTAATAAGTCTTCGCCATCTGGTATATAGATATAGCTGATATCAGGATGG
AATAATCCGCTGAATATGGATTTTAAGATTTCTAATGATTTACTTTTAAGTTCTAAACTTGGTACACCTG
CAAATGCCCAGATGATAATAGCTAAGACGACAATTGCAATAAGCCATCTTTTAATCAATTTTCGTTTGTG
TGCTTTTGTGTGAACATTATATTTTGCTATTTCCTGTGTCATGCGAGATGTGCCCTCACTTTCGTACTGA
TGTAATCAATGACGACGACGATAACTAAAGTAAATAAAATAATCGTTGCTGTTTTTGGATATTGAAATAA
ACCAAGTGTTTGATCATAAAACAATCCAATACCGCCAGCGCCGACTAATCCAAGCACAGCTGAAGCACGT
ATATTTACTTCAAATGCATATAATACGTATGACATAAATGACGATATGGCTTGTGGTACAACACCGAAAA
CAATCCATTTTATTTTATTAGCGCCAACAGCCGTCATTGCTTCCATTGGACCTGGATCTATCGTTTCCAA
TGATTCATATAATAATTTTCCAATAATACAGATAGTTAAAATAAACAGTGCTAATATCCCTGGAATTTGA
CCGATTCCAAATACAGCCACAAAGATTGCTGCTAATAACAAATCTGGAATAGTACGAACTATATTTAAAA
TAAAGCGCGAGGGTATTGAAATCCACTTTTGATGAACGATATTGCTAGCACATAATAACGCAATTGGTAT
TGAAACGATGCTACCTAATACTGTACTTACGATAGCCATTCGAATGGTATCTAACATTGGCGTTGTAATT
TGTTGTAAATACTCGAAATCAGGTGGAATCATTTGTTTGAATAGATCACCTATTTGAGGTATTCCTATCA
TTAAATCTCCAAAATTAAATCCCGTATAAATGAAGCTCCAAATGATAAGCACAATGATTAACATGAAGGT
AAAACTCGTTTTTAAAGAAACCTTTTTCTTTAAAAGGGAGTCATACTTTGTAGGTATTTCTAAAGGCATG
TTAGTTCACTCCTAGCTTTTCATCTTCTTTAATTGTACGTCCATATATTTCACTAAATACGTCATCTGTT
GCTTCAGATGCAGGACCATCATAGACAACTTCACCATCACGTAAACCAATGATGCGTGTGCCATATTCTT
TTGCCAAGTCAACAAAATGTAAATTAATTAAAATTGTGATGCCTAATTCTTGGTTGATTTTTCTTAAATC
ATCCATAACCTGTTTCGTAGTTAATGGGTCTAATGAAGCAACTGGTTCATCTGCAAGAATAATTTCAGAT
TCTTGGCATAGCGCACGTGCAATAGATATACGTTGTTGTTGGCCACCTGATAATTCATCAGAGCGTTGAT
TATATTTATCTAAGATATTGACGCGTTCTAGTGCATCCATTGCCTTAATTTTGTCTTCTTTTGGGAATAA
ACCTAATACCATTTTCCAAGTAGGGTGATAACCTACACGTCCACTTAGTACATTTCGTAATACACTTGAC
CGTTTAACTAAATTAAAATGTTGGAAAATCATACCTATATTTCGGCGCATTTCTAATAATGCTTTACCAT
GGGCTTTAGTGATTGATTTACCTTGGATGAAAATTTCACCTGACGTGATATCATGCAAACGATTTACAGA
TCTTAATAACGTGGATTTCCCAGCACCAGATAGTCCGACAATAACTGCAAATTCACCTTTTTCAATATTT
AAGTTAATATTTTTCAAGCCTACATGACCGTTAGGATAGACTTTACTGACGTTTTTAAATTCGATTTGAC
TCATGTTTGACACCTTTCTTTAATAAAGAAAAGCCAAGGTAAACTGAAATATCATAGACATTGGGTTTCT
GTGTACAAATGGTTTCAAATGTACAGCGACACCTTTCTATATCTGTTTCAGTTTAGCGCCTAGCTTTTCA
ATATTAATCTAGAATATATCTATTGAACGAAAGCTTTTAATACCAAATTCGCTAATGATTCATTTGTTAA
ATAATGATTATTTCATATCTTTAACTAATTTTTCGTACTCTCTTACAATGTCGAAATTTGAATCTTTCGT
TTCTGTGTATCCTTCATGTGAATAAACTTCGCTAATAATTTTGTGACCTTCTTTTGATTTAGCAATGTCT
ATAAAAGCTTTTTTCAATTTTTCTTGAAAATCTTTATCCATATCTGGTCTTACAGAAATTGTGTCATTCG
GAATAGCTTGTGTTAATTTTAAAATTCGTGTGTCTTTAAATACATTTGGTTGGTCTTTTTTCACAGTATT
ACGTGCATCGTTAAATACAGCCGCAGCATCTACATCTCCATTTAATAATGAGATAACTGCTTGGTCATGA
CCTTTAACATTCACAATTTTCATATCTTTAGTTGCATTAATACCTGCTTCGTTTTTTAACATCGCAAGTG
GGAATGTATATCCAGCAGTTGATGTTACATCTTGTAAGGCAATTTTCTTACCTTTTAAATCTTTCAAGCT
TTTAATTTTTGAGTCTTTTTTAACAAGAATTTCTGATTTATAACTATCTACAAGTTCTTTACTTGCTGAA
CCATCTTCTTTTACACCGAAACGTTGTGCTTGTAATAATAAATCAGCTGCTTTTTGATCATGTGCTAATG
TGTATGCCGTTGGTGGTAAGAAACCAACATCAACTTTTTTAGACTTCATAGCTTCAACAATTGTATTGTA
GTTAGTTGATACAGACACTTTAACTGGAATCCCTAATTCTTTAGATAGTAATTTTTCTAATGGTTTTGCT
TTAGCTTCTAATGTTCCAGCATTTTGCGAAGGTACAAATTGAACGGTTAATTCTTTAGGTTTGTATCCTC
CTGATTTAGAATCCGAATCATTACTAGCGTTCTTTTGATTATCTAAAGAACTTGAGTTTC
>contig_08 len=5386 cov=412x
GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT
GATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCTGGCGGAAAATGAGAAA
ATTCGACCTATCCTTGCGCAGCTCGAGAAGCTCTTACTTTGCGACCTTTCGCCATCAACTAACGATTCTG
TCAAAAACTGACGCGTTGGATGAGGAGAAGTGGCTTAATATGCTTGGCACGTTCGTCAAGGACTGGTTTA
GATATGAGTCACATTTTGTTCATGGTAGAGATTCTCTTGTTGACATTTTAAAAGAGCGTGGATTACTATC
TGAGTCCGATGCTGTTCAACCACTAATAGGTAAGAAATCATGAGTCAAGTTACTGAACAATCCGTACGTT
TCCAGACCGCTTTGGCCTCTATTAAGCTCATTCAGGCTTCTGCCGTTTTGGATTTAACCGAAGATGATTT
CGATTTTCTGACGAGTAACAAAGTTTGGATTGCTACTGACCGCTCTCGTGCTCGTCGCTGCGTTGAGGCT
TGCGTTTATGGTACGCTGGACTTTGTGGGATACCCTCGCTTTCCTGCTCCTGTTGAGTTTATTGCTGCCG
TCATTGCTTATTATGTTCATCCCGTCAACATTCAAACGGCCTGTCTCATCATGGAAGGCGCTGAATTTAC
GGAAAACATTATTAATGGCGTCGAGCGTCCGGTTAAAGCCGCTGAATTGTTCGCGTTTACCTTGCGTGTA
CGCGCAGGAAACACTGACGTTCTTACTGACGCAGAAGAAAACGTGCGTCAAAAATTACGTGCGGAAGGAG
TGATGTAATGTCTAAAGGTAAAAAACGTTCTGGCGCTCGCCCTGGTCGTCCGCAGCCGTTGCGAGGTACT
AAAGGCAAGCGTAAAGGCGCTCGTCTTTGGTATGTAGGTGGTCAACAATTTTAATTGCAGGGGCTTCGGC
CCCTTACTTGAGGATAAATTATGTCTAATATTCAAACTGGCGCCGAGCGTATGCCGCATGACCTTTCCCA
TCTTGGCTTCCTTGCTGGTCAGATTGGTCGTCTTATTACCATTTCAACTACTCCGGTTATCGCTGGCGAC
TCCTTCGAGATGGACGCCGTTGGCGCTCTCCGTCTTTCTCCATTGCGTCGTGGCCTTGCTATTGACTCTA
CTGTAGACATTTTTACTTTTTATGTCCCTCATCGTCACGTTTATGGTGAACAGTGGATTAAGTTCATGAA
GGATGGTGTTAATGCCACTCCTCTCCCGACTGTTAACACTACTGGTTATATTGACCATGCCGCTTTTCTT
GGCACGATTAACCCTGATACCAATAAAATCCCTAAGCATTTGTTTCAGGGTTATTTGAATATCTATAACA
ACTATTTTAAAGCGCCGTGGATGCCTGACCGTACCGAGGCTAACCCTAATGAGCTTAATCAAGATGATGC
TCGTTATGGTTTCCGTTGCTGCCATCTCAAAAACATTTGGACTGCTCCGCTTCCTCCTGAGACTGAGCTT
TCTCGCCAAATGACGACTTCTACCACATCTATTGACATTATGGGTCTGCAAGCTGCTTATGCTAATTTGC
ATACTGACCAAGAACGTGATTACTTCATGCAGCGTTACCATGATGTTATTTCTTCATTTGGAGGTAAAAC
CTCTTATGACGCTGACAACCGTCCTTTACTTGTCATGCGCTCTAATCTCTGGGCATCTGGCTATGATGTT
GATGGAACTGACCAAACGTCGTTAGGCCAGTTTTCTGGTCGTGTTCAACAGACCTATAAACATTCTGTGC
CGCGTTTCTTTGTTCCTGAGCATGGCACTATGTTTACTCTTGCGCTTGTTCGTTTTCCGCCTACTGCGAC
TAAAGAGATTCAGTACCTTAACGCTAAAGGTGCTTTGACTTATACCGATATTGCTGGCGACCCTGTTTTG
TATGGCAACTTGCCGCCGCGTGAAATTTCTATGAAGGATGTTTTCCGTTCTGGTGATTCGTCTAAGAAGT
TTAAGATTGCTGAGGGTCAGTGGTATCGTTATGCGCCTTCGTATGTTTCTCCTGCTTATCACCTTCTTGA
AGGCTTCCCATTCATTCAGGAACCGCCTTCTGGTGATTTGCAAGAACGCGTACTTATTCGCCACCATGAT
TATGACCAGTGTTTCCAGTCCGTTCAGTTGTTGCAGTGGAATAGTCAGGTTAAATTTAATGTGACCGTTT
ATCGCAATCTGCCGACCACTCGCGATTCAATCATGACTTCGTGATAAAAGATTGAGTGTGAGGTTATAAC
GCCGAAGCGGTAAAAATTTTAATTTTTGCCGCTGAGGGGTTGACCAAGCGAAGCGCGGTAGGTTTTCTGC
TTAGGAGTTTAATCATGTTTCAGACTTTTATTTCTCGCCATAATTCAAACTTTTTTTCTGATAAGCTGGT
TCTCACTTCTGTTACTCCAGCTTCTTCGGCACCTGTTTTACAGACACCTAAAGCTACATCGTCAACGTTA
TATTTTGATAGTTTGACGGTTAATGCTGGTAATGGTGGTTTTCTTCATTGCATTCAGATGGATACATCTG
TCAACGCCGCTAATCAGGTTGTTTCTGTTGGTGCTGATATTGCTTTTGATGCCGACCCTAAATTTTTTGC
CTGTTTGGTTCGCTTTGAGTCTTCTTCGGTTCCGACTACCCTCCCGACTGCCTATGATGTTTATCCTTTG
AATGGTCGCCATGATGGTGGTTATTATACCGTCAAGGACTGTGTGACTATTGACGTCCTTCCCCGTACGC
CGGGCAATAACGTTTATGTTGGTTTCATGGTTTGGTCTAACTTTACCGCTACTAAATGCCGCGGATTGGT
TTCGCTGAATCAGGTTATTAAAGAGATTATTTGTCTCCAGCCACTTAAGTGAGGTGATTTATGTTTGGTG
CTATTGCTGGCGGTATTGCTTCTGCTCTTGCTGGTGGCGCCATGTCTAAATTGTTTGGAGGCGGTCAAAA
AGCCGCCTCCGGTGGCATTCAAGGTGATGTGCTTGCTACCGATAACAATACTGTAGGCATGGGTGATGCT
GGTATTAAATCTGCCATTCAAGGCTCTAATGTTCCTAACCCTGATGAGGCCGCCCCTAGTTTTGTTTCTG
GTGCTATGGCTAAAGCTGGTAAAGGACTTCTTGAAGGTACGTTGCAGGCTGGCACTTCTGCCGTTTCTGA
TAAGTTGCTTGATTTGGTTGGACTTGGTGGCAAGTCTGCCGCTGATAAAGGAAAGGATACTCGTGATTAT
CTTGCTGCTGCATTTCCTGAGCTTAATGCTTGGGAGCGTGCTGGTGCTGATGCTTCCTCTGCTGGTATGG
TTGACGCCGGATTTGAGAATCAAAAAGAGCTTACTAAAATGCAACTGGACAATCAGAAAGAGATTGCCGA
GATGCAAAATGAGACTCAAAAAGAGATTGCTGGCATTCAGTCGGCGACTTCACGCCAGAATACGAAAGAC
CAGGTATATGCACAAAATGAGATGCTTGCTTATCAACAGAAGGAGTCTACTGCTCGCGTTGCGTCTATTA
TGGAAAACACCAATCTTTCCAAGCAACAGCAGGTTTCCGAGATTATGCGCCAAATGCTTACTCAAGCTCA
AACGGCTGGTCAGTATTTTACCAATGACCAAATCAAAGAAATGACTCGCAAGGTTAGTGCTGAGGTTGAC
TTAGTTCATCAGCAAACGCAGAATCAGCGGTATGGCTCTTCTCATATTGGCGCTACTGCAAAGGATATTT
CTAATGTCGTCACTGATGCTGCTTCTGGTGTGGTTGATATTTTTCATGGTATTGATAAAGCTGTTGCCGA
TACTTGGAACAATTTCTGGAAAGACGGTAAAGCTGATGGTATTGGCTCTAATTTGTCTAGGAAATAACCG
TCAGGATTGACACCCTCCCAATTGTATGTTTTCATGCCTCCAAATCTTGGAGGCTTTTTTATGGTTCGTT
CTTATTACCCTTCTGAATGTCACGCTGATTATTTTGACTTTGAGCGTATCGAGGCTCTTAAACCTGCTAT
TGAGGCTTGTGGCATTTCTACTCTTTCTCAATCCCCAATGCTTGGCTTCCATAAGCAGATGGATAACCGC
ATCAAGCTCTTGGAAGAGATTCTGTCTTTTCGTATGCAGGGCGTTGAGTTCGATAATGGTGATATGTATG
TTGACGGCCATAAGGCTGCTTCTGACGTTCGTGATGAGTTTGTATCTGTTACTGAGAAGTTAATGGATGA
ATTGGCACAATGCTACAATGTGCTCCCCCAACTTGATATTAATAACACTATAGACCACCGCCCCGAAGGG
GACGAAAAATGGTTTTTAGAGAACGAGAAGACGGTTACGCAGTTTTGCCGCAAGCTGGCTGCTGAACGCC
CTCTTAAGGATATTCGCGATGAGTATAATTACCCCAAAAAGAAAGGTATTAAGGATGAGTGTTCAAGATT
GCTGGAGGCCTCCACTATGAAATCGCGTAGAGGCTTTGCTATTCAGCGTTTGATGAATGCAATGCGACAG
GCTCATGCTGATGGTTGGTTTATCGTTTTTGACACTCTCACGTTGGCTGACGACCGATTAGAGGCGTTTT
ATGATAATCCCAATGCTTTGCGTGACTATTTTCGTGATATTGGTCGTATGGTTCTTGCTGCCGAGGGTCG
CAAGGCTAATGATTCACACGCCGACTGCTATCAGTATTTTTGTGTGCCTGAGTATGGTACAGCTAATGGC
CGTCTTCATTTCCATGCGGTGCACTTTATGCGGACACTTCCTACAGGTAGCGTTGACCCTAATTTTGGTC
GTCGGGTACGCAATCGCCGCCAGTTAAATAGCTTGCAAAATACGTGGCCTTATGGTTACAGTATGCCCAT
CGCAGTTCGCTACACGCAGGACGCTTTTTCACGTTCTGGTTGGTTGTGGCCTGTTGATGCTAAAGGTGAG
CCGCTTAAAGCTACCAGTTATATGGCTGTTGGTTTCTATGTGGCTAAATACGTTAACAAAAAGTCAGATA
TGGACCTTGCTGCTAAAGGTCTAGGAGCTAAAGAATGGAACAACTCACTAAAAACCAAGCTGTCGCTACT
TCCCAAGAAGCTGTTCAGAATCAGAATGAGCCGCAACTTCGGGATGAAAATGCTCACAATGACAAATCTG
TCCACGGAGTGCTTAATCCAACTTACCAAGCTGGGTTACGACGCGACGCCGTTCAACCAGATATTGAAGC
AGAACGCAAAAAGAGAGATGAGATTGAGGCTGGGAAAAGTTACTGTAGCCGACGTTTTGGCGGCGCAACC
TGTGACGACAAATCTGCTCAAATTTATGCGCGCTTCGATAAAAATGATTGGCGTATCCAACCTGCA
>contig_09 len=6001 cov=3x
AGGACATCCCGATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAAAGTCCTACGTGATCTGA
GTTCAGACCGGAGTAATCCAGGTCGGTTTCTATCTACNTTCAAATTCCTCCCTGTACGAAAGGACAAGAG
AAATAAGGCCTACTTCACAAAGCGCCTTCCCCCGTAAATGATATCATCTCAACTTAGTATTATACCCACA
CCCACCCAAGAACAGGGTTTGTTAAGATGGCAGAGCCCGGTAATCGCATAAAACTTAAAACTTTACAGTC
AGAGGTTCAATTCCTCTTCTTAACAACATACCCATGGCCAACCTCCTACTCCTCATTGTACCCATTCTAA
TCGCAATGGCATTCCTAATGCTTACCGAACGAAAAATTCTAGGCTATATACAACTACGCAAAGGCCCCAA
CGTTGTAGGCCCCTACGGGCTACTACAACCCTTCGCTGACGCCATAAAACTCTTCACCAAAGAGCCCCTA
AAACCCGCCACATCTACCATCACCCTCTACATCACCGCCCCGACCTTAGCTCTCACCATCGCTCTTCTAC
TATGAACCCCCCTCCCCATACCCAACCCCCTGGTCAACCTCAACCTAGGCCTCCTATTTATTCTAGCCAC
CTCTAGCCTAGCCGTTTACTCAATCCTCTGATCAGGGTGAGCATCAAACTCAAACTACGCCCTGATCGGC
GCACTGCGAGCAGTAGCCCAAACAATCTCATATGAAGTCACCCTAGCCATCATTCTACTATCAACATTAC
TAATAAGTGGCTCCTTTAACCTCTCCACCCTTATCACAACACAAGAACACCTCTGATTACTCCTGCCATC
ATGACCCTTGGCCATAATATGATTTATCTCCACACTAGCAGAGACCAACCGAACCCCCTTCGACCTTGCC
GAAGGGGAGTCCGAACTAGTCTCAGGCTTCAACATCGAATACGCCGCAGGCCCCTTCGCCCTATTCTTCA
TAGCCGAATACACAAACATTATTATAATAAACACCCTCACCACTACAATCTTCCTAGGAACAACATATGA
CGCACTCTCCCCTGAACTCTACACAACATATTTTGTCACCAAGACCCTACTTCTAACCTCCCTGTTCTTA
TGAATTCGAACAGCATACCCCCGATTCCGCTACGACCAACTCATACACCTCCTATGAAAAAACTTCCTAC
CACTCACCCTAGCATTACTTATATGATATGTCTCCATACCCATTACAATCTCCAGCATTCCCCCTCAAAC
CTAAGAAATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGGAGCTTAAACCCCCTTATTTC
TAGGACTATGAGAATCGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCT
AAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTATACCCTTCCCGTACTA
ATTAATCCCCTGGCCCAACCCGTCATCTACTCTACCATCTTTGCAGGCACACTCATCACAGCGCTAAGCT
CGCACTGATTTTTTACCTGAGTAGGCCTAGAAATAAACATGCTAGCTTTTATTCCAGTTCTAACCAAAAA
AATAAACCCTCGTTCCACAGAAGCTGCCATCAAGTATTTCCTCACGCAAGCAACCGCATCCATAATCCTT
CTAATAGCTATCCTCTTCAACAATATACTCTCCGGACAATGAACCATAACCAATACTACCAATCAATACT
CATCATTAATAATCATAATAGCTATAGCAATAAAACTAGGAATAGCCCCCTTTCACTTCTGAGTCCCAGA
GGTTACCCAAGGCACCCCTCTGACATCCGGCCTGCTTCTTCTCACATGACAAAAACTAGCCCCCATCTCA
ATCATATACCAAATCTCTCCCTCACTAAACGTAAGCCTTCTCCTCACTCTCTCAATCTTATCCATCATAG
CAGGCAGTTGAGGTGGATTAAACCAAACCCAGCTACGCAAAATCTTAGCATACTCCTCAATTACCCACAT
AGGATGAATAATAGCAGTTCTACCGTACAACCCTAACATAACCATTCTTAATTTAACTATTTATATTATC
CTAACTACTACCGCATTCCTACTACTCAACTTAAACTCCAGCACCACGACCCTACTACTATCTCGCACCT
GAAACAAGCTAACATGACTAACACCCTTAATTCCATCCACCCTCCTCTCCCTAGGAGGCCTGCCCCCGCT
AACCGGCTTTTTGCCCAAATGGGCCATTATCGAAGAATTCACAAAAAACAATAGCCTCATCATCCCCACC
ATCATAGCCACCATCACCCTCCTTAACCTCTACTTCTACCTACGCCTAATCTACTCCACCTCAATCACAC
TACTCCCCATATCTAACAACGTAAAAATAAAATGACAGTTTGAACATACAAAACCCACCCCATTCCTCCC
CACACTCATCGCCCTTACCACGCTACTCCTACCTATCTCCCCTTTTATACTAATAATCTTATAGAAATTT
AGGTTAAATACAGACCAAGAGCCTTCAAAGCCCTCAGTAAGTTGCAATACTTAATTTCTGTAACAGCTAA
GGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTACTA
GACCAATGGGACTTAAACCCACAAACACTTAGTTAACAGCTAAGCACCCTAATCAACTGGCTTCAATCTA
CTTCTCCCGCCGCCGGGAAAAAAGGCGGGAGAAGCCCCGGCAGGTTTGAAGCTGCTTCTTCGAATTTGCA
ATTCAATATGAAAATCACCTCGGAGCTGGTAAAAAGAGGCCTAACCCCTGTCTTTAGATTTACAGTCCAA
TGCTTCACTCAGCCATTTTACCTCACCCCCACTGATGTTCGCCGACCGTTGACTATTCTCTACAAACCAC
AAAGACATTGGAACACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCC
TTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCAC
AGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCATCATAATCGGAGGCTTTGGCAACTGACTA
GTTCCCCTAATAATCGGTGCCCCCGATATGGCGTTTCCCCGCATAAACAACATAAGCTTCTGACTCTTAC
CTCCCTCTCTCCTACTCCTGCTCGCATCTGCTATAGTGGAGGCCGGAGCAGGAACAGGTTGAACAGTCTA
CCCTCCCTTAGCAGGGAACTACTCCCACCCTGGAGCCTCCGTAGACCTAACCATCTTCTCCTTACACCTA
GCAGGTGTCTCCTCTATCTTAGGGGCCATCAATTTCATCACAACAATTATCAATATAAAACCCCCTGCCA
TAACCCAATACCAAACGCCCCTCTTCGTCTGATCCGTCCTAATCACAGCAGTCCTACTTCTCCTATCTCT
CCCAGTCCTAGCTGCTGGCATCACTATACTACTAACAGACCGCAACCTCAACACCACCTTCTTCGACCCC
GCCGGAGGAGGAGACCCCATTCTATACCAACACCTATTCTGATTTTTCGGTCACCCTGAAGTTTATATTC
TTATCCTACCAGGCTTCGGAATAATCTCCCATATTGTAACTTACTACTCCGGAAAAAAAGAACCATTTGG
ATACATAGGTATGGTCTGAGCTATGATATCAATTGGCTTCCTAGGGTTTATCGTGTGAGCACACCATATA
TTTACAGTAGGAATAGACGTAGACACACGAGCATATTTCACCTCCGCTACCATAATCATCGCTATCCCCA
CCGGCGTCAAAGTATTTAGCTGACTCGCCACACTCCACGGAAGCAATATGAAATGATCTGCTGCAGTGCT
CTGAGCCCTAGGATTCATCTTTCTTTTCACCGTAGGTGGCCTGACTGGCATTGTATTAGCAAACTCATCA
CTAGACATCGTACTACACGACACGTACTACGTTGTAGCCCACTTCCACTATGTCCTATCAATAGGAGCTG
TATTTGCCATCATAGGAGGCTTCATTCACTGATTTCCCCTATTCTCAGGCTACACCCTAGACCAAACCTA
CGCCAAAATCCATTTCACTATCATATTCATCGGCGTAAATCTAACTTTCTTCCCACAACACTTTCTCGGC
CTATCCGGAATGCCCCGACGTTACTCGGACTACCCCGATGCATACACCACATGAAACATCCTATCATCTG
TAGGCTCATTCATTTCTCTAACAGCAGTAATATTAATAATTTTCATGATTTGAGAAGCCTTCGCTTCGAA
GCGAAAAGTCCTAATAGTAGAAGAACCCTCCATAAACCTGGAGTGACTATATGGATGCCCCCCACCCTAC
CACACATTCGAAGAACCCGTATACATAAAATCTAGACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGG
TTTCAAGCCAACCCCATGGCCTCCATGACTTTTTCAAAAAGGTATTAGAAAAACCATTTCATAACTTTGT
CAAAGTTAAATTATAGGCTAAATCCTATATATCTTAATGGCACATGCAGCGCAAGTAGGTCTACAAGACG
CTACTTCCCCTATCATAGAAGAGCTTATCACCTTTCATGATCACGCCCTCATAATCATTTTCCTTATCTG
CTTCCTAGTCCTGTATGCCCTTTTCCTAACACTCACAACAAAACTAACTAATACTAACATCTCAGACGCT
CAGGAAATAGAAACCGTCTGAACTATCCTGCCCGCCATCATCCTAGTCCTCATCGCCCTCCCATCCCTAC
GCATCCTTTACATAACAGACGAGGTCAACGATCCCTCCCTTACCATCAAATCAATTGGCCACCAATGGTA
CTGAACCTACGAGTACACCGACTACGGCGGACTAATCTTCAACTCCTACATACTTCCCCCATTATTCCTA
GAACCAGGCGACCTGCGACTCCTTGACGTTGACAATCGAGTAGTACTCCCGATTGAAGCCCCCATTCGTA
TAATAATTACATCACAAGACGTCTTGCACTCATGAGCTGTCCCCACATTAGGCTTAAAAACAGATGCAAT
TCCCGGACGTCTAAACCAAACCACTTTCACCGCTACACGACCGGGGGTATACTACGGTCAATGCTCTGAA
ATCTGTGGAGCAAACCACAGTTTCATGCCCATCGTCCTAGAATTAATTCCCCTAAAAATCTTTGAAATAG
GGCCCGTATTTACCCTATAGCACCCCCTCTACCCCCTCTAGAGCCCACTGTAAAGCTAACTTAGCATTAA
CCTTTTAAGTTAAAGATTAAGAGAACCAACACCTCTTTACAGTGAAATGCCCCAACTAAATACTACCGTA
TGGCCCACCATAATTACCCCCATACTCCTTACACTATTCCTCATCACCCAACTAAAAATATTAAACACAA
ACTACCACCTACCTCCCTCACCAAAGCCCATAAAAATAAAAAATTATAACAAACCCTGAGAACCAAAATG
AACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCCTAGGCCTACCCGCCGCAGTACTGATCATTC
TATTTCCCCCTCTATTGATCCCCACCTCCAAATATCTCATCAACAACCGACTAATCACCACCCAACAATG
ACTAATCAAACTAACCTCAAAACAAATGATAACCATACACAACACTAAAGGACGAACCTGATCTCTTATA
CTAGTATCCTTAATCATTTTTATTGCCACAACTAACCTCCTCGGACTCCTGCCTCACTCATTTACACCAA
CCACCCAACTATCTATAAACCTAGCCATGGCCATCCCCTTATGAGCGGGCACAGTGATTATAGGCTTTCG
CTCTAAGATTAAAAATGCCCTAGCCCACTTCTTACCACAAGGCACACCTACACCCCTTATCCCCATACTA
GTTATTATCGAAACCATCAGCCTACTCATTCAACCAATAGCCCTGGCCGTA
>contig_10 len=6001 cov=4x
CGCCCAGACGATTTGTGGGCTAGGCAACCGCAAGAACGACCTGTTCAATACAGTGCTGGCCCGGGACGGC
ATCCGGCCTTACCCGGGGTCGCGACGTTGGGTGGATAGGCTGCACGAGAGTGGTATGGCGATGGCGGTGG
TCTCATCGTCACGCAATGCCGCAGCAGTGCTCAAGGCGGCGGGGATGGTCGAGGATTTCTCAGTGCTGGT
GGACGGCAACCGGTCGAAGGCTGAAAGATTACCCGGCAAGCCGGCCCCCGATACGTACCTGAGAGGAGCG
GAGCTACTCGGTGTGCCTGCTGAACAGTGCGTCGTCGTCGAGGACGCGGTGTCGGGAGTACGCGCGGGAG
CTGCCGGAGGATTTGGAATGGTGCTTGGCGTCAACCGTGGTGTCGGGGCCGACCGGCTGCGTGAGGCCGG
CGCCGACCGGGTGGTGGATGACTTAGATGAGATGGTTGAGGAGATGGCATGACGCGGACGGGCACTGACC
CGATGGATCGATGGCGCTTTCCCGACGATCCGTGGGCGATTGTTGAGACTGCCCCCAGCACAGATGATCT
AGGAACTACTGAGACGCTGTTTTCAATTGGTAATGGCTATCTCGGTATGCGTGGTAATCCGTCGGAAGGC
CGCGACTCTTTCAGTCACGGCACCTACATCAACGGGTTCCATGAAATTTGGGACATCCATCACGCTGAGA
ATGCCTACGGTTTCGCGCGTACCGGTCAGACGATTGTCAACGTTCCGGACGCCAAGCTTATGAAGCTCTA
CGTGGATGACGAACCGCTCCTGCTCAGCATCAACGAGATCCAGTCGTATAAACGGTGGATCGACTTCCGC
GAGGGTGTCCTGCGTCGGGAGCTGGTATGGCGCACTCCGGCGGGCAAGTTGCTGCGGGTGGCAACCTCGC
GGATGGTCTCCTTCACTCACCGGCACATGGCTCTAATTAGCATCGAGGTGACGATGCTCGAGGGCGACGC
CCCGATCGTCATTAGCTCGCAAATCCTTAACCGTCAGGACGGCATGGACGAGTACCACGCTCGACCTGAC
GACATCGAGAAGACTGCCGACCCTCGCCAAGCTCGTAAATTTGCCGACAAGGTACTCATTCCTGAAAAAG
ATTGGCACTCTGACAGGCGCATGATTTTGGGATTCCGCACCGCACGTTCGGGGATGACCCTGGCTGTAGG
CGCCGATCACGAGATCATCACGGAAAATGAGTACGCGGCCCTTATCGATACCGAGCCCGGGATGGGAAAA
CGGGTGTATCGTGTTGAGGCCACCCAAGGCCGACCAATTCGCATCGATAAGGCGGTTGCTTATCACACTT
CTCGCGGCGTGCCGGTACATGAACTGTTTGACCGAGTGCGCCGCAGCTTAGATCGAGTGCGTGAACAGGG
GCACAACGTCTACTACGACGAACAACGTGCATGGCTTGATGATTACTGGGCAACGGCTGATGTCGAGGTC
GAAGGTGCCCCGACCGGTATTCAGCAGGCTGTCAGGTGGAACCTTTTCCAGATTGCTCAGGCATCAGCCC
GTGCAGATCAACTTGGCATTCCGGCAAAGGGTGTAACCGGGTCAGGCTATGAAGGCCACTACTTTTGGGA
CACTGAGGTTTATGTCATCCCGATGCTGACCTACACTCATCCAAGAATCGCTGAGAATGCGCTGAGATTT
CGGGTGAATACCCTTCCGCAAGCTCGACGTCGGGCTAAGGAGTTGTCTGAACGAGGTGCCCTTTTCCCGT
GGCGAACAATCACCGGTGAGGAGGCGTCAGCTTATTACGCGGCGGGTACTGCTCAATACCACATCAATGC
CGATATTGTCCATGCGATTATGAACCACGCCCGGGCGACCGAAGATAAGACCTTCCTCTTCCGGGACGCG
GCACCGGTCCTCGTCGAGACCGCACGAATGTGGGCTGATCTTGGGTTCTGGCGCATCAACGGTGGTCGAG
AGTTCCACATCCACGGGGTTACCGGGCCCGACGAGTACACCACTGTCGTTAATAATAACCTCTACACCAA
CGTTATGGCGAGGGCCAACCTCATCGATGCCGCCGGTGTTATTAGGCGTATGCGTGACGAGGACCCACTG
TGGTACGAGCACCTGTGTTCAGAGTTGGACCTCACCGAAGATGAAGTTGGTGGCTGGGAGGAGTGTGCGG
CGGGGATGGTGATCCCGTTCGACGACACCTTCGGTATCCATCCGCAGGATGACCAATTCCTTTCACGTGA
GCTTTGGGACCTCAAGAACACCCCTGGTAACAAGCGGCCGCTGTTGCTGCATTATCACCCGCTGGTGATC
TACCGATTCCAGGTACTCAAGCAGGCAGACGTTGTGCTGGCTTTGTTCCTGCAGGGAGATCAATTCACCC
AGGCTGTCAAGCTTGCTGACTTTGAGTATTACGACCCCATCACCACCGGTGATTCGAGTCTGTCGGCGGT
AGTCCAGTCGGTGATGGCGGCAGAAGTTGGCCATCAGGAGATGGCGATGGGCTATTTTCTCGACGCGCTC
TTCTGCGACCTGGCTGATACCCATCACAACGCTTCGGACGGAGTTCACGTTGCCTCGACCGGTGGAGTGT
GGGGATGCCTCGTGCACGGTTTCGTTGGCATGCGCAATGACCGGGAGAACTTGCGTTTTGATCCGAGACT
TCCAGCCCAATGGACGTCGATCAAACACCACATGCTCATTGGCGACTCTCACATGCTCGTTTTCCTGGAA
CGTGACGCCATTACGTTCCAGATTCTGTCGGGCCATGACCGCGACGTGACAGTGCGCGGTGAGCTCTACC
ACATTGGGGTTGAGCCGGTGAGGGTGCCGTTGTCCGATCAGGGGCCGTTGCGTCCTAGCCTGCGGGTTAC
CCATCCGATCTCGGGTTTGCGTCGAGCTGACGGTTCTCTTATCACTGCAGAAGTTCCCGGCAGCATTGCT
GAGACGATTGTTTCTTCTCCGATCTCGACGTCGTGAGACTGATAACGACGGGAGGAGGCCGGCACCATTC
GCGGTGTCGGCCTCCTCCATATCGTGGAGGTACGGTGCGGTCCGGTTACCCGGTTTGCCTCGAACTCACT
TCACAACAGTGATCACCGGGGTACCTGTAGTGACGTCGTCAGCAGGAATGCCGGTCACGTTGGAGTAGTC
GAACCCGTTGGTCACTAACACCGGGGTGATCTTGGAGTATCCGGCTTCGTCAATGACGGTGGGATCAAAA
ACGACCAGGGGCTGTCCAGCACTGACGTGATCGCCAGCCTTGACCTTGACGTCAAACCCTTTGCCACCGA
GCTCGACGGTATCGATGCCGACATGGATGAGTAGCTCGACACCCGAATCCAGTTTGAGCCCGAAAGCGTG
GCCGGACTCCGGGGCGGTGACAACGGTGGCGTCGGCTGGGCAGACAACAGTGTTGCTGGTGGGCTCGATG
GCGATGCCATCACCGACGGCACCAGTGGAGAAAACCGGATCGGGGATTTCAGAGATTGGCACGGCACGTC
CCTCAAGAGGAGAGGTGACGGTGACGGCCTCGGCGGCAACATCGGCTGCAGTCGGGGCCTCGGCCCCCTC
GACTGTTGCGGTCTCGGTTATGGATTCGGTGGCGGCAGCTGGGGCAGCGGCGGTTGATTCTGCCTCGTGA
GCGCTGGCCTGACGGCGACGGGCCTTGGCTGCAGCCTTCTCGTCGTTATCACGATAGTCGCGCACGATAA
CCAGAGCCATCGCGGTGAAGAACGAAGCTGCGATCGCTACCGAGTACAGCCCGATCTTGTCGAAGGCCGG
GATGGTCAGCAGGGAAGTGAAGACAAAGGCCGGGGTGTGCAGACCGCCGCCCGCACCGATGATGACACCA
CCGACGAGACAGCCGACCAACATGGATGGGTAAATCTTTTTGTAACGCAGGTGGATGCCGTAGAGGCTGG
GCTCGGAGATACCACCGAGCAGACCGGCGGCTAGGGCGCCAATCGCGGTCTGGCGCATTGTCACGTCCTT
GTCACGCACAGACAGGAAGAGAACGCCAGCGGTAGCACCGAAACAGGCGAAGTTCCACGCACCCATCGGG
CCCTGGATGTAGTCGTAACCCAGAGTCTGGATGTTGAGCAGCATGACAGCGTTGAGCGGCCAGTGTAGGC
CCAACGGAACCATGAATGGGTACAGCAGGGGGATGAGGACCGCGAAGATGAAAGGGCTGACGTCGTTGAC
CGCCTTGAGCCCGTTGGCCAGCCCAGCACCGGCGTAGACACCGATCGGGCCGATGAGGAAGGCCGTCAGC
GGGATCATGATGAGCATGGTGAGGAAGGGTACGAAGATCATCTGGACATTGGGGGAGATGATCTTCTTGA
GCCATTTGTACAGCGGGGCCAGGATTGCAGCCATGATGAGCGGCGGGAACACCTGAGAGCCGTAATTAAA
GATCGTCAGCGGTAGGCCGAAAACGTCGACCACGTTGACTTTGAATCCGCCGAGGGTCAGTGGATGGGCA
GTCTTGCCGAGCTCGGTGAACGCAGGCAGCATGACGACAGCCATGACGGAGAACCCAACCCACGGGTCGG
CATTCAGACGCTTAGTGGCGTTGTAGGCGACCAACAGTGGCAGGAAGTAGAAGACACACTGCCAAGCCAA
GTTGACGAATTGCCAGCTCGGGGGCAGGGTGACGCGCGAGTCGGCCCAGTTTCCGATGACGCCCAGAGAA
CCCATGAGGGCCATGAACGTGATGAACAGGGAGGCGCCGAGCAGGGCACCCATGATCGGCCGGAAAGAGT
CGGACAGGTACTCGAAGAGGTTGTCGAGCCACGAGAACTTGCCTCGGGGACCCTTGGAGCGATGCTCGGC
CTTGACGTCCGCGTCCGACTGATTGCCGCTGCTGGAAGAAAGCCCTTTCATGGCGGGTAGGTCCATGATG
TCGTTGTAGACGCCGGCAACGGCGCCTCCGATGACGATTTGGTAGCGATTGCCCGACTGAGGAACGGCTC
CCATGACGCCGGGGATCGCTTCGACCTCAGCCTTGTCGATAGGCGAGGAATCGCGAAGTTGGAACCTCAG
ACGAGTGGCGCAGTGAGTGAGACCGACGATGTTGTCGGGGCCGCCCACCGCATCCACGATCTGGGACGCG
GGTGATGACATGTATGTGGTCCTCCTGTAATATTGCCGCCCAGGGCAGCGAGGCAGGATTCTAATCGTTC
CAGGATGCTGTTACCACACCGGGGTAGGTCCCGGCTTGATCCACCGATTTGGTGCTGCGGTCGCGAATCT
CCTGGGAGTTGTGAGATAGCGTAATGTATACGTCGCAAACATTATCGAGCGCGTCGCTGTGCTGAGACTG
CCGCGCGGATTCTCGGGGGCACAGTCGAGGCTCCCCGGGGCAGGAGGGTCGAACCTTCTAACGCGGGCCT
CACACAGGGAATCCGTTTTCTGTGGCGCGAGGCATGCAATCGGACAAGATGGCCTCGATGCGTCGGGCCG
ACAGTGCGGTCACGTGATCGGGTATCAGCTGCGGCTCGACCCAATGCAGGCGACTGATTTCCGAGTCGGC
TGGGTGCATGACCTTGACGAGGCTTGGCTCGGCTACACCAGCGTCATAGATGAACTCGATAGCGTCAGAC
CATCCCAACGCCGGAGGCATCCAGTCGATGAGACTCGGCGTGCCCGGAATGTCAAGGGGGGTTCCTAACT
CCTCCCGGCATTCCCGGACGGCTCCGGCATGCGGGGATTCTATAGGTTCGACGACGCCCCCGGGCAGTTC
CCAATCCGGCTTATATGTCGTCTCGCACATGAGCACCCGCCCGGCATCATCACGGAGGAGGACGTGGGAG
ATGACCCGTTTGGTGGGTAGAACAGAATCGAGCATCCCGGTGTGGCCGGTGCGGGTATAGATCTCGTCGG
TTGCCAGTCGGGAATAGATGTATTCGTCTATCCAGCCGTTTTCACGCACGCGCGCCTGGCGCAGCCGACC
CTCGCGACGATACCCAGCATGGTGCAAGGCGTGGATGACCTCATCGTCATC
>contig_11 len=15000 cov=115x
TCTCTATTTTCTAACGTTAGTATTGCATCAAACAATTCATCTAACGCTGCACCTCGTAATTTTTCAATTT
GCATATCGAAACAACCCCTATAATTTTATTCTTTATGTCATTTTAACGCATTAATCTATTAAAGTGTAGC
GATTTATATTTTATTAAATCTGAATCGTTACTTTATTTAATTTTATGCTAATCCAGCGCGTTCGAAAATA
GTGTCAACTTGATTCAAATGATGTTTAGGATCGAAACATTCATCCAATTCTTCTTTTGTTAAAACACTTG
TAATAGACTCATCTTGTTCGATTAATTCACGGAACGGTGTTTTCGTTTCCCAAGATATCATCGCTTTTGG
TTGTACTTTGTCGTATGCTTCTTCACGAACCATACCTTTATTAATTAATGCTAATAAGACACGTTGTGAG
AAAATCAGACCAAATGTTTTATCTATGTTATTACGCATATTATCTTCAAATACAGTTAAACGGTCCACAA
TATTTGTGAAACGATTCAATGCATAATCTAGTGCTATTGTAACATCTGGTAACATAATACGCTCAGCAGA
AGAATGAGAAATATCTCTTTCATGCCATAATGGCACATTCTCATAAGCTGTAGTAATATAACCACGAATG
ACTCTTGAAATACCTGTGATATTTTCAGAACCAATTGGATTTCGTTTATGAGGCATTGCAGATGAACCTT
TTTGGCCTTTTGCAAATGCTTCTTCAACTTCTCTCGTTTCGGTTTTTTGAAGGTTACGTATTTCAACGGC
AAATTTTTCTAGTGATGTCGCGATTAATGCTAATGTCGCAATATAGTATGCATGTCGATCGCGTTGCAAT
GTTTGCGTTGATACAGGCGCTGTGCCAATACCTAAATGTTTACACACATAACTTTCTATTTCAGGAGGAA
TGTTAGCAAAAGTACCTACTGCACCACTCATTTTCCCTACTTCAATTTCTTCTCTTACTTGTTTGAAACG
TTGTAAGTTACGTTGCATTTCCGTGTACCACAATGCCATTTTGACACCAAATGTAGTTGGTTCTGCATGC
ACTCCATGTGTACGTCCCATCATCAATGTATATTTATAATTTTTTGCTTTTTCAGCTAAAACGTCGATAA
ATCTTTCTAAATCTTTTTCAATAATGTCATTTGCTTGTTTAATAACGAAACTTAAAGCTGTATCTACAAC
ATCAGTAGAAGTTAAACCATAATGTACCCACTTACGTTCTTCACCTAGCGTTTCAGAAACTTGTCTAGTA
AAGGCTACAACATCATGGCGCGTTTCTTGTTCAATTTCTTGTGCACGTTCGACATTTACCTTTGCGTTTT
GACGAATTTTTTGTACGTCAGCTTTCGGTATATGTCCTAATTCACTCCATGCTTCACATGCTAAAATTTC
CACTTCTAACCATGCTTCATAGCGATTTTGATCAGTCCAAATATTAGACATTTCTTCTCTAGAATAGCGT
TCAATCATTTTTTTAACTCCTATCGTATGTATACTCGAATTTTCAAAAACAGAACATTTTTTAATAAAAA
TAAACAAAAAGTTCGTTTTTAACTCTACTCTTCACAGTTTAACAGAATTGGTTATAAATGTATAAAAAAT
TTCCGTATATTCTCATTCTATTTTATTTTTTACATTAATACAATTATCTATGGCAGATTAATACTGATTT
CTAATTGTAAATGTTCGCAATACTTTAATCTAATTGGTATTTATAATGATCCATTGAGAATATATTTTTG
CATTTTAGTTTTTCTCCCCTCAACTACACTATCCAACTTTTTATTACAATTAGGACAATCAATAGCTTTA
GTAATACCAGAATTTATTCACTTTACTTCTCCTATTGTATGGAATTTTCAACAATCTACTTTTCCTAGTT
TCATAATATTTCTTCTTACATTTCATAATATGTTTTAAAATTATTCGCTTCTAATATAACAAATTTTCTA
CTCACAAAAAGAATATATATTCTATTGCATAGGTGTGGCATTTTAAATATAAGACATTAATATTGAACTC
CAGAATCGTATAAATACGAAACTGATTATAGAAAAATACGCAGAGATTTTTTTAATAGTAAATTTTTATA
GAAAGACATAAATTGACTACGCTTCTTTCAAAATCGACTAATATATTTTCATGAAACTTTCGATCATATT
AATTATTGAATGAGATGTTTGCAATTTAGTAAATAAGCTACTAAATTATTACATAAACTGATAATGTGTA
AAAATTAATTATTTTTTAAACATAGATTAAAATTTTATATTGAAACTCTTAACTAAATTTAATATAATTA
TTTTATAATAAAATAAAGGAGTTTTTATATGAAAAGAAACTTTCCAAAATTAATAGCTTTAAGTTTAATA
TTTTCTTTAAGTGTTACCCCTATTGCAAACGCTGAGAGCAATTCAAATATCAAAGCAAAAGATAAGAAGC
ATGTGCAAGTTAATGTCGAGGACAAGAGTGTTCCTACTGACGTACGTAATTTAGCGCAAAAAGACTATCT
TTCTTATGTTACTTCACTTGATAAAATCTATAATAAAGAGAAAGCCTCGTATACACTTGGTGAACCTTTT
AAAATTTATAAATTCAACAAAAAGAGTGATGGTAATTATTATTTTCCAGTTTTGAACACTGAAGGTAATA
TCGATTATATAGTTACAATTTCACCAAAAATCACTAAGTATTCTAGTTCATCTTCTAAATACACAATCAA
TGTTTCTCCATTTTTATCTAAAGTCTTAAATCAATATAAAGATCAACAAATTACTATTTTGACAAACTCA
AAAGGTTACTATGTGGTTACTCAAAATCATAAAGCTAAACTAGTACTAAAAACACCACGTCTAGAAGATA
AAAAGTTAAAAAAGACTGAGAGTATACCTACTGGCAATAATGTAACGCAACTAAAGCAAAAAGCTTCTGT
CACTATGCCAACTAGTCAATTTAAAAGTAATAATTATACTTACAACGAGCAATATATAAATAAACTTGAG
AATTTCAAGATAAGAGAAACTCAAGGTAACAATGGTTGGTGTGCAGGCTATACGATGTCTGAATTACTAA
ATGCGACTTATAATACTAACAAATATCATGCAGAAGCTGTGATGAGATTTTTACACCCAAACTTACAAGG
ACAACGATTCCAATTTACTGGTTTAACGCCAAGAGAGATGATTTATTTTGGACAAACTCAGGGAAGAAGT
CCACAACTTCTTAATAGAATGACTACATATAATGAAGTTGACAATCTTACTAAAAATAATAAAGGAATTG
CTGTTTTAGGTAGTAGAGTAGAATCAAGAAATGGTATGCACGCAGGGCATGCCATGGCAGTTGTTGGGAA
CGCTAAATTAGATAATGGACAAGAGGTTATAATTATATGGAATCCATGGGATAATGGTTTTATGACTCAA
GATGCTAAAAACAATGTTATTCCAGTATCCAATGGCGACCATTATCGATGGTATTCATCTATTTATGGTT
ATTAAAATTTGTATTTTAGAAAGGTAAATAATATTATGGAACAAATTGAATTATTTAGTATTGATAAATT
CAAATGTAATTCAGAGGCTAAGTATTATCTTAATATTATTGAGGGAGAATGGCATCCTCAAGATTTAAAT
GATAGTCCTTTAAAATTTATTCTCAGTACCTCAGACGATTCTGATTACATTTGCAAATATATAAATACAG
AACACAAACAACTCACATTATATAATAAAAATAATAGCTCAATTGTTATTGAAATATTTATACCAAATGA
TAATAAAATACTACTAACAATTATGAACACAGAAGCTTTGGGAACTTCTCCTAGAATGACTTTCATTAAG
CATAAGTCATAAGACTATCTAGATATGATTAAGTAAACTATTAATCTTATTAAAGGTCTCTAGGAAGAAA
AAATTTCATGCGATGAGATTGTTTTAGAAAAAGCGACGATATATATAATCTAACAATAATTTTTAAAATA
TTAAGAACCCCCAGCATATGCTTTATATAGCTTGCTGGGGGTTTCTCATACAATATGGATCGAGCATCTT
GCATGTATAAACAAATAGGCTATTCAAAGTAAAATTACTCAAACTTTGAATAGCCTTTCCTCTTCAAATA
GTCTTTATTTAAATTGGATTATGTTATTGGATTTTCTTACCTACAAAAACGATTTCTTGATGTGATACTT
CAACCGTGCCATCTTTTTGCTCTTGAAACACGGGTTTTTCAATACGTTTTACCGGCATATACCCTTGCAT
TTTCATTCGTGACAAACAATCAGCAATTGTTTCATTCTCATTTACTTTGAATTTCATGTTCTACTACACT
TTCCAGTTCTTCTGGTTCATAATCTTTAAGCTTTCTACGACGAATGCCTTTAGTCCAGAATCCACCATGT
ATGTTTCGAGGTTCATACGCAATGATAAATGCTTTCGGATCTAAATTTTTTATCGTATCCATCAATTTAC
GTTCATATTTTCTTGGTGTTAAAATTTGCATCACCATACGACTACCATCTCTACCAAACGCAGCATAGTG
CGTAACGCCATATCCTAAATTTCGAAGTTCATTCGGTAAATCTAACTCATATTCTGCTGAAGTTACATTT
ACAACTGTATAACCTAATGCCAGTTTTTCTTCTATTTTCATACCAACAATGATACCTATTGAAAAACCAA
ATGCGTAGGCAATAATATTTTGAATATGGTCTAAATTAGACATAACCAAACCTAAACCAACGATATAAAC
TAATACTTCTAAAAAACTAACTGATGCAGCAATATAACGATAACCTTTCAACGTTAAAATTGTTCGCATC
GTTAAAAACGTTACATAACAAACGTTAATGATAAATATAGTTAGTACCATTAACCATGGATTTTCTGTTA
CGAATGACATAAGGCGCCCCACTCTCTTTCAATTACACTTTAAATACCAAACTGCTTCAACCTAACCTTA
ATGAATGCATATTAATATAATATAGTACGAAAACTCTCATAAATAAATAGCGAATTTGTCATAAATAAAA
ATTTATTCTCAAATGAACATACTTACTTTAGTTCATACACTTTCAATTTATACACTGAAGTAGTCTAGCT
GAATGTGGTTTAATGTTATTAATTCTAACTCATATTTAATTTAAGTCACACGTTAGAAGAAAAAATTAAA
TTAGGATTTTGGCCACGTGTATCTTGTATATGCAAGTTCACGTTTGTGTGCATTTCGTATATAATGATTT
TCAATTACTTTTTGTTCTTCTGGCGTAACTGGCTTACCTTCTAAATAATTATCAATCGCCTCATAAGTTA
CACCTAATGCATCTTCATCTGGAAGCTGTGGTTTATCATCTTCTAAATCAGCAGTTGGCGTTTTTTCATA
TAATTCCTTTGGCGCACCAAGATACGCTAATAATTGACGACCTTGTCGTTTATTCAAACCAAATATAGGT
GCGATATCTGCAGCACCATCACCGTACTTCGTATAAAACCCAGTTATATTTTCAGCTGAATGATCTGTTC
CTACTACAATACCTTGTCGGTTTGAAGCAATTGAAAATTGTACTTTCATACGTTCACGCGCTTTTTCATT
TCCTTTTTGGAAATCTGTAAGAACAATACCGGCTTCTTTTAATGATTGCACACTTTGATCAACTGCAGGC
TTAATATTGACTGTTACTATTTCATCTGGTTCAATGAATCGCAAAGCTTGCTCAACTTCATCAGCATCTT
TTTGAACTCCATAAGGTAATTTAACTGCAATAAACGTACAATCAATGCCTTCTTCACGTAATTCGTTAAC
AGACATTTGTACTAGTTTTCCAACTAATGTAGAATCCTGTCCTCCTGAAATACCTAACACTAAAGATTTT
ATAAATGAATGTGATTGTACATAATTTTTTATAAATTGCTTTAATTCCATAATTTCTTCAGCACTATCGA
TACGCTTTTTCACTTTCATTTCTTGTACAATAACGTCTTGTAATTTACTCATTATCTTCTTCCATCTCCT
TAACGTGTTCCGCAACTTCAAAAATACGTTTATGTTTATTATCCCAACATGCCTTGCTTAAATCGACTGG
ATATTCTTGTGGATTCAGGAAACGCTTATTTTCATCCCAAATAGATTGTAATCCTAGTGCTAAATATTCA
CGTGATTCATCTTCTGTTGGCATTTGATATACTAATTTACCATTTTCATAAATATTATGATGCAAATCAA
TGGCTTCGAAAGATTTTATAAATTTCATTTTATAAGTATGCACTGGATGGAATAATTTTAAAGGTTGTTC
ATCGTATGGATTTTCATTTTCCAAAGTAATATAATCGCCTTCTGCCTTACCTGTTTTCTTGTTTATAATG
CGATATACATTTTTCTTACCTGGCGTCGTAACCTTTTCAGCGTTATTTGATAATTTAATACGATCACTAT
ATGAACCATCTTCATTTTCAATAGCTACAAGTTTATATACTGCACCTAATGCTGGTTGATCGTATCCTGT
AATCAGCTTTGTACCAACGCCCCAAGAATCTACTTTTGCACCTTGTGCTTTCAAACTCGTAATCGTTTCT
TCATCCAAATCATTAGACGCGATAATTTTAGTTTCAGTAAATCCTGCTTCATCAAGCATACGTCTTGCCT
CTTTAGATAAATAAGCGATATCTCCAGAATCTAATCGAATACCTACAAAGTTAATTTTGTCACCTAATTC
TTTTGCAACTTTTATTGCATTTGGCACGCCAGATTTTAAAGTATGGAATGTATCTACTAGGAACACACAA
TTTTTATGTCTTTCAGCATATTTTTTGAAGGCAACATATTCGTCTCCATAAGTTTGGACAAATGCATGTG
CATGTGTACCAGACACAGGTATACCAAATAATTTCCCCGCCCTAACATTACTTGTAGAATCAAAGCCCCC
GATGTAAGCAGCTCTAGCGCCCCACAATGCTGCATCAATTTCTTGCGCACGACGTGTACCAAACTCCATT
AATTTATCATTTGATGCAATTTGACGAATTCTGCTAGCCTTTGTTGTAATTAATGTATGGAAATTTACAA
TGTTTAATAAAATTGTTTCTATTAATTGCGCTTGAATCAATGGTGCTTCTACGCGTAACAATGGTTCGTT
ACCAAAGCATAATTCGCCTTCTTGCATCGAACGGATGCTGCCTGTGAATTTTAAATCTTTTAAATATGAT
AAGAAATCATCCTTGTAGCCAATAGACTTTAAATATTCCAAATCAGATTCTGAAAATCCAAAATGTTCTA
TAAAATCAATGACGCGTTTTAAACCATTAAAAACAGCATAGCCACTATTAAATGGCATTTTTCTAAAATA
CAAATCAAATACAGCCATTTTTTCATGAATATTATCATTCCAATAACTTTCAGCCATATTTATTTGATAT
AAGTCATTATGTAACATTAAACTGTCGTCTTCTAATTGGTACACTTGTATCTCTCCAATCGACCTAAATA
TTTTCTTACATTTTATCATAATTCATTTTTTTATATACATAAGAGCCCCTTAATTTCCATACTTTTAATT
AAAATCAACCAACAATTTAATGACATATACATAATTTTTAAGAGTATTTTAATAATGTAGACTATAATAT
AAAGCGAGGTGTTGTTAATGTTATTTAAAGAGGCTCAAGCTTTCATAGAAAACATGTATAAAGAGTGTCA
TTATGAAACGCAAATTATCAATAAACGTTTACATGACATTGAACTAGAAATAAAAGAAACTGGGACATAT
ACACATACAGAAGAAGAACTTATTTATGGTGCTAAAATGGCTTGGCGTAATTCAAATCGTTGCATTGGTC
GTTTATTTTGGGATTCGTTAAATGTCATTGATGCAAGAGATGTTACTGACGAAGCATCGTTCTTATCATC
AATTACTTATCATATTACACAGGCTACAAATGAAGGTAAATTAAAGCCGTATATTACTATATATGCTCCA
AAGGATGGACCTAAAATTTTCAACAATCAATTAATTCGCTATGCTGGCTATGACAATTGTGGTGATCCTG
CTGAAAAAGAAGTTACACGCTTAGCAAATCACTTAGGTTGGAAAGGAAAAGGTACTAATTTTGACGTGTT
ACCACTGATTTACCAATTACCTAATGAGTCAGTTAAATTTTACGAATATCCTACTTCATTAATTAAAGAA
GTACCTATTGAACATAATCATTATCCAAAATTAAGAAAATTGAACTTAAAATGGTATGCAGTCCCTATCA
TTTCCAATATGGACTTAAAAATCGGTGGCATTGTATATCCAACTGCACCCTTTAACGGTTGGTATATGGT
AACTGAAATTGGCGTACGTAACTTTATTGATGATTACCGTTACAATTTACTAGAAAAAGTTGCAGATGCG
TTTGAATTTGATACACTTAAAAATAATTCATTTAATAAAGATCGAGCACTTGTTGAATTGAACTATGCTG
TGTATCATTCCTTTAAAAAAGAAGGCGTATCAATTGTCGATCATTTGACCGCTGCAAAGCAATTCGAACT
ATTCGAACGTAACGAAGCACAACAAGGTCGTCAAGTTACCGGAAAATGGTCTTGGCTAGCACCGCCATTA
TCTCCAACATTGACGTCAAATTATCATCACGGATATGACAATACAGTAAAAGATCCAAACTTTTTCTATA
AAAAGAAAGAATCAAATGCTAACCAGTGCCCTTTCCATCATTAAGTTAGTAGAGGTGTAGCATATGCAAT
TATATTATTTAGGCCCTAAAGGTACATTTTCATATTTGGCTTGTCGTCAGTATTTTTCTGAAAATGAAGC
AACTTTTCAACCAAAGTCAAATTTATTTGAAGTCATCAAGGCAGTTGCAGATGATGACACATCGATTGGT
GTTGTCCCTATTGAAAATTCTATAGAAGGAACAATCAACATTGTGGCAGATGCACTTGCACAGCAAGATG
TCTTTGCACACGGTGAAATACGTTTAGATATTAACTTTGCATTATATGGTAATGGAACAGATTCCATTTC
TGATATTAAAAAAGTGTATTCTATTGCACCAGCAATCAGTCAAACAACAAATTACATACATCAACATCAA
TTTGACTATGATTATGTCGACAGTACAATTCAAAGTTTGACTAAAATTGAAAATGGCGTTGCAGCAATTG
CACCTTTAGGTAGTGGCGAAGCTTATGGTTTTACACCTATTGATACACATATTGAAGACTATCCACATAA
CGTGACACGTTTTTTAGTTATAAAAAATCAGCAGCAATTTGATCAAAATGCGACATCATTGATGTTTTTA
ATTACACCAATGCATGATAAACCTGGGTTGTTAGCAAGTGTGCTTAATACATTTGCTTTGTTTAATATTA
ATTTATCCTGGATTGAGTCTCGTCCTTTAAAAACACAATTAGGTATGTATCGATTTTTCGTTCAGGCTGA
TTCTGCAATTACAACCGATATAAAGAAAGTCATTGCTATTTTAGAAACCTTAGATTTCAAAGTAGAAATG
ATTGGTGCTTTTAATTAAACTTTATTTTCTAGGTTCTACAATTTAAATAAAATAAAGAGTCAACCAAGGA
TTAAAATAGGAATCCATCCTTGGTTGACTTTCTTTTGTTTGTATTGTTTTTACACGATGCTGCGGTATAC
GTTTCGTTACCACTCGTCCAAGCAATTGCTAATGTAACCTTTTATTTTTCGTTCTAATCTGCAATTACTA
TTTCAATGGCAGTGGTTGATTTATATCAATACCTAAAACTATAGGCATTACATAATACACGACTAAAATA
ATAATTATTGCACTGATTAAGTTTACCCAGAATCCTACTGATGCCATTTGTTTGATAGATATTTTACCAG
AACCAAAGATAATTGCATTCGGTGGTGTCCCTACTGGTAACATGTATGCACAGTTAGCCGCCATAGCTGC
AGGTGCCATAAGTAGTAATGGATGCACTCCAACAGCAACAGACAACGTTGCTAAAATCGGTAAAATCATC
GTTGCAGTTGCAGTATTAGATGTCACTTCAGTTAAAAATAAGACAAAGATTGTTATGACAATTACAATAA
GAATCGGACTAACACCATTTAATGATTTCAACTGTTCGCCTAACCATTTTGCTAAACCACTTTCAGAAAT
ACCTTTCGCTAATGCTAAACCGCCACCAAATAAAATTAATACACCCCAAGGGAGCTCTTTTGCAACTTCC
CAGTCAATGATACGGCGATGTTTTTCAGTATTTTTAGCTGGAATAATAAATAATAATATTGATATAAAAA
TAGCAATCGTACCATCTGCAACAGATGACGTAACTTCCCATTTTTTCAGAAGAAACTCTCTTGTAATCCA
TAATAAGCTAGCAAGTACAAAGATAGTTTGTACTACCTTTTCTTCATACTTCATTTTGCCAAGCTCGTCT
AACTTTTGTTTAATTAACGTCTGACCACCAGGTAAATATTTCAAATCATGTCTAAACGCAACATATCTTA
AATAGAGCCAAGTAATACCTAACAAAACAATGACCGTTGGAATCCCTACAATCATCCATTTAGCAAAACT
AATTTCATGTCCAAAATGTTGCATGTATTGTCCTTTTAAAATAATTAATGGCGGGGTTCCGATTAATGTA
CCCAAGCCACCAATCGTACCTGCATAGCCAATTGCTAAAACTAGAGATTTTTCAAACTTTTGAATACTTG
TTTGATTCGTATTGGCTTCTTGTAAATCATGTGCTTCCTTAATAATTGCTAAACCAATCGGAATCATAAT
CATTACAGCTGCAGTGTTCGATACAAACATAGATAAGAATCCTGTTGCCACCATGAATCCAAGTAATATT
TTTGAAGTGCTTGCACCAATCAAATTAATAATTGTCAACGCTACTCTCGTATGTAAATTCCATCTTTCCA
TTGCAATTGCCAAAATAAATCCACCTAAAAACAAAAAGATAATATCATTGCCATATTCGGATGATACTTG
TTCTGGTGTAAGTATATGACCTAATGGTAATAACACAATTGGTAATAAGCTCGTTGCTGCAATAGGAATT
GCTTCAGTAATCCACCAAGTCGCAATCCATAAAGTAATCGCTAAAACATAGACGCCTTTCCAAGGTAAGT
CTTGTGGATGAAAGAATAATAATGTTAATAGGAAAAGTAATGGACCTAATATTAGTCCAATAAGTTGTCC
CGCACTAAACGGTTTCTTTTTCTTTGATTTTGAAAAGAATGTTAAATACCTTTTCGATATCATTGATTGA
TGTTGATTGAAATAAGCCATAAATTTCCCCCTCATGATGATTATGACAAGCTACTTCAACTCACTATATT
AGAGTATAGAAAGTTGTCTATCTTGTATTATCATTACCGTCAAAATACCCCTTTTATCAGATAATTGGGT
TTATTATAGCATATTTATACCGTTTTATGAATGCGCTTTCAAAAATTATCTACTATTCTTAATCCTAATG
TCATGCGAAACTTTTTACTTTTTTATTTCATCTACATTGATTGTAGTAACTTGAATACCCGCTTCAACAG
GACTTAAAAACACTGTATCAATATCTTCTAAATGAATTGATTGAGTCGAATTATTGAATGTATTCCCAAC
GATTAATACATCTTTATGTTTAACATTATTATAATTACGTACATGTATCGCATCTTTAGACATTGATCCT
TTAAATTCATTTCCAATTATATTCATATTTATGCCTGCTTGGGAACCTAAGTCTTTTCCTGTCATCACAT
CTGCTGCATTTTTACCATCTCTAACTCCTAAATATCTAATGCCACCCTCACAGTTAATAAACTTATTATT
AATAATGAAAGCATCTTTATACTTCAAGGGAGTTAATGCATAACCTTGTATATCTTCAAATATATTATTT
CTAATATGAATATTCTCATAGTATCGATTATGTCTACTTGCATGTGAGCCAATAGCACGATTCCAACTTC
CCATTTCGGGCAATTCTGAAGGTCCAAAATAACAATCTTCGATAATGACATTTTTCGTTATCGTACCATC
TGTCGTTCCGAATTTTGGAAAAGCACCAGGTACTTGAATGTCTAATTGTATTGCTTCAGAATAAAAGCGT
TCGCCACTATAGTCTATGAATCCTTCAAATGAACAGCTTTTAATATAGAGTCCGTTAATCCCACAAGCAT
CAATTGCATGACCACTCACTACATTCTTAATGGTCACACCTATTAATTGAATATCTTCAGCATGCCCAAT
GCACATAGCTGTATTGTTATAAGGATACGATACACCATTCATATCAAACTTGCCGCCTTTAATATGAATA
TGACTGTGTCCATTATAACCACGATAAAAACCAAAGCGACGACCATTTTTTAATAAAGGACCAGAATGAC
ATCGGCGTAAAATAGTTTCATTATCTAACAAAAGTGTTGTATTGCCATATATCGTTAATGGTTTGCAAAT
ATCATACGTCCCTTTCGGTATATAAACTGTCGTTGGTTTACATCTTCCACGATTTAACGCACGTTGTATT
GCTCGAGTATCTTTAGTTTTGTTTATGCCTGTTAACCCATAATCTTTCATATTGATATACATAAATCGTC
ATCCTTTGTCTATTTTCTAATTACTCATATATTAACCAATTATTCAAATTATAATCGTGTTAGGTTGATT
TAACAGATAGTCATAAATCTCGTGATGCCATCTATTAGCATTTTTATATTCAAAACATCAACACTTGATA
TTAGTTCTAAGTTTGTCACATTTAGTATTCTATTACAGCTATTTTAAAACTATCGCAAAAACTTCTTATA
TTATTTGTTGATTGCGTGATTTGCTCTTCCAAAAAAGTTGATTAATCTATATACTTAATCTATATTTATA
TATTAACCATTAGGGTTAAAAATTACTCTAGCATTTATGAATAGATGGGAGTTTATTTTATTATTATATA
GGAGAGATGTTGAATGACACATCGCGCACTATTAGTTGTTGACTATTCATATGACTTTATCGCAGACGAC
GGCTTACTAACATGCGGTAAACCTGGACAAAATATTGAAGATTTTATTGTTTCTCGTATCAATGACTTTA
ATTATTATCAAGACCATATATTCTTTTTAATGGATTTACATTATTTACATGACATTCATCATCCTGAAAG
TAAATTATTCCCACCACACAATATCGTAGATACAAGTGGTAGAGAATTATACGGTAAAGTAGGTAAATTA
TACGAAACAATTAAAGCGCAACCTAATGTACATTTCATTGATAAAACGCGCTATGATTCGTTCTTTGGTA
CCCCGCTTGATAGTTTATTGAGAGAAAGAAGTATTAATCAAGTCGAAATCGTTGGTGTATGTACCGATAT
TTGCGTGTTACATACAGCAATTTCTGCATACAACTTAGGTTATAAAATTTCAGTACCTGCTGAGGGAGTG
GCTTCATTTAATCAAAAAGGGCATGAATGGGCACTTGCACATTTCAAAAACTCATTAGGTGCAGAGGTAG
AACAACACGTTTAAATCGTGCTAAAATAATTATAAAGAATACAATTTACAAGGGAGATATTTGACAATGG
CTAAAACATATATTTTCGGACATAAGAATCCAGACACTGATGCAATTTCATCTGCGATTATTATGGCAGA
ATTTGAACAACTTCGAGGTAATTCAGGAGCCAAAGCATACCGTTTAGGTGATGTGAGTGCAGAAACTCAA
TTCGCGTTAGATACATTTAATGTACCTGCTCCGGAATTATTAACAGATGATTTAGATGGTCAAGATGTTA
TCTTAGTTGATCATAACGAATTCCAACAAAGTTCTGATACGATTGCCTCTGCTACAATTAAGCATGTAAT
TGATCATCACAGAATTGCAAATTTCGAAACTGCTGGTCCTTTATGTTATCGTGCTGAACCAGTTGGTTGT
ACAGCTACAATTTTATACAAAATGTTTAGAGAACGTGGCTTTGAAATTAAACCTGAAATTGCCGGTTTAA
TGTTATCAGCAATTATCTCAGATAGCTTACTTTTCAAATCACCAACATGTACACAACAAGATGTTAAAGC
AGCTGAAGAATTAAAAGATATTGCTAAAGTTGATATTCAAAAGTACGGCTTAGATATGTTAAAAGCAGGT
GCTTCAACAACTGATAAATCAGTTGAATTCTTATTAAACATGGATGCTAAATCATTTACTATGGGTGACT
ATGTGACTCGTATTGCACAAGTTAATGCTGTTGACCTTGACGAAGTGTTAAATCGTAAAGAAGATTTAGA
AAAAGAAATGTTAGCTGTAAGTGCACAAGAAAAATATGACTTATTTGTACTTGTTGTTACTGACATCATT
AATAGTGATTCTAAAATTTTAGTTGTAGGTGCTGAAAAAGATAAAGTTGGCGAAGCATTCAATGTTCAAT
TAGAAGATGACATGGCCTTCTTATCTGGTGTCGTTTCTCGAAAAAAACAAATCGTACCTCAAATCACTGA
AGCATTAACAAAATAATACTATATTACTGTCTAATTATAGACATGTTGTATTTAACTAACAGTTCATTAA
AGTAGAATTTATTTCACTTTCAATGAACTGTTTTTTATTTACGTTTGACTAATTTACAACCCTTTTTCAA
TAGTAGTTTTTATTCCTTTAGCTACCCTAACCCACAGATTAGTGATTTCTATACAATTCCCCTTTTGTCT
TTACATTTTCTTAAAATATTTGCGATGTTGAGTATAAATTTTTGTTTTCTTCCTACCTTTTTCGTTATGA
TTAAAGTTATAAATATTATTATGTACACGATTCATCGCTCTATTTTCAACTTTCAACATATATAATTCGA
AAGACCATTTAAAATTAACGGCCACAACATTCAAATCAATTAATCGCTTTTTCCAAAATAATCATATAAG
GAGGTTCTTTTCATTATGAATATCATTGAGCAAAAATTTTATGACAGTAAAGCTTTTTTCAATACACAAC
AAACTAAAGATATTAGTTTT
```
</details>

### Task
- Open the nucleotide BLAST page (blastn) and upload the file
  `draft_assembly_contigs.fasta` ("Or, upload file" button). All 11 contigs
  will be searched at once.
- Database: "Core nucleotide database (core_nt)", Program: megablast
- On the results page, use the "Results for" menu to switch between the
  contigs.

Questions:
- Fill in a table: contig, length, coverage, best hit (organism), % identity,
  query coverage, and your verdict (*S. aureus* / contaminant / chimera).
- Which contaminants did you find? For each of them, suggest how it could get
  into your sample.
- One contig is covered by two different organisms - look at the "Graphic
  Summary" for this contig. Which part of the contig belongs to which organism?
  Is it a real biological sequence or an assembly error? What would you do with
  this contig?
- Compare the coverage (`cov` in the FASTA header) of the *S. aureus* contigs
  and of the contaminant contigs. Could the coverage alone be used to detect
  contamination? Why is the coverage of one of the contaminants so high?
- Some *S. aureus* contigs also have hits to other *Staphylococcus* species with
  lower identity (~75-85 %). Why? Should these hits worry you?
- NCBI also offers a specialized service for detecting vector contamination -
  VecScreen (https://www.ncbi.nlm.nih.gov/tools/vecscreen/), which searches the
  UniVec database. Run VecScreen with the contig you identified as a vector.
  Try also the *E. coli* contig with the lac operon - can you explain the
  result?

<details>
<summary>💡 Hint</summary>

PhiX174 is a small bacteriophage whose DNA is added to every Illumina
sequencing run as a quality control ("spike-in"). Many cloning vectors (pUC,
pBR322, ...) carry parts of the *E. coli* lac operon, so a genuine *E. coli*
sequence can give a "vector" hit in VecScreen.
</details>

<details>
<summary>Backup - if NCBI BLAST is not available</summary>

**1. EMBL-EBI BLAST** (https://www.ebi.ac.uk/jdispatcher/sss/ncbiblast) -
paste or upload the contigs, select program BLASTN and several databases at
once: "ENA Sequence Prokaryote", "ENA Sequence Human", "ENA Sequence Phage"
and "ENA Sequence Synthetic" (cloning vectors). The search can take several
minutes.

**2. Local BLAST** - a small database with the *S. aureus* NCTC 8325 genome
(the reference genome from the genome assembly practicals), several possible
contaminants and a related species (*S. epidermidis*), and the UniVec_Core
vector database are in `~/Desktop/Bioinformatics/data/blast_data/contamination/backup_db/`:

```bash
mkdir -p ~/blast_backup/contamination
cd ~/blast_backup/contamination
cp ~/Desktop/Bioinformatics/data/blast_data/contamination/draft_assembly_contigs.fasta .
cp ~/Desktop/Bioinformatics/data/blast_data/contamination/backup_db/*.fasta .
# combine S. aureus genome and other reference sequences into one database
cat ~/Desktop/Bioinformatics/data/genome_assembly/GCF_000013425.1_ASM1342v1_genomic.fna \
  reference_sequences.fasta > mini_nt.fasta
grep ">" mini_nt.fasta
makeblastdb -in mini_nt.fasta -dbtype nucl
makeblastdb -in UniVec_Core.fasta -dbtype nucl
# columns: contig, subject, % identity, alignment length, query start, query end,
#          contig length, E-value, bit score, subject description
blastn -query draft_assembly_contigs.fasta -db mini_nt.fasta -max_target_seqs 5 \
  -outfmt "6 qseqid sseqid pident length qstart qend qlen evalue bitscore stitle" \
  | awk '$4 > 300' | less -S
# vector screening (local equivalent of VecScreen)
blastn -query draft_assembly_contigs.fasta -db UniVec_Core.fasta -evalue 1e-20 \
  -outfmt "6 qseqid sseqid pident length qstart qend qlen evalue" | less -S
```
Note: one of the contigs has no hit in `mini_nt.fasta` - which database finds it?

**3. Pre-computed NCBI results** - if nothing else works, the results of the
NCBI search (run in October 2026) are in
`~/Desktop/Bioinformatics/data/blast_data/contamination/precomputed/`.
</details>


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
