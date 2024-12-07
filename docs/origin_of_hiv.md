## Phylogenetic Analysis of HIV to Determine Its Origin in African Primates 

### Introduction 

Human Immunodeficiency Viruses (HIV-1 and HIV-2) are closely related to Simian Immunodeficiency Viruses (SIVs) found in various African primate species. Understanding the evolutionary relationships between HIV and SIV strains can reveal the origins of HIV and its transmission pathways from primates to humans. This guide outlines the steps to perform a phylogenetic analysis of HIV sequences using open-source tools to determine their primate origins.

### Objectives 

- Align HIV and SIV DNA sequences.
- Inspect and refine the sequence alignment.
- Construct a phylogenetic tree.
- Visualize and analyze the phylogenetic tree to infer the origins of HIV.

### Step-by-Step Guide 

#### 1. Sequence Alignment with MAFFT 
**MAFFT**  is used to create a multiple sequence alignment of HIV and SIV DNA sequences. Data for alignment are located in `~/Desktop/Bioinformatics/data/phylogenetic
/HIV_sequences.fasta`
`

```bash
cd
mkdir hiv_dentist
cd hiv_dentist
mafft --auto HIV_sequences.fasta > HIV_sequences_aligned.fasta
```

#### 2. Inspect Alignment in Jalview 

  - Examine the alignment for large gaps or misalignments.
  - Make necessary adjustments to ensure accurate alignment.

#### 3. Construct Phylogenetic Tree Using Galaxy’s PhyML-SMS 

  - Go to [Galaxy Pasteur](https://galaxy.pasteur.fr/) .
  - Click on `Upload Data`.
  - Select and upload `HIV_sequences_aligned.fasta`.
  - In the tools panel, search for **PhyML-SMS** .
  - Select the **PhyML-SMS**  tool from the search results.
  - Configure PhyML-SMS Parameters** : 
     - **Alignment file** : `HIV_sequences_aligned.fasta`
     - **Data type** : DNA
     - Other parameters - default
 
 
#### 4. Inspect output from PhyML-SMS

There are four output datasets: 
 - SMS compare models
 - PhyML newick tree
 - SMS Best Model
 - PhyML Statistics


 - Download the tree in Newick format for visualization.

#### 4. Visualize and Analyze the Tree in Dendroscope 
  - Open Dendroscope 
  - Navigate to `File` > `Open` and select the downloaded tree file.
  - Find `SIV-MON;_Mona_monkey;_AY340701` and set it as outgroup
  - Identify major groups and subtypes of HIV-1 and HIV-2. You can use commad window to select and color HIV-1 and HIV-2 leaves:
```
find searchtext='HIV-1' target=nodes regex=true;set labelcolor=#AA0000;deselect all; 
find searchtext='HIV-2' target=nodes regex=true;set labelcolor=#00AA00;deselect all;
```    
  - Examine the placement of SIV strains to infer the origins of HIV-1 and HIV-2.
  - Determine from which **African** primates HIVs originated. 
  - Do you think that transmission occurred just once or multiple times?



