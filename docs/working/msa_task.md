## Cysteine Proteases 
Cysteine proteases are enzymes that cleave peptide bonds in proteins and play vital roles in various biological processes such as protein catabolism, immune responses, and cell remodeling. Their catalytic function typically depends on a conserved **catalytic diad**  of residues—**cysteine (C)**  and **histidine (H)** —and sometimes a nearby asparagine (N) to stabilize this catalytic site. Although these residues form a functional diad in the three-dimensional enzyme structure, they are often not adjacent in the primary amino acid sequence. Instead, their positions are scattered along the linear chain but fold together in the protein’s tertiary structure to form the active site capable of nucleophilic attack on the substrate’s peptide bond.

---

## Your Task 
Identify the positions of the catalytic amino acids (C and H) in an **unknown cysteine protease** . To do this, you will create a multiple sequence alignment (MSA) using sequences of known cysteine proteases together with the unknown protein.

---
### Sequences of Known Cysteine Proteases

Below are examples of cysteine proteases. 
 
- Human Cathepsin B  (P07858)
- Human Cathepsin L  (P07711)

Full set of sequences can be downloaded from [here](https://raw.githubusercontent.com/kavonrtep/bioinformatics/refs/heads/master/docs/working/cystein_proteases.fasta).
---


### Unknown Protease 


```text
>unknown_protease
ipeYVDWRqkgavtpVKNQGsCGSCWAFSAVVTIEGIIKirtgnlnQYSEQELLDCDrrsygcnGGYPwsal
qlvaqYGIHYRnTYPYegvqrycrsrekgpyaaktDGVRQVqpynqgaLLYSIAnqPVSVVLQAagkdfqly
rggifvgpcgnkvDHAVAAVGYGpNYILIKNSWGTGWGenGYIRIKRgtgnsygvcglytSSFYpvkn
```


---


## Instructions 
- Download the protease sequences from the provided link and save them ia a file named `cysteine_proteases.fasta`.
- Using a text editor, add the unknown protein sequence(shown above) to the same FASTA file.
- Open the resulting FASTA file in Jalview program.
- Perform a Multiple Sequence Alignment (MSA) using MAFFT (use menu Web Service -> Alignment -> Mafft with default settings).
 - Use the Clustal coloring scheme and the “Colour by Conservation” function to highlight conserved amino acids.
 
- Identify the catalytic diad—cysteine (C) and histidine (H)—in the alignment.  
  - Refer to the UniProt record for **[Human Cathepsin B (CATB_HUMAN; P07858)](https://www.uniprot.org/uniprotkb/P07858/entry#function)**  to verify the correct identification of active site residues. (Active site description can be found in the “Function/Features” section of the UniProt entry.)
 
8. Record the positions of the catalytic C and H  in the unknown protein **based on the unknown protein’s own sequence numbering** .


---


## Questions 
 
1. **What are the positions of the active site amino acids (C and H) in the unknown protease?** 
(Note: Give their positions according to the unknown protein’s own sequence, **not**  the alignment.)*

2. **How many amino acids are 100% conserved in your alignment?** 
*(Consider all sequences in your alignment, including the unknown.)*
 
3. **Save the final alignment as an image**  
    - *(In Jalview: File > Export Image > PNG)*
    - Upload this image to Moodle.
 
4. **Save the alignment in FASTA format**  
    - *(In Jalview: File > Save As > FASTA (aligned))*
    - Upload this file to Moodle.


---


### Additional Tips 
 
- Look for columns in the alignment that show **C**  or **H**  strictly conserved across **all**  known cysteine proteases.
 
- Compare those positions to the known catalytic residues listed for **Human Cathepsin B (P07858)**  in UniProt.
 
- Make sure to use the **unknown protease’s own numbering**  (starting at 1 in its FASTA sequence) when reporting the final positions.