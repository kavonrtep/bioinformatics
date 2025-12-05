## Core Examples (For the Whole Practical)

### Hemoglobin Subunit Alpha (HBA1/HBA2)

- **UniProt:** P69905 – human hemoglobin alpha, length 142 aa. [UniProt](https://www.uniprot.org/uniprotkb/P69905/entry)
- **RCSB PDB:**
  - **4HHB** – human deoxyhemoglobin tetramer, 1.74 Å, no mutations. [RCSB PDB](https://www.rcsb.org/structure/4HHB)
  - **2HHB** – also human deoxyhemoglobin, classic textbook structure, 1.74 Å. [RCSB PDB](https://www.rcsb.org/structure/2HHB)
  - **1A01** – human hemoglobin mutant (Val β1→Met, Trp β37→Ala). [RCSB PDB](https://www.rcsb.org/structure/1A01)

### β2-Adrenergic Receptor (ADRB2_HUMAN)

- **UniProt:** P07550 / ADRB2_HUMAN, 413 aa multipass GPCR. [UniProt](https://www.uniprot.org/uniprotkb/ADRB2_HUMAN)
- Many partial X-ray and cryo-EM structures (e.g., 2R4R, 3P0G, 3SN6, etc.). [SWISS-MODEL](https://swissmodel.expasy.org/repository/uniprot/P07550)


## Used Databases and Resources

- **UniProt** = sequence & function. The UniProt site is: [https://www.uniprot.org/](https://www.uniprot.org/)
- **PDB** = 3D structure archive. The RCSB site is: [https://www.rcsb.org/](https://www.rcsb.org/)
- **PDBe** = PDB in Europe: [https://www.ebi.ac.uk/pdbe/](https://www.ebi.ac.uk/pdbe/)
- **AlphaFold DB** = predicted structures: [https://alphafold.ebi.ac.uk](https://alphafold.ebi.ac.uk)

---

## Linking UniProt to PDB Structures

### From UniProt P69905 to PDB Structures

1. Go to UniProt: [https://www.uniprot.org/](https://www.uniprot.org/)
2. In the search box, type `P69905` and press Enter.
3. Click the "Hemoglobin subunit alpha – H. sapiens" entry: [UniProt](https://www.uniprot.org/uniprotkb/P69905/entry)
4. In the left menu (or sections tabs), locate and click "Structure".
5. Scroll through the PDB structures section:
   - Note the long list of PDB IDs mapped to this UniProt entry (hundreds of structures).
   - Identify 2HHB and 4HHB in that list.
6. Check the sequence length in UniProt: it should be 142 amino acids.

### Chain Mapping & Coverage

7. Still in the Structure section for P69905:
   - Not all PDB structures cover the full length of the protein (often due to disordered tails being cleaved for crystallization).
   - Observe that the PDB chains often cover only residues 2–142 of the UniProt sequence (the initiator Met is typically missing in the crystallized protein).
   - Note distinct sources of structures - primary is PDB, others are modeled or predicted (e.g., AlphaFold).
   - Inspect also Resolution and Chain columns. Note identifiers used in PDB - four letter codes (e.g., 2HHB).

### ADRB2_HUMAN in UniProt

8. In a new tab, search UniProt for ADRB2_HUMAN or simply `P07550`. [UniProt](https://www.uniprot.org/uniprotkb/ADRB2_HUMAN)
9. Open the "Subcellular location" section:
   - Confirm it is a multi-pass membrane protein (GPCR) embedded in the plasma membrane.
10. Check "Length" in the entry header: it should be 413 aa.
11. Open the "Structure" section:
    - Note that many PDB structures cover only residues approx. 29–340, not the full 413 aa.

### AlphaFold Models from UniProt

12. In the UniProt `P69905` page [UniProt](https://www.uniprot.org/uniprotkb/P69905/entry), look for links to AlphaFold DB in "3D structure databases":
    - Open the AlphaFold model for P69905. You should see average pLDDT ~98 ("Very high") over 142 residues. [AlphaFold](https://alphafold.ebi.ac.uk/search/text/P69905)
13. Do the same for `P07550` (ADRB2) in AlphaFold DB:
    - Observe that transmembrane helices have high confidence (high pLDDT), whereas some termini and loops show lower confidence. [AlphaFold](https://alphafold.ebi.ac.uk/entry/P07550)
    - Check also the Predicted Aligned Error (PAE) plot to see which regions are predicted with high accuracy relative to each other.
    - Use the mouse to select blocks in the PAE plot to highlight corresponding regions in the 3D model.

---

## Getting Familiar with the PDB (RCSB)

### Searching the PDB

1. Go to the RCSB PDB homepage: [https://www.rcsb.org/](https://www.rcsb.org/)
2. In the search bar, use UniProt ID `P69905` to search for hemoglobin alpha structures.
   - Note: You will get different results if you search `P69905` in "all fields" versus "Accession codes only".
   - If you search in all fields, you will get a list of structures which you can further filter, for example, by experimental method or resolution.
   - If you click on "Accession codes only", you will get a list of PDB IDs which are mapped to this UniProt entry and alignment to this entry.
3. Explore the results in Sequence Alignment view, Structural Feature view, and Group Summary view.
4. Do the same search for ADRB2_HUMAN (UniProt ID `P07550`) and explore the results.
   - Specifically compare Sequence Alignment with Alignment view in P69905 search.
   - Note that in ADRB2_HUMAN many structures cover only part of the sequence.
5. Try a similar search but in PDBe database: [https://www.ebi.ac.uk/pdbe/](https://www.ebi.ac.uk/pdbe/)
   - Use Advanced search with UniProt ID `P69905`, compare results.
   - In PDBe it is easy to sort by **resolution** and **quality** metrics.


## Structure Validation Metrics

1. Compare validation metrics for 1A01 and 4HHB - open both PDB entries.
2. Check the following metrics:
   - Resolution
   - R-free
   - Ramachandran outliers
   - Clashscore
3. Open the Validation 3D Report for each structure.
4. Inspect validation metrics also in PDBe database: [https://www.ebi.ac.uk/pdbe/](https://www.ebi.ac.uk/pdbe/)
   - Use the "Model Quality" tab.

## Sequence Annotation Viewer in RCSB

Open PDB entry 4HHB and explore by clicking on "Sequence annotation".

The Sequence Annotations Viewer provides graphical summaries of PDB protein biological and structural features and their relationships with UniProtKB entries:

- A white circle over the start (or end) of a block indicates that the feature in its original reference does not start (or end) in that position but before (or after).
- A dashed line between blocks indicates that in the original reference the feature is connected and therefore a sequence insertion has occurred.
- Sequence Annotations Viewer will display the full range of available features (structural and biological annotations) and the alignments between Polymer Instances (chains) and UniProtKB sequences.
- Tooltips appear when the cursor is hovering over a specific feature in the Sequence Annotations Viewer.