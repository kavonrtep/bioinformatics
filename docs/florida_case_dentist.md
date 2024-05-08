# The case of the Florida dentist

This exercise is based on a real incident that occurred in Southern Florida in the 1990s, involving a dentist and a group of patients who claimed to have contracted HIV from their dental practitioner. If this assertion holds true, we would expect the HIV virus isolated from the dentist to be more closely related to the patients' virus than to other 'control' sequences from the wider population. To test this hypothesis, we will compare HIV gene sequences isolated from the dentist, his patients, and the wider population, using these alignments to construct a phylogenetic tree.

We will utilize the *env* gene, which codes for the outer coat of the HIV virus. The *env* sequences, obtained via PCR from the patients,dentist (three samples), control samples and outgroup, are available at the following link: https://raw.githubusercontent.com/kavonrtep/bioinformatics/master/data/phylogenetic/hiv_env_renamed.fna

List of DNA sequences in the file:

Control from Florida area:
- Control_1
- Control_2
- Control_3
- Control_4
- Control_5
- Control_6
- Control_7
- Control_8
- Control_9

Three isolates from the dentist:
- Dentist_s1
- Dentist_s2
- Dentist_s3

Outgroup from Africa:
- isolate_ELI

Patients:
- PatientA
- PatientB
- PatientC
- PatientD
- PatientE
- PatientF
- PatientG
- PatientH

Instructions:
1. Create a multiple sequence alignment (MSA) using the *env* sequences. Use mafft program with option  and `--reorder`:

```shell
mkdir hiv_env 
cd hiv_env
cp ~/Desktop/Bioinformatics/data/phylogenetic/hiv_env_renamed.fna .
mafft --reorder hiv_env_renamed.fna > hiv_env_aligned.fna
```

3. Calculate the genetic distances between the sequences in the MSA file. Use the distmat program with the following command:

```shell
distmat --help
distmat -sequence hiv_env_aligned.fna -nucmethod 1 -outfile hiv_env_distances.csv
```

4. Import the distance matrix to LibreOffice Calc or Excel and visualize it as a heatmap using conditional formatting. What do you observe? 

5. Construct a phylogenetic tree using the MSA file. Use the PhyML program on Galaxy server. Server is available at https://galaxy.pasteur.fr

Use Phyml-SMS program with default settings. This script runs SMS to select the substitution model which best fits the input data. It also runs PhyML with the selected model.

6. Download tree file in Newick format and visualize it using Dendroscopeprogram. Alternatively tree can be visualized using web-based Presto tool http://www.atgc-montpellier.fr/presto/

7. Use command tool to color dentist, patients, and control sequences in the tree. Use the following command:

```dendroscope
find searchtext='Dentist' target=nodes regex=true;set labelcolor=#FF0000;deselect all; 
find searchtext='Patient' target=nodes regex=true;set labelcolor=#00FF00;deselect all;
find searchtext='Control' target=nodes regex=true;set labelcolor=#0000FF;deselect all;
```
8. Select outgroup and use Edit/reroot on outgroup to root the tree.

9. Answer the following questions:
   - Are the sequences from the dentist more closely related to the patients or the control sequences?
   - What can you infer from the phylogenetic tree about the origin of the HIV sequences in the dentist and patients?
   - Which patients are most likely to have contracted the virus from the dentist if any?
   - Is divergence in control sequences higher than in the dentist and patients?
