**ChIP-Seq Analysis of CENH3 in *Pisum sativum* Using Galaxy**

In this exercise, you will analyze ChIP-seq data from Pisum sativum (garden pea) to study centromere locations. Pisum sativum has metapolycentric chromosomes, which means it has multiple centromeric regions on a single chromosome. Your task is to identify these regions associated with the centromere-specific histone H3 variant CENH3.

You will compare the identified CENH3-associated regions with existing gene annotations and tandem repeat sequences. This comparison aims to provide an understanding of the structural organization of the pea genome, particularly how centromeric regions relate to gene distribution and repetitive DNA elements.

Your task will be to map ChIP and Input sequences to a genome assembly to identify CENH3-associated regions. Explore the types of sequences associated with CENH3 regions in *Pisum sativum* and investigate their characteristics.


---

### **Materials:**

- **Galaxy Server:** Access [Galaxy Europe](https://usegalaxy.eu)
- **Data to Import:**
  - **Genome Assembly (FASTA):** Partial assembly representing a 177.6 Mb region of chromosome 6 of *Pisum sativum*, including its 81.6 Mb centromere and adjacent chromosome arms.
  - **ChIP Sequences (FASTQ):** Sequenced reads obtained using CENH3 antibody.
  - **Input Sequences (FASTQ):** Control sequenced reads without antibody enrichment.
- **Additional Resources (Available on Student Desktops):**
  - **Gene Annotation (GFF):** `genes_jon_220520.gff`
  - **Tandem Repeat Annotation (GFF):** `tandem_repeats_curated_track_220419.gff`
- **Visualization Tool:** **IGV** (Already installed on your computer)
- **Worksheet: [link](https://docs.google.com/document/d/1z8KmdVRGPSiNWWpnt_xvPKSikyK9vVgDgSf5MjABchk/edit?usp=sharing)**

---

### **Workflow Overview:**

1. **Data Import**
2. **Quality Control and Trimming**
3. **Read Mapping**
4. **ChIP vs. Input Comparison**
5. **Data Visualization in IGV**
6. **Exploration of CENH3-Enriched Regions**

---

### **Step-by-Step Instructions:**

#### **Step 1: Data Import**

1. **Login to Galaxy Europe:**

   - Go to [https://usegalaxy.eu](https://usegalaxy.eu) and log in or create an account.

2. **Import Data:**

   - Access the shared history containing the datasets: [ChIP-Seq CENH3 Data](https://usegalaxy.eu/u/petr_novak/h/chip-seq-cenh3).
   - Import the following files into your history:
     - **Genome assembly FASTA file**
     - **ChIP FASTQ file**
     - **Input FASTQ file**

---

#### **Step 2: Quality Control and Trimming**

1. **Initial Quality Assessment with FastQC:** This tool is used to assess the quality of raw sequencing data, providing metrics such as per-base quality scores and identifying any issues like adapter contamination.

   - Run **FastQC** on both the ChIP and Input FASTQ files.
     - **Tool:** *FastQC: Read Quality reports*
     - **Parameters:** Use default settings.

2. **Identify Potential Issues:** Review the FastQC report to identify common quality issues, such as adapter contamination or low-quality bases, which can affect downstream analysis.

   - Review the FastQC report to identify:
     - **Adapter contamination**
     - **Low-quality bases at the ends**
     - **Overrepresented sequences**

3. **Trim and Filter Reads with fastp:** fastp is used to remove adapter sequences and low-quality bases from the reads, ensuring that only high-quality data is used for further analysis.

   - Use **fastp** to trim adapters and filter low-quality reads.
     - **Tool:** *fastp*
     - **Parameters:**
       - **Adapter Trimming:** Enabled by default; fastp automatically detects and removes adapters.
       - **Quality Filtering:**
         - **Quality Threshold:** By default, fastp performs quality filtering using a sliding window approach with Q15 as the threshold.
       - **Global trimming:**
         - Trim front to remove first **9** bases
       - **Length Filtering:**
         - **Length Required:** Set a minimum length (set it to **140 bp**). Reads shorter than this length after trimming will be discarded.
         - **Explanation:** This parameter ensures that only reads of sufficient length, which are more likely to map accurately to the genome, are retained.
       - **Other Parameters:** Leave as default unless instructed.
     - **Run fastp on both ChIP and Input samples.**

4. **Post-Trimming Quality Assessment:** Run FastQC on the trimmed reads to verify improvements, and use MultiQC to summarize the quality of the processed data.

   - Run **FastQC** again on the trimmed FASTQ files.
   - Use **MultiQC** to summarize both the previous and the new FastQC reports.
     - **Inputs:** Select the FastQC outputs for the trimmed ChIP and Input samples.
     - **Result:** A combined report summarizing the quality of the datasets **after** processing.
   - Did the filtering by fastp provide expected results?

---

#### **Step 3: Read Mapping**

1. **Map Reads with Bowtie2:** Bowtie2 is used to align the trimmed reads to the genome assembly, producing BAM files that represent the mapped locations of the reads.
   - Map the trimmed reads to the genome assembly for both ChIP and Input samples.
     - **Tool:** *Bowtie2*
     - **Parameters:**
       - **Input Type:** Single-end
       - **FASTQ Files:** Use the trimmed FASTQ files from fastp.
       - **Reference Genome:** Select the indexed genome from Bowtie2-Build.
       - **Preset Options:** Very fast end-to-end (`--very-fast`).
       - Set **Save the bowtie2 mapping statistics to the history** to *yes*
       - **Output Format:** BAM file.
   - Inspect mapping statistics. Did the mapping work as expected?

---

#### **Step 4: ChIP vs. Input Comparison**

1. **Compare BAM Files with bamCompare:** bamCompare calculates the log2 ratio of ChIP and Input coverage, helping identify genomic regions enriched for CENH3 binding.

   - **Tool:** *bamCompare* (from deepTools)
   - **Description:** This tool is used to compare two BAM files, typically a treatment and a control, by calculating the log2 ratio of their coverage. This allows visualization of regions enriched in the treatment compared to the control.
   - **Parameters:**
     - **Treatment File:** ChIP sorted BAM file.
     - **Control File:** Input sorted BAM file.
     - **Bin Size:** **5000** bases.
     - **Operation:** **log2 ratio (ChIP/Input)**
     - **Output Format:** BigWig file.
     - **Other Parameters:** Use defaults unless specified.

2. **Convert BigWig to BedGraph:** This tool is used to convert BigWig format files to BedGraph format, which is a more accessible format for downstream analysis and visualization of coverage data.

   - **Tool:** *bigWigToBedGraph*
   - **Input:** BigWig file from bamCompare
   - **Output:** BedGraph format file for further analysis

3. **Select Enriched Regions with MACS2 Broadpeakcall:** MACS2 is used to call broad peaks in the data, identifying regions of significant enrichment that are likely associated with CENH3.

   - **Tool:** *MACS2 Broadpeakcall*
   - **Input:** BedGraph from previous step
   - **Output:** BED file with significant regions
   - **Parameters:**
     - **Cutoff for peaks:** 4
     - **Cutoff for peaks:** 3.5
     - **Minimum length of peak:** 200
     - **Maximum gap between significant peaks:** 20000
     - **Maximum linking between significant peaks:** 21000

---

#### **Step 5: Data Visualization in IGV**

1. **Prepare Files for IGV:** Prepare all necessary files for loading into IGV, which will allow for visualization of ChIP-seq results alongside annotations.

   - **Files to Load:**
     - **Genome Assembly:** FASTA file. Save it as `genome.fasta`
     - **BAM Files:** ChIP and Input sorted BAM files.
     - **BigWig File:** Output from bamCompare.
     - **BED File:** Output from MACS2 tool, save it as `enriched_regions.bed`
     - **Annotations:** Gene and tandem repeat annotation GFF files (already available in the folder described above).

2. **Set Up IGV:** Load the genome assembly and data tracks into IGV for detailed visual exploration of CENH3 enrichment across the chromosome.

   - **Open IGV:** Launch IGV on your computer.
   - **Load Genome:**
     - **Menu:** *Genomes* > *Load Genome from File...*
     - **Select:** The genome assembly FASTA file.
   - **Load Data Tracks:**
     - **Menu:** *File* > *Load from File...*
     - **Select:** ChIP and Input BAM files.
     - **Select:** BigWig file from bamCompare.
     - **Select:** Gene annotation and tandem repeat GFF files.

---

#### **Step 6: Exploration of CENH3-Enriched Regions**

1. **Identify Enriched Regions:** Use IGV to visually identify regions with high CENH3 enrichment based on the loaded BigWig track and BED from MACS2.

   - Use the BigWig track to visualize regions with high CENH3 enrichment.
   - **Observation:** Regions with higher log2(ChIP/Input) values indicate CENH3 enrichment.

2. **Investigate Sequence Features:** Explore enriched regions to determine if they overlap with genes, repetitive elements, or other sequence features.

   - **Questions to Explore:**
     - **What types of sequences are associated with CENH3 domains?**
       - Are they mainly repetitive sequences, genes, or a mix?
     - **How many CENH3-associated domains are present?**
       - Count the number of distinct enriched regions.
     - **Are there genes located within these domains?**
       - Use the gene annotation track to identify any genes.
     - **What is the nature of the satellite DNA in these regions?**
       - Examine the tandem repeat annotation track.

3. **Inspect Enriched Regions with Gepard Dotplot Program:** Use Gepard to generate dot plots for visualizing sequence similarities within the enriched regions.

   1. **Extract Sequences with Seqkit:**
      - **Tool:** *Seqkit*
      - **Command:** Use the following command to extract sequences from the regions identified by MACS2:
        ```
        seqkit subseq -f enriched_regions.bed genome.fasta -o enriched_sequences.fasta
        ```
      - **Input:** BED file from MACS2 and genome FASTA file.
      - **Output:** Extracted sequences in FASTA format.
   2. **Generate Dot Plot with Gepard:**
      - **Tool:** *Gepard*
      - **Description:** Use Gepard to create a dot plot of the extracted sequences to visualize similarities and repetitive structures within enriched regions.
      - **Command:** Open Gepard and load the `enriched_sequences.fasta` file to generate the dot plot.
   3. **Interpret the Dot Plot:** What is your interpretation of the dot plot? Is there one type of satellite sequence or multiple types? What is the monomer size of tandem repeat(s) associated with the centromere?

