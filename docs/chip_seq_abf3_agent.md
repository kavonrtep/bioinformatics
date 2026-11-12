# ChIP-seq with a coding agent - how does ABF3 respond to abscisic acid?

This is an alternative version of the
[ABF3 ChIP-seq exercise](./chip_seq_abf3.md). The data and the questions are
the same, but this time **you do not type the commands**. You will:

1. read the description of the experiment and the data,
2. write an **analysis plan** and discuss it with the teacher,
3. let a **coding agent** ([opencode](https://opencode.ai)) run the analysis
   step by step,
4. **check** what the agent did and **interpret** the results yourself.

The agent is a language model that can run commands in your terminal. It is
fast, but it makes mistakes, and it reports them with the same confidence as
correct results. We use a free model, which is weaker than the best
commercial ones, so it needs clear and detailed instructions. What you learn:
how to describe a bioinformatics task precisely, how to steer an agent, and
how to check its work.

> The free model runs on an external server. Your prompts, the commands and
> their outputs are sent there and may be used for training the model. Do not
> put any personal or confidential data into the prompts.

Contents:

- [The experiment](#the-experiment)
- [Installation and data](#installation-and-data)
- [Step 1 - Write the analysis plan](#step-1---write-the-analysis-plan)
- [Step 2 - First contact with the agent](#step-2---first-contact-with-the-agent)
- [Step 3 - Run the analysis step by step](#step-3---run-the-analysis-step-by-step)
- [Step 4 - Check and interpret the results](#step-4---check-and-interpret-the-results)
- [Optional - the final plan](#optional---the-final-plan)
- [What to hand in](#what-to-hand-in)

## The experiment

Read this part carefully - you cannot write a good plan (and you cannot
check the agent) unless you know exactly what each sample is and what is
compared with what.

### The biological question

The plant hormone **abscisic acid (ABA)** accumulates in plants exposed to
drought, salt or cold and switches on thousands of stress-response genes.
Part of this response are the **ABFs** (ABRE-binding factors), transcription
factors of the bZIP family. They bind the **ABA-responsive element (ABRE)**, a
short DNA motif with the core `ACGTG`, in the promoters of their target genes.
Without ABA the ABFs are inactive; ABA activates protein kinases that
phosphorylate them.

The question of this exercise: **where in the genome does ABF3 bind, and does
ABA change it?** Does ABF3 find new binding sites after ABA, or does it bind
more strongly to sites where it already was? And are the genes near these
sites switched on?

Data: [Song et al. 2016, *Science* 354: aag1550](https://doi.org/10.1126/science.aag1550),
who studied 21 ABA-related transcription factors this way.

### The ABF3-YPet plants - how to pull down ABF3

ChIP needs an antibody that recognises the protein of interest. Good
antibodies against single plant transcription factors are rarely available
(and the authors studied 21 of them). So they used a **tag** instead: they
made transgenic *Arabidopsis* plants with an **extra copy of the ABF3 gene**,
in which the coding sequence of **YPet** (a yellow fluorescent protein, a
variant of GFP) is attached to the end of the ABF3 coding sequence. The extra
copy is the whole genomic ABF3 gene with its own promoter, so the tagged
protein is made in the same cells and at the same time as the normal one.

![ABF3-YPet: two copies of ABF3, the anti-GFP antibody pulls down only the tagged one](./img/chipseq_abf3_tag.png)

The ABF3-YPet plant therefore has both: the normal ABF3 and the tagged
ABF3-YPet. The **anti-GFP antibody** recognises YPet, so it pulls down only
the tagged ABF3 together with the DNA it is bound to. We assume that the tag
does not change where ABF3 binds.

### The control - what does the antibody pull down without a target?

ChIP never gives pure binding sites. Some DNA sticks to the beads or the
antibody, some regions break more easily into fragments, some are present in
many copies or map badly. All of this gives reads that have nothing to do
with ABF3. To recognise this **background**, the authors did the same anti-GFP
ChIP in **wild-type plants (Col-0)**, which have no YPet at all - same
antibody, same procedure, but nothing to recognise. A region with many reads
in the ABF3-YPet ChIP **and** in the control is not a binding site; a binding
site has many more reads in the ABF3-YPet ChIP than in the control.

![The control: the same anti-GFP ChIP in plants without YPet shows only the background](./img/chipseq_abf3_control.png)

The coverage tracks in the figure are a made-up example, not the exercise data.

### The treatment - ABA and mock

Seedlings were grown for 3 days in the dark and then kept for 4 hours in a
solution with 10 µM **ABA**, or in the same solution without ABA (**mock**).
ABA is dissolved in ethanol, so the mock solution contains the same amount of
ethanol. The only difference between ABA and mock is therefore ABA itself.

Each condition has independent **biological replicates** (different
batches of plants): 3 + 3 for the ChIP, 2 + 2 for the RNA-seq. Replicates tell
us how much the signal varies by chance, so that we can say whether a
difference between ABA and mock is larger than that.

![Plants, treatment and the samples with their file names](./img/chipseq_abf3_design.png)

### What is compared with what

![What is compared with what: A = ChIP vs control, B = ChIP ABA vs mock, C = RNA-seq ABA vs mock](./img/chipseq_abf3_compare.png)

| | Compared | Question | Tool |
|---|---|---|---|
| A | each ABF3-YPet ChIP replicate vs the control | Where does ABF3 bind? (binding sites = peaks) | MACS3 |
| B | 3 ABA vs 3 mock ABF3-YPet ChIP replicates, read counts in the binding sites | Is the binding stronger or weaker after ABA? | DiffBind |
| - | sequences of the sites gained after ABA vs random genome regions | Which DNA motif does ABF3 recognise there? | HOMER |
| C | 2 ABA vs 2 mock RNA-seq replicates | Which genes are induced or repressed by ABA? | salmon + DESeq2 |
| B + C | genes with gained binding vs the other genes | Are the genes with gained binding induced? | your own R script |

Note that the control is **never compared with ABA or mock directly**. It
only tells us what background looks like, so that we can find the binding
sites (A). The effect of ABA is the comparison of ABF3-YPet ChIP with ABF3-YPet
ChIP (B), and of RNA-seq with RNA-seq (C).

To make the exercise fast, the reads were reduced to 3 million (ChIP) and 4
million (RNA-seq) per sample and shortened to 50 bp.

### The files

Sequencing data - one file per row, all in `~/chipseq_abf3/data/`:

| File | Assay | Plants | Antibody | Treatment | Replicate | Reads | Read type |
|---|---|---|---|---|---|---|---|
| `ABF3_ABA_rep1.fastq.gz` | ChIP-seq | ABF3-YPet | anti-GFP | ABA 4 h | 1 | 3 M | single-end 50 bp |
| `ABF3_ABA_rep2.fastq.gz` | ChIP-seq | ABF3-YPet | anti-GFP | ABA 4 h | 2 | 3 M | single-end 50 bp |
| `ABF3_ABA_rep3.fastq.gz` | ChIP-seq | ABF3-YPet | anti-GFP | ABA 4 h | 3 | 3 M | single-end 50 bp |
| `ABF3_mock_rep1.fastq.gz` | ChIP-seq | ABF3-YPet | anti-GFP | mock 4 h | 1 | 3 M | single-end 50 bp |
| `ABF3_mock_rep2.fastq.gz` | ChIP-seq | ABF3-YPet | anti-GFP | mock 4 h | 2 | 3 M | single-end 50 bp |
| `ABF3_mock_rep3.fastq.gz` | ChIP-seq | ABF3-YPet | anti-GFP | mock 4 h | 3 | 3 M | single-end 50 bp |
| `control_rep1.fastq.gz` | ChIP-seq control | wild type (Col-0) | anti-GFP | not stated | 1 | 2.6 M | single-end 50 bp |
| `control_rep2.fastq.gz` | ChIP-seq control | wild type (Col-0) | anti-GFP | not stated | 2 | 3 M | single-end 50 bp |
| `rnaseq_ABA_rep1.fastq.gz` | RNA-seq | wild type (Col-0) | - | ABA 4 h | 1 | 4 M | single-end 50 bp |
| `rnaseq_ABA_rep2.fastq.gz` | RNA-seq | wild type (Col-0) | - | ABA 4 h | 2 | 4 M | single-end 50 bp |
| `rnaseq_mock_rep1.fastq.gz` | RNA-seq | wild type (Col-0) | - | mock 4 h | 1 | 4 M | single-end 50 bp |
| `rnaseq_mock_rep2.fastq.gz` | RNA-seq | wild type (Col-0) | - | mock 4 h | 2 | 4 M | single-end 50 bp |

The same table is in two files, which you copy to the working directory
below:

- `sample_table.ods` - for you, open it in LibreOffice Calc,
- `sample_table.tsv` - the same table as tab-separated text, for the agent.

**This table is the description of the data you give to the agent.** In every
prompt, tell the agent which samples to use by referring to the table, not to
the file names, for example *"the rows of sample_table.tsv with assay =
ChIP-seq and treatment = ABA"*. This is how a real analysis is set up: the
sample sheet, not the file names, says what each file is. If you need groups
that are not in the table, add your own column (e.g. `group`) to
`sample_table.tsv` and tell the agent to compare the groups in it.

Reference files:

| File | Content |
|---|---|
| `data/TAIR10_genome.fa.gz` | genome sequence (chromosomes `1`-`5`, `Mt`, `Pt`) |
| `data/TAIR10_annotation.gff3.gz` | gene annotation (for IGV) |
| `data/TAIR10_transcripts.fa.gz` | transcript sequences (for RNA-seq) |
| `~/Desktop/Bioinformatics/data/chip_seq_abf3/TAIR10_genes.bed` | protein-coding genes (BED6: chromosome, start, end, gene ID, symbol, strand) |
| `~/Desktop/Bioinformatics/data/chip_seq_abf3/TAIR10_tss.bed` | transcription start sites of the same genes |
| `~/Desktop/Bioinformatics/data/chip_seq_abf3/TAIR10_tx2gene.tsv` | transcript → gene table |

The same directory also contains the R scripts for DiffBind and DESeq2 and
reference results (`precomputed/`).

**Biological questions:**

1. Where in the genome does ABF3 bind? Close to genes?
2. Does the binding change after ABA treatment?
3. Which DNA motif is in the sites with changed binding?
4. *(optional)* Are the genes with stronger ABF3 binding also more expressed
   after ABA?

**Suggested tools:** `fastp` (read QC), `minimap2` + `samtools` (mapping),
`deeptools` (coverage tracks, profiles), `macs3` (peak calling), `DiffBind`
(differential binding, ready R script), `bedtools` (intervals, closest gene,
sequences), `homer` (motifs), IGV (you look at the data), optionally `salmon` +
`DESeq2` (RNA-seq). See the [introduction slides](./slides/chipseq_abf3_agent_intro.pdf)
(the experiment, what is compared with what, and how to read the output of each method).

## Installation and data

Do this part yourself, without the agent, **before the practical** - the
download (1.9 GB) and the conda environment take a while. The environment is
the same as in the [original exercise](./chip_seq_abf3.md#installation-and-data).

### Conda environment and data

```bash
cd ~/Desktop/Bioinformatics && git pull
conda create -n chipseq -y --override-channels -c conda-forge -c bioconda \
    fastp minimap2 samtools macs3 deeptools bedtools seqkit homer salmon \
    bioconductor-diffbind bioconductor-deseq2 bioconductor-tximport bioconductor-genomeinfodb

D=~/Desktop/Bioinformatics/data/chip_seq_abf3
bash $D/download_data.sh
```

The script downloads the reads, the genome and the annotation from
[Zenodo](https://zenodo.org/records/23984926) into `~/chipseq_abf3/data` and
checks each file with `md5sum`. You can run it again at any time: complete
files are skipped, damaged ones are downloaded again. It must end with
`All files downloaded and checked.`

### Install opencode

[opencode](https://opencode.ai) is a coding agent that runs in the terminal.
Install it with the official script (no administrator rights needed, it goes
to `~/.opencode/bin`):

```bash
curl -fsSL https://opencode.ai/install | bash
```

The script adds `~/.opencode/bin` to your `PATH` in `~/.bashrc`. Open a new
terminal (or run `source ~/.bashrc`) and check that it works:

```bash
opencode --version
opencode models | grep free
```

The second command lists the free models (e.g. `opencode/big-pickle`). They
need no account and no API key.

- `opencode: command not found` - the new `PATH` is not active yet; open a
  new terminal, or run `export PATH=~/.opencode/bin:$PATH`.
- No free model listed - check the internet connection; the list of models is
  downloaded from the opencode server.

### Prepare the working directory

Copy the instructions for the agent and the plan template into the working
directory:

```bash
D=~/Desktop/Bioinformatics/data/chip_seq_abf3
cd ~/chipseq_abf3
cp $D/agent/AGENTS.md $D/agent/opencode.json .
cp $D/agent/sample_table.tsv $D/agent/sample_table.ods .
cp $D/agent/plan_template.md plan.md
```

- `AGENTS.md` - rules for the agent: where the data are, what it may and may
  not do, and tested settings for each tool. opencode reads it automatically
  at the start of every session. **Read it** - it tells you a lot about how
  the analysis is done.
- `sample_table.tsv` and `sample_table.ods` - one row per sequencing file
  (assay, plants, antibody, treatment, replicate, ...), the same as the table
  in [The files](#the-files). You refer to it in your prompts; the agent
  reads the `.tsv` to find the files.
- `opencode.json` - the model (`opencode/big-pickle`, free) and permissions:
  the agent may run commands, but asks before deleting or moving files or
  installing software.

- `plan.md` - your analysis plan, you fill it in during Steps 1-3.

If the free model has a different name, list the available ones with
`opencode models` and change `"model"` in `opencode.json`.

Open the working directory in VS Code (answer *Yes* if it asks whether you
trust the authors) and open `plan.md`:

```bash
code ~/chipseq_abf3
```

In VS Code switch on **File → Auto Save** (a check mark appears). The agent
reads `plan.md` from the disk - text you have not saved, it does not see. The
agent also writes into `plan.md` (the script names); with Auto Save VS Code
shows the change at once and you cannot overwrite it by accident.

## Step 1 - Write the analysis plan

**First inspect the sample table and the files.** Open `sample_table.ods` in
LibreOffice Calc and look at the same table as text, the way the agent sees
it:

```bash
cd ~/chipseq_abf3
column -t -s $'\t' sample_table.tsv
ls -lh data/
zcat data/ABF3_ABA_rep1.fastq.gz | head -8
```

- How many ChIP-seq, control and RNA-seq samples are there? Which columns
  tell them apart?
- Which rows (column = value) are the input of each comparison A, B and C in
  [What is compared with what](#what-is-compared-with-what)?
- Is every file of the table in `data/`? Do the reads look single-end 50 bp?

Use these column values in your plan and in your prompts. If you need to add
a column (e.g. `group`), open `sample_table.tsv` in LibreOffice Calc
(separator: Tab), add the column and save it with *Keep current format*.

Then write **Part 1 of `plan.md` - the outline**: your biological questions
and the list of steps, one line each, with the rows of the sample table that
each step uses. No details yet - you describe each step just before you run
it (Step 3). Think about the order of the steps (what needs what), which
samples go into which step, and what the control is used for.

**Discuss the outline with the teacher before you continue.**

<details>
<summary>💡 Hints for the outline</summary>

A ChIP-seq analysis usually has these steps. For each of them think: what
goes in, what comes out, what do I expect?

1. **Quality control of reads** - are the reads good? How many reads per
   sample?
2. **Mapping to the genome** - what fraction of the reads maps? Where else in
   the cell does DNA come from (`Mt`, `Pt`)? The two controls can be merged
   into one.
3. **Coverage tracks** (bigWig) - for looking at the data in IGV, normalised
   so that samples with different numbers of reads can be compared.
4. **Peak calling** - comparison A: which samples are the ChIP and which the
   control? Peaks are called in each replicate separately. Do you expect
   tens, thousands or millions of peaks?
5. **Profile around genes** (TSS) - does ABF3 bind close to genes? Before or
   after the start of the gene?
6. **Differential binding** (DiffBind) - comparison B: which samples are
   compared? Comparing the numbers of peaks is not enough to say that the
   binding changed.
7. **Closest genes** of the binding sites - which genes may ABF3 regulate?
8. **Motifs** in the sites with gained binding (HOMER) - which DNA sequence
   does ABF3 recognise? What do you expect (see the introduction slides)?
9. *(optional)* **RNA-seq** - comparison C: which genes change expression
   after ABA (salmon + DESeq2)?
10. *(optional)* **Binding vs expression** - B and C together, with your own
    specification ([see below](#optional---binding-vs-expression-your-own-specification)).

Each step should give something you can check yourself: a number, a table, a
picture.

</details>

## Step 2 - First contact with the agent

```bash
conda activate chipseq
cd ~/chipseq_abf3
opencode
```

Start opencode from the activated environment - the agent then finds all the
tools.

opencode has two agents. You switch between them with `Shift+Tab`; always
look which one is active before you send a prompt:

- **plan** - reads files and proposes what to do, but does not change files
  (and asks before it runs a command). Use it to **think**: *"how would you
  do this step?"*
- **build** - writes scripts and runs commands. Use it to **do** the step you
  agreed on.

Working with an agent on a bigger task always goes like this: first plan,
check the plan, then let it run.

| Key / command | What it does |
|---|---|
| `Shift+Tab` | switch between the **plan** and the **build** agent |
| `Esc` | stop the agent when it goes wrong |
| `/undo` | remove your last message and the agent's answer (do not rely on it to restore files) |
| `/new` | start a new conversation (`AGENTS.md` is read again) |
| `/export` | open the conversation as Markdown in a text editor - save it as `conversation.md` |

The keys are for opencode 2.x. Older versions switch agents with `Tab`; if
a key does not work, type `/help`.

Warm-up (in the **build** agent):

1. Ask the agent to read `sample_table.tsv` and tell you how many samples
   there are per assay and treatment. Does it agree with what you found?
2. Ask it to prepare the genome for mapping (decompress, index for samtools
   and minimap2), and to tell you how many chromosomes the genome has and how
   long they are.

- Watch the commands it runs. Did it follow the rules from `AGENTS.md`
  (script in `scripts/`, log in `logs/`)?
- Check the chromosome lengths yourself: `cut -f 1,2 TAIR10.fa.fai`. Do they
  agree with what the agent said?

Then switch to the **plan** agent and let it review your outline. Do not
paste the plan into the prompt - refer to the file:

```text
Read plan.md, Part 1. Is a step missing? Are the steps in the right order?
Are the right samples used in each step? Do not change the file.
```

Do you agree with its comments? If yes, change `plan.md` yourself.

## Step 3 - Run the analysis step by step

Go through your outline one step at a time. Every step has the same four
moves:

1. **Describe the step** in Part 2 of `plan.md`: goal, input (rows of the
   sample table), expected result, how you check it. The first step is
   filled in as an example.
2. **Plan agent** - ask it how it would do the step:

   ```text
   Read plan.md, Step 2. Propose how to do it: the commands and settings,
   the output files, and what you will report to me. Do not run anything.
   ```

   Read the proposal. Ask about anything you do not understand. If you want
   something different, say so. Note what you changed in *Changed in the
   agent's proposal*.
3. **Build agent** - `Do Step 2 as proposed.` The agent writes a script into
   `scripts/`, runs it, reports the results and writes the script name into
   `plan.md` (you see it in VS Code at once).
4. **Check** one number yourself (`wc -l`, `samtools flagstat`, `head` of a
   file, a picture ...) and fill in *Result* and *My check*.

The agent is a weak model and does exactly what the step description says.
A weak and a good description of the same step:

```text
## Step 2 - Mapping
Goal: do the mapping
```

```text
## Step 2 - Mapping
Goal: map the trimmed reads to TAIR10; get % mapped and the number of reads
on Mt and Pt for every sample
Input: trimmed reads of the rows with assay = ChIP-seq or ChIP-seq control
(8 samples); merge the two ChIP-seq control rows into one file
Expected result: most reads map
How I check it: samtools flagstat of one BAM file myself
```

A good description says **what** and **why**, the **input as rows of the
sample table** (not file names), and **what you want to get back** (a table,
a number, a picture). Settings of the tools are in `AGENTS.md` - you do not
need to write them, but you should understand them when the plan agent
proposes them.

Once you know the routine, you can skip the plan agent for simple steps
(e.g. coverage tracks), but not for the steps that compare samples (peaks,
DiffBind, RNA-seq).

Rules for you:

- **One step at a time.** Read the agent's report before you describe the
  next step.
- **Read the script** the agent wrote (`scripts/`). Can you say what each
  command does?
- When the agent goes in circles (the same error repeatedly, installing
  software, changing the plan on its own): press `Esc`, tell it exactly what is
  wrong and what to do instead. If it does not help, ask the teacher.
- If a step really fails, the reference results are in
  `~/Desktop/Bioinformatics/data/chip_seq_abf3/precomputed/` (peaks, DiffBind
  results) - tell the agent to use them and continue.

Mapping and peak calling of all samples take a few minutes - the agent runs
them in the background. Use the time to look at the outputs of the previous
step.

### Optional - binding vs expression: your own specification

The last question - *are the genes with stronger ABF3 binding also induced by
ABA?* - has no ready script. The agent writes it in R, from **your**
specification. Once you have the DiffBind sites with the closest genes
(`sites_closest_gene.tsv`) and the DESeq2 results (`deseq2_ABA_vs_mock.tsv`),
write the specification as a new step into `plan.md`. Decide:

- which binding sites count - all, or only those close to a TSS (how close)?
- how you define the groups of genes: "gained binding", "unchanged binding",
  all genes (which FDR / fold change thresholds)?
- what "induced" and "repressed" means (DESeq2 `padj`, `log2FoldChange`)?
- what you want to see: a table with numbers / percentages per group, a plot
  (which one?), a statistical test (which one, comparing what)?
- a table of the individual genes with gained binding and their expression.

Ask the **plan** agent to show you the R script it would write. Check: does
the script do what you specified? Then let the **build** agent save and run
it. Take two genes from the output table and find them in both input files
yourself.

## Step 4 - Check and interpret the results

The agent cannot look at pictures or at IGV - this is your part. Answer the
questions in `answers.md` in your own words, using the numbers and pictures
you produced. Do not ask the agent to answer them.

**Reads, mapping, peaks**

- How many reads were mapped in each sample? Is there a sample that differs
  from the others?
- How many peaks were found in each sample? Compare your numbers with
  `precomputed/peaks/` (`wc -l`). If they differ, find out why.
- What fragment length did MACS3 estimate?

**Binding sites in IGV** - load `TAIR10.fa`, the annotation, bigWig files of
one ABA and one mock replicate and the control, and their peaks (see
[1.6 of the original exercise](./chip_seq_abf3.md#16-binding-sites-in-igv)).
Go to RD29B (AT5G52300), RAB18 (AT5G66400), ABI1 (AT4G26080), ACT2
(AT3G18780) and ABF3 (AT4G34000).

- Where are the peaks relative to the genes? At which genes does the binding
  increase after ABA?
- The signal over the whole ABF3 gene is huge in all ChIP samples and absent
  in the control. Is it a binding site? (How were the plants made?)

**Around genes and differential binding**

- Look at the TSS profile. Where is the maximum of the signal? Why are the
  ABA and mock curves almost the same?
- How many sites change significantly after ABA, and in which direction? Do
  the samples separate by treatment in the PCA plot (`diffbind_plots.pdf`)?
- Which genes are closest to the most significant gained sites? Look up two
  of them in [TAIR](https://www.arabidopsis.org/). Do they make sense for an
  ABA response?

**Motifs**

- What is the best de novo motif in the gained sites? Does it contain
  `ACGTG`? Which known motifs are on top, and from which TF family?

**Binding vs expression** *(optional)*

- What fraction of all genes is induced by ABA? What fraction of the genes
  with gained ABF3 binding, and with unchanged binding?
- Is gained binding a good predictor of induction? Give reasons why a gene
  with gained binding may not be induced, and why an induced gene may have no
  ABF3 binding.
- How did your choices (distance to TSS, thresholds) change the result? Try
  one other setting.

**About the agent**

- Find one claim of the agent (a number or a statement) that you checked
  yourself. Was it correct?
- Where did the agent get stuck or do something you did not ask for? How did
  you get it back on track?
- Which part of your plan turned out to be too vague for the agent? How would
  you write it next time?

## Optional - the final plan

At the end you know all the steps and what each of them produces. Write
Part 3 of `plan.md`: the instructions you would give an agent to run the
whole analysis at once, without you checking every step. What must be in
it so that the agent does not go wrong? Which checks would you ask the agent
to do itself, and which results would you still check yourself?

This is how agents are used on bigger projects - but a good plan for the
whole analysis is only possible when you understand the single steps.

## What to hand in

From `~/chipseq_abf3`:

- `plan.md` - your outline and all steps with results and your checks
  (and the final plan, if you wrote it),
- `scripts/` - the scripts the agent wrote,
- the conversation exported with `/export`,
- `answers.md` - your answers to the questions in Step 4.
