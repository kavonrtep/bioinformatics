# Analysis plan - ChIP-seq of ABF3

Name:

## Biological questions

1.
2.
3.

## Steps

Copy the block below for each step. One step = one prompt to the agent.

### Step N - <short name>

- **Question / goal:** what do I want to learn from this step?
- **Input:** which files? (names, a directory, or rows of `sample_table.tsv`,
  e.g. assay = RNA-seq)
- **Tool and key settings:** which program, important parameters (see the
  tool notes in `AGENTS.md`)
- **Output:** which files, in which directory?
- **Check:** how will I know the step worked? (a number, a file, a picture
  that I look at myself)
- **Expected result:** what do I expect, and what would surprise me?

---

### Example: Step 1 - Quality control of reads

- **Question / goal:** Are the reads of good quality? Do I need to trim them?
  Do all samples have a similar number of reads?
- **Input:** the eight rows of `sample_table.tsv` with assay = ChIP-seq or
  ChIP-seq control (`ABF3_*` and `control_*` files, not the RNA-seq)
- **Tool and key settings:** fastp, single-end, 4 threads
- **Output:** trimmed reads in `trimmed/`, reports `qc/*.html` and
  `qc/*.json`; a table with reads before and after filtering for every sample
- **Check:** eight HTML reports exist; I open one report myself and compare
  its read count with the agent's table
- **Expected result:** about 3 M reads per sample (2.6 M in `control_rep1`),
  only a small fraction of reads removed; a sample that loses many more
  reads than the others would be suspicious
