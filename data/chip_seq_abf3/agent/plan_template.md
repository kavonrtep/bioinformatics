# Analysis plan - ChIP-seq of ABF3

Name:

Write plain text after the colons. Lines starting with # are headings - keep
them, the agent uses them to find the steps.


# Part 1 - Outline (before you start the agent)

Biological questions:
1.
2.
3.

Steps (one line each: what is done, with which rows of sample_table.tsv):
1. Quality control of reads - assay = ChIP-seq or ChIP-seq control
2.
3.
4.
5.
6.
7.
8.


# Part 2 - Steps (fill in one step at a time)

## Step 1 - Quality control of reads (example)

Goal: Are the reads good? Do all samples have a similar number of reads?
Input: rows of sample_table.tsv with assay = ChIP-seq or ChIP-seq control (8 files)
Expected result: about 3 M reads per sample (2.6 M in control_rep1), few reads removed
How I check it: open one fastp report (qc/) and compare its numbers with the agent's table

Changed in the agent's proposal:
Script (the agent fills in):
Result:
My check:


## Step 2 -

Goal:
Input:
Expected result:
How I check it:

Changed in the agent's proposal:
Script (the agent fills in):
Result:
My check:


(copy the block "## Step ..." for every next step)


# Part 3 - Final plan (optional, at the end)

Now that you know all the steps: write the instructions you would give an
agent to run the whole analysis at once, without you watching every step.
