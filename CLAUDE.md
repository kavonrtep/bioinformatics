# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

Teaching materials for a university bioinformatics practicals course (Galaxy training group `kmb613-25`). It is a content repo, not a software project: there is no build, lint, or test step. The deliverables are the Markdown exercise documents in `docs/` and the example data in `data/`.

Students clone this repo to `~/Desktop/Bioinformatics` on lab machines (or use the `kavonrtep/bioinformatics` Docker image via `.gitpod.yml`) and run `git pull` before each session. Exercises reference data with that absolute path, e.g. `~/Desktop/Bioinformatics/data/dotter_sequences2`. Keep any new data under `data/<topic>/` and reference it the same way.

## Layout

- `docs/Practicals.md` — the week-by-week schedule. It is a list of links to exercise anchors in the other docs (`./sequence_alignment.md/#exercise-31---...`). **Renaming an exercise heading breaks these links** — update `Practicals.md` (and `README.md`) whenever a heading changes. Anchors follow GitHub slug rules (lowercase, punctuation stripped, spaces → `-`).
- `docs/*.md` — exercise documents, one per topic (`sequence_alignment.md`, `blast_search.md`, `sequence_assembly.md`, `chip_seq_analysis.md`, `molstar.md`, `pdb_database.md`, ...). Exercises are `### Exercise N.M - Title` headings grouped under `##` topic sections, each with a task description, data location, commands to run, and questions for students.
- `docs/*.org` — older Org-mode versions of some of the same documents (last edited 2024). The `.md` files are the maintained ones; do not edit `.org` files unless asked. `README.md` still links a few `.org`-only topics (mapping, databases, transcriptomics, phylogenetics tools, protein structure).
- `docs/img/`, `fig/` — images embedded in the docs.
- `data/` — input files for exercises (FASTA, FASTQ, GTF, PDB, etc.), grouped by topic directory.
- `scripts/` — reference shell scripts shown in exercises: MetaCentrum PBS job scripts (`#PBS` headers + `module add ...`) and one-off tool install scripts. Not meant to be run from this repo.
- `bin/` — bundled AliView installer (jar + `install.sh`).
- `shell-lesson-data/` — Software Carpentry shell lesson data used by `docs/shell_introduction.md`.

## Environments referenced in exercises

Commands in the docs assume one of three environments; be consistent with the one the surrounding exercise uses:

- Lab desktop / Docker image: tools are in conda envs activated with `conda activate <env>` (`assembly`, `gepard`, `quast`, `bmge`, `bowtie`, ...), plus GUI tools (`dotter`, `aliview`, Jalview, IGV, Mol*).
- MetaCentrum cluster: PBS scripts submitted with `qsub`, software loaded with `module add <name-version>` (see `docs/using_metacentrum.md`).
- Galaxy (`https://usegalaxy.eu/`): web-based workflows, see `docs/Introduction_to_Galaxy.md`.

## Conventions

- Commit messages are short and lowercase, typically `<topic> updated` / `<topic> added`.
- Exercises are written for students: state where the data is, give the exact command, then ask questions. Match the existing tone and heading numbering within the file.
