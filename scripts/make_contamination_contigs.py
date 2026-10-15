#!/usr/bin/env python3
"""Generate a small "draft assembly" with contaminated contigs for the BLAST
contamination-check exercise (docs/blast_search.md, Exercise 1.11).

Most contigs are fragments of the Staphylococcus aureus NCTC 8325 genome (the
genome assembled in the genome assembly practicals), the rest are typical
contaminants of sequencing projects:

  - PhiX174       Illumina spike-in control
  - pUC19         cloning vector
  - human mtDNA   DNA of the person who handled the sample
  - E. coli       lab strain / reagent contamination
  - Cutibacterium acnes  skin bacterium, common "kitome" contaminant
  - one chimeric contig (S. aureus + E. coli) = misassembly

Fake k-mer coverage is written to the FASTA header; contaminants typically
have coverage very different from the target genome.

Usage (input fragments are downloaded by scripts/build_blast_backup_data.sh):
  python3 scripts/make_contamination_contigs.py <fragment_dir> <output.fasta>

The answer key (contig -> origin) is printed to stdout.
"""
import os
import random
import sys

SEED = 20261015


def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def main(frag_dir, out_path):
    rng = random.Random(SEED)
    sa = read_fasta(os.path.join(frag_dir, "saureus.fa"))
    frag = {name: read_fasta(os.path.join(frag_dir, name + ".fa"))
            for name in ["ecoli", "ecoli2", "human_mt", "phix", "puc19", "cacnes"]}

    # (origin label, sequence, coverage)
    contigs = [
        ("S. aureus 100001-110000", sa[100000:110000], 118),
        ("S. aureus 600001-612000", sa[600000:612000], 124),
        ("S. aureus 1200001-1209000 (reverse complement)", revcomp(sa[1200000:1209000]), 121),
        ("S. aureus 2000001-2015000", sa[2000000:2015000], 115),
        ("S. aureus 2600001-2607000", sa[2600000:2607000], 130),
        ("chimera: S. aureus 1500001-1504000 + E. coli NC_000913.3:2000000-2003500",
         sa[1500000:1504000] + frag["ecoli2"], 67),
        ("E. coli K-12 NC_000913.3:361000-369000 (lac operon region)", frag["ecoli"], 6),
        ("human mitochondrion NC_012920.1:3000-9000", frag["human_mt"], 3),
        ("phiX174 NC_001422.1 (Illumina spike-in)", frag["phix"], 412),
        ("pUC19 cloning vector L09137.2", frag["puc19"], 25),
        ("Cutibacterium acnes NC_006085.1:1200000-1206000", frag["cacnes"], 4),
    ]
    rng.shuffle(contigs)
    with open(out_path, "w") as out:
        for i, (origin, seq, cov) in enumerate(contigs, 1):
            name = "contig_%02d" % i
            out.write(">%s len=%d cov=%dx\n" % (name, len(seq), cov))
            for j in range(0, len(seq), 70):
                out.write(seq[j:j + 70] + "\n")
            print("%s\t%d\t%s" % (name, len(seq), origin))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
