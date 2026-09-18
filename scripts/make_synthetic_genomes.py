#!/usr/bin/env python3
"""Generate three synthetic ~30-35 kb "bacteriophage" genomes for the whole-genome dotplot
exercise (docs/sequence_alignment.md, Exercise 1.11).

All sequence material is cut from a real bacterial genome (by default the
Mycoplasma hyopneumoniae 232 genome used in Exercise 1.10) so that base
composition looks realistic. Three "syntenic" blocks (S1, S2, S3; 3-5 kb each)
are shared by all genomes, everything between them is genome-specific
(non-overlapping chunks of the source genome that share no 16-mers with any
other chunk, so no unintended diagonals appear in a dotplot).

  phage_A  (reference):  V  S1  V  S2  V  S3  V
  phage_B  (rearranged): V  S2  V  S3  V  S1  V
  phage_C  (inversion):  V  S1  V  S2  V  rc(S3)  V

Copies of the syntenic blocks in B and C carry ~3 % substitutions and a few
short indels relative to A, so they are similar but not identical.

Usage:
  wget https://zenodo.org/record/4485547/files/mycoplasma-232.fasta
  python3 scripts/make_synthetic_genomes.py mycoplasma-232.fasta data/phage_genomes

The layout (coordinates of every block in every genome) is printed to stdout;
it is the answer key for the exercise.
"""
import os
import random
import sys

SEED = 20260918
K_CHECK = 16          # k-mer size used to make sure chunks are unrelated
MAX_SHARED = 2        # tolerated number of shared k-mers between chunks
SUBST_RATE = 0.03     # substitutions in S blocks of genomes B and C
N_INDELS = 4          # short indels per S block in genomes B and C
MAX_INDEL = 12

# block name -> length (bp); S* are syntenic, V* are genome-specific
SYNTENIC = {"S1": 4000, "S2": 3500, "S3": 5000}
LAYOUT = {
    "phage_A": [("V", 4500), ("S1", 0), ("V", 5000), ("S2", 0), ("V", 4000), ("S3", 0), ("V", 5500)],
    "phage_B": [("V", 5500), ("S2", 0), ("V", 4500), ("S3", 0), ("V", 6000), ("S1", 0), ("V", 4000)],
    "phage_C": [("V", 5000), ("S1", 0), ("V", 4000), ("S2", 0), ("V", 5500), ("S3", 0), ("V", 5000)],
}
INVERTED = {("phage_C", "S3")}

COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def read_fasta(path):
    seq = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()


def revcomp(s):
    return s.translate(COMP)[::-1]


def kmers(s, k=K_CHECK):
    return {s[i:i + k] for i in range(len(s) - k + 1)}


class ChunkPicker:
    """Cut non-overlapping chunks from the source genome that share
    (almost) no k-mers with previously picked chunks on either strand."""

    def __init__(self, source, rng):
        self.source = source
        self.rng = rng
        self.used = []          # (start, end) intervals already taken
        self.seen = set()       # k-mers of all chunks picked so far

    def pick(self, length):
        for _ in range(10000):
            start = self.rng.randrange(0, len(self.source) - length)
            end = start + length
            if any(s < end and start < e for s, e in self.used):
                continue
            chunk = self.source[start:end]
            if set(chunk) - set("ACGT"):
                continue
            km = kmers(chunk) | kmers(revcomp(chunk))
            if len(km & self.seen) > MAX_SHARED:
                continue
            self.used.append((start, end))
            self.seen |= km
            return chunk
        raise RuntimeError("could not find an unrelated chunk of length %d" % length)


def mutate(seq, rng):
    """Introduce point substitutions and a few short indels."""
    s = list(seq)
    for i in range(len(s)):
        if rng.random() < SUBST_RATE:
            s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    for _ in range(N_INDELS):
        pos = rng.randrange(100, len(s) - 100)
        n = rng.randint(1, MAX_INDEL)
        if rng.random() < 0.5:
            del s[pos:pos + n]
        else:
            s[pos:pos] = rng.choices("ACGT", k=n)
    return "".join(s)


def write_fasta(path, name, seq, width=60):
    with open(path, "w") as fh:
        fh.write(">%s\n" % name)
        for i in range(0, len(seq), width):
            fh.write(seq[i:i + width] + "\n")


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    source_path, out_dir = sys.argv[1:]
    os.makedirs(out_dir, exist_ok=True)
    rng = random.Random(SEED)
    source = read_fasta(source_path)
    picker = ChunkPicker(source, rng)

    syntenic = {name: picker.pick(length) for name, length in SYNTENIC.items()}

    print("genome\tblock\tstart\tend\tlength\tstrand")
    for genome, blocks in LAYOUT.items():
        parts = []
        pos = 0
        for name, length in blocks:
            if name == "V":
                seq = picker.pick(length)
                strand = "."
            else:
                seq = syntenic[name]
                if genome != "phage_A":
                    seq = mutate(seq, rng)
                strand = "+"
                if (genome, name) in INVERTED:
                    seq = revcomp(seq)
                    strand = "-"
            parts.append(seq)
            print("%s\t%s\t%d\t%d\t%d\t%s" % (genome, name, pos + 1, pos + len(seq), len(seq), strand))
            pos += len(seq)
        write_fasta(os.path.join(out_dir, genome + ".fasta"), genome, "".join(parts))


if __name__ == "__main__":
    main()
