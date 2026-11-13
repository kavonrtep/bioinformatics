#!/usr/bin/env bash
# Prepare the working directory for the ABF3 ChIP-seq exercise with an agent.
# Run it INSIDE your (empty) working directory:
#
#   mkdir -p ~/chipseq_abf3
#   cd ~/chipseq_abf3
#   bash ~/Desktop/Bioinformatics/data/chip_seq_abf3/agent/setup.sh
#
# It copies the files for the agent, the plan template, the helper files and
# the reference results into the current directory, then downloads the data
# (1.9 GB) into data/. It can be run again: files you may have edited
# (plan.md, sample_table.*, AGENTS.md, opencode.json) are never overwritten,
# complete downloads are skipped.

set -euo pipefail

AGENT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)   # .../chip_seq_abf3/agent
TOPIC=$(dirname "$AGENT")                             # .../chip_seq_abf3
WD=$PWD

if [ "$WD" = "$HOME" ] || [[ "$WD/" == "$TOPIC/"* ]]; then
    echo "Run this script inside your working directory, e.g.:" >&2
    echo "  mkdir -p ~/chipseq_abf3 && cd ~/chipseq_abf3" >&2
    exit 1
fi
echo "Setting up $WD"

copy_new() {   # copy_new <source> <target> - never overwrite
    if [ -e "$2" ]; then
        echo "kept        $2 (already exists)"
    else
        cp "$1" "$2"
        echo "copied      $2"
    fi
}

copy_new "$AGENT/AGENTS.md"         AGENTS.md
copy_new "$AGENT/opencode.json"     opencode.json
copy_new "$AGENT/sample_table.tsv"  sample_table.tsv
copy_new "$AGENT/sample_table.ods"  sample_table.ods
copy_new "$AGENT/plan_template.md"  plan.md

# helper files and reference results (read-only for the agent)
mkdir -p help
for f in TAIR10_genes.bed TAIR10_tss.bed TAIR10_tx2gene.tsv diffbind.R deseq2.R; do
    cp "$TOPIC/$f" help/
done
echo "copied      help/ (genes, TSS, tx2gene, diffbind.R, deseq2.R)"
cp -r "$TOPIC/precomputed" .
echo "copied      precomputed/"

# the data from Zenodo, checked with md5sum
bash "$TOPIC/download_data.sh" "$WD/data"

echo
echo "Setup finished. Your working directory: $WD"
