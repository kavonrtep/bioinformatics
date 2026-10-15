#!/bin/bash
# Rebuild query files, backup databases and pre-computed NCBI results for the
# BLAST exercises 1.7-1.11 (docs/blast_search.md) and exercise 3.7
# (docs/sequence_alignment.md).
#
# Not meant to be run by students. Requires BLAST+, seqkit, python3, curl and
# internet access.
#
# The pre-computed NCBI results were produced by running the searches on the
# NCBI BLAST web page with exactly the settings given in the exercises (BLAST+
# "-remote" and scripted URL-API jobs are unreliable: they hang or wait in the
# queue for hours, and protein searches against full nr restricted to a taxon
# fail with "CPU usage limit was exceeded"). The results were then downloaded
# by their RID with blast_formatter (see the RIDS table below). RIDs expire
# after ~36 hours - to refresh the results, re-run the searches on the web,
# put the new RIDs into the table and run this script with --from-rids.
#
# Backup databases derived from BLAST hits (mcr-1 regions, MCR-like proteins,
# HERV-W loci, syncytin-like proteins, bacterial globins) are rebuilt from the
# pre-computed tables in data/blast_data/*/precomputed/.
#
# Usage: scripts/build_blast_backup_data.sh <repo_dir> [--from-rids]
set -euo pipefail

REPO=$(realpath "${1:?usage: $0 <repo_dir> [--from-rids]}")
FROM_RIDS=${2:-}
D=$REPO/data/blast_data
E=https://eutils.ncbi.nlm.nih.gov/entrez/eutils
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

efetch_nt() { # accession [start end] -> fasta on stdout
    local range=""
    [ $# -ge 3 ] && range="&seq_start=$2&seq_stop=$3"
    curl -sf "$E/efetch.fcgi?db=nuccore&id=$1&rettype=fasta$range"
    sleep 0.4
}

# columns of the pre-computed tables (identical to what is shown on the web page)
COLS="qseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs staxids stitle"

precompute() { # name dir rid n_descriptions n_alignments ; tabular + pairwise text
    local name=$1 dir=$2 rid=$3 nd=$4 na=$5
    mkdir -p "$dir/precomputed"
    blast_formatter -rid "$rid" -outfmt 11 -out "$WORK/$name.asn"
    {
        echo "# NCBI web BLAST search run on $(date +%F) (RID $rid)"
        echo "# ${COLS// /	}"
        blast_formatter -archive "$WORK/$name.asn" -outfmt "6 $COLS"
    } > "$dir/precomputed/$name.tsv"
    blast_formatter -archive "$WORK/$name.asn" -outfmt 0 \
        -num_descriptions "$nd" -num_alignments "$na" > "$dir/precomputed/$name.txt"
    sleep 5
}

precompute_web() { # name dir rid n_alignments ; for archives blast_formatter cannot read
    local name=$1 dir=$2 rid=$3 na=$4 url="https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=$rid&FORMAT_TYPE=Text"
    mkdir -p "$dir/precomputed"
    {
        echo "# NCBI web BLAST search run on $(date +%F) (RID $rid)"
        curl -sf "$url&ALIGNMENT_VIEW=Tabular&ALIGNMENTS=5000&DESCRIPTIONS=500" \
            | grep -v -e 'Status=' -e 'QBlastInfo' -e '^<' -e '^-->' -e '^$'
    } > "$dir/precomputed/$name.tsv"
    sleep 5
    curl -sf "$url&ALIGNMENTS=$na&DESCRIPTIONS=100" | sed -n '/^T\?BLAST/,$p' > "$dir/precomputed/$name.txt"
    sleep 5
}

fetch_proteins() { # file with accessions -> fasta on stdout (batches of 100)
    split -l 100 "$1" "$WORK/ids_"
    for f in "$WORK"/ids_*; do
        curl -sf "$E/efetch.fcgi?db=protein&rettype=fasta&id=$(paste -sd, "$f")"
        sleep 1
    done
    rm "$WORK"/ids_*
}

hit_ids() { # precomputed tsv -> unique subject accessions in hit order
    grep -v '^#' "$1" | cut -f2 | awk '!seen[$1]++'
}

## gyrA (exercise 1.7)
mkdir -p "$D/gyrA"
efetch_nt NZ_CP021520.1 > "$D/gyrA/NZ_CP021520.1.fasta"
awk '/^>strain_1/,/^```$/' "$REPO/docs/blast_search.md" | grep -v '```' > "$D/gyrA/gyrA_strains.fasta"

## mcr-1 (exercise 1.8)
mkdir -p "$D/mcr1/backup_db"
curl -sf "$E/efetch.fcgi?db=nuccore&id=KP347127.1&rettype=fasta_cds_na" \
    | seqkit grep -r -n -p "gene=mcr-1" \
    | sed '1s/.*/>mcr-1 KP347127.1:22413-24038 Escherichia coli SHP45 plasmid pHNSHP45, mcr-1 gene (phosphoethanolamine transferase)/' \
    > "$D/mcr1/mcr1_cds.fasta"
curl -sf "$E/efetch.fcgi?db=protein&id=AKF16168.1&rettype=fasta" \
    | sed '1s/.*/>MCR-1 AKF16168.1 phosphoethanolamine transferase MCR-1 [Escherichia coli]/' \
    > "$D/mcr1/mcr1_protein.fasta"

## SARS-CoV-2 (exercise 1.9)
mkdir -p "$D/sars_cov2/backup_db"
efetch_nt NC_045512.2 > "$D/sars_cov2/sars_cov2_genome.fasta"
# spike RBD, aa 319-541
efetch_nt NC_045512.2 22517 23185 \
    | sed '1s/.*/>SARS-CoV-2_RBD NC_045512.2:22517-23185 spike receptor-binding domain (aa 319-541)/' \
    > "$D/sars_cov2/sars_cov2_rbd.fasta"
# SARS-CoV-2 relatives (bat, pangolin), SARS-CoV-1 and other sarbecoviruses, MERS-CoV as outgroup
COV=MN996532.2,MZ937000.1,MZ937001.1,MZ937002.1,MZ937003.2,MZ937004.1,MT040333.1,MT040334.1,MT121216.1
COV=$COV,NC_004718.3,AY304486.1,KF367457.1,KC881005.1,MG772933.1,MG772934.1,NC_014470.1,KY352407.1
COV=$COV,MW251308.1,MZ081381.1,MZ190138.1,OP963576.1,OR233302.1,NC_019843.3
efetch_nt "$COV" > "$D/sars_cov2/backup_db/coronavirus_genomes.fasta"

## syncytin (exercise 1.10)
mkdir -p "$D/syncytin/backup_db"
curl -sf https://rest.uniprot.org/uniprotkb/Q9UQF0.fasta > "$D/syncytin/syncytin1_human.fasta"
curl -sf https://rest.uniprot.org/uniprotkb/Q5G5D5.fasta > "$D/syncytin/syncytinA_mouse.fasta"

## contamination (exercise 1.11)
mkdir -p "$D/contamination/backup_db"
F=$WORK/fragments
mkdir -p "$F"
cp "$REPO/data/genome_assembly/GCF_000013425.1_ASM1342v1_genomic.fna" "$F/saureus.fa"
efetch_nt NC_000913.3 361000 369000 > "$F/ecoli.fa"
efetch_nt NC_000913.3 2000000 2003500 > "$F/ecoli2.fa"
efetch_nt NC_012920.1 3000 9000 > "$F/human_mt.fa"
efetch_nt NC_001422.1 > "$F/phix.fa"
efetch_nt L09137.2 > "$F/puc19.fa"
efetch_nt NC_006085.1 1200000 1206000 > "$F/cacnes.fa"
# answer key goes to stdout - it is not stored in the repository
echo "== contamination answer key =="
python3 "$REPO/scripts/make_contamination_contigs.py" "$F" "$D/contamination/draft_assembly_contigs.fasta"
{
    efetch_nt NC_000913.3 300001 400000
    efetch_nt NC_000913.3 1950001 2050000
    efetch_nt NC_012920.1
    efetch_nt NC_001422.1
    efetch_nt NC_006085.1 1150001 1250000
    efetch_nt NC_004461.1
} | sed -E 's/^>(N[CZ]_[0-9.]+):([0-9]+)-([0-9]+) (.*), complete (genome|sequence)/>\1_\2-\3 \4, chromosome region \2-\3/' \
    > "$D/contamination/backup_db/reference_sequences.fasta"
curl -sf https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core > "$D/contamination/backup_db/UniVec_Core.fasta"

## NCBI web searches -> pre-computed results
# RIDs of the searches run on 2026-10-15 (settings as described in the exercises):
#   mcr1_blastn_core_nt                  blastn megablast, core_nt, max 500 targets
#   mcr1_blastp_clusterednr              blastp, ClusteredNR, max 250 targets
#   sars2_genome_blastn_core_nt          megablast, core_nt, exclude txid2697049 + txid32630, max 250
#   sars2_rbd_blastn_core_nt             blastn, subrange 22517-23185, same exclusions, max 250
#   syncytin1_tblastn_human_genome       tblastn, refseq_genomes, Homo sapiens
#   syncytin1_blastp_clusterednr         blastp, ClusteredNR, max 250
#   contigs_blastn_core_nt               megablast, core_nt, max 20 targets
#   neuroglobin_blastp_clusterednr_bacteria  blastp, ClusteredNR, Bacteria, max 50, E 1e-5
if [ "$FROM_RIDS" = "--from-rids" ]; then
    precompute mcr1_blastn_core_nt "$D/mcr1" D290VNY0014 100 10
    precompute mcr1_blastp_clusterednr "$D/mcr1" D2E442JK014 250 20
    precompute sars2_genome_blastn_core_nt "$D/sars_cov2" D29E4AN5014 250 0
    precompute sars2_rbd_blastn_core_nt "$D/sars_cov2" D29EJFM6014 250 30
    # blast_formatter (2.12) crashes on this archive - download from the web instead
    precompute_web syncytin1_tblastn_human_genome "$D/syncytin" D2DGVN2G014 10
    precompute syncytin1_blastp_clusterednr "$D/syncytin" D2CKM6MU014 250 20
    precompute contigs_blastn_core_nt "$D/contamination" D293S1HB014 20 1
    precompute neuroglobin_blastp_clusterednr_bacteria "$D/neuroglobin" D2E2200B014 50 10
fi

## backup databases derived from the pre-computed hits
# mcr-1: hit regions +-1000 bp, up to 8 plasmid and 3 chromosomal hits per genus
grep -v '^#' "$D/mcr1/precomputed/mcr1_blastn_core_nt.tsv" | python3 -c '
import sys, collections
seen = collections.Counter()
for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    acc, s, e, title = f[1], int(f[8]), int(f[9]), f[14]
    genus = title.split()[0]
    if genus not in ("Escherichia", "Salmonella", "Klebsiella", "Shigella", "Citrobacter"):
        continue
    if "DH5-alpha" in title:  # laboratory construct
        continue
    kind = "plasmid" if "plasmid" in title.lower() else "chrom"
    if seen[(genus, kind)] >= (8 if kind == "plasmid" else 3):
        continue
    seen[(genus, kind)] += 1
    print(acc, max(1, min(s, e) - 1000), max(s, e) + 1000, title, sep="\t")
' | while IFS=$'\t' read -r acc s e title; do
    efetch_nt "$acc" "$s" "$e" | sed "1s/.*/>${acc}_${s}-${e} ${title}, region ${s}-${e}/"
done > "$D/mcr1/backup_db/mcr1_hits.fasta"

hit_ids "$D/mcr1/precomputed/mcr1_blastp_clusterednr.tsv" > "$WORK/mcr.ids"
fetch_proteins "$WORK/mcr.ids" > "$D/mcr1/backup_db/mcr_like_proteins.fasta"

# syncytin: HERV-W env loci = tblastn HSPs (bit score >= 100, primary chromosomes)
# +-500 bp, overlapping regions merged
grep -v '^#' "$D/syncytin/precomputed/syncytin1_tblastn_human_genome.tsv" \
    | awk -F'\t' '$12 >= 100 && $2 ~ /^NC_0000/ {
        s = ($9 < $10) ? $9 : $10; e = ($9 < $10) ? $10 : $9
        print $2 "\t" s - 500 "\t" e + 500 }' \
    | sort -k1,1 -k2,2n \
    | awk -F'\t' 'BEGIN{OFS="\t"} $1 == c && $2 <= e {if ($3 > e) e = $3; next}
        {if (c != "") print c, s, e; c = $1; s = $2; e = $3} END{print c, s, e}' \
    | while IFS=$'\t' read -r acc s e; do
        n=$(echo "$acc" | sed -E 's/NC_0+([0-9]+)\..*/\1/')
        case $n in 23) chr=X ;; 24) chr=Y ;; *) chr=$n ;; esac
        efetch_nt "$acc" "$s" "$e" \
            | sed "1s/.*/>${acc}_${s}-${e} Homo sapiens chromosome ${chr} (GRCh38), region ${s}-${e}/"
    done > "$D/syncytin/backup_db/herv_w_env_loci.fasta"

hit_ids "$D/syncytin/precomputed/syncytin1_blastp_clusterednr.tsv" > "$WORK/syn.ids"
fetch_proteins "$WORK/syn.ids" > "$D/syncytin/backup_db/syncytin1_blastp_hits.fasta"

# neuroglobin (exercise 3.7): bacterial globins found by NCBI BLAST
hit_ids "$D/neuroglobin/precomputed/neuroglobin_blastp_clusterednr_bacteria.tsv" > "$WORK/ngb.ids"
fetch_proteins "$WORK/ngb.ids" > "$D/neuroglobin/bacterial_globins_50.fasta"
head -10 "$WORK/ngb.ids" > "$WORK/ngb10.ids"
fetch_proteins "$WORK/ngb10.ids" > "$D/neuroglobin/bacterial_globins_top10.fasta"
