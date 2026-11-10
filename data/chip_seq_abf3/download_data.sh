#!/usr/bin/env bash
# Download the data for the ABF3 ChIP-seq exercise from Zenodo (1.9 GB) and
# check the files with md5sum. Run it before the practical:
#
#   bash ~/Desktop/Bioinformatics/data/chip_seq_abf3/download_data.sh
#
# Files go to ~/chipseq_abf3/data (or to the directory given as the first
# argument). The script can be run again: complete files are skipped,
# interrupted downloads continue, damaged files are downloaded again.

set -uo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MD5=$HERE/md5sums.txt
OUT=${1:-$HOME/chipseq_abf3/data}
Z=https://zenodo.org/records/23984926/files

mkdir -p "$OUT"
cd "$OUT" || exit 1
echo "Downloading to $OUT"

for attempt in 1 2 3; do
    failed=0
    while read -r sum f; do
        if [ -f "$f" ] && echo "$sum  $f" | md5sum -c --quiet - 2>/dev/null; then
            echo "OK          $f"
            continue
        fi
        echo "downloading $f"
        wget -q --show-progress -c "$Z/$f"
        if ! echo "$sum  $f" | md5sum -c --quiet - 2>/dev/null; then
            # a damaged file cannot be fixed by continuing - start again
            echo "FAILED      $f (will be downloaded again)"
            rm -f "$f"
            failed=1
        fi
    done < "$MD5"
    [ $failed -eq 0 ] && break
done

if [ $failed -eq 0 ]; then
    echo "All files downloaded and checked."
else
    echo "Some files failed after 3 attempts - run the script again later." >&2
    exit 1
fi
