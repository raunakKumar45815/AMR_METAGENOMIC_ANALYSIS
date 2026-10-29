#!/bin/bash
set -euo pipefail

KRAKEN_OUT_DIR="/data/Students/Raunak/Chilka_lake_WMG/AMR_deeparg/deeparg_mapping/reports/amr_fastq/kraken_amr/pk_out"     # folder with *.kraken files
OUTDIR="kraken_read_taxid"

mkdir -p "$OUTDIR"

for kout in ${KRAKEN_OUT_DIR}/*_kraken2_output; do

    sample=$(basename "$kout" _kraken2_output)
    outfile="${OUTDIR}/${sample}_read_taxid.tsv"

    echo "▶ Processing $sample"

    echo -e "read_id\ttaxid" > "$outfile"

    awk '
    {
        read_id = $2

        if (match($0, /taxid[[:space:]]+[0-9]+/)) {
            taxid = substr($0, RSTART+6, RLENGTH-6)
            print read_id "\t" taxid
        }
    }' "$kout" >> "$outfile"

    echo "✔ $sample done"
done

echo "🎉 All samples processed successfully"

