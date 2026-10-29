#!/bin/bash
set -euo pipefail

READ_TAXID_DIR="/data/Students/Raunak/Chilka_lake_WMG/AMR_deeparg/deeparg_mapping/reports/amr_fastq/kraken_amr/pk_out/kraken_read_taxid"
KRAKEN_REPORT_DIR="/data/Students/Raunak/Chilka_lake_WMG/AMR_deeparg/deeparg_mapping/reports/amr_fastq/kraken_amr/reports/pk"
OUTDIR="kraken_read_taxoninfo"

mkdir -p "$OUTDIR"

for tsv in ${READ_TAXID_DIR}/*_read_taxid.tsv; do

    sample=$(basename "$tsv" _read_taxid.tsv)
    report="${KRAKEN_REPORT_DIR}/${sample}_kraken2_report"
    outfile="${OUTDIR}/${sample}_taxoninfo.tsv"

    if [[ ! -f "$report" ]]; then
        echo "⚠️  Kraken report missing for $sample, skipping"
        continue
    fi

    echo "▶ Processing $sample"

    echo -e "read_id\ttaxon_level\ttaxon_name" > "$outfile"

    awk '
    NR==FNR {
        # Kraken report
        taxid = $5
        level[taxid] = $4

        # remove leading spaces from column 6 (taxon name)
        name = $6
        sub(/^[ \t]+/, "", name)
        taxname[taxid] = name
        next
    }
    NR>1 {
        # read_id-taxid TSV (skip header)
        read_id = $1
        taxid = $2

        if (taxid in level) {
            print read_id "\t" level[taxid] "\t" taxname[taxid]
        }
    }
    ' "$report" "$tsv" >> "$outfile"

    echo "✔ $sample done"
done

echo "🎉 All samples processed successfully"
