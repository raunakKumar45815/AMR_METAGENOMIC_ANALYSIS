#!/bin/bash
set -euo pipefail

INDIR="/data/Students/Raunak/Chilka_lake_WMG/AMR_deeparg/deeparg_mapping/reports/amr_fastq/kraken_amr/pk_out/ARG_host_mapping"
OUTFILE="ARG_host_master.tsv"

# Write header
echo -e "sample\tread_id\tARG\tARG_class\ttaxon_level\ttaxon_name" > "$OUTFILE"

for file in ${INDIR}/*_ARG_host.tsv; do

    sample=$(basename "$file" _ARG_host.tsv)

    echo "▶ Adding $sample"

    awk -v sample="$sample" '
        NR>1 {
            print sample "\t" $0
        }
    ' "$file" >> "$OUTFILE"

done

echo "🎉 Master ARG–host table created: $OUTFILE"
