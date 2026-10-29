#!/bin/bash
set -euo pipefail

DEEPARG_DIR="/data/Students/Raunak/Chilka_lake_WMG/AMR_deeparg/deeparg_mapping/reports/reports/pk"
TAXON_DIR="/data/Students/Raunak/Chilka_lake_WMG/AMR_deeparg/deeparg_mapping/reports/amr_fastq/kraken_amr/pk_out/kraken_read_taxoninfo"
OUTDIR="ARG_host_mapping"

mkdir -p "$OUTDIR"

for deeparg in ${DEEPARG_DIR}/*.mapping.ARG; do

    sample=$(basename "$deeparg" .mapping.ARG)
    taxon="${TAXON_DIR}/${sample}_taxoninfo.tsv"
    outfile="${OUTDIR}/${sample}_ARG_host.tsv"

    if [[ ! -f "$taxon" ]]; then
        echo "⚠️  Taxonomy file missing for $sample, skipping"
        continue
    fi

    echo "▶ Processing $sample"

    echo -e "read_id\tARG\tARG_class\ttaxon_level\ttaxon_name" > "$outfile"

    awk '
    NR==FNR {
        # Kraken-derived taxonomy file
        # skip header
        if (NR>1) {
            tax_level[$1] = $2
            tax_name[$1]  = $3
        }
        next
    }
    {
        # DeepARG mapping file (no header)
        arg       = $1
        read_id   = $4
        arg_class = $5

        if (read_id in tax_level) {
            print read_id "\t" arg "\t" arg_class "\t" tax_level[read_id] "\t" tax_name[read_id]
        }
    }
    ' "$taxon" "$deeparg" >> "$outfile"

    echo "✔ $sample done"
done

echo "🎉 ARG–host mapping completed for all samples"
