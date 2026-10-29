#!/usr/bin/env bash
set -euo pipefail

# -------- paths --------
DEEPARG_DIR="/data/Students/Raunak/RF-AMR-23-STP/deeparg/all_mapping_23/rest/reports/rest"
FASTQ_DIR="/data/Students/Raunak/RF-AMR-23-STP/filtered_reads/reads/filtered/fastq_files/rest"

OUT_ID_DIR="/data/Students/Raunak/RF-AMR-23-STP/AMR-fastq-sankey/amr_read_ids"
OUT_FQ_DIR="/data/Students/Raunak/RF-AMR-23-STP/AMR-fastq-sankey/amr_fastq"

mkdir -p "${OUT_ID_DIR}" "${OUT_FQ_DIR}"

# -------- main loop --------
for rpt in ${DEEPARG_DIR}/*.mapping.ARG; do

    # Extract sample name (BH_242 from BH_242_deeparg_output.txt.mapping.ARG)
    sample=$(basename "$rpt" .mapping.ARG)

    echo "Processing sample: ${sample}"

    # 1. Extract Nanopore read IDs (column 4), skip header
    awk -F'\t' 'NR>1 && $4 != "" {print $4}' "$rpt" | \
        sort -u \
        > "${OUT_ID_DIR}/${sample}.amr.read_ids.txt"

    # 2. Subset Nanopore FASTQ using read IDs
    seqtk subseq \
        "${FASTQ_DIR}/${sample}.fastq.gz" \
        "${OUT_ID_DIR}/${sample}.amr.read_ids.txt" \
        > "${OUT_FQ_DIR}/${sample}.amr.fastq"

    echo "  AMR reads extracted: $(wc -l < ${OUT_ID_DIR}/${sample}.amr.read_ids.txt)"
done

echo "All samples processed. AMR-only FASTQs are ready for Kraken2."

