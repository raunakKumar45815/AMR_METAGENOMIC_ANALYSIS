#!/bin/bash

# Path to your folder containing all 120 fastq.gz files
INPUT_DIR="/data/Students/Raunak/RF-AMR-23-STP/test"
OUTPUT_DIR="/data/Students/Raunak/RF-AMR-23-STP/test"
HTML_DIR="/data/Students/Raunak/RF-AMR-23-STP/test"

mkdir -p "$OUTPUT_DIR"  # make output folder if it doesn't exist
mkdir -p "$HTML_DIR"
# Loop through all fastq.gz files
for file in "$INPUT_DIR"/*.fastq.gz; do
    base=$(basename "$file" .fastq.gz)

    fastplong -i "$file" \
    	      --verbose \
              --qualified_quality_phred 8 \
              --length_required 100 \
              --thread 16 \
              -o "$OUTPUT_DIR/${base}.filtered.fastq.gz" \
              -h "$HTML_DIR/${base}.html"
done
