#!/bin/bash


INPUT_DIR="path_to_input_dir"
DB="path_to_diamond_db"
OUTPUT_DIR="path_to_output"
LOG_DIR="$OUTPUT_DIR/logs"
THREADS=60

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

shopt -s nullglob
files=("$INPUT_DIR"/*.fastq)
if [ ${#files[@]} -eq 0 ]; then
    echo "ERROR: no .filtered.fastq.gz files found in $INPUT_DIR"; exit 1
fi
echo "Found ${#files[@]} samples."

for file in "${files[@]}"; do
    sample=$(basename "$file" .fastq)
    out="$OUTPUT_DIR/${sample}.tsv"

    if [ -s "$out" ]; then
        echo "SKIP  $sample (output exists)"; continue
    fi

    echo "RUN   $sample"
    diamond blastx \
        --db "$DB" \
        --query "$file" \
        --out "$out" \
        --outfmt 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore \
        --long-reads \
        --more-sensitive \
        --id 70 \
        --subject-cover 50 \
        --evalue 1e-10 \
        --threads "$THREADS" \
        > "$LOG_DIR/${sample}.log" 2>&1

    if [ $? -ne 0 ]; then
        echo "FAIL  $sample — see $LOG_DIR/${sample}.log"
        rm -f "$out"
    else
        echo "OK    $sample ($(wc -l < "$out") hits)"
    fi
done

echo "Done. Results in $OUTPUT_DIR"
