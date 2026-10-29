#!/bin/bash

INPUT_DIR="path_to_your_input_dir"

DB="path_to_your_database_dir"

OUTPUT_DIR="path_to_your_out_dir"
mkdir -p "$OUTPUT_DIR"

for file in "$INPUT_DIR"/*.fastq.gz; do
    sample=$(basename "$file" .fastq.gz)  # Extract sample name without extension
    
    echo "Running DIAMOND on $sample ..."
    
    diamond blastx \
        --db "$DB" \
        --query "$file" \
        --id 60 \
        --outfmt 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore \
        --header \
        --long-reads \
        --threads [number_of_threads] \
        --more-sensitive \
        --evalue 1e-10 \
        --out "$OUTPUT_DIR/${sample}.tsv"
done	

