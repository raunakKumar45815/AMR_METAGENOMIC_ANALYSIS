#!/bin/bash

# Directory containing FASTQ files
input_dir="/data/Students/Raunak/sdm-dharwad/wmg/new_samples/baseline/1st_part/output/trim/03_host_removed/rem"
# Directory to store Kraken 2 results
output_dir="/data/Students/Raunak/sdm-dharwad/wmg/new_samples/baseline/1st_part/output/trim/kraken2"

path="/data/Students/Raunak/sdm-dharwad/wmg/taxonomy/kraken2"

# Kraken 2 database
kraken2_db="/data/Students/Raunak/all_database/database_taxonomy/kraken_standard"

# Loop through each FASTQ file in the input directory
for fastq_file in "$input_dir"/*.clean.fastq.gz; do
    # Extracting file name without exension
    file_name=$(basename "$fastq_file" .clean.fastq.gz)
    
    # Run Kraken 2 on the current file
    "$path" --db "$kraken2_db" --threads 20 --use-names --confidence 0.1 --output "$output_dir"/"$file_name"_kraken2_output --report "$output_dir"/"$file_name"_kraken2_report "$fastq_file"
done
