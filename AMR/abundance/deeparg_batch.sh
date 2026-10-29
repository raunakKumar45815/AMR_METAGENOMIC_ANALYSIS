#!/bin/bash

input_directory="path_to_input_dir"

output_directory="path_to_output_dir"

deeparg_data_directory="path_database_dir"

for input_file in "$input_directory"/*.filtered.fastq.gz; do
    # Extract the base filename (e.g., sample1.fastq.gz → sample1)
    filename=$(basename -- "$input_file")
    filename_no_ext="${filename%%.*}"

    # Define the output file path
    output_file="$output_directory/${filename_no_ext}_deeparg_output.txt"

    # Run DeepARG prediction
    deeparg predict \
        --model LS \
        --type nucl \
        -d "$deeparg_data_directory" \
        -i "$input_file" \
        -o "$output_file"
done

echo "✅ DeepARG prediction completed for all files."

