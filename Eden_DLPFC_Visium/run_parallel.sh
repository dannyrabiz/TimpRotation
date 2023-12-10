#!/bin/bash

output_dir='/wilee/Data/drabiza1/Hope/BR8667_mid_k2_bams'
input_file="BR8667_mid_k2_cluster_assignments.csv"
bam_path="/wilee/Data/drabiza1/Hope/BAMS/br8667_mid.bam"

# Create output directory
mkdir -p "$output_dir"

# Use GNU Parallel to process each line
parallel --colsep ',' ./process_barcode.sh {1} {2} "$bam_path" "$output_dir" :::: "$input_file"
