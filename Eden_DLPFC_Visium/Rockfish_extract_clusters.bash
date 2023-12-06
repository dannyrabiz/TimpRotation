#!/bin/bash

#tail -n +2 cluster_assignments.csv | split -d -l $((4046 / 20)) - --additional-suffix=.csv split_file_

input_file="test100.csv"
output_dir="output"
bam_file="merged_output.bam"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Skip the header and read the rest of the lines
tail -n +2 "$input_file" | while IFS=, read -r barcode cluster; do
  # Remove quotes from barcode and cluster
  barcode=$(echo "$barcode" | tr -d '"')
  cluster=$(echo "$cluster" | tr -d '"')

  # Create cluster directory if it doesn't exist
  cluster_dir="$output_dir/Cluster$cluster"
  mkdir -p "$cluster_dir"

  #echo Cluster$cluster/$barcode.bam
  # Run samtools command
  sbatch -t 30:00 --wrap='ml samtools; samtools view -b -d "CB:Z:$barcode" merged_output.bam > "Cluster$cluster/$barcode.bam"'
done
