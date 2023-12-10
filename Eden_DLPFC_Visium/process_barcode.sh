#!/bin/bash

barcode="$1"
cluster="$2"
bam_path="$3"
output_dir="$4"

# Remove quotes
barcode=$(echo "$barcode" | tr -d '"')
cluster=$(echo "$cluster" | tr -d '"')

# Run samtools
samtools view -@ 24 -d "CB:$barcode" "$bam_path" >> "$output_dir/cluster.$cluster.sam"

# Echo the barcode
echo "$barcode"
