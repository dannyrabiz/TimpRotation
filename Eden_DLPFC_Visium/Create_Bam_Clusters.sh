"""
Pipeline to create BAM files for each cluster identified in BayesSapce
"""


#!/bin/bash

# Function to process each sample
process_sample() {
    output_dir="$1"
    bam="$2"
    cluster_assignments="$3"
    k="$4"

    # Create the output directory
    mkdir -p "$output_dir"

    # Loop from 1 to k and process each cluster
    for (( i=1; i<=k; i++ )); do
        # Extract barcodes for the current cluster and remove quotes
        grep "$i" "$cluster_assignments" | tr -d '"' | cut -d, -f1 > barcodes.txt

        # Run samtools to extract reads and output to the desired file
        samtools view -@ 24 -D CB:barcodes.txt "$bam" -b -o "$output_dir/Cluster$i.bam"
	echo "Sample: $bam, Cluster:$i"
    done
}

# Skip the header line and read each line from the CSV file
tail -n +2 temp_sample_paths.csv | while IFS=, read -r output_dir bam cluster_assignments k; do
    # Process each sample
    process_sample "$output_dir" "$bam" "$cluster_assignments" "$k"
done
