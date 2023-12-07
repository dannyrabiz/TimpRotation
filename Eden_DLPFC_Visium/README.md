## Long-read sequencing spatial transcriptomics with Visium to characterize the Dorsal Lateral Pre-Frontal Cortex
Project led by Hope Eden 

### General Pipeline

### Problems/Thoughs
I initally used the Rockfish cluster to perform the step where I extract specific reads from Bam file and sort them into different files based on assigned cluster. 
However, as the rest of the Timp lab uses their own servers I tried to re-write this script to be run on the servers. The script I came up with uses GNU parallel to 
process many barcodes at the same time, nevertheless, it takes much longer to run that the previous script on the HPC (where I submitted a unique job for every single cluster 
to be run simultanely on different threads). Going forward if I were to further improve this pipeline, specifically this step, I would have to create an index or 
mapping of the barcodes and split the bam file up into smaller files and have a map/dictionary on where to locate the read for each barcode. I imagine this would be a lot faster 
assuming I could index it problerly.
