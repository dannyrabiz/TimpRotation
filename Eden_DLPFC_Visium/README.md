# Long-read sequencing spatial transcriptomics with Visium to characterize the Dorsal Lateral Pre-Frontal Cortex
Project led by Hope Eden 

![](https://github.com/dannyrabiz/TimpRotation/blob/main/Eden_DLPFC_Visium/Br8667_mid.png)
![](https://github.com/dannyrabiz/TimpRotation/blob/main/Eden_DLPFC_Visium/Br6522_ant.png)
## General Pipeline
### 1. Run Sockeye
  Generate readcount matrices (gene and transcript) with annotated barcodes using Nanoporetech's [sockeye](https://github.com/nanoporetech/sockeye) pipeline. **Note:** Hope added her own formatting step, so you'll need to run those as well in order for the the next commands to work using the files in this github.
### 2. Run Bayesspace 
  2a. Create a file (control_file_gene.csv) with Sample ID's and paths to expression matrices. Here is the general format I used:
 
  ```
  ID,path,CountMat,ColData
  Sample1_ID, folder_path, matrix_path, barcode_bath
  Sample2_ID, folder_path, matrix_path, barcode_bath ` 
  ..
  ..
```
2b. Run the BayesSpace Rscript. Make sure to specify A) transcripts or genes and B) the number of clusters.  

### 3. Create BAM file for each cluster



## To Do:
1. Find out of how to compare clusteres across samples
2. How to flag false-positive novel isoforms
3. Compare different isoform callers
4. QC with Illumina data

## Problems/Thoughts
~~I initally used the Rockfish cluster to perform the step where I extract specific reads from Bam file and sort them into different files based on assigned cluster. 
However, as the rest of the Timp lab uses their own servers I tried to re-write this script to be run on the servers. The script I came up with uses GNU parallel to 
process many barcodes at the same time, nevertheless, it takes much longer to run that the previous script on the HPC (where I submitted a unique job for every single cluster 
to be run simultanely on different threads). Going forward if I were to further improve this pipeline, specifically this step, I would have to create an index or 
mapping of the barcodes and split the bam file up into smaller files and have a map/dictionary on where to locate the read for each barcode. I imagine this would be a lot faster 
assuming I could index it problerly.~~
