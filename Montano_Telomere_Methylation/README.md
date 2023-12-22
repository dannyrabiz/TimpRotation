# Measuring telomere length x sub-telomeric methylation using Nanopore sequqncing
Project led by Dr. Carolina Montano 

![](https://github.com/dannyrabiz/TimpRotation/blob/main/Montano_Telomere_Methylation/chromosomes41.png)

## General Pipeline
### 1. Generate pod5 files from fast5
`pod5 convert fast5 path/to/file.fast5 -o output/allp5`

### 2. Peform modified basecalling with Dorado 
`dorado basecaller dna_r9.4.1_e8_sup@v3.3 Sample.pod5 --modified-bases-models dna_r9.4.1_e8_sup@v3.3_5mCG@v0.1> Sample.bam`

### 3. Extract telomere reads, convert to fastq
__Telomere reads are retrived from __Reads_table/__ folder__

`grep -F -f Sample_Telomere_Reads.txt *$sample*sam > Sample.Telomere_Only.sam`

__Keep Methylation tags from Sam file__

`samtools fastq -@24 -T MM,ML Sample.Telomere_Only.sam | gzip > Sample.fastq.gz`

## Additional Tools that need to be added:
[Tidehunter:](https://github.com/Xinglab/TideHunter) - Tandem repeat detection 

[Noise Cancelling Repeat Finder:](https://github.com/makovalab-psu/NoiseCancellingRepeatFinder) - Repeat detection algorithm used in the Salk paper

[Porechop:](https://github.com/rrwick/Porechop) - Adapter trimming 

Heng Li's Bonito model trained on Telomere reads. 


## Problems to work on and think about 

 1. Heterogeneity in telomere length. Each chromosome in each sample has a lot of variance. What's the best methodology to define length.
 2.  Salk paper counted number of reads with telomere length >10kb and length < 1kb.  

 3. How to compare methylation. How to plot it?
    - Salk paper compared patterns for each arm.
   
4. Should we distinguish between the alleles?

5. Salk paper excluded reads that "did not commence with the 5'-(CCCTAA)n motif or conclude with the 5'-(TTAGGG)n motif"

