"""
Basecalling with dorado and subsequent processing steps. 
I wasnt't sure if you can run GPU steps in parallel (assume not), so I ran them one by one. 
Also included subsequent formatting steps for demulitplexing prior to alignment
"""


/home/drabiza/./dorado-0.4.3-linux-x64/bin/dorado basecaller /home/drabiza/dorado-0.4.3-linux-x64/bin/dna_r9.4.1_e8_sup@v3.3 JH101_JH102/output.pod5 --modified-bases-models dna_r9.4.1_e8_sup@v3.3_5mCG@v0.1> JH101_JH102.bam 
/home/drabiza/./dorado-0.4.3-linux-x64/bin/dorado basecaller /home/drabiza/dorado-0.4.3-linux-x64/bin/dna_r9.4.1_e8_sup@v3.3 JH103_JH104/output.pod5 --modified-bases-models dna_r9.4.1_e8_sup@v3.3_5mCG@v0.1> JH103_JH104.bam
/home/drabiza/./dorado-0.4.3-linux-x64/bin/dorado basecaller /home/drabiza/dorado-0.4.3-linux-x64/bin/dna_r9.4.1_e8_sup@v3.3 JH76_JH77/output.pod5 --modified-bases-models dna_r9.4.1_e8_sup@v3.3_5mCG@v0.1> JH76_JH77.bam
/home/drabiza/./dorado-0.4.3-linux-x64/bin/dorado basecaller /home/drabiza/dorado-0.4.3-linux-x64/bin/dna_r9.4.1_e8_sup@v3.3 JH80_JH81/output.pod5 --modified-bases-models dna_r9.4.1_e8_sup@v3.3_5mCG@v0.1> JH80_JH81.bam
/home/drabiza/./dorado-0.4.3-linux-x64/bin/dorado basecaller /home/drabiza/dorado-0.4.3-linux-x64/bin/dna_r9.4.1_e8_sup@v3.3 JH82_JH88/output.pod5 --modified-bases-models dna_r9.4.1_e8_sup@v3.3_5mCG@v0.1> JH82_JH88.bam
/home/drabiza/./dorado-0.4.3-linux-x64/bin/dorado basecaller /home/drabiza/dorado-0.4.3-linux-x64/bin/dna_r9.4.1_e8_sup@v3.3 JH89_JH90/output.pod5 --modified-bases-models dna_r9.4.1_e8_sup@v3.3_5mCG@v0.1> JH89_JH90.bam
/home/drabiza/./dorado-0.4.3-linux-x64/bin/dorado basecaller /home/drabiza/dorado-0.4.3-linux-x64/bin/dna_r9.4.1_e8_sup@v3.3 JH91_JH93/output.pod5 --modified-bases-models dna_r9.4.1_e8_sup@v3.3_5mCG@v0.1> JH91_JH93.bam
/home/drabiza/./dorado-0.4.3-linux-x64/bin/dorado basecaller /home/drabiza/dorado-0.4.3-linux-x64/bin/dna_r9.4.1_e8_sup@v3.3 JH94_JH95/output.pod5 --modified-bases-models dna_r9.4.1_e8_sup@v3.3_5mCG@v0.1> JH94_JH95.bam
/home/drabiza/./dorado-0.4.3-linux-x64/bin/dorado basecaller /home/drabiza/dorado-0.4.3-linux-x64/bin/dna_r9.4.1_e8_sup@v3.3 JH96_JH97/output.pod5 --modified-bases-models dna_r9.4.1_e8_sup@v3.3_5mCG@v0.1> JH96_JH97.bam
/home/drabiza/./dorado-0.4.3-linux-x64/bin/dorado basecaller /home/drabiza/dorado-0.4.3-linux-x64/bin/dna_r9.4.1_e8_sup@v3.3 JH98_JH100/output.pod5 --modified-bases-models dna_r9.4.1_e8_sup@v3.3_5mCG@v0.1> JH98_JH100.bam

echo "Part 1 Done"

for sample in *bam
do
	id=$(basename "$sample" .bam)
	samtools view $sample -o $id.sam -@ 24
done

echo "Part 2 Done"

#mv to to NAS
mkdir /wilee/Data/drabiza1/Carolina/BAMS_SAMS_V2/

#Go through each sample and then extract reads
for sample in JH100 JH101 JH102 JH103 JH104 JH76 JH77 JH80 JH81 JH82 JH88 JH89 JH90 JH91 JH93 JH94 JH95 JH96 JH97 JH98
do
awk '{print $1}' Reads_table/*$sample* > $sample.txt
grep -F -f $sample.txt *$sample*sam > $Sample.new.sam
done

echo "Part 3 Done"

# Get sam to bam
for sample in *new.sam
do    
	id=$(basename "$sample" .new.sam)
    	samtools view -b -o "${id}.bam" "$sample"
    	#echo $sample
done

echo "Part 4 Done"

#Move all samples to SAM
mv *sam /wilee/Data/drabiza1/Carolina/BAMS_SAMS_V2/

#need the right samtools version for next line
conda activate py2 

## strip out (of the BAM) the reads including the modification tags [MM/ML] and gzip compress, outputting to fastq_pass.mod.fastq.gz
for sample in *bam 
do  
	samtools fastq -@24 -T MM,ML $sample | gzip > $sample.fastq.gz
       	#echo $sample
done

echo "Part 5 Done"
