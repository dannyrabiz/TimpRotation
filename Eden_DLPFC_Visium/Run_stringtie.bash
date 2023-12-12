"""
Run stringtie2 for both sets of clusters (k=2 and k=9) on all samples
"""

for folder in *2clusters
do
	mkdir $folder/stringtie
	for sample in {1..2}
	do 
		stringtie $folder/Cluster$sample.bam -e -G /home/drabiza/gencode.v24.annotation.gtf -L -g 150 -p 24 -m 300 -o $folder/stringtie/Cluster$sample.gtf -A $folder/stringtie/Cluster$sample.GeneAbundance.txt
		echo "Done:$folder, $sample" 
	done
done

for folder in *9clusters
do
        mkdir $folder/stringtie
        for sample in {1..9}
        do
                stringtie $folder/Cluster$sample.bam -e -G /home/drabiza/gencode.v24.annotation.gtf -L -g 150 -p 24 -m 300 -o $folder/stringtie/Cluster$sample.gtf -A $folder/stringtie/Cluster$sample.GeneAbundance.txt
                echo "Done:$folder, $sample"
        done
done
