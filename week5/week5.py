#!/bin/bash

# bwa index c_elegans.PRJNA13758.WS283.genomic.fa

# for sample in 09 11 23 24 27 31 35 39 62 63
# do
# 	echo "Aligning sample:" ${sample}
# 	bwa mem -R "@RG\tID:${sample}\tSM:${sample}" \
# 	  sacCer3.fa \
# 	  A01_${sample}.fastq > A01_${sample}.sam
# done

# for sample in 09 11 23 24 27 31 35 39 62 63
# do 
# 	samtools view -bSh A01_${sample}.sam | samtools sort -o A01_${sample}.srt.bam

# done

# for sample in 09 11 23 24 27 31 35 39 62 63
# do 
# 	samtools index A01_${sample}.srt.bam

# done

for sample in 09 11 23 24 27 31 35 39 62 63
do 
	freebayes -f sacCer3.fa -p 2 --genotype-qualities A01_${sample}.srt.bam > A01_${sample}.vcf
done