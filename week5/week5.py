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

# freebayes -f sacCer3.fa --genotype-qualities -p 1 -b A01_09.srt.bam -b A01_11.srt.bam -b A01_23.srt.bam -b A01_24.srt.bam -b A01_27.srt.bam -b A01_31.srt.bam -b A01_35.srt.bam -b A01_39.srt.bam -b A01_62.srt.bam -b A01_63.srt.bam -v A01_sample.vcf 

# vcffilter -f "QUAL > 0.99" A01_sample.vcf > A01_sample_filtered.vcf

# vcfallelicprimitives -k -g A01_sample_filtered.vcf > A01_sample_decomposed.vcf

#snpEff download R64-1-1.86 
snpEff ann R64-1-1.86 A01_sample_decomposed.vcf > annotated.vcf


