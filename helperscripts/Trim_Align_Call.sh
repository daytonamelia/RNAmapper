#!/bin/bash




# Cut adapters

# Mutant 
/usr/bin/time -v cutadapt -a $mutAdapter1 -A $mutAdapter2 \
-o cut_Mutant_R1.fastq.gz -p cut_Mutant_R2.fastq.gz \
$mutantRead1 \
$mutantRead2

# WT 
/usr/bin/time -v cutadapt -a $wtAdapter1 -A $wtAdapter2 \
-o cut_WT_R1.fastq.gz -p cut_WT_R2.fastq.gz \
$wtRead1 \
$wtRead2




# Trim Mutant
/usr/bin/time -v trimmomatic PE cut_Mutant_R1.fastq.gz \
cut_Mutant_R2.fastq.gz \
For_Paired_Mutant_R1.fastq.gz \
For_UnPaired_Mutant_R1.fastq.gz \
Rev_Paired_Mutant_R2.fastq.gz \
Rev_UnPaired_Mutant_R2.fastq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35



# Trim WT
/usr/bin/time -v trimmomatic PE cut_WT_R1.fastq.gz \
cut_WT_R2.fastq.gz \
For_Paired_WT_R1.fastq.gz \
For_UnPaired_WT_R1.fastq.gz \
Rev_Paired_WT_R2.fastq.gz \
Rev_UnPaired_WT_R2.fastq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35




# Run FASTQC on Trimmed data
/usr/bin/time -v fastqc \
For_Paired_Mutant_R1.fastq.gz \
Rev_Paired_Mutant_R2.fastq.gz \
For_Paired_WT_R1.fastq.gz \
Rev_Paired_WT_R2.fastq.gz \
-o . -t 4




# Get splicing info for HiSat2
/usr/bin/time -v hisat2_extract_splice_sites.py $gtfFile > ./splice_sites.ss



# Create HiSat2 index build for Genome
/usr/bin/time -v hisat2-build $refGenome Genome_index




# HiSat2 Alignment 

# Mutant alignment
/usr/bin/time -v hisat2 --phred33 --dta-cufflinks -x Genome_index --known-splicesite-infile splice_sites.ss \
-1 For_Paired_Mutant_R1.fastq.gz \
-2 Rev_Paired_Mutant_R2.fastq.gz \
-S mut.sam

# WT alignment
/usr/bin/time -v hisat2 --phred33 --dta-cufflinks -x Genome_index --known-splicesite-infile splice_sites.ss \
-1 For_Paired_WT_R1.fastq.gz \
-2 Rev_Paired_WT_R2.fastq.gz \
-S wt.sam




# SAMTOOLS

#Convert into bam file
samtools view -bS mut.sam > mut.bam
samtools view -bS wt.sam > wt.bam



# Keep same output name to use igV
#Sort bam:
samtools sort mut.bam -o mut.bam
samtools sort wt.bam -o wt.bam



# Index bam file
samtools index mut.bam
samtools index wt.bam



# Call Variants with bcftools
bcftools mpileup -ABf $refGenome mut.bam | bcftools view - > mut.bcf
bcftools mpileup -ABf $refGenome wt.bam | bcftools view - > wt.bcf



# Split bcf file into separate chromosome vcf files
/usr/bin/time -v for i in {1..25}; do bcftools view mut.bcf | grep -w ^$i > mut_chr$i.vcf &;done
/usr/bin/time -v for i in {1..25}; do bcftools view wt.bcf | grep -w ^$i > wt_chr$i.vcf &;done

wait
echo "Done"
