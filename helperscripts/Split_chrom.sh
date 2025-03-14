#!/bin/bash


# Call Variants with bcftools
bcftools mpileup -ABf $refGenome mut.bam | bcftools view - > mut.bcf
bcftools mpileup -ABf $refGenome wt.bam | bcftools view - > wt.bcf



# Split bcf file into separate chromosome vcf files
/usr/bin/time -v for i in {1..25}; do bcftools view mut.bcf | grep -w ^$i > mut_chr$i.vcf &;done
/usr/bin/time -v for i in {1..25}; do bcftools view wt.bcf | grep -w ^$i > wt_chr$i.vcf &;done

