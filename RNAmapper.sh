#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=1-0
#SBATCH --job-name=RNAmappertest
#SBATCH --output=slurm_out/slurm%j_%x.out
#SBATCH --error=slurm_out/slurm%j_%x.err
#SBATCH --mail-user=adayton@uoregon.edu

wt=../testVCFs/Hox40_wt_chr1.vcf
mut=../testVCFs/Hox40_mut_chr1.vcf

/usr/bin/time -v ./RNAmapper.py \
-wt $wt \
-mut $mut


