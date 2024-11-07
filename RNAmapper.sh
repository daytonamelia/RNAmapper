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

wt=/projects/bgmp/shared/groups/2024/fw-genetics/envs/HiSat2_alignment/Split_Chrom_VCFs/Hox40_wt_chr1.vcf
mut=/projects/bgmp/shared/groups/2024/fw-genetics/envs/HiSat2_alignment/Split_Chrom_VCFs/Hox40_mut_chr1.vcf

/usr/bin/time -v ./RNAmapper.py \
-wt $wt \
-mut $mut


