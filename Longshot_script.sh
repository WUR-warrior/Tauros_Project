#!/bin/bash

#SBATCH --time=2-25:00:00

#SBATCH -N 1

#SBATCH -c 1

#SBATCH --mem=80000

#SBATCH --qos=Std

#SBATCH --output=slurm.output_longshot%j.txt

#SBATCH --error=slurm.error_longhsot%j.txt

#for the input  I need both the sorted BAM file and the .bai file in the same directory
longshot --bam /lustre/nobackup/WUR/ABGC/shared/Tauros/minimap2/Taurus.alignment_sorted.bam --ref /lustre/nobackup/WUR/ABGC/shared/Taurus/ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_genomic.fna --out /lustre/nobackup/WUR/ABGC/shared/Tauros/longshot/Taurus.alignment_variants.vcf