#!/bin/bash

#SBATCH --time=2-25:00:00

#SBATCH -N 1

#SBATCH -c 1

#SBATCH --mem=60000

#SBATCH --qos=Std

#SBATCH --output=slurm.output_norm%j.txt

#SBATCH --error=slurm.error_norm%j.txt
#SBATCH --job-name=norm
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

ml legacy 
ml bcftools/gcc/64/1.9

cd /lustre/nobackup/WUR/ABGC/shared/Tauros/VCF_merged

bcftools annotate -x INFO/AC,INFO/AN,INFO/AF merged_autosomes.vcf.gz -Oz -o merged_autosomes_clean.vcf.gz

bcftools index merged_autosomes_clean.vcf.gz

bcftools norm -m -both -f /lustre/nobackup/WUR/ABGC/shared/Taurus/ARS-UCD1.2/renamed_ARS-UCD1.2_autosomes.fna -Oz -o merged_norm.vcf.gz merged_autosomes_clean.vcf.gz

bcftools norm -d both -Oz -o merged_norm_dedup.vcf.gz merged_norm.vcf.gz

bcftools index -f merged_norm_dedup.vcf.gz