#!/bin/bash
#SBATCH --comment=773320000
#SBATCH --time=1-25:00:00
#SBATCH --mem=40000
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm.output_zip%j.txt
#SBATCH --error=slurm.error_zip%j.txt
#SBATCH --job-name=zip
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

cd /lustre/nobackup/WUR/ABGC/shared/Tauros/VCF
ml legacy
ml bcftools/gcc/64/1.9

bgzip Breed_AURG_var_fixed.vcf
bcftools index Breed_AURG_var_fixed.vcf.gz
