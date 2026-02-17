#!/bin/bash
#SBATCH --comment=773320000
#SBATCH --time=2-25:00:00
#SBATCH --mem=40000
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm.output_filter_merge%j.txt
#SBATCH --error=slurm.error_filter_merge%j.txt
#SBATCH --job-name=filter_merge
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

ml legacy
ml bcftools/gcc/64/1.9
sname=$1
#PLINKID
cd /lustre/nobackup/WUR/ABGC/shared/Tauros/VCF

bcftools annotate --rename-chrs chr_name_map.txt -o ${sname}_filtered.vcf.gz -O z ${sname}var.vcf.gz

bcftools index ${sname}_filtered.vcf.gz

bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29 -o ${sname}_autosomes_only.vcf.gz -O z ${sname}_filtered.vcf.gz