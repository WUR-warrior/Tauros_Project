#!/bin/bash
#SBATCH --time=2-25:00:00
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=40000
#SBATCH --qos=Std
#SBATCH --output=slurm.output_Annotate%j.txt
#SBATCH --error=slurm.error_Annotate%j.txt
#SBATCH --job-name=annotation
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

ml legacy
ml bcftools/gcc/64/1.9

sname=$1
#IDfile=PLINK_ID
cd /lustre/nobackup/WUR/ABGC/shared/Tauros/VCF

bcftools index ${sname}_autosomes_only.vcf.gz

bcftools view -v snps -O z -o "${sample}_autosomes_only_snps.vcf.gz" "${sample}_autosomes_only.vcf.gz
"

bcftools index -f -t "${sample}_autosomes_only_snps.vcf.gz"