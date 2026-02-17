#!/bin/bash
#SBATCH --comment=773320000
#SBATCH --time=1-25:00:00
#SBATCH --mem=40000
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm.output_merge%j.txt
#SBATCH --error=slurm.error_merge%j.txt
#SBATCH --job-name=merge
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

ml legacy
ml bcftools/gcc/64/1.9

dir=/lustre/nobackup/WUR/ABGC/shared/Tauros/VCF
outdir=/lustre/nobackup/WUR/ABGC/shared/Tauros/VCF_merged

bcftools sort ${dir}/Tauros_autosomes_only_snps.vcf.gz -O z -o ${dir}/Tauros_autosomes_only_snps_sorted.vcf.gz
bcftools sort ${dir}/Breed_AURG_autosomes_only_snps.vcf.gz -O z -o ${dir}/Breed_AURG_autosomes_only_snps_sorted.vcf.gz

bcftools index ${dir}/Tauros_autosomes_only_snps_sorted.vcf.gz
bcftools index ${dir}/Breed_AURG_autosomes_only_snps_sorted.vcf.gz

bcftools merge --force-samples ${dir}/Tauros_autosomes_only_snps_sorted.vcf.gz ${dir}/Breed_AURG_autosomes_only_snps_sorted.vcf.gz -O z -o ${outdir}/merged_autosomes.vcf.gz

bcftools index ${outdir}/merged_autosomes.vcf.gz

echo "Sample count:"
bcftools query -l ${outdir}/merged_autosomes.vcf.gz | wc -l

echo "Variant count:"
bcftools view -H ${outdir}/merged_autosomes.vcf.gz | wc -l

