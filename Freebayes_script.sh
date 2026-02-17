#!/bin/bash
#SBATCH --comment=773320000
#SBATCH --time=1-25:00:00
#SBATCH --mem=40000
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm.output_freebayes%j.txt
#SBATCH --error=slurm.error_freebayes%j.txt
#SBATCH --job-name=nanopore_assembly
#SBATCH --mail-type=ALL
#SBATCH --mail-user=https://matt.watts@wur.nl

ml legacy
ml freebayes/1.3.6

dir=/lustre/nobackup/WUR/ABGC/shared/Tauros/Sorted_BAM
outdir=/lustre/nobackup/WUR/ABGC/shared/Tauros/VCF


freebayes -f /lustre/nobackup/WUR/ABGC/shared/Taurus/ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_genomic.fna ${dir}/AURG_sorted_fixed.bam ${dir}/LM1_sorted_fixed.bam ${dir}/MA1_sorted_fixed.bam ${dir}/MA2_sorted_fixed.bam ${dir}/MA3_sorted_fixed.bam ${dir}/MN1_sorted_fixed.bam ${dir}/PA1_sorted_fixed.bam ${dir}/PA2_sorted_fixed.bam ${dir}/PO1_sorted_fixed.bam ${dir}/SA1_sorted_fixed.bam ${dir}/SA2_sorted_fixed.bam > ${outdir}/Breed_AURG_var_fixed.vcf

