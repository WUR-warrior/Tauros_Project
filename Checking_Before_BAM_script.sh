#!/bin/bash
#SBATCH --comment=773320000
#SBATCH --time=1-25:00:00
#SBATCH --mem=4000
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm.output_checks%j.txt
#SBATCH --error=slurm.error_checks%j.txt
#SBATCH --job-name=checks
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

ml legacy
ml samtools

dir=/lustre/nobackup/WUR/ABGC/shared/Tauros/MarkedBAM
outdir=/lustre/nobackup/WUR/ABGC/shared/Tauros/BAM_stats
sname=$1

samtools flagstat ${dir}/${sname}_markedup.bam > ${outdir}/${sname}_before_summary.txt
