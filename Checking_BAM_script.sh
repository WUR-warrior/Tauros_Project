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

cd /lustre/nobackup/WUR/ABGC/shared/Tauros
mkdir BAM_stats
cd /lustre/nobackup/WUR/ABGC/shared/Tauros/Sorted_BAM

ml legacy
ml samtools
ml bedtools

dir=/lustre/nobackup/WUR/ABGC/shared/Tauros/Sorted_BAM
outdir=/lustre/nobackup/WUR/ABGC/shared/Tauros/BAM_stats
sname=$1

#samtools flagstat ${dir}/${sname}_sorted.bam > ${outdir}/${sname}_summary.txt
samtools coverage -m ${dir}/${sname}_sorted.bam > ${outdir}/${sname}_coverage_histogram.txt
samtools depth ${dir}/${sname}_sorted.bam | cut -f3 | sort -n | uniq -c > ${outdir}/${sname}_depth_distribution.txt
samtools coverage ${dir}/${sname}_sorted.bam > ${outdir}/${sname}_coverage_stats.txt

