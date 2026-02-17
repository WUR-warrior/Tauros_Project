#!/bin/bash
#SBATCH --comment=773320000
#SBATCH --time=1-25:00:00
#SBATCH --mem=40000
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm.output_sorted%j.txt
#SBATCH --error=slurm.error_sorted%j.txt
#SBATCH --job-name=sorted
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

cd /lustre/nobackup/WUR/ABGC/shared/Tauros/Sorted_BAM

ml legacy
ml samtools

sname=$1

samtools index ${sname}_sorted.bam ${sname}_index.bai
