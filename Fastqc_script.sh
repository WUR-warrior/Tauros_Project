#!/bin/bash
#SBATCH --comment=773320000
#SBATCH --time=1-25:00:00
#SBATCH --mem=40000
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm.output_fastqc%j.txt
#SBATCH --error=slurm.error_fastqc%j.txt
#SBATCH --job-name=fastqc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

ml legacy
ml WUR/RIKILT/FastQC-0.11.9

dir=/lustre/nobackup/WUR/ABGC/shared/Tauros/MarkedBAM
outdir=/lustre/nobackup/WUR/ABGC/shared/Tauros/fastqc
sname=$1

fastqc ${dir}/${sname}_markedup.bam -o ${outdir}
