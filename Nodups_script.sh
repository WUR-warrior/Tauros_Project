#!/bin/bash
#SBATCH --comment=773320000
#SBATCH --time=1-25:00:00
#SBATCH --mem=4000
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm.output_nodups%j.txt
#SBATCH --error=slurm.error_nodups%j.txt
#SBATCH --job-name=nodups
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

ml legacy
ml samtools

sname=$1

cd /lustre/nobackup/WUR/ABGC/shared/Tauros/MarkedBAM
samtools view -b -F 0x400 ${sname}_markedup.bam > ${sname}_nodup.bam
