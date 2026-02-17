#!/bin/bash

#SBATCH --time=3-25:00:00

#SBATCH -N 1

#SBATCH -c 1

#SBATCH --mem=60000

#SBATCH --qos=Std

#SBATCH --output=slurm.output_BAM%j.txt

#SBATCH --error=slurm.error_BAM%j.txt
#SBATCH --job-name=BAM_script
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

ml legacy
ml samtools
ml bwa

dir=/lustre/nobackup/WUR/ABGC/shared/Tauros/Aurochs_trimming
outdir=/lustre/nobackup/WUR/ABGC/shared/Tauros/Aurochs_BAM

bwa mem -t 4 /lustre/nobackup/WUR/ABGC/shared/Tauros/Index/GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz ${dir}/AURG.fastq.gz | samtools sort -o ${outdir}/AURG_sorted.bam 
