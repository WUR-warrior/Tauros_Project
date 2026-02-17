#!/bin/bash

#SBATCH --time=2-25:00:00

#SBATCH -N 1

#SBATCH -c 1

#SBATCH --mem=80000

#SBATCH --qos=Std

#SBATCH --mail-user=matt.watts@wur.nl
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --output=slurm.output_BWA%j.txt

#SBATCH --error=slurm.error_BWA%j.txt

ml legacy
ml bwa 

dir=/lustre/nobackup/WUR/ABGC/shared/Tauros/Breed_Trimming
outdir=/lustre/nobackup/WUR/ABGC/shared/Tauros/Mapping
sname=$1

bwa mem -t 1 /lustre/nobackup/WUR/ABGC/shared/Tauros/Index/GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz ${dir}/${sname}_1_run_1.fastq.gz ${dir}/${sname}_2_run_1.fastq.gz > ${outdir}/${sname}.sam


