#!/bin/bash

#SBATCH --time=2-25:00:00

#SBATCH -N 1

#SBATCH -c 1

#SBATCH --mem=60000

#SBATCH --qos=Std

#SBATCH --output=slurm.output_Sort%j.txt

#SBATCH --error=slurm.error_sort%j.txt
#SBATCH --job-name=sorting
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

cd /home/WUR/watts001/trimmomatic/0.39
ml legacy
ml trimmomatic 

dir=/lustre/nobackup/WUR/ABGC/shared/Taurus/Aurochs_Sequences/Concatenated_Sequences
outdir=/lustre/nobackup/WUR/ABGC/shared/Tauros/Aurochs_trimming

java -jar trimmomatic-0.39.jar SE ${dir}/AURH++.fastq.gz ${outdir}/AURH+_trimmed.fastq.gz -threads 1 ILLUMINACLIP:/home/WUR/watts001/trimmomatic/0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:20:20 MINLEN:30  
