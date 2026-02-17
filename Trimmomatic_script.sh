#!/bin/bash

#SBATCH --time=2-25:00:00

#SBATCH -N 1

#SBATCH -c 1

#SBATCH --mem=60000

#SBATCH --qos=Std

#SBATCH --output=slurm.output_trimm%j.txt

#SBATCH --error=slurm.error_trimm%j.txt

cd /home/WUR/watts001/trimmomatic/0.39

ml legacy
ml trimmomatic 

dir=/lustre/nobackup/WUR/ABGC/shared/Taurus/Founder_Breed_Sequences
outdir=/lustre/nobackup/WUR/ABGC/shared/Tauros/Breed_trimming
sname=$1

#with loop
java -jar trimmomatic-0.39.jar PE ${dir}/${sname}_1.fastq.gz ${dir}/${sname}_2.fastq.gz ${outdir}/${sname}_1_run_1.fastq.gz ${outdir}/${sname}_1_run_1_unpaired.fastq.gz ${outdir}/${sname}_2_run_1.fastq.gz ${outdir}/${sname}_2_run_1_unpaired.fastq.gz  -threads 1 ILLUMINACLIP:/home/WUR/watts001/trimmomatic/0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:20:20 MINLEN:30

