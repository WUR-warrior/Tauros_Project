#!/bin/bash
#SBATCH --comment=773320000
#SBATCH --time=1-25:00:00
#SBATCH --mem=4000
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm.output_freebayes%j.txt
#SBATCH --error=slurm.error_freebayes%j.txt
#SBATCH --job-name=nanopore_assembly
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

ml legacy
ml samtools

dir=/lustre/nobackup/WUR/ABGC/shared/Tauros/BAM
outdir=/lustre/nobackup/WUR/ABGC/shared/Tauros/MarkedBAM
sname=$1

samtools sort -n -o ${outdir}/${sname}_name_sorted.bam ${dir}/${sname}.bam
samtools fixmate -m ${outdir}/${sname}_name_sorted.bam ${outdir}/${sname}_fixmate.bam
samtools sort -o ${outdir}/${sname}_fixmate_sorted.bam ${outdir}/${sname}_fixmate.bam
samtools markdup ${outdir}/${sname}_fixmate_sorted.bam ${outdir}/${sname}_markedup.bam

