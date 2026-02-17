#!/bin/bash
#SBATCH --comment=773320000
#SBATCH --time=1-25:00:00
#SBATCH --mem=4000
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm.output_quast%j.txt
#SBATCH --error=slurm.error_quast%j.txt
#SBATCH --job-name=quast
#SBATCH --mail-type=ALL
#SBATCH --mail-user=https://matt.watts@wur.nl

ml legacy
ml groups
ml WUR/RIKILT/QUAST/5.2.0

quast /lustre/nobackup/WUR/ABGC/shared/Tauros/Assembly_Run_4/assembly/Tauros_Test/assembly.fasta -r /lustre/nobackup/WUR/ABGC/shared/Taurus/ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_genomic.fna -o /lustre/nobackup/WUR/ABGC/shared/Tauros/QUAST -t 8 --eukaryote --large

