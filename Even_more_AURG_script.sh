#!/bin/bash
#SBATCH --job-name=haplo_AURG_SAMPLE
#SBATCH --time=4:00:00
#SBATCH --mem=64000
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm.output_AURG_SAMPLE_%j.txt
#SBATCH --error=slurm.error_AURG_SAMPLE_%j.txt
#SBATCH --mail-user=matt.watts@wur.nl
#SBATCH --mail-type=FAIL,END

ml legacy
module load python/3.9.4  

cd /lustre/nobackup/WUR/ABGC/shared/Tauros/Beagle_two_sample_phasing_No_LD/final_output
mkdir results

python compare_haplotypes_with_background.py --vcf two_sample_phased_AURG_SAMPLE.vcf.gz --admixed SAMPLE --references AURG --window_size 500 --background_samples 1000 --z_threshold 0.5 --out results/SAMPLE_vs_AURG_zscores.csv
