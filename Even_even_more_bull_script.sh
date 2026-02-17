#!/bin/bash
#SBATCH --job-name=compare_haplotypes
#SBATCH --time=8:00:00
#SBATCH --mem=400000
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm.output_haplotype_%j.txt
#SBATCH --error=slurm.error_haplotype_%j.txt
#SBATCH --mail-user=matt.watts@wur.nl
#SBATCH --mail-type=FAIL,END

ml legacy
module load python/3.9.4  

cd /lustre/nobackup/WUR/ABGC/shared/Tauros/Beagle_phasing_no_LD/final_output

python compare_haplotypes_with_background.py --vcf cattle_phased_no_LD_filter.vcf.gz --admixed SAMPLE --references PO1 PA1 PA2 MA1 MA2 MA3 SA1 SA2 LM1 MN1 --window_size 500 --background_samples 1000 --z_threshold 0.5 --out results/haplotype_matches.csv