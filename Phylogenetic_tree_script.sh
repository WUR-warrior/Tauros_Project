#!/bin/bash
#SBATCH --comment=773320001
#SBATCH --time=2-12:00:00
#SBATCH --mem=20000
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm.output_phylo%j.txt
#SBATCH --error=slurm.error_phylo%j.txt
#SBATCH --job-name=PhylogenetictreewithLD
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

ml legacy
ml plink/2.0
ml RAxML/gcc/64/8.2.9

cd /lustre/nobackup/WUR/ABGC/shared/Tauros/PLINK/Plink_with_LD_trimming
plink2 --bfile PLINK_pruned_final_relaxed --allow-extra-chr --chr-set 29 no-xy --recode vcf-iid --out PLINK_pruned_final_relaxed

wget https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py

python vcf2phylip.py -i PLINK_pruned_final_relaxed.vcf -o phy

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s PLINK_pruned_final_relaxed.min4.phy -n SNP_tree --asc-corr=lewis

plink2 --bfile PLINK_pruned_final_relaxed_no_ancestral --allow-extra-chr --chr-set 29 no-xy --recode vcf-iid --out PLINK_pruned_final_relaxed_no_ancestral

python vcf2phylip.py -i PLINK_pruned_final_relaxed_no_ancestral.vcf -o phy

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s PLINK_pruned_final_relaxed_no_ancestral.min4.phy -n SNP_tree_no_ancestor --asc-corr=lewis
