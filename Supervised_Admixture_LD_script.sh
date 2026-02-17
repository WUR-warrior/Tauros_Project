#!/bin/bash
#SBATCH --time=2-25:00:00

#SBATCH -N 1

#SBATCH -c 1

#SBATCH --mem=80000

#SBATCH --qos=Std

#SBATCH --output=slurm.output_PLINK%j.txt

#SBATCH --error=slurm.error_PLINK%j.txt
#SBATCH --job-name=PLink
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

cd /lustre/nobackup/WUR/ABGC/shared/Tauros/PLINK/Plink_with_LD_trimming
echo "AURG AURG" > remove_ancestral.txt
ml legacy
ml plink/2.0
plink2 --bfile PLINK_pruned_final_relaxed --remove remove_ancestral.txt --make-bed --out PLINK_pruned_final_relaxed_no_ancestral --allow-extra-chr --geno 0.99 --chr-set 29 no-xy
cat > PLINK_pruned_final_relaxed_no_ancestral.pop << 'EOF'
-
SA
PO
PA
PA
MN
SA
MA
MA
MA
LM
EOF
admixture --supervised PLINK_pruned_final_relaxed_no_ancestral.bed 6