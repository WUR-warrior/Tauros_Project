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

ml legacy 
ml plink/2.0
ml bcftools/gcc/64/1.9

PLINK_BASE="/PLINK_filtered_no_ancestral"
PLINK_DIR="/lustre/nobackup/WUR/ABGC/shared/Tauros/PLINK/Plink_no_LD_trimming"
OUTDIR="/lustre/nobackup/WUR/ABGC/shared/Tauros/PLINK"
VCF_OUT="${OUTDIR}/No_LD_ancestor_trimming_for_beagle.vcf.gz"

echo "Converting PLINK binary files to VCF format for BEAGLE..."

plink2 --bfile ${PLINK_DIR}/${PLINK_BASE} --recode vcf-iid --allow-extra-chr --chr-set 29 no-xy --out ${OUTDIR}/No_LD_ancestor_trimming_for_beagle
bgzip -f ${OUTDIR}/No_LD_ancestor_trimming_for_beagle.vcf
tabix -f -p vcf ${OUTDIR}/No_LD_ancestor_trimming_for_beagle.vcf.gz

echo "VCF file created: ${OUTDIR}/No_LD_ancestor_trimming_for_beagle.vcf.gz"

echo "Checking VCF file..."
echo "Number of variants: $(bcftools view -H ${OUTDIR}/No_LD_ancestor_trimming_for_beagle.vcf.gz | wc -l)"
echo "Number of samples: $(bcftools query -l ${OUTDIR}/No_LD_ancestor_trimming_for_beagle.vcf.gz | wc -l)"

echo ""
echo "BEAGLE input file ready: ${OUTDIR}/No_LD_ancestor_trimming_for_beagle.vcf.gz"
echo ""
