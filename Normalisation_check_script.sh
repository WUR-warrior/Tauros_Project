#!/bin/bash

#SBATCH --time=2-25:00:00

#SBATCH -N 1

#SBATCH -c 1

#SBATCH --mem=4000

#SBATCH --qos=Std

#SBATCH --output=slurm.output_norm%j.txt

#SBATCH --error=slurm.error_norm%j.txt
#SBATCH --job-name=norm
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

dir=/lustre/nobackup/WUR/ABGC/shared/Tauros/VCF_merged

VCF_FILE="${dir}/merged_norm_dedup.vcf.gz"

FIELDS=("GQ" "SAP" "MQM" "AB")

for field in "${FIELDS[@]}"; do
    echo "Checking field: $field"
    
    if bcftools view -h "$VCF_FILE" | grep -q "ID=${field},"; then
        echo "  $field is present in VCF."
        
        bcftools query -f "%CHROM\t%POS[\t%${field}]\n" "$VCF_FILE" | \
        awk -v field=$field '
        {
            for(i=3; i<=NF; i++) {
                if($i != "." && $i != "0") {
                    sample[i]++;
                }
            }
        }
        END {
            if (length(sample) == 0) {
                print "  No non-zero values found for", field;
            } else {
                for (i in sample) {
                    printf("  Sample %d has %d non-zero %s entries\n", i-2, sample[i], field);
                }
            }
        }'
    else
        echo "  $field is NOT present in the VCF."
    fi

    echo ""
done
