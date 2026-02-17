#!/bin/bash

#SBATCH --time=2-25:00:00

#SBATCH -N 1

#SBATCH -c 1

#SBATCH --mem=60000

#SBATCH --qos=Std

#SBATCH --output=slurm.output_filter%j.txt

#SBATCH --error=slurm.error_filter%j.txt
#SBATCH --job-name=filter
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

dir=/lustre/nobackup/WUR/ABGC/shared/Tauros/VCF_merged
ml legacy
ml bcftools/gcc/64/1.9

summary_file="${dir}/filtering_summary8.txt"
echo "Cattle SNP Filtering Pipeline Summary" > ${summary_file}
echo "Date: $(date)" >> ${summary_file}
echo "======================================" >> ${summary_file}
echo "" >> ${summary_file}

count_vcf_stats() {
    local vcf_file=$1
    local step_name=$2
    
    if [[ -f ${vcf_file} ]]; then
        variants=$(bcftools view -H ${vcf_file} | wc -l)
        samples=$(bcftools query -l ${vcf_file} | wc -l)
        echo "${step_name}: ${variants} variants, ${samples} samples" >> ${summary_file}
    else
        echo "${step_name}: File not found" >> ${summary_file}
    fi
}

count_vcf_stats "${dir}/merged_norm_dedup.vcf.gz" "Initial dataset"

bcftools view --include 'QUAL >= 10' -Oz -o ${dir}/qual_filtered.vcf.gz ${dir}/merged_norm_dedup.vcf.gz
count_vcf_stats "${dir}/qual_filtered.vcf.gz" "After quality filtering (QUAL >= 10)"

bcftools view --include 'FORMAT/DP >= 4' -Oz -o ${dir}/depth_filtered.vcf.gz ${dir}/qual_filtered.vcf.gz
count_vcf_stats "${dir}/depth_filtered.vcf.gz" "After per-sample depth filtering (FORMAT/DP >= 4)"

bcftools view --include 'SAP < 80 && MQM > 25' -Oz -o ${dir}/final_filtered.vcf.gz ${dir}/depth_filtered.vcf.gz
count_vcf_stats "${dir}/final_filtered.vcf.gz" "After allele filtering (SAP < 80, MQM > 25)"

rm -f ${dir}/final_filtered.vcf.gz.csi
rm -f ${dir}/final_filtered.vcf.gz.tbi

bcftools index ${dir}/final_filtered.vcf.gz
count_vcf_stats "${dir}/final_filtered.vcf.gz" "Final filtered dataset (MAF >= 0.00)"

echo "DIAGNOSTIC - Direct count after creation:" >> ${summary_file}
direct_count=$(bcftools view -H ${dir}/final_filtered.vcf.gz | wc -l)
echo "Direct bcftools count: ${direct_count}" >> ${summary_file}

initial_variants=$(bcftools view -H ${dir}/merged_norm_dedup.vcf.gz | wc -l)
final_variants=$(bcftools view -H ${dir}/final_filtered.vcf.gz | wc -l)
retention_rate=$(echo "scale=2; ${final_variants} * 100 / ${initial_variants}" | bc)

echo "" >> ${summary_file}
echo "Filtering Summary:" >> ${summary_file}
echo "Retention rate: ${retention_rate}% (${final_variants}/${initial_variants} variants retained)" >> ${summary_file}
echo "" >> ${summary_file}
rm ${dir}/snps_only.vcf.gz
rm ${dir}/qual_filtered.vcf.gz
rm ${dir}/depth_filtered.vcf.gz

echo "Filtering complete. Summary saved to: ${summary_file}"

sync
sleep 5

echo "Final verification before PLINK processing:" >> ${summary_file}
final_check=$(bcftools view -H ${dir}/final_filtered.vcf.gz | wc -l)
echo "Final check count: ${final_check}" >> ${summary_file}

sync ${dir}/final_filtered.vcf.gz*
sleep 2

echo "File ready for PLINK processing" >> ${summary_file}