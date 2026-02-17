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

dir=/lustre/nobackup/WUR/ABGC/shared/Tauros/VCF_merged
outdir=/lustre/nobackup/WUR/ABGC/shared/Tauros/PLINK

plink_summary="${outdir}/plink_analysis_summary8.txt"
echo "PLINK Analysis Pipeline Summary" > ${plink_summary}
echo "Date: $(date)" >> ${plink_summary}
echo "======================================" >> ${plink_summary}
echo "" >> ${plink_summary}

extract_plink_stats() {
    local log_file=$1
    local step_name=$2
    
    if [[ -f ${log_file} ]]; then
        # Extract variant and sample counts from PLINK log
        variants=$(grep "variants loaded from" ${log_file} | tail -1 | awk '{print $1}')
        samples=$(grep "samples (" ${log_file} | tail -1 | awk '{print $1}')
        
        if [[ -z "$variants" ]]; then
            variants=$(grep "variants remaining after" ${log_file} | tail -1 | awk '{print $1}')
        fi
        
        if [[ -z "$samples" ]]; then
            samples=$(grep "people remaining after" ${log_file} | tail -1 | awk '{print $1}')
        fi
        
        echo "${step_name}: ${variants:-'N/A'} variants, ${samples:-'N/A'} samples" >> ${plink_summary}
    else
        echo "${step_name}: Log file not found" >> ${plink_summary}
    fi
}

sleep 10

while [[ $(bcftools view -H ${dir}/final_filtered.vcf.gz | wc -l) -eq 0 ]]; do
    echo "Waiting for file to be accessible..." >> ${plink_summary}
    sleep 5
done

echo "DIAGNOSTIC - PLINK input verification:" >> ${plink_summary}
plink_input_count=$(bcftools view -H ${dir}/final_filtered.vcf.gz | wc -l)
echo "bcftools count of PLINK input: ${plink_input_count}" >> ${plink_summary}

echo "File timestamp: $(stat -c %y ${dir}/final_filtered.vcf.gz)" >> ${plink_summary}

echo "Chromosomes in VCF:" >> ${plink_summary}
bcftools view -H ${dir}/final_filtered.vcf.gz | cut -f1 | sort | uniq -c >> ${plink_summary}

echo "Step 1: Calculating sample missingness..." >> ${plink_summary}
plink2 --vcf ${dir}/final_filtered.vcf.gz --missing --allow-extra-chr --chr-set 29 no-xy --out ${outdir}/sample_missingness
extract_plink_stats "${outdir}/sample_missingness.log" "Sample missingness calculation"

echo "" >> ${plink_summary}
echo "Step 2: Creating filtered BED file..." >> ${plink_summary}

pre_plink_count=$(bcftools view -H ${dir}/final_filtered.vcf.gz | wc -l)
echo "Variants before PLINK processing: ${pre_plink_count}" >> ${plink_summary}

plink2 --vcf ${dir}/final_filtered.vcf.gz --allow-extra-chr --chr-set 29 no-xy --set-missing-var-ids @:# --make-bed --out ${outdir}/PLINK_filtered

extract_plink_stats "${outdir}/PLINK_filtered.log" "After PLINK conversion"

post_plink_count=$(wc -l < ${outdir}/PLINK_filtered.bim)
echo "Variants in BIM file: ${post_plink_count}" >> ${plink_summary}

echo "PLINK conversion log (last 20 lines):" >> ${plink_summary}
tail -20 ${outdir}/PLINK_filtered.log >> ${plink_summary}
echo "" >> ${plink_summary}

echo "PLINK warnings/errors:" >> ${plink_summary}
grep -i "warning\|error\|skipped\|duplicate" ${outdir}/PLINK_filtered.log >> ${plink_summary}

echo "" >> ${plink_summary}
echo "Step 3: LD pruning..." >> ${plink_summary}
plink2 --bfile ${outdir}/PLINK_filtered --allow-extra-chr --chr-set 29 no-xy --indep-pairwise 50 10 0.2 --out ${outdir}/LD_pruning
extract_plink_stats "${outdir}/LD_pruning.log" "LD pruning analysis"

if [[ -f ${outdir}/LD_pruning.prune.in ]]; then
    pruned_in=$(wc -l < ${outdir}/LD_pruning.prune.in)
    echo "Variants retained after LD pruning: ${pruned_in}" >> ${plink_summary}
fi

if [[ -f ${outdir}/LD_pruning.prune.out ]]; then
    pruned_out=$(wc -l < ${outdir}/LD_pruning.prune.out)
    echo "Variants removed by LD pruning: ${pruned_out}" >> ${plink_summary}
fi

echo "" >> ${plink_summary}
echo "Step 4: Creating final pruned dataset..." >> ${plink_summary}

sort -u ${outdir}/LD_pruning.prune.in > ${outdir}/LD_pruning.prune.in.clean

echo "Original prune.in file lines: $(wc -l < ${outdir}/LD_pruning.prune.in)" >> ${plink_summary}
echo "Cleaned prune.in file lines: $(wc -l < ${outdir}/LD_pruning.prune.in.clean)" >> ${plink_summary}

if [[ $(wc -l < ${outdir}/LD_pruning.prune.in.clean) -lt 1000 ]]; then
    echo "Cleaning seems to have failed, using original prune.in file" >> ${plink_summary}
    cp ${outdir}/LD_pruning.prune.in ${outdir}/LD_pruning.prune.in.clean
fi

variants_to_extract=$(wc -l < ${outdir}/LD_pruning.prune.in.clean)
echo "Variants to extract: ${variants_to_extract}" >> ${plink_summary}

plink2 --bfile ${outdir}/PLINK_filtered --allow-extra-chr --chr-set 29 no-xy --extract ${outdir}/LD_pruning.prune.in.clean --make-bed --out ${outdir}/PLINK_pruned_final_relaxed

if [[ -f ${outdir}/PLINK_pruned_final.bim ]]; then
    final_variants=$(wc -l < ${outdir}/PLINK_pruned_final_relaxed.bim)
    final_samples=$(wc -l < ${outdir}/PLINK_pruned_final_relaxed.fam)
    echo "Final pruned dataset: ${final_variants} variants, ${final_samples} samples" >> ${plink_summary}
else
    echo "ERROR: Final pruned dataset creation failed" >> ${plink_summary}
fi

echo "" >> ${plink_summary}
echo "Step 5: Principal Component Analysis..." >> ${plink_summary}
plink2 --bfile ${outdir}/PLINK_pruned_final_relaxed --allow-extra-chr --chr-set 29 no-xy --pca 4 --out ${outdir}/PLINK_PCA_relaxed
echo "PCA completed: 4 principal components calculated" >> ${plink_summary}

echo "" >> ${plink_summary}
echo "Step 6: Creating similarity matrix..." >> ${plink_summary}
plink2 --bfile ${outdir}/PLINK_pruned_final_relaxed --make-rel square --allow-extra-chr --chr-set 29 no-xy --out ${outdir}/similarity_matrix_relaxed
echo "Similarity matrix created" >> ${plink_summary}

echo "" >> ${plink_summary}
echo "Final Analysis Summary:" >> ${plink_summary}
echo "======================================" >> ${plink_summary}

if [[ -f ${outdir}/PLINK_pruned_final_relaxed.fam ]]; then
    final_samples=$(wc -l < ${outdir}/PLINK_pruned_final_relaxed.fam)
    echo "Final sample count: ${final_samples}" >> ${plink_summary}
fi

if [[ -f ${outdir}/PLINK_pruned_final_relaxed.bim ]]; then
    final_variants=$(wc -l < ${outdir}/PLINK_pruned_final_relaxed.bim)
    echo "Final variant count: ${final_variants}" >> ${plink_summary}
fi

echo "" >> ${plink_summary}
echo "Output files generated:" >> ${plink_summary}
echo "- Sample missingness: ${outdir}/sample_missingness.smiss" >> ${plink_summary}
echo "- Final dataset: ${outdir}/PLINK_pruned_final_relaxed.{bed,bim,fam}" >> ${plink_summary}
echo "- PCA results: ${outdir}/PLINK_PCA_relaxed.eigenvec, ${outdir}/PLINK_PCA_relaxed.eigenval" >> ${plink_summary}
echo "- Similarity matrix: ${outdir}/similarity_matrix_relaxed.rel" >> ${plink_summary}

echo ""
echo "PLINK analysis complete. Summary saved to: ${plink_summary}"
echo "Filtering summary available at: ${summary_file}"