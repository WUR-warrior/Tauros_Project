#!/bin/bash

#SBATCH --time=2-25:00:00

#SBATCH -N 1

#SBATCH -c 1

#SBATCH --mem=80000

#SBATCH --qos=Std

#SBATCH --output=slurm.output_rename%j.txt

#SBATCH --error=slurm.error_rename%j.txt
#SBATCH --job-name=rename
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

ml legacy
ml samtools
ml bcftools/gcc/64/1.9

bcftools reheader -s <(echo "SAMPLE Tauros") Taurosvar.vcf.gz > Taurosvar_renamed.vcf.gz

bcftools index Taurosvar_renamed.vcf.gz

for bam in AURG_sorted.bam LM1_sorted.bam MA1_sorted.bam MA2_sorted.bam MA3_sorted.bam MN1_sorted.bam PA1_sorted.bam PA2_sorted.bam PO1_sorted.bam SA1_sorted.bam SA2_sorted.bam; do
  {
    sample=$(basename "$bam" _sorted.bam)
    echo "Adding RG to $bam with sample name: $sample"
    samtools addreplacerg -r "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:${sample}_lib" \
      -o "${sample}_sorted_fixed.bam" "$bam"
    samtools index "${sample}_sorted_fixed.bam"
    echo "Finished processing $sample"
  } &
done

wait
echo "All samples processed!"