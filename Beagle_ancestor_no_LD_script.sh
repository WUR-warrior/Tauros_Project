#!/bin/bash
#SBATCH --comment=773320001
#SBATCH --time=2-12:00:00
#SBATCH --mem=64000
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm.output_beagle%j.txt
#SBATCH --error=slurm.error_beagle%j.txt
#SBATCH --job-name=beagle_phase_two_samples
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

echo "=== BEAGLE PHASING PIPELINE FOR TWO SAMPLES ==="
echo "Starting phasing analysis: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"

ml legacy
ml bcftools/1.9
ml plink/2.0

INPUT_DIR="/lustre/nobackup/WUR/ABGC/shared/Tauros/PLINK"
OUTPUT_DIR="/lustre/nobackup/WUR/ABGC/shared/Tauros/Beagle_two_sample_phasing_No_LD/"

SAMPLE1="AURG"
SAMPLE2="SAMPLE"

echo "Target samples: $SAMPLE1 and $SAMPLE2"

mkdir -p "$OUTPUT_DIR"

INPUT_VCF="$INPUT_DIR/No_LD_trimming_for_beagle.vcf.gz"
OUTPUT_PREFIX="two_sample_phased_${SAMPLE1}_${SAMPLE2}"

echo "=== STEP 1: PREPARE INPUT VCF ==="
echo "Checking and preparing input VCF..."

if [ ! -f "$INPUT_VCF" ]; then
    echo "ERROR: Input VCF file not found: $INPUT_VCF"
    echo "Please check the input file path"
    exit 1
fi

if [ ! -f "$INPUT_VCF.tbi" ] && [ ! -f "$INPUT_VCF.csi" ]; then
    echo "Indexing input VCF..."
    bcftools index -t "$INPUT_VCF"
fi

echo "Extracting sample information..."
bcftools query -l "$INPUT_VCF" > "$OUTPUT_DIR/all_samples.txt"
TOTAL_SAMPLES=$(wc -l < "$OUTPUT_DIR/all_samples.txt")
echo "Total samples in VCF: $TOTAL_SAMPLES"

echo "Checking for target samples..."
if ! grep -q "^${SAMPLE1}$" "$OUTPUT_DIR/all_samples.txt"; then
    echo "ERROR: Sample '$SAMPLE1' not found in VCF"
    echo "Available samples (first 10):"
    head -10 "$OUTPUT_DIR/all_samples.txt"
    exit 1
fi

if ! grep -q "^${SAMPLE2}$" "$OUTPUT_DIR/all_samples.txt"; then
    echo "ERROR: Sample '$SAMPLE2' not found in VCF"
    echo "Available samples (first 10):"
    head -10 "$OUTPUT_DIR/all_samples.txt"
    exit 1
fi

echo "? Found both target samples: $SAMPLE1 and $SAMPLE2"

echo "$SAMPLE1" > "$OUTPUT_DIR/target_samples.txt"
echo "$SAMPLE2" >> "$OUTPUT_DIR/target_samples.txt"

echo "=== STEP 2: EXTRACT TARGET SAMPLES ==="
echo "Extracting only the two target samples..."

bcftools view \
    -S "$OUTPUT_DIR/target_samples.txt" \
    -Oz \
    -o "$OUTPUT_DIR/two_samples.vcf.gz" \
    "$INPUT_VCF"

bcftools index -t "$OUTPUT_DIR/two_samples.vcf.gz"

echo "Extracted VCF statistics:"
bcftools stats "$OUTPUT_DIR/two_samples.vcf.gz" | grep -E "^SN" | head -10

echo "Samples in extracted VCF:"
bcftools query -l "$OUTPUT_DIR/two_samples.vcf.gz"

echo "=== STEP 3: SPLIT BY CHROMOSOME ==="
echo "Splitting VCF by chromosome for parallel processing..."

mkdir -p "$OUTPUT_DIR/chromosome_splits"
mkdir -p "$OUTPUT_DIR/phased_chromosomes"
mkdir -p "$OUTPUT_DIR/logs"
mkdir -p "$OUTPUT_DIR/final_output"

bcftools query -f '%CHROM\n' "$OUTPUT_DIR/two_samples.vcf.gz" | sort -u > "$OUTPUT_DIR/chromosome_list.txt"

echo "Found chromosomes:"
cat "$OUTPUT_DIR/chromosome_list.txt"

while read -r CHROM; do
    echo "Extracting chromosome $CHROM..."
    bcftools view \
        -r "$CHROM" \
        -Oz \
        -o "$OUTPUT_DIR/chromosome_splits/chr${CHROM}.vcf.gz" \
        "$OUTPUT_DIR/two_samples.vcf.gz"
    
    if [ -f "$OUTPUT_DIR/chromosome_splits/chr${CHROM}.vcf.gz" ]; then
        # Index chromosome VCF
        bcftools index -t "$OUTPUT_DIR/chromosome_splits/chr${CHROM}.vcf.gz"
        
        VARIANTS=$(bcftools view -H "$OUTPUT_DIR/chromosome_splits/chr${CHROM}.vcf.gz" | wc -l)
        echo "Successfully created chr${CHROM}.vcf.gz with $VARIANTS variants"
    else
        echo "ERROR: Failed to create chr${CHROM}.vcf.gz"
    fi
done < "$OUTPUT_DIR/chromosome_list.txt"

echo "=== STEP 4: PHASE EACH CHROMOSOME WITH BEAGLE ==="
echo "Running BEAGLE phasing on each chromosome..."

phase_chromosome() {
    local CHROM=$1
    local INPUT_CHR="$OUTPUT_DIR/chromosome_splits/chr${CHROM}.vcf.gz"
    local OUTPUT_CHR="$OUTPUT_DIR/phased_chromosomes/chr${CHROM}_phased"
    local LOG_FILE="$OUTPUT_DIR/logs/chr${CHROM}_beagle.log"
    
    echo "Phasing chromosome $CHROM..."
    
    if [ ! -f "$INPUT_CHR" ]; then
        echo "ERROR: Input file for chromosome $CHROM not found: $INPUT_CHR"
        return 1
    fi

    VARIANTS=$(bcftools view -H "$INPUT_CHR" | wc -l)
    if [ "$VARIANTS" -eq 0 ]; then
        echo "WARNING: No variants found for chromosome $CHROM, skipping..."
        return 1
    fi
    
    echo "Processing $VARIANTS variants on chromosome $CHROM..."
    
    java -Xmx8g -jar "$OUTPUT_DIR/beagle.jar" \
        gt="$INPUT_CHR" \
        out="$OUTPUT_CHR" \
        > "$LOG_FILE" 2>&1
    
    if [ $? -eq 0 ] && [ -f "${OUTPUT_CHR}.vcf.gz" ]; then
        echo "Successfully phased chromosome $CHROM"
        
        bcftools index -t "${OUTPUT_CHR}.vcf.gz"
        
        PHASED_VARIANTS=$(bcftools view -H "${OUTPUT_CHR}.vcf.gz" | wc -l)
        echo "Chromosome $CHROM: $PHASED_VARIANTS variants phased"
        
        return 0
    else
        echo "ERROR: Phasing failed for chromosome $CHROM"
        echo "Check log file: $LOG_FILE"
        return 1
    fi
}

SUCCESS_COUNT=0
TOTAL_CHROM=0

while read -r CHROM; do
    TOTAL_CHROM=$((TOTAL_CHROM + 1))
    if phase_chromosome "$CHROM"; then
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    fi
done < "$OUTPUT_DIR/chromosome_list.txt"

echo "Phasing summary: $SUCCESS_COUNT/$TOTAL_CHROM chromosomes successfully phased"

if [ $SUCCESS_COUNT -eq 0 ]; then
    echo "ERROR: No chromosomes were successfully phased. Check the log files."
    exit 1
fi

echo "=== STEP 5: MERGE PHASED CHROMOSOMES ==="
echo "Merging phased chromosomes into single VCF..."

ls "$OUTPUT_DIR/phased_chromosomes/chr"*"_phased.vcf.gz" 2>/dev/null | sort -V > "$OUTPUT_DIR/phased_files_list.txt"

if [ ! -s "$OUTPUT_DIR/phased_files_list.txt" ]; then
    echo "ERROR: No phased chromosome files found!"
    exit 1
fi

echo "Phased files to merge:"
cat "$OUTPUT_DIR/phased_files_list.txt"

bcftools concat \
    -f "$OUTPUT_DIR/phased_files_list.txt" \
    -Oz \
    -o "$OUTPUT_DIR/final_output/${OUTPUT_PREFIX}.vcf.gz"

bcftools index -t "$OUTPUT_DIR/final_output/${OUTPUT_PREFIX}.vcf.gz"

echo "=== STEP 6: QUALITY ASSESSMENT ==="
echo "Assessing phasing quality..."

TOTAL_VARIANTS=$(bcftools view -H "$OUTPUT_DIR/final_output/${OUTPUT_PREFIX}.vcf.gz" | wc -l)
PHASED_VARIANTS=$(bcftools view -H "$OUTPUT_DIR/final_output/${OUTPUT_PREFIX}.vcf.gz" | grep -c "|")

echo "Total variants in final VCF: $TOTAL_VARIANTS"
echo "Phased variants: $PHASED_VARIANTS"

echo "Calculating phasing rates per sample..."
bcftools query -f '[%SAMPLE\t%GT\n]' "$OUTPUT_DIR/final_output/${OUTPUT_PREFIX}.vcf.gz" | \
    awk '{
        total[$1]++; 
        if($2 ~ /\|/) phased[$1]++
    } END {
        for(sample in total) {
            rate = (phased[sample] ? phased[sample] : 0) / total[sample] * 100
            print sample "\t" rate"%"
        }
    }' > "$OUTPUT_DIR/final_output/phasing_rates_per_sample.txt"

echo "Phasing rates:"
cat "$OUTPUT_DIR/final_output/phasing_rates_per_sample.txt"

echo "Phasing rate summary:"
if [ -f "$OUTPUT_DIR/final_output/phasing_rates_per_sample.txt" ]; then
    awk '{gsub("%", "", $2); print $2}' "$OUTPUT_DIR/final_output/phasing_rates_per_sample.txt" | \
        awk '{sum+=$1; count++} END {
            if(count>0) {
                mean=sum/count; 
                print "Mean phasing rate: " mean "%"
            }
        }'
fi

echo "=== STEP 7: PREPARE FINAL OUTPUT ==="
echo "Copying final files to main directory..."

cp "$OUTPUT_DIR/final_output/${OUTPUT_PREFIX}.vcf.gz" "$INPUT_DIR/two_sample_phased.vcf.gz"
cp "$OUTPUT_DIR/final_output/${OUTPUT_PREFIX}.vcf.gz.tbi" "$INPUT_DIR/two_sample_phased.vcf.gz.tbi"

echo "$SAMPLE1" > "$INPUT_DIR/two_sample_list.txt"
echo "$SAMPLE2" >> "$INPUT_DIR/two_sample_list.txt"

echo "=== STEP 8: VALIDATE PHASED VCF ==="
echo "Validating phased VCF..."

echo "Samples in final VCF:"
bcftools query -l "$INPUT_DIR/two_sample_phased.vcf.gz"

echo "Checking phasing for each sample:"
bcftools query -f '[%SAMPLE: %GT\n]' "$INPUT_DIR/two_sample_phased.vcf.gz" | head -10

echo "Final VCF validation:"
bcftools stats "$INPUT_DIR/two_sample_phased.vcf.gz" | grep -E "^SN" | head -5

echo "=== BEAGLE TWO-SAMPLE PHASING COMPLETE ==="
echo "Phased VCF created: $INPUT_DIR/two_sample_phased.vcf.gz"
echo "Sample list created: $INPUT_DIR/two_sample_list.txt"
echo ""
echo "Files created:"
echo "  - two_sample_phased.vcf.gz (phased VCF with $SAMPLE1 and $SAMPLE2)"
echo "  - two_sample_list.txt (sample list)"
echo "  - $OUTPUT_DIR/final_output/phasing_rates_per_sample.txt"
echo ""
echo "Analysis completed: $(date)"

echo "Cleaning up intermediate files..."
echo "You can manually remove intermediate files from $OUTPUT_DIR if needed"
echo "Directories: chromosome_splits/, phased_chromosomes/, logs/"