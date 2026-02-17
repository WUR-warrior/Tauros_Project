#!/bin/bash
#SBATCH --comment=773320001
#SBATCH --time=2-12:00:00
#SBATCH --mem=64000
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm.output_beagle%j.txt
#SBATCH --error=slurm.error_beagle%j.txt
#SBATCH --job-name=beagle_phase
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

echo "=== BEAGLE PHASING PIPELINE FOR CATTLE BREEDS ==="
echo "Starting phasing analysis: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"

ml legacy
ml bcftools/1.9
ml plink/2.0

INPUT_DIR="/lustre/nobackup/WUR/ABGC/shared/Tauros/PLINK"
OUTPUT_DIR="/lustre/nobackup/WUR/ABGC/shared/Tauros/Beagle_phasing_no_LD"

mkdir -p "$OUTPUT_DIR"

INPUT_VCF="$INPUT_DIR/No_LD_ancestor_trimming_for_beagle.vcf.gz"
OUTPUT_PREFIX="cattle_phased_no_LD_filter"
REFERENCE_PANEL="" 

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

echo "Input VCF statistics:"
bcftools stats "$INPUT_VCF" | grep -E "^SN" | head -10

echo "Extracting sample information..."
bcftools query -l "$INPUT_VCF" > "$OUTPUT_DIR/sample_list.txt"
TOTAL_SAMPLES=$(wc -l < "$OUTPUT_DIR/sample_list.txt")
echo "Total samples: $TOTAL_SAMPLES"

echo "=== STEP 2: SPLIT BY CHROMOSOME ==="
echo "Splitting VCF by chromosome for parallel processing..."

mkdir -p "$OUTPUT_DIR/chromosome_splits"
mkdir -p "$OUTPUT_DIR/phased_chromosomes"
mkdir -p "$OUTPUT_DIR/logs"
mkdir -p "$OUTPUT_DIR/final_output"

bcftools query -f '%CHROM\n' "$INPUT_VCF" | sort -u > "$OUTPUT_DIR/chromosome_list.txt"

echo "Found chromosomes:"
cat "$OUTPUT_DIR/chromosome_list.txt"

while read -r CHROM; do
    echo "Extracting chromosome $CHROM..."
    bcftools view \
        -r "$CHROM" \
        -Oz \
        -o "$OUTPUT_DIR/chromosome_splits/chr${CHROM}.vcf.gz" \
        "$INPUT_VCF"
    
    if [ -f "$OUTPUT_DIR/chromosome_splits/chr${CHROM}.vcf.gz" ]; then
        # Index chromosome VCF
        bcftools index -t "$OUTPUT_DIR/chromosome_splits/chr${CHROM}.vcf.gz"
        echo "Successfully created $OUTPUT_DIR/chromosome_splits/chr${CHROM}.vcf.gz"
    else
        echo "ERROR: Failed to create $OUTPUT_DIR/chromosome_splits/chr${CHROM}.vcf.gz"
    fi
done < "$OUTPUT_DIR/chromosome_list.txt"

echo "=== STEP 3: PHASE EACH CHROMOSOME WITH BEAGLE ==="
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
    
    java -jar "$OUTPUT_DIR/beagle.jar" \
        gt="$INPUT_CHR" \
        out="$OUTPUT_CHR" \
        > "$LOG_FILE" 2>&1
    
    if [ $? -eq 0 ] && [ -f "${OUTPUT_CHR}.vcf.gz" ]; then
        echo "Successfully phased chromosome $CHROM"
        
        bcftools index -t "${OUTPUT_CHR}.vcf.gz"
        
        VARIANTS=$(bcftools view -H "${OUTPUT_CHR}.vcf.gz" | wc -l)
        echo "Chromosome $CHROM: $VARIANTS variants phased"
        
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

echo "=== STEP 4: MERGE PHASED CHROMOSOMES ==="
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

echo "=== STEP 5: QUALITY ASSESSMENT ==="
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

echo "Phasing rate summary:"
if [ -f "$OUTPUT_DIR/final_output/phasing_rates_per_sample.txt" ]; then
    awk '{gsub("%", "", $2); print $2}' "$OUTPUT_DIR/final_output/phasing_rates_per_sample.txt" | \
        awk '{sum+=$1; sumsq+=$1*$1} END {
            if(NR>0) {
                mean=sum/NR; 
                print "Mean phasing rate: " mean "%"
                print "Std deviation: " sqrt(sumsq/NR - mean*mean) "%"
            }
        }'
fi

echo "=== STEP 6: CREATE SAMPLE LISTS FOR ANCESTRY ANALYSIS ==="
echo "Creating breed-specific sample lists..."

echo "Creating sample lists by breed..."
echo "NOTE: Update the grep patterns below to match your sample naming convention"

grep "^LM" "$OUTPUT_DIR/sample_list.txt" > "$OUTPUT_DIR/final_output/LM_samples.txt" 2>/dev/null || echo "No LM samples found"
grep "^MA" "$OUTPUT_DIR/sample_list.txt" > "$OUTPUT_DIR/final_output/MA_samples.txt" 2>/dev/null || echo "No MA samples found"
grep "^MN" "$OUTPUT_DIR/sample_list.txt" > "$OUTPUT_DIR/final_output/MN_samples.txt" 2>/dev/null || echo "No MN samples found"
grep "^PA" "$OUTPUT_DIR/sample_list.txt" > "$OUTPUT_DIR/final_output/PA_samples.txt" 2>/dev/null || echo "No PA samples found"
grep "^PO" "$OUTPUT_DIR/sample_list.txt" > "$OUTPUT_DIR/final_output/PO_samples.txt" 2>/dev/null || echo "No PO samples found"
grep "^SA" "$OUTPUT_DIR/sample_list.txt" > "$OUTPUT_DIR/final_output/SA_samples.txt" 2>/dev/null || echo "No SA samples found"
grep -E "^TAUROS|^Tauros" "$OUTPUT_DIR/sample_list.txt" > "$OUTPUT_DIR/final_output/admixed_samples.txt" 2>/dev/null || echo "No Tauros samples found"

echo "Sample counts by breed:"
for breed in LM MA MN PA PO SA; do
    if [ -f "$OUTPUT_DIR/final_output/${breed}_samples.txt" ] && [ -s "$OUTPUT_DIR/final_output/${breed}_samples.txt" ]; then
        count=$(wc -l < "$OUTPUT_DIR/final_output/${breed}_samples.txt")
        echo "  $breed: $count samples"
    else
        echo "  $breed: 0 samples"
    fi
done

if [ -f "$OUTPUT_DIR/final_output/admixed_samples.txt" ] && [ -s "$OUTPUT_DIR/final_output/admixed_samples.txt" ]; then
    count=$(wc -l < "$OUTPUT_DIR/final_output/admixed_samples.txt")
    echo "  Tauros (admixed): $count samples"
else
    echo "  Tauros (admixed): 0 samples"
fi

echo "=== STEP 7: PREPARE FOR ANCESTRY ANALYSIS ==="
echo "Copying final files to main directory..."

cp "$OUTPUT_DIR/final_output/${OUTPUT_PREFIX}.vcf.gz" "$INPUT_DIR/cattle_phased.vcf.gz"
cp "$OUTPUT_DIR/final_output/${OUTPUT_PREFIX}.vcf.gz.tbi" "$INPUT_DIR/cattle_phased.vcf.gz.tbi"

cp "$OUTPUT_DIR/final_output/"*"_samples.txt" "$INPUT_DIR/" 2>/dev/null || echo "No sample lists to copy"

echo "=== STEP 8: VALIDATE PHASED VCF ==="
echo "Validating phased VCF for ancestry analysis..."

echo "Checking for required samples..."
required_samples=("AURG" "Tauros")
for sample in "${required_samples[@]}"; do
    if bcftools query -l "$INPUT_DIR/cattle_phased.vcf.gz" | grep -q "$sample"; then
        echo "? Found required sample: $sample"
    else
        echo "? WARNING: Required sample '$sample' not found in VCF"
        echo "  Available samples (first 10):"
        bcftools query -l "$INPUT_DIR/cattle_phased.vcf.gz" | head -10
    fi
done

echo "Final VCF validation:"
bcftools stats "$INPUT_DIR/cattle_phased.vcf.gz" | grep -E "^SN" | head -5

echo "=== BEAGLE PHASING COMPLETE ==="
echo "Phased VCF created: $INPUT_DIR/cattle_phased.vcf.gz"
echo "Sample lists created for ancestry analysis"
echo "Ready to run the ancestry analysis script"
echo ""
echo "Next steps:"
echo "1. Review and edit the sample list files if needed"
echo "2. Ensure AURG and Tauros samples are correctly named"
echo "3. Run the ancestry analysis script"
echo ""
echo "Files created:"
echo "  - cattle_phased.vcf.gz (main phased VCF)"
echo "  - *_samples.txt (breed-specific sample lists)"
echo "  - $OUTPUT_DIR/final_output/phasing_rates_per_sample.txt"
echo ""
echo "Analysis completed: $(date)"

echo "Cleaning up intermediate files..."
echo "You can manually remove intermediate files from $OUTPUT_DIR if needed"
echo "Directories: chromosome_splits/, phased_chromosomes/, logs/"