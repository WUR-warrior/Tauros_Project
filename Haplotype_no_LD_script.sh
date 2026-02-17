#!/bin/bash
#SBATCH --comment=773320000
#SBATCH --time=3-25:00:00
#SBATCH --mem=80000
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm.output_haplotype%j.txt
#SBATCH --error=slurm.error_haplotype%j.txt
#SBATCH --job-name=local_ancestry
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

echo "=== SIMPLIFIED TAUROS LOCAL ANCESTRY ANALYSIS ==="
echo "Starting analysis: $(date)"
echo "Job ID: $SLURM_JOB_ID"

ml legacy 
ml plink/2.0
ml bcftools/gcc/64/1.9

cd /lustre/nobackup/WUR/ABGC/shared/Tauros/Beagle_phasing_no_LD/final_output

mkdir -p ancestry_analysis
mkdir -p results

TAUROS_SAMPLE="SAMPLE"  # Change this to match your actual sample name

AURG_VCF="/lustre/nobackup/WUR/ABGC/shared/Tauros/PLINK/No_LD_trimming_for_beagle.vcf.gz" 
echo "=== STEP 1: DIAGNOSE SAMPLE NAMES ==="
echo "Checking samples in VCF file..."
bcftools query -l cattle_phased_no_LD_filter.vcf.gz > ancestry_analysis/vcf_samples.txt
echo "Found $(wc -l < ancestry_analysis/vcf_samples.txt) samples in VCF"

if grep -q "^${TAUROS_SAMPLE}$" ancestry_analysis/vcf_samples.txt; then
    echo "? Tauros sample '$TAUROS_SAMPLE' found in VCF"
else
    echo "? Tauros sample '$TAUROS_SAMPLE' NOT found in VCF"
    echo "Available samples containing 'SAMPLE' or similar:"
    grep -i "sample" ancestry_analysis/vcf_samples.txt || echo "None found"
    echo "Please check the sample name and update TAUROS_SAMPLE variable"
    exit 1
fi

if grep -q "^AURG$" ancestry_analysis/vcf_samples.txt; then
    echo "? AURG sample found in main VCF"
    AURG_AVAILABLE=true
    AURG_VCF_PATH="cattle_phased_no_LD_filter.vcf.gz"
elif [ -n "$AURG_VCF" ] && [ -f "$AURG_VCF" ]; then
    echo "? AURG VCF specified: $AURG_VCF"
    if bcftools query -l "$AURG_VCF" | grep -q "^AURG$"; then
        echo "? AURG sample found in specified VCF"
        AURG_AVAILABLE=true
        AURG_VCF_PATH="$AURG_VCF"
    else
        echo "? AURG sample not found in specified VCF"
        AURG_AVAILABLE=false
    fi
else
    echo "? AURG sample not found in main VCF and no alternative VCF specified"
    AURG_AVAILABLE=false
fi

echo "Checking sample files..."
for breed in LM MA MN PA PO SA; do
    if [ -f ${breed}_samples.txt ]; then
        echo "  ${breed}_samples.txt: $(wc -l < ${breed}_samples.txt) samples"
        if [ -s ${breed}_samples.txt ]; then
            found_samples=$(grep -F -f ${breed}_samples.txt ancestry_analysis/vcf_samples.txt | wc -l)
            total_samples=$(wc -l < ${breed}_samples.txt)
            echo "    $found_samples/$total_samples samples found in VCF"
            
            echo "    First 3 samples in file:"
            head -3 ${breed}_samples.txt | sed 's/^/      /'
            echo "    Matching samples in VCF:"
            grep -F -f ${breed}_samples.txt ancestry_analysis/vcf_samples.txt | head -3 | sed 's/^/      /'
        fi
    else
        echo "  ${breed}_samples.txt: FILE NOT FOUND"
    fi
done

echo ""

echo "=== STEP 2: EXTRACT BREED PANELS (with validation) ==="
echo "Extracting breed reference panels..."

for breed in LM MA MN PA PO SA; do
    echo "Processing breed: $breed"
    
    if [ ! -f ${breed}_samples.txt ]; then
        echo "  ERROR: ${breed}_samples.txt not found, skipping..."
        continue
    fi
    
    if [ ! -s ${breed}_samples.txt ]; then
        echo "  ERROR: ${breed}_samples.txt is empty, skipping..."
        continue
    fi
    
    > ancestry_analysis/${breed}_found_samples.txt
    
    while read sample; do
        # Skip empty lines
        if [ -z "$sample" ]; then
            continue
        fi
        
        if grep -q "^${sample}$" ancestry_analysis/vcf_samples.txt; then
            echo "$sample" >> ancestry_analysis/${breed}_found_samples.txt
        fi
    done < ${breed}_samples.txt
    
    found_samples=$(wc -l < ancestry_analysis/${breed}_found_samples.txt)
    
    if [ $found_samples -eq 0 ]; then
        echo "  ERROR: No samples from ${breed}_samples.txt found in VCF, skipping..."
        continue
    fi
    
    echo "  Extracting $found_samples samples for breed $breed"
    echo "  Samples found:"
    cat ancestry_analysis/${breed}_found_samples.txt | sed 's/^/    /'
    
    bcftools view -S ancestry_analysis/${breed}_found_samples.txt cattle_phased_no_LD_filter.vcf.gz \
        -O z -o ancestry_analysis/${breed}_reference.vcf.gz
    
    if [ $? -eq 0 ]; then
        bcftools index ancestry_analysis/${breed}_reference.vcf.gz
        extracted_samples=$(bcftools query -l ancestry_analysis/${breed}_reference.vcf.gz | wc -l)
        echo "  Successfully extracted $extracted_samples samples"
    else
        echo "  ERROR: Failed to extract samples for breed $breed"
    fi
done

echo "Extracting Tauros target sample..."
if grep -q "^${TAUROS_SAMPLE}$" ancestry_analysis/vcf_samples.txt; then
    echo "$TAUROS_SAMPLE" > ancestry_analysis/tauros_sample.txt
    bcftools view -S ancestry_analysis/tauros_sample.txt cattle_phased_no_LD_filter.vcf.gz \
        -O z -o ancestry_analysis/tauros_target.vcf.gz
    
    if [ $? -eq 0 ]; then
        bcftools index ancestry_analysis/tauros_target.vcf.gz
        echo "  Successfully extracted Tauros sample: $TAUROS_SAMPLE"
    else
        echo "  ERROR: Failed to extract Tauros sample"
        exit 1
    fi
else
    echo "  ERROR: Tauros sample '$TAUROS_SAMPLE' not found in VCF"
    exit 1
fi

echo "=== DEBUGGING: VERIFY EXTRACTED SAMPLES ==="
for breed in LM MA MN PA PO SA; do
    if [ -f ancestry_analysis/${breed}_reference.vcf.gz ]; then
        echo "Breed $breed extracted samples:"
        bcftools query -l ancestry_analysis/${breed}_reference.vcf.gz | sed 's/^/  /'
        echo "  Total: $(bcftools query -l ancestry_analysis/${breed}_reference.vcf.gz | wc -l)"
    fi
done

echo "Tauros sample extracted:"
if [ -f ancestry_analysis/tauros_target.vcf.gz ]; then
    bcftools query -l ancestry_analysis/tauros_target.vcf.gz | sed 's/^/  /'
fi

echo "=== STEP 3: VALIDATE EXTRACTED FILES ==="
if [ ! -f ancestry_analysis/tauros_target.vcf.gz ]; then
    echo "ERROR: Tauros target file not created. Cannot proceed."
    exit 1
fi

available_breeds=()
for breed in LM MA MN PA PO SA; do
    if [ -f ancestry_analysis/${breed}_reference.vcf.gz ]; then
        available_breeds+=($breed)
        echo "  Breed $breed: Available"
    else
        echo "  Breed $breed: Not available (skipped)"
    fi
done

if [ ${#available_breeds[@]} -eq 0 ]; then
    echo "ERROR: No breed reference files were created. Cannot proceed."
    exit 1
fi

echo "Proceeding with ${#available_breeds[@]} breeds: ${available_breeds[*]}"

calculate_similarity() {
    local tauros_vcf=$1
    local ref_vcf=$2
    local breed=$3
    local region=$4
    local output_prefix=$5

    bcftools query -f '%CHROM\t%POS[\t%GT]\n' -r $region $tauros_vcf > ${output_prefix}_tauros.txt
    bcftools query -f '%CHROM\t%POS[\t%GT]\n' -r $region $ref_vcf > ${output_prefix}_${breed}.txt
    
    awk '
    BEGIN {matches=0; total=0}
    NR==FNR {tauros[$1":"$2]=$3; next}
    {
        pos=$1":"$2
        if(pos in tauros && tauros[pos]!="./." && $3!="./." && tauros[pos]==$3) {
            matches++
        }
        if(pos in tauros && tauros[pos]!="./." && $3!="./.")
            total++
    }
    END {
        if(total>0) print matches/total
        else print 0
    }' ${output_prefix}_tauros.txt ${output_prefix}_${breed}.txt
    
    rm -f ${output_prefix}_tauros.txt ${output_prefix}_${breed}.txt
}

bcftools query -f '%CHROM\n' ancestry_analysis/tauros_target.vcf.gz | sort -u > ancestry_analysis/chromosomes.txt

echo "=== STEP 4: SLIDING WINDOW ANALYSIS ==="
echo "Running sliding window analysis..."
echo "chromosome,start,end,breed,similarity,window_size" > results/ancestry_windows.csv

while read chrom; do
    echo "Processing chromosome: $chrom"
    
    chrom_length=$(bcftools query -f '%POS\n' -r $chrom ancestry_analysis/tauros_target.vcf.gz | tail -1)

    if [ -z "$chrom_length" ]; then
        echo "Skipping empty chromosome $chrom"
        continue
    fi
    
    window_size=5000000  
    step_size=2500000    
    
    for start in $(seq 1 $step_size $chrom_length); do
        end=$((start + window_size))
        if [ $end -gt $chrom_length ]; then
            end=$chrom_length
        fi
        
        region="${chrom}:${start}-${end}"
        echo "  Analyzing window: $region"
        
        best_breed=""
        best_similarity=0
        
        for breed in "${available_breeds[@]}"; do
            if [ -f ancestry_analysis/${breed}_reference.vcf.gz ]; then
                similarity=$(calculate_similarity \
                    ancestry_analysis/tauros_target.vcf.gz \
                    ancestry_analysis/${breed}_reference.vcf.gz \
                    $breed \
                    $region \
                    ancestry_analysis/tmp_${chrom}_${start})
                
                if [ $(echo "$similarity > $best_similarity" | bc -l) -eq 1 ]; then
                    best_similarity=$similarity
                    best_breed=$breed
                fi
            fi
        done
        
        if [ $(echo "$best_similarity > 0.5" | bc -l) -eq 1 ]; then
            echo "$chrom,$start,$end,$best_breed,$best_similarity,$window_size" >> results/ancestry_windows.csv
        fi
    done
done < ancestry_analysis/chromosomes.txt

echo "=== STEP 5: AURG COMPARISON (FIXED) ==="
if [ "$AURG_AVAILABLE" = true ]; then
    echo "Comparing with AURG (handling phased/unphased format differences)..."
    echo "chromosome,start,end,breed,breed_similarity,aurg_similarity,common_sites,total_sites" > results/aurg_comparison.csv
    
    tail -n +2 results/ancestry_windows.csv | while IFS=',' read chrom start end breed similarity window_size variant_count; do
        region="${chrom}:${start}-${end}"
        echo "  Comparing region $region with AURG"
        
        bcftools query -f '%CHROM\t%POS[\t%GT]\n' -r $region -s $TAUROS_SAMPLE ancestry_analysis/tauros_target.vcf.gz > tmp_tauros.txt
        bcftools query -f '%CHROM\t%POS[\t%GT]\n' -r $region -s AURG "$AURG_VCF_PATH" > tmp_aurg.txt
        
        aurg_result=$(awk '
        BEGIN {matches=0; total=0; common=0}
        NR==FNR {
            # Normalize Tauros genotype (convert phased to unphased)
            gt = $3
            gsub(/\|/, "/", gt)
            tauros[$1":"$2] = gt
            next
        }
        {
            pos = $1":"$2
            if(pos in tauros) {
                common++
                aurg_gt = $3
                gsub(/\|/, "/", aurg_gt)
                tauros_gt = tauros[pos]
                
                if(tauros_gt != "./." && tauros_gt != "." && aurg_gt != "./." && aurg_gt != ".") {
                    total++
                    
                    if(tauros_gt == aurg_gt) {
                        matches++
                    }
                    else if((tauros_gt == "0/1" && aurg_gt == "1/0") || (tauros_gt == "1/0" && aurg_gt == "0/1")) {
                        matches++
                    }
                }
            }
        }
        END {
            if(total > 0) {
                similarity = matches/total
                print similarity "," common "," total
            } else {
                print "0,0,0"
            }
        }' tmp_tauros.txt tmp_aurg.txt)
        
        aurg_similarity=$(echo $aurg_result | cut -d',' -f1)
        common_sites=$(echo $aurg_result | cut -d',' -f2)
        total_sites=$(echo $aurg_result | cut -d',' -f3)
        
        echo "    AURG similarity: $aurg_similarity (${total_sites} comparable sites from ${common_sites} common positions)"
        
        echo "$chrom,$start,$end,$breed,$similarity,$aurg_similarity,$common_sites,$total_sites" >> results/aurg_comparison.csv
        
        rm -f tmp_tauros.txt tmp_aurg.txt
    done
    
    echo "AURG comparison completed with format normalization"
    echo "Summary statistics:"
    echo "  Total windows: $(tail -n +2 results/aurg_comparison.csv | wc -l)"
    echo "  Windows with AURG data: $(tail -n +2 results/aurg_comparison.csv | awk -F',' '$8 > 0' | wc -l)"
    echo "  Average AURG similarity: $(tail -n +2 results/aurg_comparison.csv | awk -F',' '$6 > 0 {sum+=$6; count++} END {if(count>0) print sum/count; else print 0}')"
    
else
    echo "AURG comparison skipped (AURG sample not available)"
fi

echo "=== STEP 6: BASIC SUMMARY ==="

echo "Creating breed composition summary..."
echo "breed,segment_count,total_length_mb,avg_similarity" > results/breed_summary.csv

for breed in "${available_breeds[@]}"; do
    count=$(awk -v b="$breed" -F',' '$4==b' results/ancestry_windows.csv | wc -l)
    if [ $count -gt 0 ]; then
        total_length=$(awk -v b="$breed" -F',' '$4==b {sum+=$3-$2} END {print sum/1000000}' results/ancestry_windows.csv)
        avg_similarity=$(awk -v b="$breed" -F',' '$4==b {sum+=$5; count++} END {if(count>0) print sum/count; else print 0}' results/ancestry_windows.csv)
        echo "$breed,$count,$total_length,$avg_similarity" >> results/breed_summary.csv
    fi
done

echo "Creating chromosome summary..."
echo "chromosome,segment_count,total_length_mb,dominant_breed" > results/chromosome_summary.csv

for chrom in $(tail -n +2 results/ancestry_windows.csv | cut -d',' -f1 | sort -u); do
    count=$(awk -v c="$chrom" -F',' '$1==c' results/ancestry_windows.csv | wc -l)
    total_length=$(awk -v c="$chrom" -F',' '$1==c {sum+=$3-$2} END {print sum/1000000}' results/ancestry_windows.csv)
    dominant_breed=$(awk -v c="$chrom" -F',' '$1==c {breed[$4]++} END {max=0; for(b in breed) if(breed[b]>max) {max=breed[b]; dominant=b}} END {print dominant}' results/ancestry_windows.csv)
    echo "$chrom,$count,$total_length,$dominant_breed" >> results/chromosome_summary.csv
done

echo "=== STEP 7: SIMPLE VISUALIZATION DATA ==="

echo "Creating BED file for genome visualization..."
echo "track name=\"Tauros_Ancestry\" description=\"Tauros Local Ancestry\" visibility=2 itemRgb=On" > results/ancestry_segments.bed

declare -A breed_colors
breed_colors[LM]="255,0,0"      
breed_colors[MA]="0,255,0"     
breed_colors[MN]="0,0,255"     
breed_colors[PA]="255,255,0"   
breed_colors[PO]="255,0,255"   
breed_colors[SA]="0,255,255"    

tail -n +2 results/ancestry_windows.csv | while IFS=',' read chrom start end breed similarity window_size; do
    score=$(echo "$similarity * 1000" | bc -l | cut -d'.' -f1)
    color=${breed_colors[$breed]}
    echo "$chrom	$start	$end	${breed}_${score}	$score	.	$start	$end	$color" >> results/ancestry_segments.bed
done

echo "=== FINAL SUMMARY ==="
echo "Analysis completed: $(date)"
echo "Tauros sample used: $TAUROS_SAMPLE"
echo "AURG comparison: $([ "$AURG_AVAILABLE" = true ] && echo "Available" || echo "Not available")"
echo ""
echo "BREED COMPOSITION:"
echo "==================="
cat results/breed_summary.csv
echo ""
echo "CHROMOSOME SUMMARY:"
echo "==================="
cat results/chromosome_summary.csv
echo ""

if [ "$AURG_AVAILABLE" = true ]; then
    echo "AURG COMPARISON SUMMARY:"
    echo "========================"
    echo "Average AURG similarity by breed:"
    for breed in "${available_breeds[@]}"; do
        avg=$(awk -v b="$breed" -F',' '$4==b {sum+=$6; count++} END {if(count>0) print sum/count; else print 0}' results/aurg_comparison.csv)
        count=$(awk -v b="$breed" -F',' '$4==b' results/aurg_comparison.csv | wc -l)
        if [ $count -gt 0 ]; then
            echo "  $breed: $avg (n=$count)"
        fi
    done
    echo ""
fi

echo "OUTPUT FILES:"
echo "============="
echo "- results/ancestry_windows.csv (main results)"
if [ "$AURG_AVAILABLE" = true ]; then
    echo "- results/aurg_comparison.csv (AURG comparison)"
fi
echo "- results/breed_summary.csv (breed composition)"
echo "- results/chromosome_summary.csv (chromosome summary)"
echo "- results/ancestry_segments.bed (genome browser format)"
echo ""
echo "To visualize results, load CSV files into R, Excel, or any plotting software."
echo "Use the BED file to view results in IGV or UCSC Genome Browser."