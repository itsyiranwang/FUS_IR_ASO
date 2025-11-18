#!/bin/bash

# Load required modules
ml BEDTools

# Usage function
usage() {
    echo "Usage: $0 -g GENE_NAME [-d DATASET1,DATASET2,...] [-b BASE_DIR] [-r REGION_FILE]"
    echo "  -g GENE_NAME      Name of the gene to analyze"
    echo "  -d DATASETS       Comma-separated list of datasets to process (default: all available)"
    echo "  -b BASE_DIR       Base directory (default: /camp/lab/patanir/home/shared/for_GT/meta)"
    echo "  -r REGION_FILE    Region file name (default: exon_intron.bed)"
    echo "  -c CPU            CPUs per task (default: 4)"
    echo "  -m MEM            Memory allocation (default: 32G)"
    echo "  -t TIME           Time allocation (default: 12:00:00)"
    echo "  -h                Display this help message"
    exit 1
}

# Default values
BASE_DIR="/camp/lab/patanir/home/shared/for_GT/meta"
REGION_FILE="exon_intron.bed"
CPUS=4
MEM="32G"
TIME="12:00:00"
ALL_DATASETS=("neuronlincs","answerals","nygc")

# Parse command line arguments
while getopts "g:d:b:r:c:m:t:h" opt; do
    case $opt in
        g) GENE_NAME="$OPTARG" ;;
        d) IFS=',' read -r -a DATASETS <<< "$OPTARG" ;;
        b) BASE_DIR="$OPTARG" ;;
        r) REGION_FILE="$OPTARG" ;;
        c) CPUS="$OPTARG" ;;
        m) MEM="$OPTARG" ;;
        t) TIME="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check if gene name is provided
if [ -z "$GENE_NAME" ]; then
    echo "Error: Gene name is required"
    usage
fi

# If no datasets specified, use all
if [ -z "${DATASETS[*]}" ]; then
    DATASETS=("${ALL_DATASETS[@]}")
fi

echo "Processing gene: $GENE_NAME"
echo "Using datasets: ${DATASETS[*]}"

# Create output directory for gene if it doesn't exist
mkdir -p "$BASE_DIR/$GENE_NAME"

# Submit job for each dataset
for dataset in "${DATASETS[@]}"; do
    bed_dir="$BASE_DIR/public_data_bed/${dataset}_bam2bed/"
    out_dir="$BASE_DIR/$GENE_NAME/${dataset}_bam2bed"
    region_dir="$BASE_DIR/regions"
    
    # Ensure output directory exists
    mkdir -p "$out_dir"
    
    # Create the full path to the region file
    FULL_REGION_PATH="${region_dir}/${REGION_FILE}"
    
    # Check if region file exists
    if [ ! -f "$FULL_REGION_PATH" ]; then
        echo "Warning: Region file $FULL_REGION_PATH does not exist for dataset $dataset"
        continue
    fi
    
    # Check if bed directory exists
    if [ ! -d "$bed_dir" ]; then
        echo "Warning: Bed directory $bed_dir does not exist for dataset $dataset"
        continue
    fi
    
    # Check if any bed files exist
    if [ -z "$(ls -A "$bed_dir"/*.bed.gz 2>/dev/null)" ]; then
        echo "Warning: No bed.gz files found in $bed_dir for dataset $dataset"
        continue
    fi
    
    echo "Submitting job for dataset: $dataset"
    
    sbatch --job-name="${GENE_NAME}_${dataset}" \
           --cpus-per-task="$CPUS" \
           --mem="$MEM" \
           --time="$TIME" \
           --output="$out_dir/${GENE_NAME}_${dataset}_coverage.log" \
           --wrap="
    echo 'Processing dataset: $dataset for gene: $GENE_NAME'
    echo 'Region file path: $FULL_REGION_PATH'
    output_file='$out_dir/coverage_${GENE_NAME}_${dataset}.csv'
    
    # Write header
    echo 'sample,chr,start,end,region,dot,strand,read_count,base_coverage,region_size,coverage_fraction' > \"\$output_file\"
    
    # Process each bed file
    for bed in \"$bed_dir\"/*.bed.gz; do
        if [ -f \"\$bed\" ]; then
            sample=\$(basename \"\$bed\" .markdup.sorted.bed.gz)
            echo \"Processing sample: \$sample\"
            
            bedtools coverage -s -a \"$FULL_REGION_PATH\" -b \"\$bed\" | \\
            awk -v sample=\"\$sample\" '{print sample \",\" \$0}' | \\
            tr \"\t\" \",\" >> \"\$output_file\"
        fi
    done
    
    echo 'Finished processing dataset: $dataset for gene: $GENE_NAME'
    "
done

echo "All jobs submitted. Results will be saved in the respective dataset folders."