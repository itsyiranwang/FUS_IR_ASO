#!/bin/bash

#SBATCH --job-name=bam_to_bed_conversion
#SBATCH --output=bam_to_bed_%A_%a.out
#SBATCH --error=bam_to_bed_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Function to safely write to stats file
write_to_stats() {
    local content="$1"
    local stats_file="$2"
    local max_retries=5
    local retry_count=0
    
    while [[ $retry_count -lt $max_retries ]]; do
        if (
            flock -w 10 200
            echo -e "$content" >> "$stats_file"
        ) 200>"${stats_file}.lock"; then
            break
        else
            log_message "Failed to write to stats file, attempt $((retry_count + 1))"
            sleep 1
            ((retry_count++))
        fi
    done
    
    if [[ $retry_count -eq $max_retries ]]; then
        log_message "ERROR: Failed to write to stats file after $max_retries attempts"
        exit 1
    fi
}

# Load required modules
log_message "Loading required modules..."
if ! module load SAMtools/1.3.1-foss-2016b; then
    log_message "ERROR: Failed to load SAMtools module"
    exit 1
fi
if ! module load BEDTools; then
    log_message "ERROR: Failed to load BEDTools module"
    exit 1
fi

# Set up input/output directories
BAM_DIR=${1:-"/camp/lab/patanir/home/shared/patani-collab/proj-luscombn-patani/working/public-data/neurolincs/nfcore/star_salmon"}
BED_DIR=${2:-"/camp/lab/patanir/home/shared/for_GT/meta/public_data_bed/neurolincs_bam2bed"}
FILTERED_DIR="${BED_DIR}/filtered_bams"
STATS_FILE="${BED_DIR}/mapping_statistics.txt"
TEMP_STATS_DIR="${BED_DIR}/temp_stats"
COMPLETION_DIR="${BED_DIR}/completion_markers"

# Create all required directories
for dir in "$BED_DIR" "$FILTERED_DIR" "$TEMP_STATS_DIR" "$COMPLETION_DIR"; do
    mkdir -p "$dir" || {
        log_message "ERROR: Failed to create directory: $dir"
        exit 1
    }
done

# Get list of BAM files
mapfile -t BAM_FILES < <(find "$BAM_DIR" -name "*.bam" -type f) || {
    log_message "ERROR: Failed to find BAM files in $BAM_DIR"
    exit 1
}

TOTAL_FILES=${#BAM_FILES[@]}
if [ "$TOTAL_FILES" -eq 0 ]; then
    log_message "ERROR: No BAM files found in $BAM_DIR"
    exit 1
fi

# Initialize main stats file only once
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    echo -e "Sample\tTotal_Reads\tUniquely_Mapped_Reads\tPercentage_Unique\tDate_Processed" > "$STATS_FILE"
fi

FILE_INDEX=$((SLURM_ARRAY_TASK_ID-1))
if [ "$FILE_INDEX" -ge "$TOTAL_FILES" ]; then
    log_message "ERROR: Array task ID $SLURM_ARRAY_TASK_ID is beyond number of files ($TOTAL_FILES)"
    exit 1
fi

BAM_FILE=${BAM_FILES[$FILE_INDEX]}
BASE_NAME=$(basename "$BAM_FILE" .bam)
FILTERED_BAM="${FILTERED_DIR}/${BASE_NAME}_unique.bam"
BED_FILE="${BED_DIR}/${BASE_NAME}_unique.bed"
TEMP_STATS="${TEMP_STATS_DIR}/${BASE_NAME}_stats.txt"
COMPLETION_MARKER="${COMPLETION_DIR}/${BASE_NAME}_complete"

log_message "Processing $BAM_FILE..."

# Process BAM file
if ! samtools view -@ 4 -h -F 4 -q 20 "$BAM_FILE" | \
   grep -P '^@|NH:i:1\b' | \
   samtools sort -@ 4 -m 6G - | \
   bedtools bamtobed > "$BED_FILE"; then
    log_message "ERROR: BAM processing failed"
    exit 1
fi

# Count reads
TOTAL_READS=$(samtools view -@ 4 -c "$BAM_FILE") || {
    log_message "ERROR: Failed to count total reads"
    exit 1
}

UNIQUE_READS=$(wc -l < "$BED_FILE") || {
    log_message "ERROR: Failed to count unique reads"
    exit 1
}

PERCENTAGE=$(awk "BEGIN {printf \"%.2f\", ($UNIQUE_READS/$TOTAL_READS)*100}")
CURRENT_DATE=$(date '+%Y-%m-%d')

# Write statistics to temporary file
echo -e "${BASE_NAME}\t${TOTAL_READS}\t${UNIQUE_READS}\t${PERCENTAGE}\t${CURRENT_DATE}" > "$TEMP_STATS"

# Append to main stats file
write_to_stats "$(cat "$TEMP_STATS")" "$STATS_FILE"

# Create completion marker
touch "$COMPLETION_MARKER"

log_message "Completed processing $BAM_FILE"
log_message "Statistics for ${BASE_NAME}:"
log_message "  Total reads: $TOTAL_READS"
log_message "  Uniquely mapped reads: $UNIQUE_READS"
log_message "  Percentage unique: $PERCENTAGE%"

# If this is the last job in the array, wait for all other jobs to complete
if [ "$SLURM_ARRAY_TASK_ID" -eq "$TOTAL_FILES" ]; then
    log_message "Waiting for all jobs to complete..."
    
    # Wait for all completion markers
    while [ "$(ls -1 "$COMPLETION_DIR" | wc -l)" -lt "$TOTAL_FILES" ]; do
        sleep 10
    done
    
    log_message "All jobs completed. Generating summary statistics..."
    
    # Generate summary statistics
    {
        echo -e "\nSummary Statistics:"
        echo "Date of analysis: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "Number of files processed: $TOTAL_FILES"
        
        awk 'BEGIN {FS="\t"}
             NR>1 && $3!="" {
                total+=$3;
                count++;
                if($4!="") sum_pct+=$4
             }
             END {
                if(count>0) {
                    printf "Average unique reads per file: %.2f\n", total/count;
                    printf "Average percentage unique: %.2f%%\n", sum_pct/count
                }
             }' "$STATS_FILE"
    } >> "$STATS_FILE"
    
    # Cleanup
    rm -rf "$TEMP_STATS_DIR" "$COMPLETION_DIR" "$FILTERED_DIR"
fi