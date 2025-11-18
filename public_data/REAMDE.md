README

Public data analysis 

Methods: Coverage Calculation

For each dataset, coverage was computed using a standardized workflow based on BEDTools. 

Input Data:

1. Genome-aligned reads previously converted to BED format (*.markdup.sorted.bed.gz) using convert_bam_to_bed.sh

2. gene-specific region annotation file (FUS_exon_intron.bed) containing genomic intervals for FUS exons and introns of interest.


Coverage Computation

Coverage was calculated using BEDTools coverage with the strand-specific mode enabled.
For each sample, the following command is executed:

bedtools coverage -s -a <region_file> -b <sample.bed.gz>

Sample Labeling and Output Formatting
For every BED file (each biological sample), the script: 
extracts the sample ID from the filename.
Runs the coverage command.
Prepends the sample ID to each coverage line using AWK.
Converts the TAB-delimited BEDTools output into CSV format: coverage_<GENE>_<DATASET>.csv
