Detailed information on the Riboseq-flow Nextflow DSL2 pipeline used for the analysis can be found in the publication Iosub IA, Wilkins OG, Ule J. Riboseq-flow: A streamlined, reliable pipeline for ribosome profiling data analysis and quality control. Wellcome Open Res. 2024 Apr 11;9:179. doi: 10.12688/wellcomeopenres.21000.1. PMID: 38846930; PMCID: PMC11153996 (https://wellcomeopenresearch.org/articles/9-179).

The following is the summary of steps that were undertaken to use the pipeline to analyse the ribo-seq datasets provided in the paper.

1. Load the required modules

module purge
ml Java/21.0.2
ml Nextflow/24.04.2
ml Singularity/3.6.4
ml Graphviz/2.47.2-GCCcore-10.3.0

2. Specify the directory where Nextflow should cache Singularity container images and store its configuration and log files.

export NXF_SINGULARITY_CACHEDIR=/path/to/singularity
export NXF_HOME=/path/to/.nextflow

3. Pull the pipeline from the GitHub repository:

nextflow pull iraiosub/riboseq-flow -r v1.1.1

4. Run the pipeline on the provided test dataset to check if it works:

nextflow run iraiosub/riboseq-flow -r v1.1.1 -profile test,singularity

5. Create a riboseq_HEK_samplesheet.csv sample sheet with information about the samples as a comma-separated file with 2 columns, and a header row as below:

sample,fastq
riboseq_HEK_1,/path/to/PET7366A2_S333_L004_R1_001.fastq.gz
riboseq_HEK_2,/path/to/PET7366A4_S335_L004_R1_001.fastq.gz
riboseq_HEK_3,/path/to/PET7560A8_R1.fastq.gz

*riboseq_HEK_3 raw FASTQ files PET7560A8_S615_L008_R1_001.fastq.gz and PET7560A8_S76_L006_R1_001.fastq.gz, deposited to the Gene Expression Omnibus (GEO), were concatenated to PET7560A8_R1.fastq.gz as below:

cat PET7560A8_S615_L008_R1_001.fastq.gz PET7560A8_S76_L006_R1_001.fastq.gz > PET7560A8_R1.fastq.gz

PET7560A8_R1.fastq.gz was then used as input for the pipeline.

6. Run the pipeline on the data using the following code:

nextflow run iraiosub/riboseq-flow -r v1.1.1 \
-profile singularity,crick \
-resume \
--input /path/to/riboseq_HEK_samplesheet.csv \
--skip_umi_extract \
--adapter_threeprime AAAAAAAAAA \
--cut_end 3 \
--strandedness forward \
--expected_length "26:34" \
--skip_ribocutter \
--org GRCh38 \
--outdir /path/to/results

Pipeline parameters explained:
-profile: specifies a configuration profile
-resume: specified when restarting the pipeline, Nextflow would then use cached data and continue from where it ended previously
--input: specifies the input sample sheet
--skip_umi_extract: skips UMI extraction as UMIs were not used
--adapter_threeprime: sequence of 3' adapter - specified as instructed in the SMARTer smRNA-Seq Kit for Illumina manual (AAAAAAAAAA in this case as the library prep employs polyadenylation and template switching method)
--cut_end: number of nucleotides to be trimmed from the 5' end of reads - specified as instructed in the SMARTer smRNA-Seq Kit for Illumina manual
--strandedness: specifies library strandedness
--expected_length: expected read lengths range
--skip_ribocutter: skips Ribocutter, a tool for designing sgRNA templates to deplete unwanted abundant contaminants; the Human Ribo-Seq riboPOOL (siTOOLs Biotech) rRNA depletion approach was used instead
--org: automatically downloads and sets up reference files for the human genome
--outdir: specifies the output results directory