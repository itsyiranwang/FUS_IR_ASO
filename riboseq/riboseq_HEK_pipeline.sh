#!/bin/sh
#SBATCH --job-name=nf-riboseq
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=ncpu

module purge
ml Java/21.0.2
ml Nextflow/24.04.2
ml Singularity/3.6.4
ml Graphviz/2.47.2-GCCcore-10.3.0

export NXF_SINGULARITY_CACHEDIR=/path/to/singularity
export NXF_HOME=/path/to/.nextflow

nextflow pull iraiosub/riboseq-flow -r v1.1.1

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