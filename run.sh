#!/bin/sh

# TDP-43 iiCLIP on hiPSC-derived - d7 MN (post terminal differentiation) and iPSCs- GM samples first trial with l7 adapters 
# G. Manferrari
# 17th January 2022

#SBATCH --job-name="mn_l7_iiclip"
#SBATCH --mail-type="END,FAIL"
#SBATCH --mail-user="g.manferrari@ucl.ac.uk"
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --output=giulia-%A.out
#SBATCH --mem=8G

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/20.10.0
ml Singularity/3.6.4
ml Graphviz/2.38.0-foss-2016b

export NXF_SINGULARITY_CACHEDIR=/camp/lab/luscomben/home/shared/singularity
export ICOUNT_TMP_ROOT=/camp/lab/luscomben/home/users/manferg/tmp/icount

## UPDATE PIPLINE
nextflow pull nf-core/clipseq

## RUN PIPELINE
nextflow run nf-core/clipseq \
-profile crick \
-r master \
-resume \
--input /camp/home/manferg/lusgm/projects/nf/clip/tdp43/motor_neurons/mn_l7/files/mn_sample_sheet.csv \
--smrna_org human \
--fasta /camp/home/manferg/lusgm/ref/GRCh38.release34.primary_assembly.genome.fa.gz \
--gtf /camp/home/manferg/lusgm/ref/gencode.v38.annotation.gtf.gz \
--save_index true \