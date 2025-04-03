#!/bin/sh
# Script to run clip-viz
# A. M. Chakrabarti
# 26th November 2021
​# edited on 23/10/24

ml purge
ml Nextflow/21.10.3
ml Singularity/3.6.4
ml Graphviz/2.38.0-foss-2016b
​
export NXF_SINGULARITY_CACHEDIR=/camp/apps/misc/stp/babs/nf-core/singularity/clipseq/1.0.0/
​
REFDIR=/camp/lab/patanir/home/users/manferg/manferg2/ref/
​
cd /camp/lab/patanir/home/users/manferg/manferg2/projects/nf/clip/tdp43/tdp43_metanalysis/results_no_mouse/results/
​
nextflow pull amchakra/clip-viz -r main
​
nextflow run amchakra/clip-viz -r main \
-resume \
-profile crick \
--input samplesheet.csv \
--outdir igv_tracks \
--fai $REFDIR/GRCh38.release34.primary_assembly.genome.fai \
--bin_size 1 \
--smooth_length 25