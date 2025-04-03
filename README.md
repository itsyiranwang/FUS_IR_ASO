# FUS_IR_ASO

Code and data relative to: 
"Targeting pathogenic FUS intron retention for amyotrophic lateral sclerosis treatment", Y.Wang et al. 2025 

** Project summary/abstract**: 

** TDP-43 CLIP data **: 

TDP-43 iCLIP data generated in this study in hiPSCs Motor Neurons can be access here: GEO Accession to original fastq files 

(pre-demultiplexing) samples run across two different lanes:
`/camp/lab/patanir/home/users/manferg/manferg2/projects/nf/clip/tdp43/motor_neurons/mn_l7/files/demux/MAN4617A1_S4_L001_R1_001.fastq.gz`
`/camp/lab/patanir/home/users/manferg/manferg2/projects/nf/clip/tdp43/motor_neurons/mn_l7/files/demux/MAN4617A1_S4_L002_R1_001.fastq.gz`

De-multiplexed fastq: 
`mn_d7_ctrl_6,/camp/home/manferg/lusgm/projects/nf/clip/tdp43/motor_neurons/mn_l7/files/demux/u_plexed_merged/ultraplex_demux_5bc_NNNNGAATANNNNNN_3bc_NNNCCAGNNN.fastq.gz`

nfcore-RNA-seq pipeline: mn_sample_sheet.csv: sample sheet provided to nfcore-clipseq pipeline run.sh: bash script to run nfcore-clipseq pipeline

Clipplot-viz: (CHECK WITH NOBBY IF CODE ALREADY PUBLIC )
samplesheet.csv: sample sheet provided to run clip-viz
clip-viz.sh: script to run nfcore-clipseq clip-viz

Visualisation:
bigwig_F.tsv: input tsv
trackplot_fus.sh: code to run trackplot visualisation

See relative study accession code to access previously published TDP-43 iCLUP raw data (fastq files).