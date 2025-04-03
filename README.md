# FUS_IR_ASO

Code and data relative to: 

"Targeting pathogenic FUS intron retention for amyotrophic lateral sclerosis treatment"  

 Y.Wang et al. 2025

**Data Sets:**

1. TDP-43 iiCLIP:  [tdp43_clip](https://github.com/itsyiranwang/FUS_IR_ASO/tree/main/tdp43_clip)    
2. Riboseq: 
3. ALS progression:  [postmortem-analysis](https://github.com/itsyiranwang/FUS_IR_ASO/tree/main/postmortem-analysis)



**TDP-43 iCLIP**: 
TDP-43 iCLIP data generated in this study in hiPSC-Motor Neurons can be access here:   
GEO Accession to original fastq files   

De-multiplexed fastq: 
`mn_d7_ctrl_6,/camp/home/manferg/lusgm/projects/nf/clip/tdp43/motor_neurons/mn_l7/files/demux/u_plexed_merged/ultraplex_demux_5bc_NNNNGAATANNNNNN_3bc_NNNCCAGNNN.fastq.gz`

TDP-43 iCLIP data were processed and analysed using the [nf-core-clipseq pipeline v1.0.0](https://nf-co.re/clipseq/1.0.0/)  
- sample sheet provided to nfcore-clipseq: [mn_sample_sheet.csv](https://github.com/itsyiranwang/FUS_IR_ASO/blob/main/tdp43_clip/nfcore_clipseq/mn_sample_sheet.csv)    
- script to run nfcore-clipseq: [run.sh](https://github.com/itsyiranwang/FUS_IR_ASO/blob/main/tdp43_clip/nfcore_clipseq/run.sh)

**iiCLIP TDP-43 data visualisation**:  To reduce noise and enhance signal clarity, the identified crosslinks were subsequently processed using [clip-viz](https://github.com/amchakra/clip-viz)  (_developed by Dr.Anob Chakrabarti_ ). The smoothed crosslinks were then visualised using [trackplot (v0.4.0)](https://trackplot.readthedocs.io/en/latest/) to generate genome browser-style representations of TDP-43 binding profiles across FUS intronic regions.

**clip-viz:**   
- sample sheet provided to run clip-viz: [samplesheet.csv](https://github.com/itsyiranwang/FUS_IR_ASO/blob/main/tdp43_clip/clip_viz/samplesheet.csv)     
- script to run clip-viz: [clip-viz.sh](https://github.com/itsyiranwang/FUS_IR_ASO/blob/main/tdp43_clip/clip_viz/clip_viz.sh)   

**trackplot:**    
- trackplot input (bigwigs generated by clip_viz) : [bigwig_F.tsv](https://github.com/itsyiranwang/FUS_IR_ASO/blob/main/tdp43_clip/trackplot/bigwig_F.tsv)    
- script to run trackplot: [trackplot_fus.sh](https://github.com/itsyiranwang/FUS_IR_ASO/blob/main/tdp43_clip/trackplot/trackplot_fus.sh)   

_See relative study accession codes to access previously published TDP-43 iCLIP data (fastq files)._ 
