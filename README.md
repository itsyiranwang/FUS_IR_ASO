# FUS_IR_ASO

Code and data relative to: 

"Targeting pathogenic FUS intron retention for amyotrophic lateral sclerosis treatment"  

 Y.Wang et al. 2025

**Data Sets:**

1. TDP-43 iiCLIP: [tdp43_clip](https://github.com/itsyiranwang/FUS_IR_ASO/tree/main/tdp43_clip)    
2. Riboseq: [riboseq](https://github.com/itsyiranwang/FUS_IR_ASO/tree/main/riboseq)
3. NYGC-ALS consortium data - ALS progresison:  [postmortem-analysis](https://github.com/itsyiranwang/FUS_IR_ASO/tree/main/postmortem-analysis)


**GEO Accession**:   
Raw (fastq) and processed (bigwig) data used in the study are available through the GEO accession records GSE293827 and GSE293828.


**TDP-43 iCLIP**: 
TDP-43 iCLIP data generated in this study in hiPSC-Motor Neurons can be accessed here:   
 
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

**ENCODE eCLIP:** 
- ENCODE input: [fe_clip_IR.RData](https://github.com/itsyiranwang/FUS_IR_ASO/blob/main/encode_eclip/input_files/fe_clip_IR.RData)  
- ENCODE analysis script: [FUS_analysis.R](https://github.com/itsyiranwang/FUS_IR_ASO/blob/main/encode_eclip/script/FUS_analysis.R)  

**Riboseq:** 
- riboseq input: [riboseq_HEK_samplesheet.csv](https://github.com/itsyiranwang/FUS_IR_ASO/blob/main/riboseq/riboseq_HEK_samplesheet.csv)  
- riboseq pipeline script: [riboseq_HEK_pipeline.sh](https://github.com/itsyiranwang/FUS_IR_ASO/blob/main/riboseq/riboseq_HEK_pipeline.sh)  

**ALS progression:** 
- NYGC-ALS consortium postmortem data analysis input: [fus_7_ratio_thoracic.csv](https://github.com/itsyiranwang/FUS_IR_ASO/blob/main/postmortem-analysis/input_files/fus_7_ratio_thoracic.csv)  
- NYGC-ALS consortium postmortem data analysis scripts: [scripts](https://github.com/itsyiranwang/FUS_IR_ASO/tree/main/postmortem-analysis/scripts)  
