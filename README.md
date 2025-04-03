# FUS_IR_ASO

Code and data relative to: 

"Targeting pathogenic FUS intron retention for amyotrophic lateral sclerosis treatment"  

 Y.Wang et al. 2025


**Project summary/abstract:**

Amyotrophic lateral sclerosis (ALS) is a devastating condition characterised by the progressive degeneration of motor neurons. Evidence suggests that aberrant RNA splicing, including intron retention (IR), is a common feature across ALS models and may contribute to pathogenesis. We previously identified an abundantly elevated IR event in FUS across diverse ALS gene mutations. Here, we reveal that increasing aberrant FUS IR by a splice-modifying antisense oligonucleotide (ASO) is sufficient to disturb cellular homeostasis in both healthy and ALS-VCP mutant human induced pluripotent stem cell (hiPSC)-derived motor neurons, including oxidative stress, mitochondrial dysfunction, lysosomal dysfunction, and apoptosis. Moreover, FUS intron-retaining transgene expression in healthy mice recapitulates key ALS features, including motor deficits and lifespan reduction. In contrast, reducing aberrant FUS IR by downregulating whole FUS transcripts using ASO approaches, including Jacifusen/Ulefnersen, mitigates cellular defects in ALS-VCP mutant motor neurons. Finally, we identify TDP-43 as the most enriched FUS IR-binding protein and show that increased FUS IR disrupts TDP-43 function, leading to RNA missplicing in motor neurons. This is the first demonstration that aberrant IR plays a pathogenic role in neurodegeneration and provides a potential therapeutic opportunity in targeting FUS for forms of ALS beyond those caused by FUS mutations.


**Data Sets:**


1. TDP-43 iiCLIP:  [tdp43_clip](https://github.com/itsyiranwang/FUS_IR_ASO/tree/main/tdp43_clip)    
2. Riboseq: 



**TDP-43 iCLIP**: 
TDP-43 iCLIP data generated in this study in hiPSCs Motor Neurons can be access here:   
GEO Accession to original fastq files   

(pre-demultiplexing) samples run across two different lanes:
`/camp/lab/patanir/home/users/manferg/manferg2/projects/nf/clip/tdp43/motor_neurons/mn_l7/files/demux/MAN4617A1_S4_L001_R1_001.fastq.gz`
`/camp/lab/patanir/home/users/manferg/manferg2/projects/nf/clip/tdp43/motor_neurons/mn_l7/files/demux/MAN4617A1_S4_L002_R1_001.fastq.gz`

De-multiplexed fastq: 
`mn_d7_ctrl_6,/camp/home/manferg/lusgm/projects/nf/clip/tdp43/motor_neurons/mn_l7/files/demux/u_plexed_merged/ultraplex_demux_5bc_NNNNGAATANNNNNN_3bc_NNNCCAGNNN.fastq.gz`

TDP-43 iCLIP data were processed and analysed using the [**nf-core-clipseq pipeline v1.0.0**]
- sample sheet provided to nfcore-clipseq pipeline: [mn_sample_sheet.csv](https://github.com/itsyiranwang/FUS_IR_ASO/blob/main/tdp43_clip/nfcore_clipseq/mn_sample_sheet.csv)    
- script to run nfcore-clipseq pipeline: [run.sh](https://github.com/itsyiranwang/FUS_IR_ASO/blob/main/tdp43_clip/nfcore_clipseq/run.sh)

**iiCLIP TDP-43 data visualisation**:   
- To reduce noise and enhance signal clarity, the identified crosslinks were subsequently processed using [clip-viz](https://github.com/amchakra/clip-viz) _developed by Dr.Anob Chakrabarti _
sample sheet provided to run clip-viz: samplesheet.csv
script to run clip-viz: clip-viz.sh

- Visualisation with [trackplot (v0.4.0)](https://trackplot.readthedocs.io/en/latest/):    
bigwig_F.tsv: input tsv
trackplot_fus.sh: code to run trackplot visualisation

See relative study accession codes to access previously published TDP-43 iCLIP data (fastq files). 
