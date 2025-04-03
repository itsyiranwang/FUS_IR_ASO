
#script to generate sashimi plots with trackplot 
# G. Manferrari
# 23/10/24

#SBATCH --job-name="trackplot"
#SBATCH --mail-type="END,FAIL"
#SBATCH --mail-user="g.manferrari@ucl.ac.uk"
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --output=trackplot-%A.out
#SBATCH --mem=72GB

ml Python-bundle-PyPI #Python 3.11.3

#source /camp/lab/patanir/home/users/manferg/manferg2/bin/trackplot/

#Annotation reference: Homo_sapiens.GRCh38.99


#FUS intron 6-7
#TDP-43 iiCLIP binding (bedgraphs) from public data:
#grot et al. HEK 293
#Hallegger HEK293
#tollervey brain
#tollervey SHSY5

#FUS 	chr16:31184939-31188357:+
# introns not scaled
sbatch --cpus-per-task=4 --mem=32G --time=12:00:00 --wrap="trackplot \
-e chr16:31184939-31188357:+ \
-r /camp/lab/patanir/home/shared/for_GT/sashimi/Homo_sapiens.GRCh38.96.chr_patch_hapl_scaff.gtf \
--density /camp/lab/patanir/home/users/manferg/manferg2/projects/nf/clip/tdp43/tdp43_metanalysis/results_no_mouse/results/sashimi/bigwig_F.tsv \
-o /camp/lab/patanir/home/users/manferg/manferg2/projects/nf/clip/tdp43/tdp43_metanalysis/results_no_mouse/results/sashimi/plots/fus_intron_6_7_no_scaling.pdf \
--dpi 300 \
--width 10 \
--height 1 \
--normalize-format cpm \
-t 10"