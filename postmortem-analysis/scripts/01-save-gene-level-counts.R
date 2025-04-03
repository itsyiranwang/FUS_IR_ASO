data_dir <- normalizePath("/camp/lab/patanir/home/users/palk/")
project_dir <- file.path(data_dir, "Projects", "025-yiran-wang-postmortem-collab")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "01-postmortem-analysis")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plots_dir <- file.path(analysis_dir, "plots")
temp_dir <- file.path(analysis_dir, "tmp")
alsod_genes_df <- read.csv(file.path(project_dir, "data_files/alsod_genelist.csv"))
annotation_gtf <- "/nemo/lab/patanir/home/users/palk/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf"
tmp_dir <- file.path(analysis_dir, "slurm")
clinical_metadata_df <- read.csv("/nemo/lab/patanir/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons//analysis/00-manuscript/0103-NYGC-consortium/input_files/samplesheet.csv")
sampleinfo_df <- readRDS(file.path(analysis_save_dir, "00-nygc-patient-sampleinfo-datasheet.rds"))
postmortem_list <- readRDS("/nemo/lab/patanir/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons//analysis/00-manuscript/13-nygc-consortium/rdata_files/03-02-list_of_counts_stat_samplesheet_from_nygc.rds")
feature_ranges <- readRDS("/nemo/lab/patanir/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons//analysis/00-manuscript/05-simulated-reads-differential-profiling/rdata_files/01-all-merged-exonic-and-intronic-segments-in-genes.rds")
# ==========================================================================
# Imports
# ==========================================================================
library("stringr")
library("GenomicRanges")
library("DESeq2")
library("rtracklayer")
library("ggplot2")
library("GenomicFeatures")
library("GenomeInfoDb")
require("Rsamtools")
require("STRINGdb")
# ==========================================================================
# Analysis
# ==========================================================================


# feature_ranges_with_repeats <- significant_de_event_list$repeat_containing_ranges
counts_matrix <- postmortem_list$counts


gene_level_counts_matrix <- do.call(rbind, lapply(split(feature_ranges$hash_id[feature_ranges$type == "exon"], feature_ranges$gene_id[feature_ranges$type == "exon"]), function(hash_id){
	if(length(hash_id) == 1){
		return(counts_matrix[hash_id,])
	}
	return(colSums(counts_matrix[hash_id,]))
}))

genic_count_list <- list(counts = gene_level_counts_matrix, coldata = sampleinfo_df)

saveRDS(object = genic_count_list, file = file.path(analysis_save_dir, "01-gene-level-counts-of-postmortem-samples.rds"))