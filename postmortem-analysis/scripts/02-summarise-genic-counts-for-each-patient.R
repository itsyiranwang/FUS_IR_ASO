data_dir <- normalizePath("/nemo/lab/patanir/home/users/palk/")
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
postmortem_list <- readRDS("/nemo/lab/patanir/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons//analysis/00-manuscript/13-nygc-consortium/rdata_files/03-02-list_of_counts_stat_samplesheet_from_nygc.rds")
genic_count_list <- readRDS(file.path(analysis_save_dir, "01-gene-level-counts-of-postmortem-samples.rds"))
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
# source("scripts/src/functions.R")
# source("../src/rmats-analysis/src/source.R")
# ==========================================================================
# Analysis
# ==========================================================================
annotation_ranges <- import.gff2(annotation_gtf)
annotation_gene_ranges <- annotation_ranges[annotation_ranges$type == "gene"]


# feature_ranges_with_repeats <- significant_de_event_list$repeat_containing_ranges
counts_matrix <- genic_count_list$counts 
tpm_matrix <- counts_matrix / (width(annotation_gene_ranges[match(rownames(counts_matrix), annotation_gene_ranges$gene_id)]) / 1e3)
tpm_matrix <- t(t(tpm_matrix) / (colSums(tpm_matrix) / 1e6))

cpm_matrix <- t(t(counts_matrix) / (colSums(counts_matrix)/1e6))



coldata_df <- genic_count_list$coldata
coldata_df$stmn2_expr <- tpm_matrix[annotation_gene_ranges$gene_id[annotation_gene_ranges$gene_name == "STMN2" & !is.na(annotation_gene_ranges$gene_name)],
    coldata_df$column_id]
coldata_df$unc13a_expr <- tpm_matrix[annotation_gene_ranges$gene_id[annotation_gene_ranges$gene_name == "UNC13A" & !is.na(annotation_gene_ranges$gene_name)],
    coldata_df$column_id]
coldata_df$fus_expr <- tpm_matrix[annotation_gene_ranges$gene_id[annotation_gene_ranges$gene_name == "FUS" & !is.na(annotation_gene_ranges$gene_name)],
    coldata_df$column_id]
coldata_df$kcnq2_expr <- tpm_matrix[annotation_gene_ranges$gene_id[annotation_gene_ranges$gene_name == "KCNQ2" & !is.na(annotation_gene_ranges$gene_name)][1],
    coldata_df$column_id]
coldata_df$sort1_expr <- tpm_matrix[annotation_gene_ranges$gene_id[annotation_gene_ranges$gene_name == "SORT1" & !is.na(annotation_gene_ranges$gene_name)],
    coldata_df$column_id]
coldata_df$neat1_expr <- tpm_matrix[annotation_gene_ranges$gene_id[annotation_gene_ranges$gene_name == "NEAT1" & !is.na(annotation_gene_ranges$gene_name)],
    coldata_df$column_id]


feature_counts_matrix <- postmortem_list$counts
feature_counts_cpm_matrix <- t(t(feature_counts_matrix) / (colSums(feature_counts_matrix) / 1e6))
feature_counts_tpm_matrix <- feature_counts_matrix / (width(feature_ranges[match(rownames(feature_counts_matrix), feature_ranges$hash_id)]) / 1e3)
feature_counts_tpm_matrix <- t(t(feature_counts_tpm_matrix) / (colSums(feature_counts_tpm_matrix) / 1e6))
# feature_counts_cpm_matrix <- t(t(feature_counts_tpm_matrix) / (colSums(feature_counts_tpm_matrix) / 1e6))


unc13a_ir_ranges <- GRanges(
    seqnames = c("19", "19"),
    IRanges(start = c(17642414, 17642414),
        end = c(17642541, 17642591)),
    strand = c("-", "-"))
unc13a_ir_ranges$gene_name <- "UNC13A"
unc13a_ir_ranges$intron_number <- c("128bp", "178bp")


unc13a_ol_object <- findOverlaps(unc13a_ir_ranges, feature_ranges, type = "within")
current_feature_ranges <- feature_ranges[subjectHits(unc13a_ol_object)]
current_feature_ranges$intron_number <- unc13a_ir_ranges$intron_number[queryHits(unc13a_ol_object)]
current_feature_ranges$gene_name <- unc13a_ir_ranges$gene_name[queryHits(unc13a_ol_object)]
unc13a_exon_skipping_ranges <- unique(current_feature_ranges)


# current_feature_ranges <- feature_ranges[feature_ranges$gene_id == "ENSG00000104435"]
# current_feature_ranges$intron_number <- unc13a_ir_ranges$intron_number[queryHits(unc13a_ol_object)]
# current_feature_ranges$gene_name <- unc13a_ir_ranges$gene_name[queryHits(unc13a_ol_object)]
# unc13a_exon_skipping_ranges <- unique(current_feature_ranges)





kcnq2_ir_ranges <- GRanges(
    seqnames = "20",
    IRanges(start = 63442406, end = 63442531),
    strand = "-")
kcnq2_ir_ranges$gene_name <- "KCNQ2"
kcnq2_ir_ranges$intron_number <- "5"

kcnq2_ol_object <- findOverlaps(feature_ranges, kcnq2_ir_ranges, type = "within")
current_feature_ranges <- feature_ranges[queryHits(kcnq2_ol_object)]
current_feature_ranges$intron_number <- kcnq2_ir_ranges$intron_number[subjectHits(kcnq2_ol_object)]
current_feature_ranges$gene_name <- kcnq2_ir_ranges$gene_name[subjectHits(kcnq2_ol_object)]
kcnq_exon_skipping_ranges <- current_feature_ranges



sort1_ir_ranges <- GRanges(
    seqnames = "1",
    IRanges(start = 109315299, end = 109315397),
    strand = "-")
sort1_ir_ranges$gene_name <- "SORT1"
sort1_ir_ranges$intron_number <- "17b"

sort1_ol_object <- findOverlaps(sort1_ir_ranges, feature_ranges, type = "within")
current_feature_ranges <- feature_ranges[subjectHits(sort1_ol_object)]
current_feature_ranges$intron_number <- sort1_ir_ranges$intron_number[queryHits(sort1_ol_object)]
current_feature_ranges$gene_name <- sort1_ir_ranges$gene_name[queryHits(sort1_ol_object)]
sort1_exon_skipping_ranges <- current_feature_ranges


neat1_ir_ranges <- GRanges(
    seqnames = "11",
    IRanges(start = 65426533, end = 65445540),
    strand = "+")
neat1_ir_ranges$gene_name <- "NEAT1"
neat1_ir_ranges$intron_number <- "long"

neat1_ol_object <- findOverlaps(feature_ranges, neat1_ir_ranges, type = "within")
current_feature_ranges <- feature_ranges[queryHits(neat1_ol_object)]
current_feature_ranges$intron_number <- neat1_ir_ranges$intron_number[subjectHits(neat1_ol_object)]
current_feature_ranges$gene_name <- neat1_ir_ranges$gene_name[subjectHits(neat1_ol_object)]
neat1_exon_skipping_ranges <- current_feature_ranges


fus_ir_ranges <- GRanges(
    seqnames = c("16", "16"),
    IRanges(start = c(31185180, 31186837), end = c(31186801, 31188324)),
    strand = c("+", "+"))
fus_ir_ranges$gene_name <- "FUS"
fus_ir_ranges$intron_number <- c("6","7")

fus_ol_object <- findOverlaps(feature_ranges, fus_ir_ranges, type = "within")
current_feature_ranges <- feature_ranges[queryHits(fus_ol_object)]
current_feature_ranges$intron_number <- fus_ir_ranges$intron_number[subjectHits(fus_ol_object)]
current_feature_ranges$gene_name <- fus_ir_ranges$gene_name[subjectHits(fus_ol_object)]
fus_intron_ranges <- current_feature_ranges


all_intron_ranges <- c(unc13a_exon_skipping_ranges, kcnq_exon_skipping_ranges, sort1_exon_skipping_ranges, neat1_exon_skipping_ranges, fus_intron_ranges)


intron_ranges_split <- split(all_intron_ranges, paste(all_intron_ranges$gene_name, all_intron_ranges$intron_number, sep = "_"))


summarised_counts_matrix <- do.call(rbind, lapply(intron_ranges_split, function(x){
    all_ranges <- feature_ranges[feature_ranges$gene_id == unique(x$gene_id)]
    if(length(x) == 1){
        return(feature_counts_cpm_matrix[x$hash_id,])
    }
    return(colSums(feature_counts_cpm_matrix[x$hash_id,]))
}))

coldata_df$unc13a_128bp <- summarised_counts_matrix["UNC13A_128bp",coldata_df$column_id]
coldata_df$fus_6 <- summarised_counts_matrix["FUS_6",coldata_df$column_id]
coldata_df$fus_7 <- summarised_counts_matrix["FUS_7",coldata_df$column_id]
coldata_df$kcnq2_5 <- summarised_counts_matrix["KCNQ2_5",coldata_df$column_id]
coldata_df$neat1_long <- summarised_counts_matrix["NEAT1_long",coldata_df$column_id]
coldata_df$sort1_17b <- summarised_counts_matrix["SORT1_17b",coldata_df$column_id]



# width_split <- vapply(split(width(current_intron_ranges), paste(current_intron_ranges$gene_name, current_intron_ranges$intron_number, sep = "_")), )
summarised_counts_ratio_matrix <- do.call(rbind, lapply(intron_ranges_split, function(x){
    all_ranges <- feature_ranges[feature_ranges$gene_id == unique(x$gene_id)]
    all_genic_counts <- colSums(feature_counts_matrix[all_ranges$hash_id,]) 
    if(length(x) == 1){
        current_counts <- feature_counts_matrix[x$hash_id,]
    }else{
        current_counts <- colSums(feature_counts_matrix[x$hash_id,])
    }
    return(current_counts / all_genic_counts)
}))

coldata_df$unc13a_128bp_ratio <- summarised_counts_ratio_matrix["UNC13A_128bp",coldata_df$column_id]
coldata_df$fus_6_ratio <- summarised_counts_ratio_matrix["FUS_6",coldata_df$column_id]
coldata_df$fus_7_ratio <- summarised_counts_ratio_matrix["FUS_7",coldata_df$column_id]
coldata_df$kcnq2_5_ratio <- summarised_counts_ratio_matrix["KCNQ2_5",coldata_df$column_id]
coldata_df$neat1_long_ratio <- summarised_counts_ratio_matrix["NEAT1_long",coldata_df$column_id]
coldata_df$sort1_17b_ratio <- summarised_counts_ratio_matrix["SORT1_17b",coldata_df$column_id]



coldata_df_split <- split(coldata_df, paste(coldata_df$subject_id, coldata_df$tissue))
patient_plot_df <- do.call(rbind, lapply(coldata_df_split, function(temp_df){
    message(unique(temp_df$subject_id))
    return_df <- unique(temp_df[,!(colnames(temp_df) %in% 
        c("X", "sample_name", "alignment_file", "alignment_file_dir",
          "isPairedEnd", "strand_specificity", "samplesheet_path", 
          "column_id", "sample_id", "mutation", "site_of_motor_onset",
          "instrument", "rin", "stmn2_expr", "unc13a_expr", "fus_expr",
          "kcnq2_expr", "sort1_expr", "neat1_expr",
          "unc13a_128bp", "fus_6", "fus_7", "kcnq2_5", "neat1_long", "sort1_17b",
          "unc13a_128bp_ratio", "fus_6_ratio", "fus_7_ratio", "kcnq2_5_ratio", "neat1_long_ratio", "sort1_17b_ratio"))])
    if(nrow(return_df) > 1){
        # message(paste(unique(return_df$subject_id), unique(return_df$tissue)))
        if("N/A" %in% return_df$mnd_with_ftd & "Yes" %in% return_df$mnd_with_ftd){
            return_df <- return_df[return_df$mnd_with_ftd == "Yes",]
        }else if("N/A" %in% return_df$mnd_with_ftd & "No" %in% return_df$mnd_with_ftd){
            return_df <- return_df[return_df$mnd_with_ftd == "No",]
        }else{
            stop(paste(return_df$subject_id, return_df$tissue))
        }
        # message(nrow(return_df))
    }

    return_df$stmn2_expr <- mean(temp_df$stmn2_expr)
    return_df$unc13a_expr <- mean(temp_df$unc13a_expr)
    return_df$fus_expr <- mean(temp_df$fus_expr)
    return_df$kcnq2_expr <- mean(temp_df$kcnq2_expr)
    return_df$sort1_expr <- mean(temp_df$sort1_expr)
    return_df$neat1_expr <- mean(temp_df$neat1_expr)
    return_df$unc13a_128bp <- mean(temp_df$unc13a_128bp)
    return_df$fus_6 <- mean(temp_df$fus_6)
    return_df$fus_7 <- mean(temp_df$fus_7)
    return_df$kcnq2_5 <- mean(temp_df$kcnq2_5)
    return_df$neat1_long <- mean(temp_df$neat1_long)
    return_df$sort1_17b <- mean(temp_df$sort1_17b)
    return_df$unc13a_128bp_ratio <- mean(temp_df$unc13a_128bp_ratio)
    return_df$fus_6_ratio <- mean(temp_df$fus_6_ratio)
    return_df$fus_7_ratio <- mean(temp_df$fus_7_ratio)
    return_df$kcnq2_5_ratio <- mean(temp_df$kcnq2_5_ratio)
    return_df$neat1_long_ratio <- mean(temp_df$neat1_long_ratio)
    return_df$sort1_17b_ratio <- mean(temp_df$sort1_17b_ratio)
    return(return_df)
}))
# patient_plot_df$tissue <- droplevels(patient_plot_df$tissue)
saveRDS(object = patient_plot_df, file = file.path(analysis_save_dir, "02-summarised-counts-for-patients.rds"))