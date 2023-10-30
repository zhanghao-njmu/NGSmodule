#!/usr/bin/env Rscript

library(refGenome)
library(data.table)
library(dplyr)
library(SCP)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]
gtfFile <- args[2]
aligner <- args[3]
species <- args[4]
database <- args[5]

######### example 1 ############
# work_dir <- "/data/lab/HeXi/Rpl39l/RNC-seq/NGSmodule_work/"
# gtfFile <- "/data/database/iGenomes/mouse/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf"
# aligner <- "hisat2"
# species <- "Homo_sapiens"
# database <- "Ensembl"
##############################

if (aligner == "kallisto") {
  files <- list.files(work_dir, recursive = T, full.names = T) %>%
    grep(x = ., pattern = paste0("Alignment-kallisto/abundance.tsv$"), perl = T, value = T) %>%
    sort()
  df_list <- lapply(1:length(files), function(x) {
    df <- read.table(file = files[x], header = T, sep = "\t", stringsAsFactors = F, quote = "", fill = T, comment.char = "", check.names = F)
    df <- df[, c("target_id", "tpm")]
    colnames(df) <- c("TranscriptID", paste0(basename(dirname(dirname(files[x]))), ".kallisto.tpm"))
    return(df)
  })
  input <- Reduce(function(x, y) merge(x, y, by = 1, all = TRUE), df_list)
  write.table(x = input, file = paste("Quantification", ".", aligner, ".", "tpm", ".tab", sep = ""), sep = "\t", row.names = F)
} else {
  for (type in c("count", "rpkm", "fpkm", "tpm", "log2CPM")) {
    files <- list.files(work_dir, recursive = T, full.names = T) %>%
      grep(x = ., pattern = paste0("Quantification/.*", aligner, ".", type, "$"), perl = T, value = T) %>%
      sort()
    if (length(files) != 0) {
      df_list <- lapply(1:length(files), function(x) {
        read.table(file = files[x], header = T, sep = "\t", stringsAsFactors = F, quote = "", fill = T, comment.char = "", check.names = F)
      })
      input <- Reduce(function(x, y) merge(x, y, by = 1, all = TRUE), df_list)
      write.table(x = input, file = paste("Quantification", ".", aligner, ".", type, ".tab", sep = ""), sep = "\t", row.names = F)
    } else {
      cat("Warning: No .", type, " file in ", work_dir, "\n")
      next
    }
  }
}
