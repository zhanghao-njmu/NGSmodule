#!/usr/bin/env Rscript

library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]
gtf <- args[2]
aligner <- args[3]
species <- args[4]
database <- args[5]

######### example 1 ############
# work_dir <- "/storage/shihongjunLab/lichenxi/zhanghao/linjiayi/HFpEF/RNAseq_hongxu/NGSmodule_work/"
# gtf <- "/storage/shihongjunLab/lichenxi/zhanghao/reference/iGenomes/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf"
# aligner <- "hisat2"
# species <- "Mus_musculus"
# database <- "Ensembl"
##############################

merge_gtf_by = "gene_id"
columns = c(
  "seqname", "feature", "start", "end", "strand",
  "gene_id", "gene_name", "gene_type"
)

gtf_all <- suppressWarnings(fread(gtf, sep = "\t"))
gtf_all <- gtf_all[, 1:9]
colnames(gtf_all) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
for (type in c("gene", "transcript", "exon", "CDS")) {
  if (type %in% gtf_all[["feature"]]) {
    gtf_all <- gtf_all[gtf_all[["feature"]] == type, ]
    break
  }
}
columns1 <- intersect(colnames(gtf_all), columns)

gtf_attribute <- gtf_all[["attribute"]]
gtf_attribute <- gsub(pattern = "\"", replacement = "", x = gtf_attribute)
gtf_attribute <- strsplit(gtf_attribute, split = "; *")
gene_attr <- lapply(gtf_attribute, function(x) {
  detail <- strsplit(x, " ")
  out <- lapply(detail, function(x) x[2:length(x)])
  names(out) <- sapply(detail, function(x) x[1])
  out <- out[intersect(columns, names(out))]
  return(out)
})
gene_attr_df <- rbindlist(gene_attr, fill = TRUE)
gtf_columns <- cbind(gtf_all[, intersect(colnames(gtf_all), columns), with = FALSE], gene_attr_df)
colnames(gtf_columns) <- make.unique(colnames(gtf_columns))
gtf_columns_collapse <- aggregate(gtf_columns, by = list(rowid = gtf_columns[[merge_gtf_by]]), FUN = function(x) {
  paste0(unique(x), collapse = ";")
})
rownames(gtf_columns_collapse) <- gtf_columns_collapse[["rowid"]]
gtf_columns_collapse[["rowid"]] <- NULL


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
  output <- Reduce(function(x, y) merge(x, y, by = 1, all = TRUE), df_list)
  write.table(x = output, file = paste("Quantification", ".", aligner, ".", "tpm", ".tsv", sep = ""), sep = "\t", row.names = F)
} else {
  for (type in c("count", "rpkm", "fpkm", "tpm","cpm", "log2cpm")) {
    files <- list.files(work_dir, recursive = T, full.names = T) %>%
      grep(x = ., pattern = paste0("Quantification/.*", aligner, ".", type, "$"), perl = T, value = T) %>%
      sort()
    if (length(files) != 0) {
      df_list <- lapply(1:length(files), function(x) {
        read.table(file = files[x], header = T, sep = "\t", stringsAsFactors = F, quote = "", fill = T, comment.char = "", check.names = F)
      })
      output <- Reduce(function(x, y) merge(x, y, by = 1, all = TRUE), df_list)
      if (exists("gtf_columns_collapse")) {
        output <- cbind(output,"Annotation.GTF")
        output <- merge(x = output, by.x = "GeneID", y = gtf_columns_collapse, by.y = "row.names", all.x = TRUE)
      }
      write.table(x = output, file = paste("Quantification", ".", aligner, ".", type, ".tsv", sep = ""), sep = "\t", row.names = F)
    } else {
      cat("Warning: No .", type, " file in ", work_dir, "\n")
      next
    }
  }
}
