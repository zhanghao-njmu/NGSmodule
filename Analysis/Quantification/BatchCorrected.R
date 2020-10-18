#!/usr/bin/env Rscript
suppressWarnings(suppressPackageStartupMessages(invisible(lapply(
  c(
    "dplyr", "stringr", "ggplot2", "ggsci", "ggtree", "RColorBrewer", "cowplot",
    "aplot", "ggplotify", "edgeR", "sva", "limma", "patchwork"
  ),
  require,
  character.only = TRUE
))))

args <- commandArgs(trailingOnly = TRUE)
maindir <- args[1]
aligner <- args[2]
SampleInfoFile <- args[3]
script_path <- as.character(args[4])

# maindir <- "/data/database/SRR_collection/human/germ_cell_specification//"
# aligner <- "hisat2"
# SampleInfoFile <- "/data/database/SRR_collection/human/germ_cell_specification/temp_20200803161047.Sample_info.csv"
# script_path <- "/home/zhanghao/Program/NGS/UniversalTools/NGSmodule/Analysis/Quantification/BatchCorrected.R"

library(dplyr)
library(stringr)
library(ggplot2)
library(ggsci)
library(ggtree)
library(RColorBrewer)
library(cowplot)
library(aplot)
library(ggplotify)
library(edgeR)
library(sva)
library(limma)
set.seed(11)

script_dir <- gsub(x = script_path, pattern = "BatchCorrected.R", replacement = "")
source(paste0(script_dir, "/BatchCorrected_function.R"))

count_file <- paste0(maindir, "/NGSmodule_analysis/Quantification/Quantification.", aligner, ".count.tab")

if (!file.exists(count_file)) {
  print(paste0("Can not find count file: ", count_file, "\nPlease run NGSmodule Quantification first!\n"))
  quit(status = 1)
}
if (!file.exists(SampleInfoFile)) {
  print(paste0("Can not find SampleInfoFile: ", SampleInfoFile, "\nPlease check the config parameters!\n"))
  quit(status = 1)
}

##### Load data #####
sample_info <- read.csv(SampleInfoFile, stringsAsFactors = F, header = T)
sample_info <- as.data.frame(sapply(sample_info, as.character))
sample_info <- sample_info %>%
  group_by(SampleID) %>%
  mutate(RunID = paste0(RunID, collapse = ",")) %>%
  distinct() %>%
  as.data.frame()
rownames(sample_info) <- sample_info[, "SampleID"]
# sample_info[, "Group"] <- factor(sample_info[, "Group"],
#   levels = c(
#     "GV_oocyte", "zygote", "2-cell", "4-cell", "8-cell",
#     "morula", "early_blastocyst", "middle_blastocyst", "late_blastocyst", "ICM", "TE", "ESC"
#   )
# )
sample_info[, "Group"] <- factor(sample_info[, "Group"], unique(sample_info[, "Group"]))
sample_info <- sample_info[!is.na(sample_info[["Group"]]), ]

count_matrix_raw <- read.table(file = count_file, header = T, sep = "\t", row.names = 1, stringsAsFactors = F, check.names = F)
count_matrix <- count_matrix_raw[, which(str_detect(colnames(count_matrix_raw), pattern = paste0(".", aligner, ".count"))), drop = FALSE]
colnames(count_matrix) <- gsub(x = colnames(count_matrix), pattern = paste0(".", aligner, ".count"), replacement = "")
count_matrix <- count_matrix[, sample_info[["SampleID"]]]
count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]

logcpm <- cpm(count_matrix, log = TRUE, prior.count = 2)
logcpm_adj_raw <- logcpm

start <- which(str_detect(colnames(count_matrix_raw), ">>"))[1]
annotation_matrix <- count_matrix_raw[, start:ncol(count_matrix_raw)]


# select covariates for Combat and removeBatchEffect ----------------------
covar <- c()
for (col in c("Layout.PE.SE.")) {
  col_batch <- unique(sample_info[, c("BatchID", col)])
  if (min(table(col_batch[, col])) != 1 | length(table(col_batch[, col])) >= 3) {
    covar <- c(covar, col)
  }
}
# if (length(covar) >= 1) {
#   covar_only_fmla <- as.formula(paste("~", paste(covar, collapse = "+")))
#   covar_only <- model.matrix(covar_only_fmla, data = sample_info)
#   covar_mod_fmla <- as.formula(paste("~", paste(covar, collapse = "+"), "+ Group"))
#   covar_mod <- model.matrix(covar_mod_fmla, data = sample_info)
# } else {
covar_only <- NULL
covar_mod <- model.matrix(~Group, data = sample_info)
# }

# create design model -----------------------------------------------------
mod1 <- model.matrix(~ as.factor(Group), data = sample_info)
mod0 <- model.matrix(~1, data = sample_info)
n.sv <- num.sv(dat = logcpm, mod = mod1, method = "leek", seed = 11)

# apply raw ComBat funtion ------------------------------------------------
cat("\n>>> Apply raw ComBat funtion\n")
logcpm_adj_rawComBat <- ComBat(
  dat = logcpm,
  batch = as.factor(sample_info[["BatchID"]]),
  mod = covar_mod,
  prior.plots = FALSE,
  BPPARAM = bpparam("MulticoreParam")
)

# apply ComBat-seq funtion ------------------------------------------------
cat("\n>>> Apply ComBat-seq funtion\n")
count_adj <- ComBat_seq_custom(
  counts = count_matrix,
  batch = as.factor(sample_info[["BatchID"]]),
  covar_mod = covar_mod
)
logcpm_adj_ComBatSeq <- cpm(count_adj, log = TRUE, prior.count = 2)

# apply SVA function ------------------------------------------------------
cat("\n>>> Apply SVA funtion\n")
svobj <- sva(dat = logcpm, mod = mod1, mod0 = mod0, n.sv = n.sv)
logcpm_adj_SVA <- cleanY(y = logcpm, mod = mod1, svs = svobj$sv)

# apply removeBatchEffect function ----------------------------------------
cat("\n>>> Apply removeBatchEffect funtion\n")
logcpm_adj_removeBatchEffect <- removeBatchEffect(
  x = logcpm,
  batch = as.factor(sample_info[["BatchID"]]),
  covariates = covar_only,
  design = mod1
)


###################
methods <- c("raw", "rawComBat", "ComBatSeq", "SVA", "removeBatchEffect")
pl <- lapply(setNames(methods, methods), function(method) {
  cat("+++", method, "+++\n")
  logcpm_adj <- get(paste0("logcpm_adj_", method))
  write.table(
    x = cbind(data.frame(GeneID = rownames(logcpm_adj)), logcpm_adj),
    file = paste0("BatchCorrected_", method, ".log2CPM.tab"),
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
  )
  logcpm_adj_scale <- t(scale(t(logcpm_adj)))

  group_num <- length(unique(sample_info[, "Group"]))
  if (group_num <= 5) {
    col_color <- nord("mountain_forms")[1:group_num]
  }
  if (group_num > 5 & group_num <= 8) {
    col_color <- pal_nejm("default")(group_num)
  }
  if (group_num > 8 & group_num <= 20) {
    col_color <- pal_d3("category20")(group_num)
  }
  if (group_num > 20 & group_num <= 51) {
    col_color <- pal_igv("default")(group_num)
  }
  names(col_color) <- unique(sample_info[, "Group"])
  batch_color <- colorRampPalette(brewer.pal(n = 11, name = "Spectral")[1:11])(length(unique(sample_info[["BatchID"]])))
  names(batch_color) <- unique(sample_info[["BatchID"]])

  plot_list <- list()
  ##### Dendrogram #####
  cat(">>> Hierarchical Clustering (colored by Groups)\n")
  hc <- hclust(dist(t(logcpm_adj_scale), method = "euclidean")) ## sqrt(sum((x_i - y_i)^2)). Large expressed genes -large vars.
  p1 <- ggtree(tr = hc) + layout_dendrogram()
  p2 <- ggplot(sample_info) +
    geom_tile(aes(x = SampleID, y = 1, fill = Group), color = ifelse(nrow(sample_info) <= 30, "black", "transparent")) +
    scale_fill_manual(values = col_color, name = "Group") +
    theme_void()
  p <- p2 %>% insert_top(p1, height = 8)
  plot_list[["Hierarchical_Clustering_colored_by_group"]] <- as.ggplot(aplotGrob(p))

  cat(">>> Hierarchical Clustering (colored by BatchID)\n")
  p2 <- ggplot(sample_info) +
    geom_tile(aes(x = SampleID, y = 1, fill = BatchID), color = ifelse(nrow(sample_info) <= 30, "black", "transparent")) +
    scale_fill_manual(values = batch_color, name = "Group") +
    theme_void()
  p <- p2 %>% insert_top(p1, height = 8)
  plot_list[["Hierarchical_Clustering_colored_by_batch"]] <- as.ggplot(aplotGrob(p))

  ##### PCA #####
  cat(">>> Principal Components Analysis\n")
  df_pca <- prcomp(t(logcpm_adj_scale), center = F, scale. = F)
  df_pca <- summary(df_pca)
  PoV <- round(df_pca$importance[2, ] * 100, 2)
  x <- df_pca$x[, "PC1"]
  y <- df_pca$x[, "PC2"]
  m <- max(c(abs(x), abs(y)))

  df <- data.frame(
    x = x, y = y, sample = names(x),
    Group = sample_info[names(x), "Group"],
    Batch = sample_info[names(x), "BatchID"],
    stringsAsFactors = FALSE
  )

  p <- ggplot(data = df, aes(x = x, y = y, Group = Group, fill = Group, label = sample)) +
    geom_point(shape = 21, alpha = 0.8) +
    geom_rug(aes(color = Group), show.legend = FALSE) +
    labs(title = "Principal Components Analysis", x = paste0("PC1(", PoV[1], "%)"), y = paste0("PC2(", PoV[2], "%)")) +
    scale_fill_manual(values = col_color) +
    scale_color_manual(values = col_color) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    xlim(-m, m) +
    ylim(-m, m) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      panel.grid.major = element_line()
    )
  if (nrow(sample_info) <= 30) {
    p <- p + geom_label_repel(
      size = 2.5, color = "white",
      min.segment.length = 0, segment.color = "black", segment.alpha = 0.8
    )
  }
  plot_list[["PCA_colored_by_group"]] <- p

  p <- ggplot(data = df, aes(x = x, y = y, Batch = Batch, fill = Batch, label = sample)) +
    geom_point(shape = 21, alpha = 0.8) +
    geom_rug(aes(color = Batch), show.legend = FALSE) +
    labs(title = "Principal Components Analysis", x = paste0("PC1(", PoV[1], "%)"), y = paste0("PC2(", PoV[2], "%)")) +
    scale_fill_manual(values = batch_color) +
    scale_color_manual(values = batch_color) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    xlim(-m, m) +
    ylim(-m, m) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      panel.grid.major = element_line()
    )
  if (nrow(sample_info) <= 30) {
    p <- p + geom_label_repel(
      size = 2.5, color = "white",
      min.segment.length = 0, segment.color = "black", segment.alpha = 0.8
    )
  }
  plot_list[["PCA_colored_by_batch"]] <- p

  return(cowplot::plot_grid(plotlist = plot_list))
})

pdf("BatchCorrected.pdf", width = 12, height = 8)
invisible(lapply(pl, print))
invisible(dev.off())


##### check whether the unwanted file exists and remove it #####
if (file.exists("Rplots.pdf")) {
  invisible(file.remove("Rplots.pdf"))
}
