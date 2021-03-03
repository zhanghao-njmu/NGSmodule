#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
SCP_path <- as.character(args[1])
SCPwork_dir <- as.character(args[2])
# ##### test #####
# setwd("/ssd/lab/ZhangHao/test/SCP/NGSmodule_SCP_analysis/CellQC")
# # parameters: global settings ---------------------------------------------
# SCP_path <- "/home/zhanghao/Program/NGS/UniversalTools/NGSmodule/SingleCellPipe/"
# SCPwork_dir <- "/ssd/lab/ZhangHao/test/SCP/NGSmodule_SCP_work/"

# ##### true data #####
# setwd("/ssd/lab/HuangMingQian/scRNAseq/iPSC-ESC-iMeLC-PGCd6-CellLine-Testis/NGSmodule_SCP_analysis/CellQC/")
# # parameters: global settings ---------------------------------------------
# SCP_path <- "/home/zhanghao/Program/NGS/UniversalTools/NGSmodule/SingleCellPipe/"
# SCPwork_dir <- "/ssd/lab/HuangMingQian/scRNAseq/iPSC-ESC-iMeLC-PGCd6-CellLine-Testis/NGSmodule_SCP_work/"

# Library -----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(invisible(lapply(
  c(
    "Seurat", "SeuratDisk", "SeuratWrappers", "dplyr", "ggplot2", "ggplotify", "aplot", "cowplot",
    "reshape2", "stringr", "scuttle"
  ),
  require,
  character.only = TRUE
))))
set.seed(11)

source(paste0(SCP_path, "/SCP-workflow-funtcion.R"))
samples <- list.dirs(path = SCPwork_dir, recursive = FALSE, full.names = FALSE)


# Preprocessing: Create Seurat object ------------------------------------------------
srtList <- list()
for (i in 1:length(samples)) {
  cat("++++++", samples[i], "++++++", "\n")
  cat("Loading single cell expression data...\n")
  cell_upset <- as.data.frame(readRDS(file = paste0(SCPwork_dir, "/", samples[i], "/Alignment-Cellranger/", samples[i], "/CellCalling/cell_upset.rds")))
  rownames(cell_upset) <- cell_upset[, "Barcode"]
  cell_upset[["sample"]] <- samples[i]

  srt_matrix <- Read10X(data.dir = paste0(SCPwork_dir, "/", samples[i], "/Alignment-Cellranger/", samples[i], "/outs/raw_feature_bc_matrix/"))[, cell_upset[["Barcode"]]]
  srt <- CreateSeuratObject(counts = srt_matrix, project = samples[i])
  srt[["percent.mt"]] <- PercentageFeatureSet(object = srt, pattern = "(^MT-)|(^Mt-)|(^mt-)")
  srt[["percent.ribo"]] <- PercentageFeatureSet(object = srt, pattern = "(^RP[SL]\\d+(\\w|)$)|(^Rp[sl]\\d+(\\w|)$)|(^rp[sl]\\d+(\\w|)$)")
  srt[["cellcalling_method"]] <- cell_upset[Cells(srt), "Method_comb"]
  srt[["cellcalling_methodNum"]] <- cell_upset[Cells(srt), "Method_num"]
  srt <- RenameCells(object = srt, add.cell.id = samples[i])

  sce <- as.SingleCellExperiment(srt)
  sce <- addPerCellQC(sce, percent_top = c(20))
  srt[["percent.top_20"]] <- colData(sce)$percent.top_20
  srt[["log10_nCount_RNA"]] <- log10(srt[["nCount_RNA", drop = TRUE]])
  srt[["log10_nFeature_RNA"]] <- log10(srt[["nFeature_RNA", drop = TRUE]])

  cat("Loading single cell RNA-velocity data...\n")
  velocity <- ReadVelocity(file = paste0(SCPwork_dir, "/", samples[i], "/Alignment-Cellranger/", samples[i], "/velocyto/", samples[i], ".loom"), verbose = FALSE)
  velocity <- as.Seurat(x = velocity, project = samples[i], verbose = FALSE)
  velocity <- RenameCells(
    object = velocity,
    new.names = gsub(x = colnames(velocity), pattern = ".*:", replacement = "", perl = TRUE) %>%
      gsub(x = ., pattern = "x$", replacement = "-1", perl = TRUE) %>%
      paste0(samples[i], "_", .)
  )
  if (!identical(dim(srt), dim(velocity))) {
    stop(paste0(samples[i], ": RNA velocity have different dimensions with RNA count matrix."))
  }

  cat("Combine expression, RNA-velocity, and other annotation to an Seurat object...\n")
  srt@assays$spliced <- velocity@assays$spliced
  srt@assays$unspliced <- velocity@assays$unspliced
  srt@assays$ambiguous <- velocity@assays$ambiguous
  srt$nCount_spliced <- velocity$nCount_spliced
  srt$nFeature_spliced <- velocity$nFeature_spliced
  srt$nCount_unspliced <- velocity$nCount_unspliced
  srt$nFeature_unspliced <- velocity$nFeature_unspliced
  srt$nCount_ambiguous <- velocity$nCount_ambiguous
  srt$nFeature_ambiguous <- velocity$nFeature_ambiguous
  
  cat("Save the Seurat object to h5Seurat...\n")
  SaveH5Seurat(
    object = srt,
    filename = paste0(samples[i], ".h5Seurat"),
    overwrite = TRUE,
    verbose = FALSE
  )

  srtList[[samples[i]]] <- srt
}


# if (length(srtList) >= 2) {
#   srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
# } else {
#   srtMerge <- srtList[[1]]
# }
#
# meta <- srtMerge@meta.data
#
# p <- ggplot(meta, aes(x = orig.ident, y = scDblFinder_score, fill = orig.ident)) +
#   geom_point(position = position_jitter(), alpha = 0.5) +
#   geom_boxplot(outlier.shape = NA) +
#   stat_summary(fun = median, geom = "point", size = 4, shape = 21, color = "black", fill = "white") +
#   scale_fill_igv() +
#   theme_classic() +
#   theme(aspect.ratio = 0.8, panel.grid.major = element_line())
#
# p1 <- meta %>%
#   group_by(orig.ident, cellcalling_methodNum) %>%
#   summarise(
#     orig.ident = orig.ident,
#     cellcalling_methodNum = cellcalling_methodNum,
#     count = n()
#   ) %>%
#   distinct() %>%
#   ggplot(aes(x = cellcalling_methodNum, y = count)) +
#   geom_col(aes(fill = orig.ident), color = "black", position = position_dodge2(width = 1)) +
#   scale_fill_igv() +
#   geom_text(aes(label = count, group = orig.ident), vjust = -1, position = position_dodge2(width = 1)) +
#   scale_y_continuous(expand = expansion(0.05)) +
#   theme_classic() +
#   theme(aspect.ratio = 0.8, panel.grid.major = element_line())
#
# p1 <- ggplot(meta, aes(x = log10_nCount_RNA, y = log10_nFeature_RNA)) +
#   geom_point(colour = "steelblue") +
#   geom_smooth(method = "loess", color = "black") +
#   labs(x = "log10_nCount_RNA", y = "log10_nFeature_RNA") +
#   theme_classic()
# p2 <- ggplot(meta, aes(x = log10_nCount_RNA)) +
#   geom_histogram(fill = "steelblue", color = "black") +
#   theme_void()
# p3 <- ggplot(meta, aes(y = log10_nFeature_RNA)) +
#   geom_histogram(fill = "steelblue", color = "black") +
#   theme_void()
# p <- p1 %>%
#   insert_top(p2, height = .3) %>%
#   insert_right(p3, width = .1)
# p <- as.ggplot(aplotGrob(p)) +
#   labs(title = "nCount vs nFeature") +
#   theme(aspect.ratio = 0.8)

# df <- data.frame(cell = rownames(srt@meta.data))
# df[db_out,"db_out"] <- "db_out"
# df[qc_out,"qc_out"] <- "qc_out"
# df[umi_out,"umi_out"] <- "umi_out"
# df[gene_out,"gene_out"] <- "gene_out"
# df[mt_out,"mt_out"] <- "mt_out"
# index_upset <- reshape2::melt(df,
#   measure.vars = c("db_out", "qc_out", "umi_out", "gene_out", "mt_out"),
#   variable.name = "index",
#   value.name = "value"
# ) %>%
#   dplyr::filter(!is.na(value)) %>%
#   group_by(cell) %>%
#   summarize(
#     index_list = list(value),
#     index_comb = paste(value, collapse = ","),
#     index_num = n()
#   )
# y_max <- max(table(pull(index_upset, "index_comb")))
#
# p <- index_upset %>%
#   group_by(index_comb) %>%
#   mutate(count = n()) %>%
#   ungroup() %>%
#   filter(count %in% head(sort(unique(count), decreasing = T), 10)) %>%
#   ggplot(aes(x = index_list)) +
#   geom_bar(aes(fill = ..count..), color = "black", width = 0.5) +
#   geom_text(aes(label = ..count..), stat = "count", vjust = -0.5, hjust = 0, angle = 45) +
#   labs(title = "Cell intersection among differnent methods", x = "", y = "Cell number") +
#   scale_x_upset(sets = c("db_out", "qc_out", "umi_out", "gene_out", "mt_out")) +
#   scale_y_continuous(limits = c(0, 1.2 * y_max)) +
#   scale_fill_material(name = "Count", palette = "blue-grey") +
#   theme_combmatrix(
#     combmatrix.label.text = element_text(size = 10, color = "black"),
#     combmatrix.label.extra_spacing = 6
#   ) +
#   theme_classic() +
#   theme(
#     aspect.ratio = 0.6,
#     legend.position = "none"
#  )
