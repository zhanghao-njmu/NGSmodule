#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
SCP_path <- as.character(args[1])
SCPwork_dir <- as.character(args[2])
# ##### test #####
# setwd("/ssd/lab/ZhangHao/test/SCP/NGSmodule_SCP_analysis/CellQC")
# # parameters: global settings ---------------------------------------------
# SCP_path <- "/home/zhanghao/Program/NGS/UniversalTools/NGSmodule/SingleCellPipe/"
# SCPwork_dir <- "/ssd/lab/ZhangHao/test/SCP/NGSmodule_SCP_work/"
 
# Library -----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(invisible(lapply(
  c(
    "Seurat", "SeuratDisk", "SeuratWrappers", "sctransform", "intrinsicDimension", "scater", "Matrix", "BiocParallel",
    "future", "reticulate", "harmony", "liger", "simspec", "scMerge", "BiocSingular", "zinbwave", "plyr", "dplyr", "RColorBrewer", "scales", "gtools",
    "ggsci", "ggpubr", "ggplot2", "ggplotify", "aplot", "cowplot", "reshape2", "stringr", "scDblFinder",
    "velocyto.R", "biomaRt", "rvest", "xml2"
  ),
  require,
  character.only = TRUE
))))
set.seed(11)

source(paste0(SCP_path, "/SCP-workflow-funtcion.R"))
samples <- list.dirs(path = SCPwork_dir, recursive = FALSE, full.names = FALSE)


# Preprocessing: Create Seurat object ------------------------------------------------
cellcalling_list <- list()
srt_list <- list()
velocity_list <- list()
for (i in 1:length(samples)) {
  cat("++++++", samples[i], "(Preprocessing-LoadingData)", "++++++", "\n")
  cell_upset <- as.data.frame(readRDS(file = paste0(SCPwork_dir, "/", samples[i], "/Alignment-Cellranger/", samples[i], "/CellCalling/cell_upset.rds")))
  rownames(cell_upset) <- cell_upset[, "Barcode"]
  cell_upset[["sample"]] <- samples[i]
  cellcalling_list[[samples[i]]] <- cell_upset
  
  srt_matrix <- Read10X(data.dir = paste0(SCPwork_dir, "/", samples[i], "/Alignment-Cellranger/", samples[i], "/outs/raw_feature_bc_matrix/"))[,cell_upset[["Barcode"]]]
  srt <- CreateSeuratObject(counts = srt_matrix, project = samples[i])
  srt[["percent.mt"]] <- PercentageFeatureSet(object = srt, pattern = "(^MT-)|(^Mt-)|(^mt-)")
  srt[["percent.ribo"]] <- PercentageFeatureSet(object = srt, pattern = "(^RP[SL]\\d+$)|(^Rp[sl]\\d+$)|(^rp[sl]\\d+$)")
  srt[["cellcalling_method"]] <- get(paste0(samples[i], "_cellcalling"))[Cells(srt), "Method_comb"]
  srt[["cellcalling_methodNum"]] <- get(paste0(samples[i], "_cellcalling"))[Cells(srt), "Method_num"]
  srt <- RenameCells(object = srt, add.cell.id = samples[i])
  srt_list[[samples[i]]] <- srt
  
  velocity <- ReadVelocity(file = paste0(SCPwork_dir, "/", samples[i], "/Alignment-Cellranger/", samples[i], "/velocyto/", samples[i], ".loom"))
  velocity <- as.Seurat(x = velocity, project = samples[i])
  velocity <- RenameCells(
    object = velocity,
    new.names = gsub(x = colnames(velocity), pattern = ".*:", replacement = "", perl = TRUE) %>%
      gsub(x = ., pattern = "x$", replacement = "-1", perl = TRUE) %>%
      paste0(samples[i], "_", .)
  )
  velocity[["orig.ident"]] <- samples[i]
  velocity_list[[samples[i]]] <- velocity
}


# Preprocessing: Cell QC -----------------------------------
for (i in srt_list) {
  srt <- srt_list[[i]]
  ntotal <- ncol(srt)
  sce <- as.SingleCellExperiment(srt)
  sce <- scDblFinder(sce, verbose = FALSE)
  srt[["scDblFinder_score"]] <- sce[["scDblFinder.score"]]
  srt[["scDblFinder_class"]] <- sce[["scDblFinder.class"]]
  db_out <- colnames(srt)[srt[["scDblFinder_class"]] == "doublet"]
  
  sce <- addPerCellQC(sce, percent_top = c(20))
  srt[["percent.top_20"]] <- pct_counts_in_top_20_features <- colData(sce)$percent.top_20
  log10_total_counts <- log10(srt[["nCount_RNA", drop = TRUE]])
  log10_total_features <- log10(srt[["nFeature_RNA", drop = TRUE]])
  srt[["log10_nCount_RNA"]] <- log10_total_counts
  srt[["log10_nFeature_RNA"]] <- log10_total_features
  
  mod <- loess(log10_total_features ~ log10_total_counts)
  pred <- predict(mod, newdata = data.frame(log10_total_counts = log10_total_counts))
  featcount_dist <- log10_total_features - pred
  srt[["featcount_dist"]] <- featcount_dist

  srt_list[[i]] <- srt
}



# p1 <- ggplot(srt@meta.data, aes(x = scDblFinder_score, y = "orig.ident")) +
#   geom_point(position = position_jitter()) +
#   geom_boxplot(outlier.shape = NA) +
#   stat_summary(fun = median, geom = "point", size = 2, shape = 21, color = "black", fill = "white")
# p2 <- ggplot(srt@meta.data) +
#   geom_histogram(aes(x = scDblFinder_score))
# p1 %>% insert_top(p2)


# p1 <- srt@meta.data %>%
#   group_by(cellcalling_methodNum) %>%
#   summarise(
#     cellcalling_methodNum = cellcalling_methodNum,
#     count = n()
#   ) %>%
#   distinct() %>%
#   ggplot(aes(x = cellcalling_methodNum, y = count)) +
#   geom_col(aes(fill = count), color = "black") +
#   scale_fill_material(palette = "blue-grey") +
#   geom_text(aes(label = count), vjust = -1) +
#   scale_y_continuous(expand = expansion(0.08)) +
#   theme_classic() +
#   theme(aspect.ratio = 0.8, panel.grid.major = element_line())
#
# df <- data.frame(x = log10_total_counts, y = log10_total_features, featcount_dist = featcount_dist)
# p1 <- ggplot(df, aes(x = x, y = y)) +
#   stat_density2d(geom = "tile", aes(fill = ..density..^0.25, alpha = 1), contour = FALSE, show.legend = F) +
#   stat_density2d(geom = "tile", aes(fill = ..density..^0.25, alpha = ifelse(..density..^0.25 < 0.4, 0, 1)), contour = FALSE, show.legend = F) +
#   geom_point(aes(colour = featcount_dist)) +
#   geom_smooth(method = "loess", color = "black") +
#   scale_fill_gradientn(colours = colorRampPalette(c("white", "black"))(256)) +
#   scale_color_material(palette = "blue-grey", reverse = T, guide = FALSE) +
#   labs(x = "log10_nCount_RNA", y = "log10_nFeature_RNA") +
#   theme_classic()
# p2 <- ggplot(df, aes(x = x)) +
#   geom_histogram(fill = "steelblue") +
#   theme_void()
# p3 <- ggplot(df, aes(y = y)) +
#   geom_histogram(fill = "firebrick") +
#   theme_void()
# p <- p1 %>%
#   insert_top(p2, height = .3) %>%
#   insert_right(p3, width = .1)
# p <- as.ggplot(aplotGrob(p)) +
#   labs(title = "nCount vs nFeature") +
#   theme(aspect.ratio = 1)

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
