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


# Preprocessing: Load data ------------------------------------------------
srt_list <- list()
cellcalling_list <- list()
for (i in 1:length(samples)) {
  cat("++++++", samples[i], "(Preprocessing-LoadingData)", "++++++", "\n")
  cell_upset <- as.data.frame(readRDS(file = paste0(SCPwork_dir, "/", samples[i], "/Alignment-Cellranger/", samples[i], "/CellCalling/cell_upset.rds")))
  rownames(cell_upset) <- cell_upset[, "Barcode"]
  cell_upset[["sample"]] <- samples[i]
  cellcalling_list[[samples[i]]] <- cell_upset
  
  srt <- Read10X(data.dir = paste0(SCPwork_dir, "/", samples[i], "/Alignment-Cellranger/", samples[i], "/outs/raw_feature_bc_matrix/"))[,cell_upset[["Barcode"]]]
  srt_list[[i]] <- srt
}

# Preprocessing: Create Seurat object -------------------------------------

for (i in 1:length(samples)) {
  cat("++++++", samples[i], "(Preprocessing-CreateSeuratObject)", "++++++", "\n")
  srt <- CreateSeuratObject(counts = get(paste0(samples[i], "_10X")), project = samples[i])
  srt[["orig.ident"]] <- samples[i]
  srt[["percent.mt"]] <- PercentageFeatureSet(object = srt, pattern = "(^MT-)|(^Mt-)|(^mt-)")
  srt[["percent.ribo"]] <- PercentageFeatureSet(object = srt, pattern = "(^RP[SL]\\d+$)|(^Rp[sl]\\d+$)|(^rp[sl]\\d+$)")
  srt[["cellcalling_method"]] <- get(paste0(samples[i], "_cellcalling"))[Cells(srt), "Method_comb"]
  srt[["cellcalling_methodNum"]] <- get(paste0(samples[i], "_cellcalling"))[Cells(srt), "Method_num"]
  srt <- RenameCells(object = srt, add.cell.id = samples[i])
  srt_list[[samples[i]]] <- srt
  
  velocity <- as.Seurat(x = get(paste0(samples[i], "_velocity")), project = samples[i])
  velocity <- RenameCells(
    object = velocity,
    new.names = gsub(x = colnames(velocity), pattern = ".*:", replacement = "", perl = TRUE) %>%
      gsub(x = ., pattern = "x$", replacement = "-1", perl = TRUE) %>%
      paste0(samples[i], "_", .)
  )
  velocity[["orig.ident"]] <- samples[i]
  velocity_list[[samples[i]]] <- velocity
  
  assign(x = paste0(samples[i], "_10X"), value = NULL)
  assign(x = paste0(samples[i], "_velocity"), value = NULL)
}

# Preprocessing: Create Seurat object -------------------------------------
cat("++++++", sample, "(Preprocessing-CreateSeuratObject)", "++++++", "\n")
srt <- CreateSeuratObject(counts = get(paste0("10X_", sample)), project = sample)
srt[["orig.ident"]] <- sample
srt[["percent.mt"]] <- PercentageFeatureSet(object = srt, pattern = "(^MT-)|(^Mt-)|(^mt-)")
srt[["percent.ribo"]] <- PercentageFeatureSet(object = srt, pattern = "(^RP[SL]\\d+$)|(^Rp[sl]\\d+$)|(^rp[sl]\\d+$)")
srt[["cellcalling_method"]] <- get(paste0("Cellcalling_", sample))[Cells(srt), "Method_comb"]
srt[["cellcalling_methodNum"]] <- get(paste0("Cellcalling_", sample))[Cells(srt), "Method_num"]
srt <- RenameCells(object = srt, add.cell.id = sample)

# Preprocessing: Cell filtering -----------------------------------
ntotal <- ncol(srt)
sce <- as.SingleCellExperiment(srt)
sce <- scDblFinder(sce, verbose = FALSE)
srt[["scDblFinder_score"]] <- sce[["scDblFinder.score"]]
srt[["scDblFinder_class"]] <- sce[["scDblFinder.class"]]
db_out <- colnames(srt)[srt[["scDblFinder_class"]] == "doublet"]
# p1 <- ggplot(srt@meta.data, aes(x = scDblFinder_score, y = "orig.ident")) +
#   geom_point(position = position_jitter()) +
#   geom_boxplot(outlier.shape = NA) +
#   stat_summary(fun = median, geom = "point", size = 2, shape = 21, color = "black", fill = "white")
# p2 <- ggplot(srt@meta.data) +
#   geom_histogram(aes(x = scDblFinder_score))
# p1 %>% insert_top(p2)

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

metaCol <- c(
  "percent.mt", "percent.ribo", "cellcalling_methodNum",
  "log10_nCount_RNA", "log10_nFeature_RNA", "percent.top_20"
)
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

filters <- c(
  "log10_total_counts:higher:2.5",
  "log10_total_counts:lower:5",
  "log10_total_features:higher:2.5",
  "log10_total_features:lower:5",
  "pct_counts_in_top_20_features:both:5",
  "featcount_dist:both:5"
)
qc_out <- lapply(strsplit(filters, ":"), function(f) {
  colnames(srt)[isOutlier(get(f[1]),
                          log = FALSE,
                          nmads = as.numeric(f[3]),
                          type = f[2]
  )]
})
qc_out <- table(unlist(qc_out))
qc_out <- names(qc_out)[qc_out >= 2]

umi_out <- colnames(srt)[srt[["nCount_RNA", drop = TRUE]] <= UMI_threshold]
gene_out <- colnames(srt)[srt[["nFeature_RNA", drop = TRUE]] <= gene_threshold]

pct_counts_Mt <- srt[["percent.mt", drop = TRUE]]
if (all(pct_counts_Mt > 0 & pct_counts_Mt < 1)) {
  pct_counts_Mt <- pct_counts_Mt * 100
}
if (mito_threshold > 0 & mito_threshold < 1) {
  mito_threshold <- mito_threshold * 100
}
mt_out <- colnames(srt)[isOutlier(pct_counts_Mt, nmads = 3, type = "lower") |
                          (isOutlier(pct_counts_Mt, nmads = 2.5, type = "higher") & pct_counts_Mt > 10) |
                          (pct_counts_Mt > mito_threshold)]

total_out <- unique(c(db_out, qc_out, umi_out, gene_out, mt_out))
srt[["CellFilterng"]] <- rep(x = FALSE, ncol(srt))
srt[["CellFilterng"]][total_out, ] <- TRUE

cat(">>>", "Total cells:", ntotal, "\n")
cat(">>>", "Cells which are filtered out:", length(total_out), "\n")
cat("...", length(db_out), "potential doublets", "\n")
cat("...", length(qc_out), "unqualified cells", "\n")
cat("...", length(umi_out), "low-UMI cells", "\n")
cat("...", length(gene_out), "low-Gene cells", "\n")
cat("...", length(mt_out), "high-Mito cells", "\n")
cat(">>>", "Remained cells after filtering :", ntotal - length(total_out), "\n")

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

srt_list_filter <- lapply(setNames(samples, samples), function(sc_set) {
  cat("++++++", sc_set, "(Preprocessing-CellFiltering)", "++++++", "\n")
  srt <- srt_list_QC[[sc_set]]
  srt_filter <- subset(srt, CellFilterng == FALSE)
  return(srt_filter)
})



# Standard_SCP -----------------------------------
for (nm in normalization_method) {
  dir.create(paste0("Normalization-", nm), recursive = T, showWarnings = FALSE)
  if (file.exists(paste0("srt_list_filter_", nm, ".rds"))) {
    cat("Loading the", paste0("srt_list_filter_", nm), "from the file....\n")
    assign(
      x = paste0("srt_list_filter_", nm),
      value = readRDS(paste0("srt_list_filter_", nm, ".rds"))
    )
  } else {
    assign(
      x = paste0("srt_list_filter_", nm),
      value = lapply(setNames(samples, samples), function(sc_set) {
        cat("++++++", sc_set, paste0("Standard_SCP-", nm), "++++++", "\n")
        srt <- srt_list_filter[[sc_set]]
        srt <- Standard_SCP(
          srt = srt, normalization_method = nm, nHVF = nHVF,
          maxPC = maxPC, resolution = resolution, reduction = reduction,
          cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
          exogenous_genes = exogenous_genes
        )
        return(srt)
      })
    )
    
    invisible(lapply(samples, function(x) {
      saveRDS(get(paste0("srt_list_filter_", nm))[[x]], paste0("Normalization-", nm, "/", x, ".rds"))
    }))
    saveRDS(get(paste0("srt_list_filter_", nm)), file = paste0("srt_list_filter_", nm, ".rds"))
  }
}

future:::ClusterRegistry("stop")
