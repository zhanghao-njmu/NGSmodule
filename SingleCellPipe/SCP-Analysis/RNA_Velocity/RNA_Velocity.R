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
sc <- import("scanpy")
scv <- import("scvelo")
adata <- scv$datasets$pancreas()
adata
# scv$pl$scatter(adata, legend_loc='lower left', size=60)
## get embedding
emb <- adata$obsm["X_umap"]
clusters <- adata$obs$clusters
rownames(emb) <- names(clusters) <- adata$obs_names$values

## get clusters, convert to colors
col <- rainbow(length(levels(clusters)), s = 0.8, v = 0.8)
cell.cols <- col[clusters]
names(cell.cols) <- names(clusters)
plot(emb,
  col = cell.cols, pch = 16,
  xlab = "UMAP X", ylab = "UMAP Y"
)
legend(
  x = -13, y = 0,
  legend = levels(clusters),
  col = col,
  pch = 16
)
scv$pp$filter_genes(adata) ## filter
scv$pp$moments(adata) ## normalize and compute moments
scv$tl$recover_dynamics(adata) ## model

# adata$write('pancreas.h5ad', compression='gzip')
# adata = scv$read('pancreas.h5ad')

scv$tl$velocity(adata, mode = "dynamical")
scv$tl$velocity_graph(adata)
scv$pl$velocity_embedding_stream(adata, basis = "umap")
 

## top dynamic genes
topgenes <- adata$var["fit_likelihood"]
topgenes_vals <- topgenes[, 1]
names(topgenes_vals) <- rownames(topgenes)
topgenes_vals <- sort(topgenes_vals, decreasing = TRUE)
head(topgenes_vals)


scv$pl$scatter(adata, basis = names(topgenes_vals)[1:5], ncols = 5, frameon = FALSE)










srt_logCPM <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/iPSC-ESC-iMeLC-PGCd6-CellLine-Testis/NGSmodule_SCP_analysis/Integration/Normalization-logCPM/separate_HVF/Testis-d10,Testis-d20,Testis-d30,Testis-d40,Testis-d50,Testis-d70.rds")
velocity_list <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/iPSC-ESC-iMeLC-PGCd6-CellLine-Testis/NGSmodule_SCP_analysis/Integration/velocity_list.rds")
ldat <- Reduce(function(x, y) merge(x, y), velocity_list)
ldat2 <- subset(ldat, cell = colnames(srt_logCPM))
srt_logCPM[["spliced"]] <- ldat[["spliced"]]
srt_logCPM[["unspliced"]] <- ldat[["unspliced"]]

for (i in 1:length(samples)) {
  cat("++++++", samples[i], "(Preprocessing-LoadingData)", "++++++", "\n")
  cell_upset <- as.data.frame(readRDS(file = paste0(SCPwork_dir, "/", samples[i], "/Alignment/Cellranger/", samples[i], "/CellCalling/cell_upset.rds")))
  rownames(cell_upset) <- cell_upset[, "Barcode"]
  cells <- cell_upset %>%
    filter(Method_num >= cell_calling_methodNum) %>%
    pull("Barcode")
  cat(length(cells), "cells with calling methods >=", cell_calling_methodNum, "\n", sep = " ")
  assign(
    x = paste0(samples[i], "_cellcalling"),
    value = cell_upset
  )

  assign(
    x = paste0(samples[i], "_10X"),
    value = Read10X(data.dir = paste0(SCPwork_dir, "/", samples[i], "/Alignment/Cellranger/", samples[i], "/outs/raw_feature_bc_matrix/"))[, cells]
  )
  assign(
    x = paste0(samples[i], "_velocity"),
    value = ReadVelocity(file = paste0(SCPwork_dir, "/", samples[i], "/Alignment/Cellranger/", samples[i], "/velocyto/", samples[i], ".loom"))
  )
}





# Preprocessing: Create Seurat object -------------------------------------
srt_list <- list()
velocity_list <- list()
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




# !/usr/bin/env Rscript
setwd("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/NGSmodule_SCP_analysis/Prepare/")
# parameters: global settings ---------------------------------------------
SCP_path <- "/home/zhanghao/Program/NGS/UniversalTools/NGSmodule/SingleCellPipe/"
SCPwork_dir <- "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/NGSmodule_SCP_work/"

# Library -----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(invisible(lapply(
  c(
    "Seurat", "SeuratDisk", "SeuratWrappers", "dplyr", "ggplot2", "ggplotify", "aplot", "cowplot",
    "reshape2", "stringr", "scuttle", "RColorBrewer", "velocyto.R"
  ),
  require,
  character.only = TRUE
))))
set.seed(11)

source(paste0(SCP_path, "/SCP-workflow-funtcion.R"))
samples <- list.dirs(path = SCPwork_dir, recursive = FALSE, full.names = FALSE)

samples <- c(
  "22_ESC", "22_PGCLC_6d", "22_Cellline",
  "22_Cellline_injection_10d",
  "22_Cellline_injection_20d",
  "22_Cellline_injection_30d",
  "22_Cellline_injection_40d",
  "22_Cellline_injection_50d",
  "22_Cellline_injection_70d",
  "22_Cellline_injection_90d",
  "22_Cellline_injection_120d"
)
color_sample <- brewer.pal(11, "Paired")[c(1:4, 7:8, 5:6, 9:11)]
names(color_sample) <- samples
theme_zh <- function() {
  theme_classic() + theme(
    text = element_text(color = "black"),
    title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    strip.text = element_text(color = "black"),
    strip.text.x = element_text(color = "black"),
    strip.text.y = element_text(color = "black")
  )
}

velocityList <- list()
for (i in 1:length(samples)) {
  cat("Loading", samples[i], "RNA-velocity data...\n")
  velocity <- ReadVelocity(file = paste0(SCPwork_dir, "/", samples[i], "/Alignment-Cellranger/", samples[i], "/velocyto/", samples[i], ".loom"), verbose = FALSE)
  velocity <- as.Seurat(x = velocity, project = samples[i], verbose = FALSE)
  velocity <- RenameCells(
    object = velocity,
    new.names = gsub(x = colnames(velocity), pattern = ".*:", replacement = "", perl = TRUE) %>%
      gsub(x = ., pattern = "x$", replacement = "-1", perl = TRUE) %>%
      paste0(samples[i], "_", .)
  )

  cat("Combine expression, RNA-velocity, and other annotation to an Seurat object...\n")
  # srt@assays$spliced <- velocity@assays$spliced
  # srt@assays$unspliced <- velocity@assays$unspliced
  # srt@assays$ambiguous <- velocity@assays$ambiguous
  # srt$nCount_spliced <- velocity$nCount_spliced
  # srt$nFeature_spliced <- velocity$nFeature_spliced
  # srt$nCount_unspliced <- velocity$nCount_unspliced
  # srt$nFeature_unspliced <- velocity$nFeature_unspliced
  # srt$nCount_ambiguous <- velocity$nCount_ambiguous
  # srt$nFeature_ambiguous <- velocity$nFeature_ambiguous
  velocityList[[samples[i]]] <- velocity
}
RenameFeatures <- function(srt, newnames = NULL) {
  for (i in Seurat::Assays(srt)) {
    assay <- GetAssay(srt, i)
    for (d in c("counts", "data", "scale.data", "meta.features")) {
      if (dim(slot(assay, d))[1] == length(newnames)) {
        rownames(slot(assay, d)) <- newnames
      } else {
        message(paste0("Slot ", d, " have a different number of features."))
      }
    }
    srt@assays[[i]] <- assay
  }
  return(srt)
}
velocityMerge <- Reduce(merge, x = velocityList)

srt <- LoadH5Seurat("~/Normalization-logCPM/global/Integration/22_Cellline,22_Cellline_injection_10d,22_Cellline_injection_20d,22_Cellline_injection_30d,22_Cellline_injection_40d,22_Cellline_injection_50d,22_Cellline_injection_70d,22_Cellline_injection_90d.h5Seurat")
human_gene <- rownames(velocityMerge)[grep(pattern = "^GRCh38-", x = rownames(velocityMerge))]
velocityMerge <- subset(x = velocityMerge, cells = colnames(srt), features = human_gene)
velocityMerge <- RenameFeatures(velocityMerge, newnames = gsub(x = rownames(velocityMerge), pattern = "GRCh38-", replacement = ""))
emat <- velocityMerge$spliced
nmat <- velocityMerge$unspliced
emb <- srt@reductions$UncorrectedUMAP2d@cell.embeddings
cell.dist <- as.dist(1 - armaCor(t(srt@reductions$UncorrectedUMAP2d@cell.embeddings)))
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat, nmat,
  deltaT = 2,
  kCells = 10,
  cell.dist = cell.dist,
  fit.quantile = fit.quantile,
  n.cores = 1
)
labels <- levels(srt[["Uncorrectedclusters", drop = TRUE]])
labels_tb <- table(srt[["Uncorrectedclusters", drop = TRUE]])
gg <- DimPlot(srt, reduction = "UncorrectedUMAP2d", group.by = "Uncorrectedclusters") +
  scale_color_manual(
    name = paste0("Total cells: ", ncol(srt)),
    values = colorRampPalette(brewer.pal(12, "Paired"))(nlevels(srt[["Uncorrectedclusters", drop = TRUE]])),
    labels = paste0(labels, "(", labels_tb[labels], ")")
  ) +
  labs(
    title = "Uncorrected UMAP",
    x = "UMAP_1",
    y = "UMAP_2"
  ) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3))) +
  theme_zh() +
  theme(aspect.ratio = 1)
ggplot_build(gg)$data
colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
names(colors) <- rownames(emb)
p1 <- show.velocity.on.embedding.cor(emb, rvel.cd,
  n = 30, scale = "sqrt",
  cell.colors = ac(colors, alpha = 0.5),
  cex = 0.8, arrow.scale = 2, show.grid.flow = T,
  min.grid.cell.mass = 1.0, grid.n = 50, arrow.lwd = 1,
  do.par = F, cell.border.alpha = 0.1,
  n.cores = 120, main = "Cell Velocity"
)
