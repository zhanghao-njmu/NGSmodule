#!/usr/bin/env Rscript

# parameters: global settings ---------------------------------------------
work_dir <- "/data/lab/LiLaiHua/scRNA-seq/Gonadal_ridge/analysis_zh/"
NGSmodule_SCP_dir <- "/data/lab/LiLaiHua/scRNA-seq/Gonadal_ridge/NGSmodule_SCP_work/"
threads <- 80
samples <- c("d0", "d1","d3","d5","d7")
compare_list <- list(
  c("d0", "d1","d3","d5","d7")
)

cc_S_genes <- cc.genes.updated.2019$s.genes
cc_G2M_genes <- cc.genes.updated.2019$g2m.genes
exogenous_gene <- NULL

# parameters: cell filtering ----------------------------------------------
cell_calling <- 4
minGene <- 1000
maxGene <- 10000
minUMI <- 1000
maxUMI <- 40000
maxMitoPercent <- 20

# parameters: integration -------------------------------------------------
HVF_source <- "separate"
nHVF <- 3000
anchor_dims <- 1:30
integrate_dims <- 1:30

# parameters: clustering --------------------------------------------------
maxPC <- 100
resolution <- 0.8




########################### Start the workflow ############################
library(sctransform)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(intrinsicDimension)
library(Matrix)
library(BiocParallel)
library(future)
library(reticulate)
library(harmony)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(scales)
library(gtools)
library(ggsci)
library(ggpubr)
library(ggplot2)
library(ggtree)
library(cowplot)
library(reshape2)
library(stringr)
library(velocyto.R)

setwd(work_dir)
options(expressions = 5e5)
options(future.globals.maxSize = 754 * 1000 * 1024^2)
options(future.fork.enable = TRUE)
plan(multiprocess, workers = threads, gc = TRUE) # stop with the command 'future:::ClusterRegistry("stop")'
plan()

source("/data/lab/LiLaiHua/scRNA-seq/Gonadal_ridge/analysis_zh/scRNA-SeuratWorkflow-function.R")
source("/home/zhanghao/Documents/pipeline/Single_cell/customize_Seurat_FeaturePlot.R")
plotlist <- list()

# Preprocessing: load data ------------------------------------------------
for (i in 1:length(samples)) {
  cat("[", i, "]", "samples:", samples[i], "\n", sep = "")
  assign(
    x = paste0(samples[i], "_10X"),
    value = Read10X(data.dir = paste0(NGSmodule_SCP_dir,"/",samples[i], "/", "Alignment/Cellranger/",samples[i], "/outs/filtered_feature_bc_matrix/"))
  )
  assign(
    x = paste0(samples[i], "_velocity"),
    value = ReadVelocity(file = paste0(NGSmodule_SCP_dir,"/",samples[i], "/", "Alignment/Cellranger/", samples[i], "/velocyto/", samples[i], ".loom"))
  )
}

# Preprocessing: create Seurat object -------------------------------------
sc_list <- list()
velocity_list <- list()
for (i in 1:length(samples)) {
  cat(samples[i], "\n")
  srt <- CreateSeuratObject(counts = get(paste0(samples[i], "_10X")), project = samples[i])
  srt <- RenameCells(object = srt, add.cell.id = samples[i])
  srt[["orig.ident"]] <- samples[i]
  srt[["percent.mt"]] <- PercentageFeatureSet(object = srt, pattern = "^(MT-)|(mt-)")
  sc_list[[samples[i]]] <- srt

  velocity <- as.Seurat(x = get(paste0(samples[i], "_velocity")), project = samples[i])
  velocity <- RenameCells(
    object = velocity,
    new.names = gsub(x = colnames(velocity), pattern = ".*:", replacement = "",perl = TRUE) %>%
      gsub(x = ., pattern = "x$", replacement = "-1", perl = TRUE) %>% 
      paste0(samples[i],"_",.)
  )
  velocity[["orig.ident"]] <- samples[i]
  velocity_list[[samples[i]]] <- velocity
}
if (length(sc_list) == 1) {
  sc_merge <- sc_list[[1]]
  velocity_merge <- velocity_list[[1]]
} else {
  sc_merge <- Reduce(function(x, y) merge(x, y), sc_list)
  Idents(sc_merge) <- factor(sc_merge[["orig.ident"]], levels = samples)
  # p <- VlnPlot(object = sc_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, combine = F) %>%
  #   custom_Vln_plot(
  #     alpha = 0, y_cutoff = c(minGene, minUMI, maxMitoPercent), pass = c("h", "h", "l"), combine = T,
  #     ncol = 3, AddBox = T, palette = "nejm"
  #   )

  velocity_merge <- Reduce(function(x, y) merge(x, y), velocity_list)
  Idents(velocity_merge) <- factor(velocity_merge[["orig.ident"]], levels = samples)
}

# Preprocessing: cell filtering -----------------------------------
sc_list_filter <- bplapply(setNames(samples, samples), function(sc_set) {
  srt <- sc_list[[sc_set]]
  srt <- subset(
    x = srt,
    subset = nFeature_RNA > minGene &
      nFeature_RNA < maxGene &
      nCount_RNA > minUMI &
      nCount_RNA < maxUMI &
      percent.mt < maxMitoPercent
  )
  return(srt)
}, BPPARAM = MulticoreParam())

# Integration: Standard workflow ------------------------------------------
sc_list_filter_Standard <- bplapply(setNames(samples, samples), function(sc_set) {
  srt <- sc_list_filter[[sc_set]]
  srt <- Standard_SCP(sc = srt,nHVF = nHVF,maxPC = maxPC,resolution = resolution,
                      cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
                      exogenous_genes = NULL, assay = "RNA")

  return(srt)
}, BPPARAM = MulticoreParam())

srt_list_Standard <- list()
for (compare in compare_list) {
  cat(paste0(compare, collapse = "-"), "\n")
  if (length(compare) == 0)  {
    next
  }
  if (length(compare) == 1)  {
    srt_integrated <-  Standard_SCP(sc = sc_list_filter_Standard[[compare]],nHVF = nHVF,maxPC = maxPC,resolution = resolution,
                                    cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
                                    exogenous_genes = NULL, assay = "RNA")
  }
  if (length(compare) > 1) {
    srt_integrated <- Standard_integrate(sc_list = sc_list_filter_Standard[compare], nHVF = nHVF, anchor_dims = anchor_dims, integrate_dims = anchor_dims, maxPC = maxPC, resolution = resolution,
                                         HVF_source = HVF_source,
                                         cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
                                         exogenous_genes = NULL)
  } 

  srt_list_Standard[[srt_integrated@project.name]] <- srt_integrated
}
if (!file.exists("srt_list_Standard.rds")) {
  saveRDS(object = srt_list_Standard, file = "srt_list_Standard.rds")
}


# Integration: SCTransform workflow  --------------------------------------
sc_list_filter_SCT <- bplapply(setNames(samples, samples), function(sc_set) {
  srt <- sc_list_filter[[sc_set]]
  srt <- SCTransform(
    object = srt,
    variable.features.n = nHVF,
    return.only.var.genes = FALSE
  )
  return(srt)
}, BPPARAM = MulticoreParam())

srt_list_SCT <- list()
 for (compare in compare_list) {
  cat(paste0(compare, collapse = "-"), "\n")
   if (length(compare) == 0)  {
     next
   }
   if (length(compare) == 1)  {
     srt_integrated <-  SCTransform_SCP(sc = sc_list_filter_SCT[[compare]],nHVF = nHVF,maxPC = maxPC,resolution = resolution,
                                        cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
                                        exogenous_genes = NULL, assay = "RNA")
   }
   if (length(compare) > 1) {
     srt_integrated <- SCTransform_integrate(sc_list = sc_list_filter_SCT[compare], nHVF = nHVF, anchor_dims = anchor_dims, integrate_dims = anchor_dims, maxPC = maxPC, resolution = resolution,
                                             HVF_source = HVF_source,
                                             cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
                                             exogenous_genes = NULL)
   }

  srt_list_SCT[[srt_integrated@project.name]] <- srt_integrated
}
if (!file.exists("srt_list_SCT.rds")) {
  saveRDS(object = srt_list_SCT, file = "srt_list_SCT.rds")
}


# Integration: fastMNN workflow -------------------------------------------
srt_list_fastMNN <- list()
for (compare in compare_list) {
  cat(paste0(compare, collapse = "-"), "\n")
  if (length(compare) == 0)  {
    next
  }
  if (length(compare) == 1)  {
    srt_integrated <-  Standard_SCP(sc = sc_list_filter_Standard[[compare]],nHVF = nHVF,maxPC = maxPC,resolution = resolution,
                                    cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
                                    exogenous_genes = NULL, assay = "RNA")
  }
  if (length(compare) > 1) {
    srt_integrated <- fastMNN_integrate(sc_list = sc_list_filter_Standard[compare], nHVF = nHVF, anchor_dims = anchor_dims, integrate_dims = anchor_dims, maxPC = maxPC, resolution = resolution,
                                         HVF_source = HVF_source,
                                         cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
                                         exogenous_genes = NULL)
  }
  srt_list_fastMNN[[srt_integrated@project.name]] <- srt_integrated
}
if (!file.exists("srt_list_fastMNN.rds")) {
  saveRDS(object = srt_list_fastMNN, file = "srt_list_fastMNN.rds")
}


# Integration: Harmony workflow -------------------------------------------
srt_list_Harmony <- list()
for (compare in compare_list) {
  cat(paste0(compare, collapse = "-"), "\n")
  if (length(compare) == 0)  {
    next
  }
  if (length(compare) == 1)  {
    srt_integrated <-  Standard_SCP(sc = sc_list_filter_Standard[[compare]],nHVF = nHVF,maxPC = maxPC,resolution = resolution,
                                    cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
                                    exogenous_genes = NULL, assay = "RNA")
  }
  if (length(compare) > 1) {
    srt_integrated <- Harmony_integrate(sc_list = sc_list_filter_Standard[compare], nHVF = nHVF, anchor_dims = anchor_dims, integrate_dims = anchor_dims, maxPC = maxPC, resolution = resolution,
                                        HVF_source = HVF_source,
                                        cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
                                        exogenous_genes = NULL)
  }
  srt_list_Harmony[[srt_integrated@project.name]] <- srt_integrated
}
if (!file.exists("srt_list_Harmony.rds")) {
  saveRDS(object = srt_list_Harmony, file = "srt_list_Harmony.rds")
}


# Run velocity ------------------------------------------------------------
# srt_list_Standard <- readRDS(file = "srt_list_Standard.rds")
# srt_list_SCT <- readRDS(file = "srt_list_SCT.rds")
# srt_list_fastMNN <- readRDS(file = "srt_list_fastMNN.rds")
# srt_list_Harmony <- readRDS(file = "srt_list_Harmony.rds")

for (method in c("Standard", "SCT", "fastMNN", "Harmony")) {
  srt_integrated <- get(paste0("srt_list_", method))[[1]]

  velocity_merge <- subset(velocity_merge, cells = Cells(srt_integrated))
  srt_integrated[["spliced"]] <- velocity_merge[["spliced"]]
  srt_integrated[["unspliced"]] <- velocity_merge[["unspliced"]]
  srt_integrated[["ambiguous"]] <- velocity_merge[["ambiguous"]]

  srt_integrated <- RunVelocity(
    object = srt_integrated,
    spliced = "spliced",
    unspliced = "unspliced",
    ambiguous = "ambiguous",
    reduction = "umap",
    deltaT = 1,
    kCells = 25,
    fit.quantile = 0.02,
    ncores = threads
  )

  ident_colors <- scales::hue_pal()(length(levels(Idents(srt_integrated))))
  names(ident_colors) <- levels(Idents(srt_integrated))
  cell_colors <- ident_colors[Idents(object = srt_integrated)]
  names(cell_colors) <- colnames(srt_integrated)

  vel_out <- show.velocity.on.embedding.cor(
    emb = Embeddings(object = srt_integrated, reduction = "umap"),
    vel = Tool(object = srt_integrated, slot = "RunVelocity"),
    n = 200,
    scale = "sqrt",
    cell.colors = ac(x = cell_colors, alpha = 0.5),
    cex = 0.8,
    arrow.scale = 3,
    show.grid.flow = TRUE,
    min.grid.cell.mass = 0.5,
    grid.n = 40,
    arrow.lwd = 1,
    do.par = FALSE,
    cell.border.alpha = 0.1,
    n.cores = threads,
    return.details = T
  )
  
  if (!file.exists(paste0("srt_list_", method, ".velocity.rds"))) {
    saveRDS(object = vel_out, file = paste0("srt_list_", method, ".velocity.rds"))
  }
  
  
}
