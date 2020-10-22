#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
script_path <- as.character(args[1])
SCPanalysis_dir <- as.character(args[2])
NGSmodule_SCP_dir <- as.character(args[3])
threads <- as.numeric(args[4])
datasets_raw <- as.character(args[5])
species <- as.character(args[6])
exogenous_genes <- as.character(args[7])
cell_calling_methodNum <- as.numeric(args[8])
HVF_source <- as.character(args[9])
nHVF <- as.numeric(args[10])
anchor_dims <- 1:as.numeric(args[11])
integrate_dims <- 1:as.numeric(args[12])
maxPC <- as.numeric(args[13])
resolution <- as.numeric(args[14])
Ensembl_version <- 101


# ##### test #####
# # parameters: global settings ---------------------------------------------
# SCPanalysis_dir <- "/data/lab/HuangMingQian/scRNA-seq/ESC-PGC-GSCLC-new/NGSmodule_SCP_analysis/Integration/"
# NGSmodule_SCP_dir <- "/data/lab/HuangMingQian/scRNA-seq/ESC-PGC-GSCLC-new/NGSmodule_SCP_work/"
# threads <- 120
# datasets_raw <- "ESC,iPSC,PGCd6,CellLine1,CellLine2"
# 
# species <- "Homo_sapiens"
# exogenous_genes <- "GFP"
# 
# # parameters: cell filtering ----------------------------------------------
# cell_calling_methodNum <- 3
# 
# # parameters: integration -------------------------------------------------
# HVF_source <- "global"
# nHVF <- 4000
# anchor_dims <- 1:30
# integrate_dims <- 1:30
# 
# # parameters: clustering --------------------------------------------------
# maxPC <- 100
# resolution <- 1

# Library -----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(invisible(lapply(
  c(
    "sctransform", "Seurat", "SeuratWrappers", "intrinsicDimension", "scater", "Matrix", "BiocParallel",
    "future", "reticulate", "harmony", "plyr", "dplyr", "RColorBrewer", "scales", "gtools",
    "ggsci", "ggpubr", "ggplot2", "ggtree", "cowplot", "reshape2", "stringr", "scDblFinder", 
    "velocyto.R", "biomaRt", "rvest", "xml2"
  ),
  require,
  character.only = TRUE
))))
set.seed(11)

datasets <- strsplit(datasets_raw, split = ";") %>%
  unlist() %>%
  lapply(., function(x) {
    strsplit(x, split = ",") %>% unlist()
  })
samples <- datasets %>%
  unlist() %>%
  unique()

if (species == "Homo_sapiens") {
  cc_S_genes <- Seurat::cc.genes.updated.2019$s.genes
  cc_G2M_genes <- Seurat::cc.genes.updated.2019$g2m.genes
} else {
  species_split <- unlist(strsplit(species, split = "_"))
  species_homolog <- paste0(tolower(substring(species_split[1], 1, 1)), species_split[2], "_homolog_associated_gene_name")

  archives <- listEnsemblArchives()
  # web <- read_html(httr::RETRY("GET", "http://www.ensembl.org/info/website/archives/index.html?redirect=no", times = 1000, timeout(1000)))
  # urls <- web %>% html_nodes("ul") %>% html_nodes("strong") %>% html_nodes("a") %>% html_attr("href")
  # version <- web %>% html_nodes("ul") %>% html_nodes("strong") %>% html_nodes("a") %>% html_text(trim = TRUE) %>%
  #   gsub(pattern = "(Ensembl )|(:.*)",replacement = "",x = .,perl = T)
  # archives <- data.frame(version=version,url=urls,stringsAsFactors = F)
  url <- archives[which(archives$version == Ensembl_version), "url"]

  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = url)
  homolog <- listAttributes(mart)$name

  if (species_homolog %in% homolog) {
    cc_S_genes <- getBM(
      mart = mart,
      attributes = c(species_homolog),
      filters = c("external_gene_name"),
      values = list(Seurat::cc.genes.updated.2019$s.genes)
    )[[1]]
    cc_G2M_genes <- getBM(
      mart = mart,
      attributes = c(species_homolog),
      filters = c("external_gene_name"),
      values = list(Seurat::cc.genes.updated.2019$g2m.genes)
    )[[1]]

    if (length(cc_S_genes) < 3 | length(cc_G2M_genes) < 3) {
      warning(paste0("number of cell-cycle homolog genes is too small. CellCycleScoring will not performed."))
    }
  } else {
    warning(paste0("Can not find the homolog attributes for the species: ", species, " (", species_homolog, ")"))
  }
}


########################### Start the workflow ############################
setwd(SCPanalysis_dir)
options(expressions = 5e5)
options(future.globals.maxSize = 754 * 1000 * 1024^2)
# options(future.fork.enable = TRUE)
if (threads >= 125) {
  cat("Threads number is too large. Re-set it to the 125.\n")
  threads <- 125
}
plan(multiprocess, workers = threads, gc = TRUE) # stop with the command 'future:::ClusterRegistry("stop")'
plan()

script_dir <- gsub(x = script_path, pattern = "Integration.R", replacement = "")
source(paste0(script_dir, "/SCP-workflow-funtcion.R"))

# source("/data/lab/LiLaiHua/scRNA-seq/Gonadal_ridge/analysis_zh/scRNA-SeuratWorkflow-function.R")
# source("/home/zhanghao/Documents/pipeline/Single_cell/customize_Seurat_FeaturePlot.R")

# Preprocessing: load data ------------------------------------------------
for (i in 1:length(samples)) {
  cat("++++++", samples[i], "++++++", "\n")
  cell_upset <- as.data.frame(readRDS(file = paste0(NGSmodule_SCP_dir, "/", samples[i], "/Alignment/Cellranger/", samples[i], "/CellCalling/cell_upset.rds")))
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
    value = Read10X(data.dir = paste0(NGSmodule_SCP_dir, "/", samples[i], "/Alignment/Cellranger/", samples[i], "/outs/raw_feature_bc_matrix/"))[, cells]
  )
  assign(
    x = paste0(samples[i], "_velocity"),
    value = ReadVelocity(file = paste0(NGSmodule_SCP_dir, "/", samples[i], "/Alignment/Cellranger/", samples[i], "/velocyto/", samples[i], ".loom"))
  )
}

# Preprocessing: create Seurat object -------------------------------------
sc_list <- list()
velocity_list <- list()
for (i in 1:length(samples)) {
  cat("++++++", samples[i], "++++++", "\n")
  srt <- CreateSeuratObject(counts = get(paste0(samples[i], "_10X")), project = samples[i])
  srt[["orig.ident"]] <- samples[i]
  srt[["percent.mt"]] <- PercentageFeatureSet(object = srt, pattern = "^(MT-)|(mt-)|(Mt-)")
  srt[["cellcalling_method"]] <- get(paste0(samples[i], "_cellcalling"))[Cells(srt), "Method_comb"]
  srt[["cellcalling_methodNum"]] <- get(paste0(samples[i], "_cellcalling"))[Cells(srt), "Method_num"]
  srt <- RenameCells(object = srt, add.cell.id = samples[i])
  sc_list[[samples[i]]] <- srt

  velocity <- as.Seurat(x = get(paste0(samples[i], "_velocity")), project = samples[i])
  velocity <- RenameCells(
    object = velocity,
    new.names = gsub(x = colnames(velocity), pattern = ".*:", replacement = "", perl = TRUE) %>%
      gsub(x = ., pattern = "x$", replacement = "-1", perl = TRUE) %>%
      paste0(samples[i], "_", .)
  )
  velocity[["orig.ident"]] <- samples[i]
  velocity_list[[samples[i]]] <- velocity
}

if (length(sc_list) == 1) {
  sc_merge <- sc_list[[1]]
  velocity_merge <- velocity_list[[1]]
} else {
  sc_merge <- Reduce(function(x, y) merge(x, y), sc_list)
  Idents(sc_merge) <- sc_merge[["orig.ident"]] <- factor(sc_merge[["orig.ident", drop = TRUE]], levels = samples)
  velocity_merge <- Reduce(function(x, y) merge(x, y), velocity_list)
  Idents(velocity_merge) <- velocity_merge[["orig.ident"]] <- factor(velocity_merge[["orig.ident", drop = TRUE]], levels = samples)
}
if (!file.exists("sc_merge.rds")) {
  saveRDS(sc_merge,file = "sc_merge.rds")
}
if (!file.exists("velocity_merge.rds")) {
  saveRDS(velocity_merge,file = "velocity_merge.rds")
}

# Preprocessing: cell filtering -----------------------------------
sc_list_filter <- lapply(setNames(samples, samples), function(sc_set) {
  cat("++++++", sc_set, "++++++", "\n")
  srt <- sc_list[[sc_set]]
  ntotal <- ncol(srt)
  sce <- as.SingleCellExperiment(srt)
  sce <- scDblFinder(sce, verbose = FALSE)
  ndoublets <- sum(sce$scDblFinder.class == "doublet")
  sce <- subset(sce, , scDblFinder.class == "singlet")
  srt <- subset(x = srt, cells = colnames(sce))

  log10_total_counts <- log10(srt[["nCount_RNA", drop = TRUE]])
  log10_total_features <- log10(srt[["nFeature_RNA", drop = TRUE]])

  mod <- loess(log10_total_features ~ log10_total_counts)
  pred <- predict(mod, newdata = data.frame(log10_total_counts = log10_total_counts))
  featcount_dist <- log10_total_features - pred
  # ggplot(data.frame(x=log10_total_counts,y=log10_total_features,diff=diff),
  #        aes(x=x, y=y, colour=diff)) +
  #   geom_point() + geom_smooth(method = "loess", col="black")

  sce <- addPerCellQC(sce, percent_top = c(20))
  pct_counts_in_top_20_features <- colData(sce)$percent_top_20

  pct_counts_Mt <- srt[["percent.mt", drop = TRUE]]

  filters <- c(
    "log10_total_counts:higher:2.5",
    "log10_total_counts:lower:5",
    "log10_total_features:higher:2.5",
    "log10_total_features:lower:5",
    "pct_counts_in_top_20_features:both:5",
    "featcount_dist:both:5"
  )
  out <- lapply(strsplit(filters, ":"), function(f) {
    which(isOutlier(get(f[1]),
      log = FALSE,
      nmads = as.numeric(f[3]),
      type = f[2]
    ))
  })
  mtout <- isOutlier(pct_counts_Mt, nmads = 3, type = "lower") |
    (isOutlier(pct_counts_Mt, nmads = 2.5, type = "higher") & pct_counts_Mt > 0.08)
  out <- c(out, list(mt = which(mtout)))
  out <- table(unlist(out))
  out <- as.numeric(names(out)[which(out >= 2)])

  cat(">>>", "Total cells:", ntotal, "\n")
  cat(">>>", "Filter out ", ndoublets + length(out), " cells (potential doublets: ", ndoublets, " and ", " unqualified cells: ", length(out), ")", "\n",sep = "")
  cat(">>>", "Filtered cells:", ntotal - ndoublets - length(out), "\n")

  if (length(out) > 0) {
    srt <- subset(srt, cell = colnames(srt)[-out])
  }

  return(srt)
})
if (!file.exists("sc_list_filter.rds")) {
  saveRDS(sc_list_filter,file = "sc_list_filter.rds")
}

# Integration: Standard workflow ------------------------------------------
sc_list_filter_Standard <- lapply(setNames(samples, samples), function(sc_set) {
  cat("++++++", sc_set, "++++++", "\n")
  srt <- sc_list_filter[[sc_set]]
  srt <- Standard_SCP(
    sc = srt, nHVF = nHVF, maxPC = maxPC, resolution = resolution,
    cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
    exogenous_genes = exogenous_genes, assay = "RNA", reduction = NULL
  )
  return(srt)
})

srt_list_Standard <- lapply(setNames(datasets,sapply(datasets, function(x) paste0(x,collapse = ","))),function(dataset){
  cat("++++++", paste0(dataset, collapse = "-"), "++++++", "\n")
  if (length(dataset) == 0) {
    srt_integrated <- NULL
  }
  if (length(dataset) == 1) {
    srt_integrated <- Standard_SCP(
      sc = sc_list_filter_Standard[[dataset]], nHVF = nHVF, maxPC = maxPC, resolution = resolution,
      cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
      exogenous_genes = exogenous_genes, assay = "RNA"
    )
  }
  if (length(dataset) > 1) {
    srt_integrated <- Standard_integrate(
      sc_list = sc_list_filter_Standard[dataset], nHVF = nHVF, anchor_dims = anchor_dims, integrate_dims = integrate_dims, maxPC = maxPC, resolution = resolution,
      HVF_source = HVF_source,
      cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
      exogenous_genes = exogenous_genes
    )
  }
  return(srt_integrated)
})

if (!file.exists("srt_list_Standard.rds")) {
  saveRDS(object = srt_list_Standard, file = "srt_list_Standard.rds")
}


# Integration: SCTransform workflow  --------------------------------------
sc_list_filter_SCT <- lapply(setNames(samples, samples), function(sc_set) {
  cat("++++++", paste0(dataset, collapse = "-"), "++++++", "\n")
  srt <- sc_list_filter[[sc_set]]
  srt <- SCTransform(
    object = srt,
    variable.features.n = nHVF,
    return.only.var.genes = TRUE,
    assay = "RNA"
  )
  VariableFeatures(srt) <- HVFInfo(srt) %>%
    filter((!rownames(.) %in% exogenous_genes)) %>%
    arrange(desc(residual_variance)) %>%
    rownames(.) %>%
    head(n = nHVF)
  return(srt)
})

srt_list_SCT <- lapply(setNames(datasets,sapply(datasets, function(x) paste0(x,collapse = ","))),function(dataset){
  cat("++++++", paste0(dataset, collapse = "-"), "++++++", "\n")
  if (length(dataset) == 0) {
    srt_integrated <- NULL
  }
  if (length(dataset) == 1) {
    srt_integrated <- SCTransform_SCP(
      sc = sc_list_filter_SCT[[dataset]], nHVF = nHVF, maxPC = maxPC, resolution = resolution,
      cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
      exogenous_genes = exogenous_genes, assay = "RNA"
    )
  }
  if (length(dataset) > 1) {
    srt_integrated <- SCTransform_integrate(
      sc_list = sc_list_filter_SCT[dataset], nHVF = nHVF, anchor_dims = anchor_dims, integrate_dims = integrate_dims, maxPC = maxPC, resolution = resolution,
      HVF_source = HVF_source,
      cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
      exogenous_genes = exogenous_genes
    )
  }
  return(srt_integrated)
})

if (!file.exists("srt_list_SCT.rds")) {
  saveRDS(object = srt_list_SCT, file = "srt_list_SCT.rds")
}


# Integration: fastMNN workflow -------------------------------------------
srt_list_fastMNN <- lapply(setNames(datasets,sapply(datasets, function(x) paste0(x,collapse = ","))),function(dataset){
  cat("++++++", paste0(dataset, collapse = "-"), "++++++", "\n")
  if (length(dataset) == 0) {
    srt_integrated <- NULL
  }
  if (length(dataset) == 1) {
    srt_integrated <- Standard_SCP(
      sc = sc_list_filter_Standard[[dataset]], nHVF = nHVF, maxPC = maxPC, resolution = resolution,
      cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
      exogenous_genes = exogenous_genes, assay = "RNA"
    )
  }
  if (length(dataset) > 1) {
    srt_integrated <- fastMNN_integrate(
      sc_list = sc_list_filter_Standard[dataset], nHVF = nHVF, maxPC = maxPC, resolution = resolution,
      HVF_source = HVF_source,
      cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
      exogenous_genes = exogenous_genes
    )
  }
  return(srt_integrated)
})

if (!file.exists("srt_list_fastMNN.rds")) {
  saveRDS(object = srt_list_fastMNN, file = "srt_list_fastMNN.rds")
}


# Integration: Harmony workflow -------------------------------------------
srt_list_Harmony <- lapply(setNames(datasets,sapply(datasets, function(x) paste0(x,collapse = ","))),function(dataset){
  cat("++++++", paste0(dataset, collapse = "-"), "++++++", "\n")
  if (length(dataset) == 0) {
    srt_integrated <- NULL
  }
  if (length(dataset) == 1) {
    srt_integrated <- Standard_SCP(
      sc = sc_list_filter_Standard[[dataset]], nHVF = nHVF, maxPC = maxPC, resolution = resolution,
      cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
      exogenous_genes = exogenous_genes, assay = "RNA"
    )
  }
  if (length(dataset) > 1) {
    srt_integrated <- Harmony_integrate(
      sc_list = sc_list_filter_Standard[dataset], nHVF = nHVF, maxPC = maxPC, resolution = resolution,
      HVF_source = HVF_source,
      cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
      exogenous_genes = exogenous_genes
    )
  }
 return(srt_integrated)
})

if (!file.exists("srt_list_Harmony.rds")) {
  saveRDS(object = srt_list_Harmony, file = "srt_list_Harmony.rds")
}

# for (srt_name in c("srt_list_Standard","srt_list_SCT","srt_list_fastMNN","srt_list_Harmony")) {
#   srt_use <- get(srt_name)[[1]]
#   DefaultAssay(srt_use) <- "RNA"
#   srt_use$batch <- srt_use$orig.ident
#   project_name <- srt_use@project.name
#   p <- summary_plot(
#     srt = srt_use, return_list = F,features = "DDX4", color_by = "seurat_clusters", reduction = "umap", split_by = "batch", palette = "nejm",
#     do_save = T, file_save = paste0(srt_name, ".summary.png")
#   )
# }


future:::ClusterRegistry("stop")
