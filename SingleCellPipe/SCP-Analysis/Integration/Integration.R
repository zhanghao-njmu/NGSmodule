#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
script_path <- as.character(args[1])
SCPanalysis_dir <- as.character(args[2])
SCPwork_dir <- as.character(args[3])
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
# SCPwork_dir <- "/data/lab/HuangMingQian/scRNA-seq/ESC-PGC-GSCLC-new/NGSmodule_SCP_work/"
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
datasets <- datasets[sapply(datasets, length) > 1]

# samples <- datasets %>%
#   unlist() %>%
#   unique()
samples <- list.dirs(path = SCPwork_dir, recursive = FALSE, full.names = FALSE)

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
if (file.exists("sc_list.rds") & file.exists("velocity_list.rds")) {
  cat("Loading the sc_list and velocity_list from the existing file....\n")
  sc_list <- readRDS("sc_list.rds")
  velocity_list <- readRDS("velocity_list.rds")
} else {
  for (i in 1:length(samples)) {
    cat("++++++", samples[i], "++++++", "\n")
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

  # Preprocessing: create Seurat object -------------------------------------
  sc_list <- list()
  velocity_list <- list()
  for (i in 1:length(samples)) {
    cat("++++++", samples[i], "++++++", "\n")
    srt <- CreateSeuratObject(counts = get(paste0(samples[i], "_10X")), project = samples[i])
    srt[["orig.ident"]] <- samples[i]
    srt[["percent.mt"]] <- PercentageFeatureSet(object = srt, pattern = "^(MT-)|(Mt-)|(mt-)")
    srt[["percent.ribo"]] <- PercentageFeatureSet(object = srt, pattern = "^(RP[SL]\\d+)|((Rp[sl]\\d+))|(rp[sl]\\d+)")
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

    assign(x = paste0(samples[i], "_10X"), value = NULL)
    assign(x = paste0(samples[i], "_velocity"), value = NULL)
  }
  saveRDS(sc_list, file = "sc_list.rds")
  saveRDS(velocity_list, file = "velocity_list.rds")
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


# Preprocessing: cell filtering -----------------------------------
if (file.exists("sc_list_filter.rds")) {
  cat("Loading the sc_list_filter from the existing file....\n")
  sc_list_filter <- readRDS("sc_list_filter.rds")
} else {
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
    out <- table(unlist(out))
    out <- as.numeric(names(out)[which(out >= 2)])

    pct_counts_Mt <- srt[["percent.mt", drop = TRUE]]
    if (all(pct_counts_Mt > 0 & pct_counts_Mt < 1)) {
      pct_counts_Mt <- pct_counts_Mt * 100
    }
    mt_out <- isOutlier(pct_counts_Mt, nmads = 3, type = "lower") |
      (isOutlier(pct_counts_Mt, nmads = 2.5, type = "higher") & pct_counts_Mt > 10) |
      (pct_counts_Mt > 20)
    total_out <- unique(c(out,as.numeric(which(mt_out))))

    cat(">>>", "Total cells:", ntotal, "\n")
    cat(">>>", "Filter out", ndoublets + length(total_out), "cells (potential doublets: ", ndoublets, "and unqualified cells:", length(total_out), ")", "\n")
    cat(">>>", "Filtered cells:", ntotal - ndoublets - length(total_out), "\n")

    if (length(total_out) > 0) {
      srt <- subset(srt, cell = colnames(srt)[-total_out])
    }

    return(srt)
  })

  saveRDS(sc_list_filter, file = "sc_list_filter.rds")
}


# Preprocessing: basic srt normalization -----------------------------------
if (file.exists("sc_list_filter_Standard.rds")) {
  cat("Loading the sc_list_filter_Standard from the existing file....\n")
  sc_list_filter_Standard <- readRDS("sc_list_filter_Standard.rds")
} else {
  sc_list_filter_Standard <- lapply(setNames(samples, samples), function(sc_set) {
    cat("++++++", sc_set, "++++++", "\n")
    srt <- sc_list_filter[[sc_set]]
    srt <- Standard_SCP(
      sc = srt, nHVF = nHVF, maxPC = maxPC, resolution = resolution,
      cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
      exogenous_genes = exogenous_genes, assay = "RNA"
    )
    return(srt)
  })
  saveRDS(sc_list_filter_Standard, file = "sc_list_filter_Standard.rds")
}

if (file.exists("sc_list_filter_SCT.rds")) {
  cat("Loading the sc_list_filter_SCT from the existing file....\n")
  sc_list_filter_SCT <- readRDS("sc_list_filter_SCT.rds")
} else {
  sc_list_filter_SCT <- lapply(setNames(samples, samples), function(sc_set) {
    cat("++++++", sc_set, "++++++", "\n")
    srt <- sc_list_filter[[sc_set]]
    srt <- SCTransform_SCP(
      sc = srt, nHVF = nHVF, maxPC = maxPC, resolution = resolution,
      cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
      exogenous_genes = exogenous_genes, assay = "RNA"
    )
    return(srt)
  })
  saveRDS(object = sc_list_filter_SCT, file = "sc_list_filter_SCT.rds")
}


# Individual: Standard normalization ------------------------------
dir.create("Individual-Standard", recursive = T, showWarnings = FALSE)
for (sample in samples) {
  cat("++++++", sample, "(Individual-Standard)", "++++++", "\n")
  if (file.exists(paste0("Individual-Standard/", sample, ".rds"))) {
    cat(">>> Individual-Standard process for the", sample, "has finished. Skip to the next step.\n")
    next
  } else {
    srt <- sc_list_filter_Standard[[sample]]
    saveRDS(srt, paste0("Individual-Standard/", sample, ".rds"))
    cat(">>> Individual-Standard process for the", sample, "completed successfully.\n")
  }
}


# Individual: SCTransform normalization ------------------------------
dir.create("Individual-SCTransform", recursive = T, showWarnings = FALSE)
for (sample in samples) {
  cat("++++++", sample, "(Individual-SCTransform)", "++++++", "\n")
  if (file.exists(paste0("Individual-SCTransform/", sample, ".rds"))) {
    cat(">>> Individual-SCTransform process for the", sample, "has finished. Skip to the next step.\n")
    next
  } else {
    cat("++++++", sample, "++++++", "\n")
    srt <- sc_list_filter_SCT[[sample]]
    saveRDS(srt, paste0("Individual-SCTransform/", sample, ".rds"))
    cat(">>> Individual-SCTransform process for the", sample, "completed successfully.\n")
  }
}


if (length(datasets) != 0) {
  # Integration: Simple merge ----------------------------------------------
  dir.create("Integration-SimpleMerge-Standard", recursive = T, showWarnings = FALSE)
  for (dataset in datasets) {
    cat("++++++", paste0(dataset, collapse = ","), "++++++", "\n")
    if (file.exists(paste0("Integration-SimpleMerge-Standard/", paste0(dataset, collapse = ","), ".rds"))) {
      cat(">>> Integration-SimpleMerge-Standard process for the", paste0(dataset, collapse = ","), "has finished. Skip to the next step.\n")
      next
    } else {
      sc_list_filter_merge <- Reduce(function(x, y) merge(x, y), sc_list_filter[dataset])
      Idents(sc_list_filter_merge) <- sc_list_filter_merge[["orig.ident"]] <- factor(sc_list_filter_merge[["orig.ident", drop = TRUE]], levels = dataset)
      srt <- Standard_SCP(
        sc = sc_list_filter_merge, nHVF = nHVF, maxPC = maxPC, resolution = resolution,
        cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
        exogenous_genes = exogenous_genes, assay = "RNA",
      )
      saveRDS(srt, file = paste0("Integration-SimpleMerge-Standard/", paste0(dataset, collapse = ","), ".rds"))
      cat(">>> Integration-SimpleMerge-Standard process for the", paste0(dataset, collapse = ","), "completed successfully.\n")
    }
  }

  dir.create("Integration-SimpleMerge-SCTransform", recursive = T, showWarnings = FALSE)
  for (dataset in datasets) {
    cat("++++++", paste0(dataset, collapse = ","), "++++++", "\n")
    if (file.exists(paste0("Integration-SimpleMerge-SCTransform/", paste0(dataset, collapse = ","), ".rds"))) {
      cat(">>> Integration-SimpleMerge-SCTransform process for the", paste0(dataset, collapse = ","), "has finished. Skip to the next step.\n")
      next
    } else {
      sc_list_filter_merge <- Reduce(function(x, y) merge(x, y), sc_list_filter[dataset])
      Idents(sc_list_filter_merge) <- sc_list_filter_merge[["orig.ident"]] <- factor(sc_list_filter_merge[["orig.ident", drop = TRUE]], levels = dataset)
      srt <- SCTransform_SCP(
        sc = sc_list_filter_merge, nHVF = nHVF, maxPC = maxPC, resolution = resolution,
        cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
        exogenous_genes = exogenous_genes, assay = "RNA"
      )
      saveRDS(srt, file = paste0("Integration-SimpleMerge-SCTransform/", paste0(dataset, collapse = ","), ".rds"))
      cat(">>> Integration-SimpleMerge-SCTransform process for the", paste0(dataset, collapse = ","), "completed successfully.\n")
    }
  }


  # Integration: Standard workflow ------------------------------------------
  dir.create("Integration-Standard", recursive = T, showWarnings = FALSE)
  for (dataset in datasets) {
    cat("++++++", paste0(dataset, collapse = ","), "(Integration-Standard)", "++++++", "\n")
    if (file.exists(paste0("Integration-Standard/", paste0(dataset, collapse = ","), ".rds"))) {
      cat(">>> Integration-Standard process for the", paste0(dataset, collapse = ","), "has finished. Skip to the next step.\n")
      next
    } else {
      srt_integrated <- Standard_integrate(
        sc_list = sc_list_filter_Standard[dataset], nHVF = nHVF, anchor_dims = anchor_dims, integrate_dims = integrate_dims, maxPC = maxPC, resolution = resolution,
        HVF_source = HVF_source,
        cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
        exogenous_genes = exogenous_genes
      )
      saveRDS(srt_integrated, file = paste0("Integration-Standard/", paste0(dataset, collapse = ","), ".rds"))
      cat(">>> Integration-Standard process for the", paste0(dataset, collapse = ","), "completed successfully.\n")
    }
  }


  # Integration: SCTransform workflow  --------------------------------------
  dir.create("Integration-SCTransform", recursive = T, showWarnings = FALSE)
  for (dataset in datasets) {
    cat("++++++", paste0(dataset, collapse = ","), "(Integration-SCTransform)", "++++++", "\n")
    if (file.exists(paste0("Integration-SCTransform/", paste0(dataset, collapse = ","), ".rds"))) {
      cat(">>> Integration-SCTransform process for the", paste0(dataset, collapse = ","), "has finished. Skip to the next step.\n")
      next
    } else {
      srt_integrated <- SCTransform_integrate(
        sc_list = sc_list_filter_SCT[dataset], nHVF = nHVF, anchor_dims = anchor_dims, integrate_dims = integrate_dims, maxPC = maxPC, resolution = resolution,
        HVF_source = HVF_source,
        cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
        exogenous_genes = exogenous_genes
      )
      saveRDS(srt_integrated, file = paste0("Integration-SCTransform/", paste0(dataset, collapse = ","), ".rds"))
      cat(">>> Integration-SCTransform process for the", paste0(dataset, collapse = ","), "completed successfully.\n")
    }
  }


  # Integration: fastMNN workflow -------------------------------------------
  dir.create("Integration-fastMNN", recursive = T, showWarnings = FALSE)
  for (dataset in datasets) {
    cat("++++++", paste0(dataset, collapse = ","), "(Integration-fastMNN)", "++++++", "\n")
    if (file.exists(paste0("Integration-fastMNN/", paste0(dataset, collapse = ","), ".rds"))) {
      cat(">>> Integration-fastMNN process for the", paste0(dataset, collapse = ","), "has finished. Skip to the next step.\n")
      next
    } else {
      srt_integrated <- fastMNN_integrate(
        sc_list = sc_list_filter_Standard[dataset], nHVF = nHVF, maxPC = maxPC, resolution = resolution,
        HVF_source = HVF_source,
        cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
        exogenous_genes = exogenous_genes
      )
      saveRDS(srt_integrated, file = paste0("Integration-fastMNN/", paste0(dataset, collapse = ","), ".rds"))
      cat(">>> Integration-fastMNN process for the", paste0(dataset, collapse = ","), "completed successfully.\n")
    }
  }


  # Integration: Harmony workflow -------------------------------------------
  dir.create("Integration-Harmony", recursive = T, showWarnings = FALSE)
  for (dataset in datasets) {
    cat("++++++", paste0(dataset, collapse = ","), "(Integration-Harmony)", "++++++", "\n")
    if (file.exists(paste0("Integration-Harmony/", paste0(dataset, collapse = ","), ".rds"))) {
      cat(">>> Integration-Harmony process for the", paste0(dataset, collapse = ","), "has finished. Skip to the next step.\n")
      next
    } else {
      srt_integrated <- Harmony_integrate(
        sc_list = sc_list_filter_Standard[dataset], nHVF = nHVF, maxPC = maxPC, resolution = resolution,
        HVF_source = HVF_source,
        cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
        exogenous_genes = exogenous_genes
      )
      saveRDS(srt_integrated, file = paste0("Integration-Harmony/", paste0(dataset, collapse = ","), ".rds"))
      cat(">>> Integration-Harmony process for the", paste0(dataset, collapse = ","), "completed successfully.\n")
    }
  }
} else {
  cat("Integration skipped.")
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
