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
mito_threshold <- as.numeric(args[9])
gene_threshold <- as.numeric(args[10])
UMI_threshold <- as.numeric(args[11])
normalization_method <- as.character(args[12])
nHVF <- as.numeric(args[13])
maxPC <- as.numeric(args[14])
resolution <- as.numeric(args[15])
reduction <- as.character(args[16])
HVF_source <- as.character(args[17])
integration_method <- as.character(args[18])
Ensembl_version <- 101


# ##### test #####
# # parameters: global settings ---------------------------------------------
# SCPanalysis_dir <- "/ssd/lab/HuangMingQian/scRNAseq/iPSC-ESC-iMeLC-PGCd6-CellLine-Testis/NGSmodule_SCP_analysis/Integration/"
# SCPwork_dir <- "/ssd/lab/HuangMingQian/scRNAseq/iPSC-ESC-iMeLC-PGCd6-CellLine-Testis/NGSmodule_SCP_work/"
# threads <- 120
# datasets_raw <- "CellLine1,CellLine2,Testis-d10,Testis-d20,Testis-d30,Testis-d40,Testis-d50"
# # datasets_raw <- "ESC,PGCLC-d6,CellLine1,CellLine2;ESC,PGCLC-d6,CellLine1,CellLine2,Testis-d10,Testis-d20,Testis-d30,Testis-d40,Testis-d50;CellLine1,CellLine2,Testis-d10,Testis-d20,Testis-d30,Testis-d40,Testis-d50;Testis-d10,Testis-d20,Testis-d30,Testis-d40,Testis-d50"
# species <- "Homo_sapiens"
# exogenous_genes <- "GFP"
#
# # parameters: cell filtering ----------------------------------------------
# cell_calling_methodNum <- 3
# mito_threshold <- 0.2
# gene_threshold <- 1000
# UMI_threshold <- 3000
#
#
# # parameters: basic --------------------------------------------------
# normalization_method <- "logCPM,SCT" # logCPM,SCT
# nHVF <- 4000
# maxPC <- 100
# resolution <- 1
# reduction <- "umap" # umap,tsne
#
# # parameters: integration -------------------------------------------------
# HVF_source <- "global" # global,separate
# integration_method <- "Uncorrected,Seurat,fastMNN,Harmony,Scanorama,BBKNN,CSS,LIGER" # Seurat,fastMNN,Harmony,Scanorama,BBKNN,CSS,LIGER


# Library -----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(invisible(lapply(
  c(
    "Seurat", "SeuratDisk", "SeuratWrappers", "sctransform", "intrinsicDimension", "scater", "Matrix", "BiocParallel",
    "future", "reticulate", "harmony", "liger", "simspec", "scMerge", "BiocSingular", "iSMNN", "plyr", "dplyr", "RColorBrewer", "scales", "gtools",
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

reduction <- strsplit(reduction, split = ",") %>% unlist()
normalization_method <- strsplit(normalization_method, split = ",") %>% unlist()
integration_method <- strsplit(integration_method, split = ",") %>% unlist()

CCgenes<- CCgene_prefetch(species=species)
cc_S_genes <- CCgenes[["cc_S_genes"]]
cc_G2M_genes <- CCgenes[["cc_G2M_genes"]]

########################### Start the workflow ############################
setwd(SCPanalysis_dir)
options(expressions = 5e5)
options(future.globals.maxSize = 754 * 1000 * 1024^2)
options(future.fork.enable = TRUE)
if (threads >= 125) {
  cat("Threads number is too large. Re-set it to the 125.\n")
  threads <- 125
}
plan(multiprocess, workers = threads, gc = TRUE) # stop with the command 'future:::ClusterRegistry("stop")'
plan()

script_dir <- gsub(x = script_path, pattern = "Integration.R", replacement = "")
source(paste0(script_dir, "/SCP-workflow-funtcion.R"))

# source("/home/zhanghao/Program/NGS/UniversalTools/NGSmodule/SingleCellPipe/SCP-Analysis/Integration/SCP-workflow-funtcion.R")
# source("/home/zhanghao/Documents/pipeline/Single_cell/customize_Seurat_FeaturePlot.R")
save.image(file = "base_env.Rdata")

# Preprocessing: Load data ------------------------------------------------
if (file.exists("sc_list.rds") & file.exists("velocity_list.rds")) {
  if (file.exists("sc_list_filter.rds")) {
    cat("Skip the sc_list loading....\n")
  } else {
    cat("Loading the sc_list from the existing file....\n")
    sc_list <- readRDS("sc_list.rds")
  }
} else {
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
  sc_list <- list()
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
# if (length(sc_list) == 1) {
#   sc_merge <- sc_list[[1]]
#   velocity_merge <- velocity_list[[1]]
# } else {
#   sc_merge <- Reduce(function(x, y) merge(x, y), sc_list)
#   Idents(sc_merge) <- sc_merge[["orig.ident"]] <- factor(sc_merge[["orig.ident", drop = TRUE]], levels = samples)
#   velocity_merge <- Reduce(function(x, y) merge(x, y), velocity_list)
#   Idents(velocity_merge) <- velocity_merge[["orig.ident"]] <- factor(velocity_merge[["orig.ident", drop = TRUE]], levels = samples)
# }


# Preprocessing: Cell filtering -----------------------------------
if (file.exists("sc_list_filter.rds")) {
  if (all(file.exists(paste0("sc_list_filter_", normalization_method, ".rds")))) {
    cat("Skip the sc_list_filter loading....\n")
  } else {
    cat("Loading the sc_list_filter from the existing file....\n")
    sc_list_filter <- readRDS("sc_list_filter.rds")
  }
} else {
  sc_list_filter <- lapply(setNames(samples, samples), function(sc_set) {
    cat("++++++", sc_set, "(Preprocessing-CellFiltering)", "++++++", "\n")
    srt <- sc_list[[sc_set]]
    ntotal <- ncol(srt)
    sce <- as.SingleCellExperiment(srt)
    sce <- scDblFinder(sce, verbose = FALSE)
    ndoublets <- sum(sce[["scDblFinder.class"]] == "doublet")
    srt[["scDblFinder.score"]] <- sce[["scDblFinder.score"]]
    srt[["scDblFinder.class"]] <- sce[["scDblFinder.class"]]

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
    pct_counts_in_top_20_features <- colData(sce)$percent.top_20

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

    umi_out <- which(srt[["nCount_RNA", drop = TRUE]] <= UMI_threshold)
    gene_out <- which(srt[["nFeature_RNA", drop = TRUE]] <= gene_threshold)

    pct_counts_Mt <- srt[["percent.mt", drop = TRUE]]
    if (all(pct_counts_Mt > 0 & pct_counts_Mt < 1)) {
      pct_counts_Mt <- pct_counts_Mt * 100
    }
    if (mito_threshold > 0 & mito_threshold < 1) {
      mito_threshold <- mito_threshold * 100
    }
    mt_out <- which(isOutlier(pct_counts_Mt, nmads = 3, type = "lower") |
      (isOutlier(pct_counts_Mt, nmads = 2.5, type = "higher") & pct_counts_Mt > 10) |
      (pct_counts_Mt > mito_threshold))

    total_out <- unique(c(out, as.numeric(umi_out), as.numeric(gene_out), as.numeric(mt_out)))
    cat(">>>", "Total cells:", ntotal, "\n")
    cat(">>>", "Cells which are filtered out:", ndoublets + length(total_out), "\n")
    cat("...", ndoublets, "potential doublets", "\n")
    cat("...", length(out), "unqualified cells", "\n")
    cat("...", length(umi_out), "low-UMI cells", "\n")
    cat("...", length(gene_out), "low-Gene cells", "\n")
    cat("...", length(mt_out), "high-Mito cells", "\n")
    cat(">>>", "Remained cells after filtering :", ntotal - ndoublets - length(total_out), "\n")

    if (length(total_out) > 0) {
      srt <- subset(srt, cell = colnames(srt)[-total_out])
    }

    return(srt)
  })

  saveRDS(sc_list_filter, file = "sc_list_filter.rds")
}


# Preprocessing: Normalization -----------------------------------
for (nm in normalization_method) {
  dir.create(paste0("Normalization-", nm), recursive = T, showWarnings = FALSE)
  if (file.exists(paste0("sc_list_filter_", nm, ".rds"))) {
    cat("Loading the", paste0("sc_list_filter_", nm), "from the existing file....\n")
    assign(
      x = paste0("sc_list_filter_", nm),
      value = readRDS(paste0("sc_list_filter_", nm, ".rds"))
    )
  } else {
    assign(
      x = paste0("sc_list_filter_", nm),
      value = lapply(setNames(samples, samples), function(sc_set) {
        cat("++++++", sc_set, paste0("Normalization-", nm), "++++++", "\n")
        srt <- sc_list_filter[[sc_set]]
        srt <- Standard_SCP(
          sc = srt, normalization_method = nm, nHVF = nHVF,
          maxPC = maxPC, resolution = resolution, reduction = reduction,
          cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
          exogenous_genes = exogenous_genes
        )
        return(srt)
      })
    )

    invisible(lapply(samples, function(x) {
      saveRDS(get(paste0("sc_list_filter_", nm))[[x]], paste0("Normalization-", nm, "/", x, ".rds"))
    }))
    saveRDS(get(paste0("sc_list_filter_", nm)), file = paste0("sc_list_filter_", nm, ".rds"))
  }
}

# Integration ----------------------------------------------
if (length(datasets) != 0) {
  for (dataset in datasets) {
    for (nm in normalization_method) {
      cat("++++++ Use", nm, "normalized data to do integration ++++++", "\n")

      checked <- Check_scList(
        sc_list = get(paste0("sc_list_filter_", nm))[dataset], normalization_method = nm,
        HVF_source = HVF_source, nHVF = nHVF, hvf = NULL,
        exogenous_genes = exogenous_genes
      )
      sc_list <- checked[["sc_list"]]
      hvf <- checked[["hvf"]]

      for (im in integration_method) {
        if (im %in% c("Uncorrected", "Seurat", "fastMNN", "Harmony", "Scanorama", "BBKNN", "CSS", "LIGER", "scMerge")) {
          dir_path <- paste0("Normalization-", nm, "/", HVF_source, "_HVF/", "Integration-", im)
          dir.create(dir_path, recursive = T, showWarnings = FALSE)
          cat("++++++", paste0(dataset, collapse = ","), paste0("(Integration-", im, ")"), "++++++", "\n")
          if (file.exists(paste0(dir_path, "/", paste0(dataset, collapse = ","), ".rds"))) {
            cat(">>> Integration:", im, "process for the", paste0(dataset, collapse = ","), "has finished. Skip to the next step.\n")
          } else {
            srt_integrated <- Integration_SCP(
              sc_list = sc_list, integration_method = im,
              normalization_method = nm,
              HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
              maxPC = maxPC, resolution = resolution, reduction = reduction,
              cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
              exogenous_genes = exogenous_genes
            )
            p1 <- DimPlot(srt_integrated, group.by = c("orig.ident", "seurat_clusters", "cellcalling_method"), combine = FALSE)
            p1 <- lapply(p1, function(p) {
              p + theme(aspect.ratio = 1)
            })
            p2 <- FeaturePlot(object = srt_integrated, features = c(exogenous_genes, "percent.mt", "percent.ribo", "nCount_RNA", "nFeature_RNA", "scDblFinder.score", "cellcalling_methodNum"), combine = FALSE, order = TRUE)
            p2 <- lapply(p2, function(p) {
              p + theme(aspect.ratio = 1)
            })
            p <- cowplot::plot_grid(plotlist = c(p1, p2), align = "hv", axis = "tblr", ncol = 2, byrow = T)
            ggsave(
              plot = p, filename = paste0(dir_path, "/", paste0(dataset, collapse = ","), ".png"),
              units = "mm", width = 300, height = 70 * ceiling(length(c(p1, p2)) / 2),
              scale = 1.5, limitsize = FALSE
            )
            saveRDS(srt_integrated, file = paste0(dir_path, "/", paste0(dataset, collapse = ","), ".rds"))
            cat(">>> Integration:", im, "process for the", paste0(dataset, collapse = ","), "completed successfully.\n")
          }
        } else {
          cat("Warning!", im, "process skipped because it is not a supported integration method.\n")
        }
      }
    }
  }
} else {
  cat("Less than 2 datasets. Integration skipped.")
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
