# !/usr/bin/env Rscript

# Set parameters ---------------------------------------------------------
## parameters: global settings ---------------------------------------------
SCPwork_dir <- "/ssd/lab/HuangMingQian/scRNAseq/H0CellLine-project/NGSmodule_SCP_work/"
SCP_path <- "/home/zhanghao/Program/NGS/UniversalTools/NGSmodule/SingleCellPipe/"
SampleInfoFile <- "/ssd/lab/HuangMingQian/scRNAseq/H0CellLine-project/temp_20210929225849.Sample_info.csv"
Rscript_threads <- 120
Rscript_memory <- 500

## parameters: integration --------------------------------------------------
species <- "Homo_sapiens,Mus_musculus"  ## "Homo_sapiens,Mus_musculus", comma-separated.
normalization_method <- "logCPM"  ## "logCPM,SCT", comma-separated.
vars_to_regress <- NULL  ## "nCount_RNA,nFeature_RNA,percent.mt,percent.ribo" or NULL, comma-separated.
nHVF <- 3000  ## 3000. Number of high-variable features to use.
liner_reduction <- "pca"  ## One of "pca","ica","nmf","mds". Used to compute neighbors and nonlinear dimensionality reduction.
liner_reduction_dims <- 50  ## 50. Number of dimensions to calculate in nonlinear dimensionality reduction.
liner_reduction_dims_use <- NULL  ## A number or NULL (automatic). The top N dimensions used to compute neighbors and nonlinear dimensionality reduction.
liner_reduction_distance <- "cosine"  ## "cosine" or other method in parallelDist::parDist. The distance measure to be used in the MDS or other liner dimensionality reduction.
nonliner_reduction <- "umap"  ## "umap,tsne,dm", comma-separated. Methods for nonlinear dimensionality reduction.
nonliner_reduction_dims <- "2,3"  ## "2,3", comma-separated. Dimensional numbers in descending space.
nonliner_reduction_distance <- "euclidean"  ## "euclidean" or other method in parallelDist::parDist. The distance measure to be used in the DiffusionMap or other nonliner dimensionality reduction.
cluster_algorithm <- "louvain"  ## One of "louvain","slm" and "leiden". The algorithm used to identify clusters of cells.
cluster_resolution <- 0.8  ## The resolution used to identify clusters of cells.
exogenous_genes <- "GFP"  ## Gene symbol or NULL, comma-separated.
features_inspect <- "POU5F1,SOX2,SOX17,SOX9,EOMES,PAX6,PRDM1,SOX1,OTX2,nCount_RNA,nFeature_RNA,percent.mito,percent.ribo" ## Gene symbol or other qc features, comma-separated.

HVF_source <- "separate"  ## "separate,global", comma-separated. Source of high variable featues.
datasets_integration <- "all;"  ## "all;LLH_A,LLH_B;", comma-separated. Comma separates datasets in one integration, semicolon separates integrations.
datasets_names <- "all;"  ## "all;LLH_sets;", semicolon-separated. Name for each integration.
integration_method <- "Uncorrected,Seurat,fastMNN,Harmony,Scanorama,BBKNN,CSS,LIGER"  ## "Uncorrected,Seurat,fastMNN,Harmony,Scanorama,BBKNN,CSS,LIGER", comma-separated. Methods used for integration.

# Library -----------------------------------------------------------------
reticulate::py_run_string("import matplotlib.pyplot as plt")
reticulate::import("scanpy")
suppressWarnings(suppressPackageStartupMessages(invisible(lapply(
  c(
    "Seurat", "SeuratDisk", "reticulate", "dplyr", "future", "htmlwidgets", "ggplot2"
  ),
  require,
  character.only = TRUE
))))
source(paste0(SCP_path, "/SCP-workflow-funtcion.R"))
source(paste0(SCP_path, "/SCP-plot.R"))
source_python(paste0(SCP_path, "/SCP-workflow-funtcion.py"))

main_dir <- paste0(SCPwork_dir, "/../NGSmodule_SCP_analysis/")
dir.create(main_dir, recursive = T, showWarnings = FALSE)
setwd(main_dir)

# Prepare data -----------------------------------------------------------------
species <- strsplit(species, split = ",") %>% unlist()
normalization_method <- strsplit(normalization_method, split = ",") %>% unlist()
if (!is.null(vars_to_regress)) {
  vars_to_regress <- strsplit(vars_to_regress, split = ",") %>% unlist()
}
liner_reduction <- strsplit(liner_reduction, split = ",") %>% unlist()
if (length(liner_reduction) > 1) {
  warning("Only the first method in the 'liner_reduction' will be used.")
  liner_reduction <- liner_reduction[1]
}
nonliner_reduction <- strsplit(nonliner_reduction, split = ",") %>% unlist()
nonliner_reduction_dims <- strsplit(nonliner_reduction_dims, split = ",") %>%
  unlist() %>%
  as.numeric()
if (!is.null(exogenous_genes)) {
  exogenous_genes <- strsplit(exogenous_genes, split = ",") %>% unlist()
}
features_inspect <- c(exogenous_genes, strsplit(features_inspect, split = ",") %>% unlist())

HVF_source <- strsplit(HVF_source, split = ",") %>% unlist()
datasets <- strsplit(datasets_integration, split = ";") %>%
  unlist() %>%
  lapply(., function(x) {
    if (x == "all") {
      out <- as.character(list.dirs(path = SCPwork_dir, recursive = FALSE, full.names = FALSE))
    } else {
      out <- strsplit(x, split = ",") %>% unlist()
    }
    return(out)
  })
datasets_names <- strsplit(datasets_names, split = ";") %>% unlist()
if (length(datasets) != length(datasets_names)) {
  stop("length of datasets_integration is not same as datasets_names")
}
integration_method <- strsplit(integration_method, split = ",") %>% unlist()

samples <- unlist(datasets) %>% unique()
samples <- factor(samples, levels = unique(samples))
samleinfo <- read.csv(SampleInfoFile, header = T, stringsAsFactors = F)
samleinfo <- unique(samleinfo[, c("SampleID", "Group")])
rownames(samleinfo) <- samleinfo$SampleID
samleinfo <- samleinfo[as.character(samples), ]
groups <- factor(samleinfo[["Group"]], unique(samleinfo$Group))
names(groups) <- samleinfo[["SampleID"]]

# Prepare data -----------------------------------------------------------------
CCgenes <- CC_GenePrefetch(species = species[1])
cc_S_genes <- CCgenes[["cc_S_genes"]]
cc_G2M_genes <- CCgenes[["cc_G2M_genes"]]

# Global options ------------------------------------------------
set.seed(11)
options(expressions = 5e5)
options(future.globals.maxSize = Rscript_memory * 1024^3)
options(future.fork.enable = TRUE)
if (Rscript_threads >= 125) {
  cat("Rscript_threads number is too large. Re-set it to the 125.\n")
  Rscript_threads <- 125
}
plan("multicore", workers = Rscript_threads, gc = TRUE)

# Data Loading ------------------------------------------------
srtList <- list()
for (i in 1:length(samples)) {
  cat("++++++", as.character(samples[i]), "(Data Loading)", "++++++", "\n")
  for (nm in normalization_method) {
    srtList[[paste0(as.character(samples[i]), "-", nm)]] <- LoadH5Seurat(paste0(main_dir, "/Standard/Normalization-", nm, "/", as.character(samples[i]), ".h5Seurat"))
  }
}

# Integration SCP----------------------------------------------------------------
for (nm in normalization_method) {
  cat(">>> Normalization method:", nm, "\n")
  for (hvfsc in HVF_source) {
    cat(">>> HVF_source:", hvfsc, "\n")
    for (i in 1:length(datasets)) {
      dataset <- datasets[[i]]
      dataset_name <- datasets_names[[i]]
      cat("++++++ ", dataset_name, ": ", paste0(dataset, collapse = ","), " ++++++\n", sep = "")
      dir <- paste0("Integration/Normalization-", nm, "/", hvfsc, "HVF/")
      dir.create(paste0(dir, "/Plot"), recursive = T, showWarnings = FALSE)
      srtList <- srtList[paste0(as.character(dataset), "-", nm)]

      if (file.exists(paste0(dir, "/", dataset_name, ".h5Seurat"))) {
        cat("Loading the existed integration data... ...\n")
        srtMerge <- LoadH5Seurat(paste0(dir, "/", dataset_name, ".h5Seurat"))
        hvf <- srtMerge@misc$integration_HVF
      } else {
        cat("Preparing the integration data... ...\n")
        checked <- Check_srtList(
          srtList = srtList, batch = "orig.ident",
          do_normalization = TRUE, do_HVF_finding = TRUE,
          normalization_method = nm, vars_to_regress = vars_to_regress,
          HVF_source = hvfsc, nHVF = nHVF, hvf = NULL,
          exogenous_genes = exogenous_genes
        )
        srtList <- checked[["srtList"]]
        srtMerge <- Reduce(function(x, y) merge(x, y), checked[["srtList"]])
        srtMerge@misc$integration_HVF <- hvf <- checked[["hvf"]]
        srtMerge$orig.ident <- factor(srtMerge$orig.ident, levels = as.character(samples)[as.character(samples) %in% as.character(srtMerge$orig.ident)])
        srtMerge$Samples <- srtMerge$orig.ident
        srtMerge$Groups <- factor(groups[as.character(srtMerge$orig.ident)], levels = unique(groups[names(groups) %in% as.character(srtMerge$orig.ident)]))
        srtMerge <- RunDEtest(srt = srtMerge, group_by = "Groups")
        srtMerge <- Standard_SCP(srtMerge,hvf = hvf)
        SaveH5Seurat(
          object = srtMerge,
          filename = paste0(dir, "/", dataset_name, ".h5Seurat"),
          overwrite = TRUE,
          verbose = FALSE
        )
      }

      for (im in integration_method) {
        cat(">>>", dataset_name, paste0("(", nm, " normalized", "|", hvfsc, " HVF", "|", im, " integration)"), "\n")
        assays_current <- Seurat::Assays(srtMerge)
        reductions_current <- Seurat::Reductions(srtMerge)
        if (!grepl(pattern = im, x = names(srtMerge@reductions))) {
          srtMerge <- Integration_SCP(
            srtList = srtList, srtMerge = srtMerge, append = TRUE, batch = "orig.ident",
            integration_method = im,
            do_normalization = NULL, normalization_method = nm,
            do_HVF_finding = TRUE, HVF_source = hvfsc, nHVF = nHVF, hvf = hvf,
            do_scaling = TRUE, vars_to_regress = vars_to_regress,
            liner_reduction = liner_reduction, liner_reduction_dims = liner_reduction_dims, liner_reduction_distance = liner_reduction_distance, liner_reduction_dims_use = liner_reduction_dims_use,
            nonliner_reduction = nonliner_reduction, nonliner_reduction_dims = nonliner_reduction_dims, nonliner_reduction_distance = nonliner_reduction_distance,
            cluster_algorithm = cluster_algorithm, cluster_resolution = cluster_resolution, cluster_reorder = TRUE,
            exogenous_genes = exogenous_genes, seed = 11
          )
          srtMerge$orig.ident <- factor(srtMerge$orig.ident, levels = samples[as.character(samples) %in% as.character(srtMerge$orig.ident)])
          srtMerge$Samples <- srtMerge$orig.ident
          srtMerge$Groups <- factor(groups[as.character(srtMerge$orig.ident)], levels = unique(groups[names(groups) %in% as.character(srtMerge$orig.ident)]))
        }

        for (nr in nonliner_reduction) {
          nr <- paste0(toupper(nr))
          features <- features_inspect[features_inspect %in% c(colnames(srtMerge@meta.data), rownames(srtMerge))]
          p <- SummaryPlot(
            srt = srtMerge, clusters = paste0(im, "clusters"), groups = "Groups",
            reduction = paste0(im, nr, "2d"),
            features = features
          )
          ggsave(
            plot = p, filename = paste0(dir, "/Plot/", dataset_name, ".", im, ".", nr, "2d.png"),
            units = "mm", height = 200, width = 80 * max(c(2 + ceiling(length(features) / 2), length(unique(srtMerge[["Groups", drop = T]])))),
            scale = 1.5, limitsize = FALSE
          )
          if (3 %in% nonliner_reduction_dims) {
            p <- ClassDimPlot3D(srt = srtMerge, group.by = paste0(im, "clusters"), reduction = paste0(im, nr, "3d"))
            filepath <- paste0(dir, "/Plot/", dataset_name, ".", im, ".", nr, "3d.html")
            saveWidget(
              widget = plotly::as_widget(p),
              file = filepath
            )
            unlink(stringr::str_replace(filepath, "\\.html", "_files"), recursive = TRUE)
          }

          if (length(cc_S_genes) >= 3 & length(cc_G2M_genes) >= 3) {
            srtMerge <- CellCycleScoring(
              object = srtMerge,
              s.features = cc_S_genes,
              g2m.features = cc_G2M_genes,
              set.ident = FALSE
            )
            srtMerge[["CC.Difference"]] <- srtMerge[["S.Score"]] - srtMerge[["G2M.Score"]]
            srtMerge[["Phase"]] <- factor(srtMerge[["Phase", drop = TRUE]], levels = c("G1", "S", "G2M"))
            p1 <- ClassDimPlot(srtMerge,
              group.by = "Phase",
              reduction = paste0(im, nr, "2d"),
              palette = "Set1", combine = FALSE
            )
            p2 <- ExpDimPlot(srtMerge,
              features = c("S.Score", "G2M.Score"),
              reduction = paste0(im, nr, "2d"),
              combine = FALSE
            )
            p <- cowplot::plot_grid(plotlist = c(p1, p2), nrow = 1, align = "hv", axis = "tblr")
            ggsave(
              plot = p, width = 40, height = 20, units = "cm",
              filename = paste0(dir, "/Plot/", dataset_name, ".", im, ".", nr, "2d.CellCycle.png"),
              limitsize = FALSE
            )
          }
        }

        srtMerge <- RunDEtest(srt = srtMerge, group_by = paste0(im, "clusters"))
        de <- srtMerge@misc[[paste0("DEtest_", im, "clusters")]][["AllMarkers_Wilcoxon"]]
        df1 <- de %>%
          filter(p_val_adj < 0.05) %>%
          group_by(gene) %>%
          top_n(1, avg_log2FC) %>%
          group_by(group1) %>%
          top_n(3, avg_log2FC)
        p <- ExpDotPlot(srt = srtMerge, genes = df1[["gene"]], gene_groups = df1[["group1"]], columns = paste0(im, "clusters"))
        pngunit <- length(unique(srtMerge[[paste0(im, "clusters"), drop = T]]))
        ggsave(
          plot = p, width = 5 + pngunit, height = 5 + pngunit * 3, units = "cm",
          filename = paste0(dir, "/Plot/", dataset_name, ".DEtest_", im, "clusters", ".dotplot.png"), limitsize = FALSE
        )
        df2 <- de %>%
          filter(p_val_adj < 0.05)
        p <- ExpHeatmapPlot(srt = srtMerge, genes = df2[["gene"]], gene_groups = df2[["group1"]], columns = paste0(im, "clusters"))
        ggsave(
          plot = p, width = 30, height = 20, units = "cm",
          filename = paste0(dir, "/Plot/", dataset_name, ".DEtest_", im, "clusters", ".heatmap.png"), limitsize = FALSE
        )

        if (!all(c("spliced", "unspliced") %in% Seurat::Assays(srtMerge))) {
          cat("Add the RNA velocity data...\n")
          velocityList <- list()
          for (i in 1:length(samples)) {
            velocityList[[as.character(samples[i])]] <- LoadH5Seurat(paste0("CellQC/Velocity/", as.character(samples[i]), ".filtered.velocity.h5Seurat"), verbose = FALSE)
          }
          velocityMerge <- Reduce(merge, x = velocityList)
          if (!(identical(rownames(srtMerge), rownames(velocityMerge)) &
            identical(colnames(srtMerge), colnames(velocityMerge)))) {
            for (assay in Seurat::Assays(velocityMerge)) {
              velocityMerge[[assay]] <- Seurat::CreateAssayObject(counts = velocityMerge[[assay]][rownames(srtMerge), colnames(srtMerge)])
              velocityMerge[[assay]]@key <- paste0(assay, "_")
            }
          }
          srtMerge@assays$spliced <- velocityMerge@assays$spliced
          srtMerge@assays$unspliced <- velocityMerge@assays$unspliced
          srtMerge$nCount_spliced <- velocityMerge$nCount_spliced
          srtMerge$nFeature_spliced <- velocityMerge$nFeature_spliced
          srtMerge$nCount_unspliced <- velocityMerge$nCount_unspliced
          srtMerge$nFeature_unspliced <- velocityMerge$nFeature_unspliced
        }

        adata <- srt_to_adata(srtMerge)
        for (nr in nonliner_reduction) {
          nr <- paste0(toupper(nr))
          adata <- RunSCVELO(
            adata = adata,
            group_by = paste0(im, "clusters"),
            liner_reduction = "Standardpca",
            nonliner_reduction = paste0(im, nr, "2d"),
            recover_dynamics = F,
            dirpath = paste0(dir, "/Plot/"),
            fileprefix = paste0(dataset_name, ".", im, ".", nr, "2d")
          )
          adata <- RunPAGA(
            adata = adata,
            group_by = paste0(im, "clusters"),
            liner_reduction = "Standardpca",
            nonliner_reduction = paste0(im, nr, "2d"),
            dirpath = paste0(dir, "/Plot/"),
            fileprefix = paste0(dataset_name, ".", im, ".", nr, "2d")
          )
        }

        cat("Save the final Seurat object to h5Seurat...\n")
        SaveH5Seurat(
          object = srtMerge,
          filename = paste0(dir, "/", dataset_name, ".tmp.h5Seurat"),
          overwrite = TRUE,
          verbose = FALSE
        )
        file.rename(
          from = paste0(dir, "/", dataset_name, ".tmp.h5Seurat"),
          to = paste0(dir, "/", dataset_name, ".h5Seurat")
        )

        # ## SeuratDisk is not completed.
        # assays_append <- setdiff(Seurat::Assays(srtMerge), assays_current)
        # reductions_append <- setdiff(Seurat::Reductions(srtMerge), reductions_current)
        # AppendData(
        #   file = paste0(dir, "/", paste0(dataset, collapse = ","), ".h5Seurat"),
        #   object = srtMerge,
        #   assays = assays_append,
        #   reductions  = reductions_append,
        #   verbose = FALSE
        # )
      }
    }
  }
}

# The end -----------------------------------------------------------------
future:::ClusterRegistry("stop")
