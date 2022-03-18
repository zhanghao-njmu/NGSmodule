# !/usr/bin/env Rscript

# Set parameters ---------------------------------------------------------
## parameters: global settings ---------------------------------------------
SCPwork_dir <- "/ssd/lab/HuangMingQian/scRNAseq/H0CellLine-project/NGSmodule_SCP_work/"
SCP_path <- "/home/zhanghao/Program/NGS/UniversalTools/NGSmodule/SingleCellPipe/"
Rscript_threads <- 120
Rscript_memory <- 500

## parameters: standard --------------------------------------------------
species <- "Homo_sapiens,Mus_musculus" ## "Homo_sapiens,Mus_musculus", comma-separated.
normalization_method <- "logCPM" ## "logCPM,SCT", comma-separated.
vars_to_regress <- NULL ## "nCount_RNA,nFeature_RNA,percent.mt,percent.ribo" or NULL, comma-separated.
nHVF <- 3000 ## 3000. Number of high-variable features to use.
liner_reduction <- "pca" ## "pca,ica,nmf,mds", comma-separated. Used to compute neighbors and nonlinear dimensionality reduction.
liner_reduction_dims <- 50 ## 50. Number of dimensions to calculate in nonlinear dimensionality reduction.
liner_reduction_dims_use <- NULL ## A number or NULL (automatic). The top N dimensions used to compute neighbors and nonlinear dimensionality reduction.
liner_reduction_distance <- "cosine" ## "cosine" or other method in parallelDist::parDist. The distance measure to be used in the MDS or other liner dimensionality reduction.
nonliner_reduction <- "umap" ## "umap,tsne,dm", comma-separated. Methods for nonlinear dimensionality reduction.
nonliner_reduction_dims <- "2,3" ## "2,3", comma-separated. Dimensional numbers in descending space.
nonliner_reduction_distance <- "euclidean" ## "euclidean" or other method in parallelDist::parDist. The distance measure to be used in the DiffusionMap or other nonliner dimensionality reduction.
cluster_algorithm <- "louvain" ## One of "louvain","slm" and "leiden". The algorithm used to identify clusters of cells.
cluster_resolution <- 0.8 ## The resolution used to identify clusters of cells.
exogenous_genes <- "GFP" ## Gene symbol or NULL, comma-separated.
features_inspect <- "POU5F1,SOX2,SOX17,SOX9,EOMES,PAX6,PRDM1,SOX1,OTX2,nCount_RNA,nFeature_RNA,percent.mito,percent.ribo" ## Gene symbol or other qc features, comma-separated.

# Library -----------------------------------------------------------------
reticulate::py_run_string("import matplotlib.pyplot as plt")
reticulate::import("scanpy")
suppressWarnings(suppressPackageStartupMessages(invisible(lapply(
  c(
    "Seurat", "SeuratDisk", "future", "reticulate", "dplyr", "htmlwidgets"
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

species <- strsplit(species, split = ",") %>% unlist()
normalization_method <- strsplit(normalization_method, split = ",") %>% unlist()
if (!is.null(vars_to_regress)) {
  vars_to_regress <- strsplit(vars_to_regress, split = ",") %>% unlist()
}
liner_reduction <- strsplit(liner_reduction, split = ",") %>% unlist()
nonliner_reduction <- strsplit(nonliner_reduction, split = ",") %>% unlist()
if (!is.null(exogenous_genes)) {
  exogenous_genes <- strsplit(exogenous_genes, split = ",") %>% unlist()
}
nonliner_reduction_dims <- strsplit(nonliner_reduction_dims, split = ",") %>%
  unlist() %>%
  as.numeric()
features_inspect <- c(exogenous_genes, strsplit(features_inspect, split = ",") %>% unlist())

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
samples <- list.dirs(path = SCPwork_dir, recursive = FALSE, full.names = FALSE)
srtList <- list()
for (i in 1:length(samples)) {
  cat("++++++", samples[i], "(Data Loading)", "++++++", "\n")
  srtList[[samples[i]]] <- LoadH5Seurat(paste0(main_dir, "/CellQC/RNA/", samples[i], ".filtered.h5Seurat"))
}

# Standard SCP ------------------------------------------------------------
for (nm in normalization_method) {
  dir <- paste0(main_dir, "/Standard/Normalization-", nm, "/")
  dir.create(paste0(dir, "/Plot"), recursive = T, showWarnings = FALSE)
  for (i in 1:length(samples)) {
    sample_current <- samples[i]
    if (file.exists(paste0(dir, sample_current, ".", nm, ".h5Seurat"))) {
      cat("Loading", sample_current, "from the h5Seurat....\n")
      srt <- LoadH5Seurat(paste0(dir, sample_current, ".", nm, ".h5Seurat"))
    } else {
      cat("++++++", sample_current, nm, "Normalization", "++++++", "\n")
      srt <- srtList[[sample_current]]
      srt <- Standard_SCP(
        srt = srt, prefix = "Standard",
        do_normalization = TRUE, normalization_method = normalization_method,
        do_HVF_finding = TRUE, nHVF = nHVF, hvf = NULL,
        do_scaling = TRUE, vars_to_regress = vars_to_regress,
        liner_reduction = liner_reduction, liner_reduction_dims = liner_reduction_dims, liner_reduction_distance = liner_reduction_distance, liner_reduction_dims_use = liner_reduction_dims_use,
        nonliner_reduction = nonliner_reduction, nonliner_reduction_dims = nonliner_reduction_dims, nonliner_reduction_distance = nonliner_reduction_distance,
        cluster_algorithm = cluster_algorithm, cluster_resolution = cluster_resolution, cluster_reorder = TRUE,
        exogenous_genes = exogenous_genes, seed = 11
      )
      for (lr in liner_reduction) {
        for (nr in nonliner_reduction) {
          nr <- paste0(toupper(nr))
          features <- features_inspect[features_inspect %in% c(colnames(srt@meta.data), rownames(srt))]
          p <- SummaryPlot(
            srt = srt, clusters = paste0("Standard", lr, "clusters"), groups = "orig.ident",
            reduction = paste0("Standard", lr, nr, "2d"),
            features = features
          )
          ggsave(
            plot = p, filename = paste0(dir, "/Plot/", sample_current, ".Standard.", lr, nr, "2d.png"),
            units = "mm", height = 200, width = 80 * max(c(2 + ceiling(length(features) / 2), length(unique(srt[["orig.ident", drop = T]])))),
            scale = 1.5, limitsize = FALSE
          )
          if (3 %in% nonliner_reduction_dims) {
            p <- ClassDimPlot3D(srt = srt, group.by = paste0("Standard", lr, "clusters"), reduction = paste0("Standard", lr, nr, "3d"))
            filepath <- paste0(dir, "/Plot/", sample_current, ".Standard.", lr, nr, "3d.html")
            saveWidget(
              widget = plotly::as_widget(p),
              file = filepath
            )
            unlink(stringr::str_replace(filepath, "\\.html", "_files"), recursive = TRUE)
          }

          if (length(cc_S_genes) >= 3 & length(cc_G2M_genes) >= 3) {
            srt <- CellCycleScoring(
              object = srt,
              s.features = cc_S_genes,
              g2m.features = cc_G2M_genes,
              set.ident = FALSE
            )
            srt[["CC.Difference"]] <- srt[["S.Score"]] - srt[["G2M.Score"]]
            srt[["Phase"]] <- factor(srt[["Phase", drop = TRUE]], levels = c("G1", "S", "G2M"))
            p1 <- ClassDimPlot(srt,
              group.by = "Phase",
              reduction = paste0("Standard", lr, nr, "2d"),
              palette = "Set1", combine = FALSE
            )
            p2 <- ExpDimPlot(srt,
              features = c("S.Score", "G2M.Score"),
              reduction = paste0("Standard", lr, nr, "2d"),
              combine = FALSE
            )
            p <- cowplot::plot_grid(plotlist = c(p1, p2), nrow = 1, align = "hv", axis = "tblr")
            ggsave(
              plot = p, width = 40, height = 20, units = "cm",
              filename = paste0(dir, "/Plot/", sample_current, ".Standard.", lr, nr, "2d.CellCycle.png"),
              limitsize = FALSE
            )
          }
        }

        srt <- RunDEtest(srt = srt, group_by = paste0("Standard", lr, "clusters"))
        de <- srt@misc[[paste0("DEtest_", paste0("Standard", lr, "clusters"))]][["AllMarkers_Wilcoxon"]]
        df1 <- de %>%
          filter(p_val_adj < 0.05) %>%
          group_by(gene) %>%
          top_n(1, avg_log2FC) %>%
          group_by(group1) %>%
          top_n(3, avg_log2FC)
        p <- ExpDotPlot(srt = srt, genes = df1[["gene"]], gene_groups = df1[["group1"]], columns = paste0("Standard", lr, "clusters"))
        pngunit <- length(unique(srt[[paste0("Standard", lr, "clusters"), drop = T]]))
        ggsave(
          plot = p, width = 5 + pngunit, height = 5 + pngunit * 3, units = "cm",
          filename = paste0(dir, "/Plot/", sample_current, ".DEtest_Standard", lr, "clusters.dotplot.png"), limitsize = FALSE
        )
        df2 <- de %>%
          filter(p_val_adj < 0.05)
        p <- ExpHeatmapPlot(srt = srt, genes = df2[["gene"]], gene_groups = df2[["group1"]], columns = paste0("Standard", lr, "clusters"))
        ggsave(
          plot = p, width = 30, height = 20, units = "cm",
          filename = paste0(dir, "/Plot/", sample_current, ".DEtest_Standard", lr, "clusters.heatmap.png"), limitsize = FALSE
        )
      }

      cat("Add the RNA velocity data...\n")
      velocity <- LoadH5Seurat(paste0("CellQC/Velocity/", sample_current, ".filtered.velocity.h5Seurat"), verbose = FALSE)
      if (!(identical(rownames(srt), rownames(velocity)) &
        identical(colnames(srt), colnames(velocity)))) {
        for (assay in Seurat::Assays(velocity)) {
          velocity[[assay]] <- Seurat::CreateAssayObject(counts = velocity[[assay]][rownames(srt), colnames(srt)])
          velocity[[assay]]@key <- paste0(assay, "_")
        }
      }
      srt@assays$spliced <- velocity@assays$spliced
      srt@assays$unspliced <- velocity@assays$unspliced
      srt$nCount_spliced <- velocity$nCount_spliced
      srt$nFeature_spliced <- velocity$nFeature_spliced
      srt$nCount_unspliced <- velocity$nCount_unspliced
      srt$nFeature_unspliced <- velocity$nFeature_unspliced

      adata <- srt_to_adata(srt)
      for (lr in liner_reduction) {
        for (nr in nonliner_reduction) {
          nr <- paste0(toupper(nr))
          adata <- RunSCVELO(
            adata = adata,
            group_by = paste0("Standard", lr, "clusters"),
            liner_reduction = paste0("Standard", lr),
            nonliner_reduction = paste0("Standard", lr, nr, "2d"),
            recover_dynamics = F,
            dirpath = paste0(dir, "/Plot/"),
            fileprefix = paste0(sample_current, ".Standard.", lr, nr, "2d")
          )
          adata <- RunPAGA(
            adata = adata,
            group_by = paste0("Standard", lr, "clusters"),
            liner_reduction = paste0("Standard", lr),
            nonliner_reduction = paste0("Standard", lr, nr, "2d"),
            dirpath = paste0(dir, "/Plot/"),
            fileprefix = paste0(sample_current, ".Standard.", lr, nr, "2d")
          )
        }
      }

      cat("Save the final Seurat object to h5Seurat...\n")
      SaveH5Seurat(
        object = srt,
        filename = paste0(dir, "/", sample_current, ".h5Seurat"),
        overwrite = TRUE,
        verbose = FALSE
      )
    }
  }
}

# The end -----------------------------------------------------------------
future:::ClusterRegistry("stop")
