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


scv$pl$scatter(adata, basis=names(topgenes_vals)[1:5], ncols=5, frameon=FALSE)











srt_logCPM <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/iPSC-ESC-iMeLC-PGCd6-CellLine-Testis/NGSmodule_SCP_analysis/Integration/Normalization-logCPM/separate_HVF/Testis-d10,Testis-d20,Testis-d30,Testis-d40,Testis-d50,Testis-d70.rds")
velocity_list <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/iPSC-ESC-iMeLC-PGCd6-CellLine-Testis/NGSmodule_SCP_analysis/Integration/velocity_list.rds")
ldat <- Reduce(function(x, y) merge(x, y), velocity_list)
ldat2<- subset(ldat,cell=colnames(srt_logCPM))
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





