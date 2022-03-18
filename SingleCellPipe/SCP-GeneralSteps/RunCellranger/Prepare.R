#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
SCPwork_dir <- as.character(args[1])

# Set parameters ---------------------------------------------------------
## parameters: global settings ---------------------------------------------
# SCPwork_dir <- "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/NGSmodule_SCP_work/"

# Library -----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(invisible(lapply(c(
  "SCP", "Seurat", "SeuratDisk", "SeuratWrappers", "dplyr"
),
require,
character.only = TRUE
))))

main_dir <- paste0(SCPwork_dir, "/../NGSmodule_SCP_analysis/Prepare")
dir.create(main_dir, recursive = T, showWarnings = FALSE)
setwd(main_dir)

# Create Seurat object ------------------------------------------------
samples <- list.dirs(path = SCPwork_dir, recursive = FALSE, full.names = FALSE)
for (i in 1:length(samples)) {
  cat("++++++", samples[i], "++++++", "\n")
  if (!file.exists(paste0("./", samples[i], ".h5Seurat"))) {
    cat("Load the RNA expression data...\n")
    srt_matrix <- Read10X(data.dir = paste0(SCPwork_dir, "/", samples[i], "/Alignment-Cellranger/", samples[i], "/outs/filtered_feature_bc_matrix/"))
    srt <- CreateSeuratObject(counts = srt_matrix, project = samples[i])
    srt <- RenameCells(object = srt, add.cell.id = samples[i])

    cat("Load the RNA velocity data...\n")
    velocity <- ReadVelocity(file = paste0(SCPwork_dir, "/", samples[i], "/Alignment-Cellranger/", samples[i], "/velocyto/", samples[i], ".loom"), verbose = FALSE)
    velocity <- as.Seurat(x = velocity, project = samples[i], verbose = FALSE)
    velocity <- RenameCells(
      object = velocity,
      new.names = gsub(x = colnames(velocity), pattern = ".*:", replacement = "", perl = TRUE) %>%
        gsub(x = ., pattern = "x$", replacement = "-1", perl = TRUE) %>%
        paste0(samples[i], "_", .)
    )

    if (!(identical(rownames(srt), rownames(velocity)) &
      identical(colnames(srt), colnames(velocity)))) {
      message("Rearrange cell or gene order of the velocity data according to the RNA slot")
      for (assay in Assays(velocity)) {
        velocity[[assay]] <- CreateAssayObject(counts = velocity[[assay]][rownames(srt), colnames(srt)])
        velocity[[assay]]@key <- paste0(assay, "_")
      }
    }
    srt[["spliced"]] <- velocity[["spliced"]]
    srt[["unspliced"]] <- velocity[["unspliced"]]
    srt[["nCount_spliced"]] <- velocity[["nCount_spliced"]]
    srt[["nFeature_spliced"]] <- velocity[["nFeature_spliced"]]
    srt[["nCount_unspliced"]] <- velocity[["nCount_unspliced"]]
    srt[["nFeature_unspliced"]] <- velocity[["nFeature_unspliced"]]

    cat("Save the Seurat object as h5Seurat...\n")
    SaveH5Seurat(
      object = srt,
      filename = paste0("./", samples[i], ".h5Seurat"),
      overwrite = TRUE,
      verbose = FALSE
    )
  } else {
    cat(samples[i], "Seurat object is ready.\n")
  }
}