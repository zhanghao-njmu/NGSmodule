#!/usr/bin/env Rscript
# Set parameters ---------------------------------------------------------
## parameters: global settings ---------------------------------------------
SCPwork_dir <- "/ssd/lab/HuangMingQian/scRNAseq/H0CellLine-project/NGSmodule_SCP_work/"
species <- "Homo_sapiens,Mus_musculus"
species_gene_prefix <- "GRCh38-,mm10-"

# Library -----------------------------------------------------------------
reticulate::py_run_string("import matplotlib.pyplot as plt")
reticulate::import("scanpy")
suppressWarnings(suppressPackageStartupMessages(invisible(lapply(c(
  "Seurat", "SeuratDisk", "dplyr", "SeuratWrappers"
),
require,
character.only = TRUE
))))

main_dir <- paste0(SCPwork_dir, "/../NGSmodule_SCP_analysis/Prepare")
dir.create(main_dir, recursive = T, showWarnings = FALSE)
setwd(main_dir)

dir.create(paste0(main_dir, "/RNA"), recursive = T, showWarnings = FALSE)
dir.create(paste0(main_dir, "/Velocity"), recursive = T, showWarnings = FALSE)

species <- strsplit(species, split = ",") %>% unlist()
species_gene_prefix <- strsplit(species_gene_prefix, split = ",") %>% unlist()

# Create Seurat object ------------------------------------------------
samples <- list.dirs(path = SCPwork_dir, recursive = FALSE, full.names = FALSE)
for (i in 1:length(samples)) {
  cat("++++++", samples[i], "++++++", "\n")
  cat("Load the RNA expression data...\n")
  srt_matrix <- Read10X(data.dir = paste0(SCPwork_dir, "/", samples[i], "/Alignment-Cellranger/", samples[i], "/outs/filtered_feature_bc_matrix/"))
  srt <- CreateSeuratObject(counts = srt_matrix, project = samples[i])
  srt <- RenameCells(object = srt, add.cell.id = samples[i])
  for (n in 1:length(species)) {
    sp <- species[n]
    prefix <- species_gene_prefix[n]
    sp_genes <- rownames(srt)[grep(pattern = paste0("^", prefix), x = rownames(srt))]
    srt[[paste0(c("nCount_RNA", sp), collapse = ".")]] <- colSums(srt@assays$RNA@counts[sp_genes, ])
    srt[[paste0(c("nFeature_RNA", sp), collapse = ".")]] <- colSums(srt@assays$RNA@counts[sp_genes, ] > 0)
    srt[[paste0(c("percent.mito", sp), collapse = ".")]] <- PercentageFeatureSet(object = srt, pattern = paste0("(", paste0("^", prefix, "-*", c("MT-", "Mt-", "mt-")), ")", collapse = "|"))
    srt[[paste0(c("percent.ribo", sp), collapse = ".")]] <- PercentageFeatureSet(object = srt, pattern = paste0("(", paste0("^", prefix, "-*", c("RP[SL]\\d+(\\w+)$", "Rp[sl]\\d+(\\w+)$", "rp[sl]\\d+(\\w+)$")), ")", collapse = "|"))
    srt[[paste0(c("percent.genome", sp), collapse = ".")]] <- PercentageFeatureSet(object = srt, pattern = paste0("^", prefix))
  }

  cat("Save the RNA expression Seurat object to h5Seurat...\n")
  SaveH5Seurat(
    object = srt,
    filename = paste0("RNA/", samples[i], ".h5Seurat"),
    overwrite = TRUE,
    verbose = FALSE
  )

  cat("Load the RNA velocity data...\n")
  velocity <- ReadVelocity(file = paste0(SCPwork_dir, "/", samples[i], "/Alignment-Cellranger/", samples[i], "/velocyto/", samples[i], ".loom"), verbose = FALSE)
  velocity <- as.Seurat(x = velocity, project = samples[i], verbose = FALSE)
  velocity <- RenameCells(
    object = velocity,
    new.names = gsub(x = colnames(velocity), pattern = ".*:", replacement = "", perl = TRUE) %>%
      gsub(x = ., pattern = "x$", replacement = "-1", perl = TRUE) %>%
      paste0(samples[i], "_", .)
  )

  cat("Save the RNA velocity Seurat object to h5Seurat...\n")
  SaveH5Seurat(
    object = velocity,
    filename = paste0("Velocity/", samples[i], ".h5Seurat"),
    overwrite = TRUE,
    verbose = FALSE
  )
}
