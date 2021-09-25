#!/usr/bin/env Rscript

setwd("/ssd/lab/HuangMingQian/scRNAseq/22ESC-project/NGSmodule_SCP_analysis/Prepare/")
# parameters: global settings ---------------------------------------------
SCP_path <- "/home/zhanghao/Program/NGS/UniversalTools/NGSmodule/SingleCellPipe/"
SCPwork_dir <- "/ssd/lab/HuangMingQian/scRNAseq/22ESC-project/NGSmodule_SCP_work/"

# Library -----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(invisible(lapply(
  c(
    "Seurat", "SeuratDisk", "SeuratWrappers", "dplyr", "ggplot2", "ggplotify", "aplot", "cowplot",
    "reshape2", "stringr", "scuttle","RColorBrewer","velocyto.R"
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
color_sample <- brewer.pal(11, "Paired")[c(1:4,7:8,5:6,9:11)]
names(color_sample) <- samples
color_cluster <- colorRampPalette(brewer.pal(12, "Paired"))(nlevels(srt[["Uncorrectedclusters", drop = TRUE]]))
names(color_cluster) <- 1:nlevels(srt[["Uncorrectedclusters", drop = TRUE]])

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
  cat("Loading",samples[i],"RNA-velocity data...\n")
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
RenameFeatures <- function(srt,newnames=NULL){
  for (i in Seurat::Assays(srt)) {
    assay <- GetAssay(srt,i)
    for (d in c("counts","data","scale.data","meta.features")) {
      if (dim(slot(assay,d))[1]==length(newnames)) {
        rownames(slot(assay,d)) <- newnames
      }else{
        message(paste0("Slot ",d," have a different number of features."))
      }
    }
    srt@assays[[i]] <- assay
  }
  return(srt)
}
velocityMerge <- Reduce(merge,x = velocityList)

srt <- LoadH5Seurat("/ssd/lab/HuangMingQian/scRNAseq/22ESC-project/NGSmodule_SCP_analysis/Seurat/srtMerge2.h5Seurat")
human_gene <- rownames(velocityMerge)[grep(pattern = "^GRCh38-", x = rownames(velocityMerge))]
velocityMerge <- subset(x=velocityMerge, cells=colnames(srt), features=human_gene)
velocityMerge <- RenameFeatures(velocityMerge,newnames = gsub(x = rownames(velocityMerge), pattern = "GRCh38-", replacement = ""))
velocityMerge@reductions <- srt@reductions
velocityMerge@meta.data <- srt@meta.data
velocityMerge[["Uncorrectedclusters_colors"]]<- color_cluster[srt[["Uncorrectedclusters",drop=T]]]
velocityMerge[["Sample_colors"]]<- color_sample[srt[["orig.ident",drop=T]]]
velocityMerge[["RNA"]] <- srt[["RNA"]]
DefaultAssay(velocityMerge) <- "RNA"
dir <- "/ssd/lab/HuangMingQian/scRNAseq/22ESC-project/NGSmodule_SCP_analysis/Seurat/"
dir.create(dir,showWarnings = F,recursive = T)
SaveH5Seurat(velocityMerge, filename = paste0(dir,"/velocity.h5Seurat"),overwrite = T)
Convert(paste0(dir,"/velocity.h5Seurat"), dest = "h5ad",overwrite = T)


















