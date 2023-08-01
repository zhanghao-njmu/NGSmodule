#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
SCPwork_dir <- as.character(args[1])
mode <- as.character(args[2])
genome_file <- as.character(args[3])
gtf_file <- as.character(args[4])

# Set parameters ---------------------------------------------------------
## parameters: global settings ---------------------------------------------
# SCPwork_dir <- "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/NGSmodule_SCP_work/"
# mode <- "rna"
#
# SCPwork_dir <- "/ssd/lab/HuangMingQian/scATACseq/NGSmodule_SCP_work/"
# mode <- "atac"
# genome_file <- "/data/reference/CellRanger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
# gtf_file <- "/data/reference/CellRanger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf"

min_peakwidths <- 100
max_peakwidths <- 10000

# Set working directory -----------------------------------------------------------------------------------
main_dir <- paste0(SCPwork_dir, "/../NGSmodule_SCP_analysis/Preparation")
dir.create(main_dir, recursive = T, showWarnings = FALSE)
setwd(main_dir)
samples <- list.dirs(path = SCPwork_dir, recursive = FALSE, full.names = FALSE)

# Library -----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(invisible(lapply(
  c(
    "SCP", "Seurat", "SeuratWrappers", "future",
    "Signac", "rtracklayer", "Rsamtools", "GenomicRanges"
  ),
  require,
  character.only = TRUE
))))
plan(multisession, workers = min(availableCores(), length(samples), 64))
print(plan())

# Create Seurat object ------------------------------------------------
f <- list()
for (i in seq_along(samples)) {
  nm <- samples[i]
  cat("++++++", nm, "++++++", "\n")
  f[[nm]] <- future(
    {
      if (!file.exists(paste0("./", nm, ".Seurat.rds"))) {
        if (mode == "rna") {
          cat("Load the RNA expression data...\n")
          srt_matrix <- Read10X(data.dir = paste0(SCPwork_dir, "/", nm, "/Alignment-Cellranger/", nm, "/outs/filtered_feature_bc_matrix/"))
          srt <- CreateSeuratObject(counts = srt_matrix, project = paste0(nm, "-RNA"))
          srt <- RenameCells(object = srt, add.cell.id = nm)

          cat("Load the RNA velocity data...\n")
          velocity <- ReadVelocity(file = paste0(SCPwork_dir, "/", nm, "/Alignment-Cellranger/", nm, "/velocyto/", nm, ".loom"), verbose = FALSE)
          velocity <- as.Seurat(x = velocity, project = nm, verbose = FALSE)
          velocity <- RenameCells(
            object = velocity,
            new.names = gsub(x = colnames(velocity), pattern = ".*:", replacement = "", perl = TRUE) %>%
              gsub(x = ., pattern = "x$", replacement = "-1", perl = TRUE) %>%
              paste0(nm, "_", .)
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
        } else if (mode == "atac") {
          cat("Load the ATAC data...\n")
          srt_matrix <- Read10X_h5(filename = paste0(SCPwork_dir, "/", nm, "/Alignment-Cellranger-atac/", nm, "/outs/filtered_peak_bc_matrix.h5"))
          chr <- unique(gsub(":.*", "", rownames(srt_matrix)))
          metadata <- read.csv(file = paste0(SCPwork_dir, "/", nm, "/Alignment-Cellranger-atac/", nm, "/outs/singlecell.csv"), header = TRUE, row.names = 1)
          gtf <- rtracklayer::import(gtf_file)
          # seqlevels(gtf, pruning.mode = "coarse") <- intersect(seqlevels(gtf), chr)
          mcols_df <- mcols(gtf)
          if ("gene_type" %in% colnames(mcols_df)) {
            colnames(mcols_df)[colnames(mcols_df) == "gene_type"] <- "gene_biotype"
          }
          if ("transcript_id" %in% colnames(mcols_df)) {
            colnames(mcols_df)[colnames(mcols_df) == "transcript_id"] <- "tx_id"
          }
          mcols(gtf) <- mcols_df
          fa <- Rsamtools::FaFile(genome_file)
          genome <- seqinfo(fa)
          genome(genome) <- genome_file
          seqlevels(genome) <- seqlevels(gtf)
          seqinfo(gtf) <- genome
          srt_assay <- CreateChromatinAssay(
            counts = srt_matrix,
            sep = c(":", "-"),
            genome = genome,
            annotation = gtf,
            min.cells = 1
          )
          Key(srt_assay) <- "atac_"
          srt <- CreateSeuratObject(
            counts = srt_assay,
            assay = "ATAC",
            meta.data = metadata[colnames(srt_assay), , drop = FALSE],
            project = paste0(nm, "-ATAC")
          )
          srt <- RenameCells(object = srt, add.cell.id = nm)
          fpath <- paste0(SCPwork_dir, "/", nm, "/Alignment-Cellranger-atac/", nm, "/outs/fragments.tsv.gz")
          frags <- CreateFragmentObject(path = fpath, cells = setNames(colnames(srt_matrix), paste0(nm, "_", colnames(srt_matrix))), verbose = FALSE)
          Fragments(srt[["ATAC"]]) <- frags
        }
      } else {
        cat(nm, "Seurat object is ready.\n")
        srt <- readRDS(paste0("./", nm, ".Seurat.rds"))
      }
      srt
    },
    seed = TRUE
  )
}
srtList <- lapply(f, FUN = value)

if (!file.exists("./srtMerge.Seurat.rds")) {
  if (length(srtList) > 1) {
    if (mode == "rna") {
      f <- list()
      for (i in seq_along(srtList)) {
        nm <- names(srtList)[i]
        srt <- srtList[[i]]
        f[[nm]] <- future({
          saveRDS(srt, file = paste0("./", nm, ".Seurat.rds"))
          srt
        })
      }
      srtList <- lapply(f, FUN = value)

      cat("Merge Seurat objects...\n")
      srtMerge <- merge(srtList[[1]], srtList[2:length(srtList)])
      saveRDS(srtMerge, file = paste0("./srtMerge.Seurat.rds"))
    }
    if (mode == "atac") {
      # combined.peaks <- suppressWarnings(reduce(do.call(base::c, unname(lapply(srtList, function(srt) granges(srt[["ATAC"]]))))))
      combined.peaks <- suppressWarnings(UnifyPeaks(object.list = lapply(srtList, function(srt) srt[["ATAC"]]), mode = "reduce"))
      peakwidths <- width(combined.peaks)
      cat("Distribution of peak width :\n")
      print(summary(peakwidths))
      message("Filter out bad peaks with widths < ", min_peakwidths, " : ", sum(peakwidths < min_peakwidths))
      message("Filter out bad peaks with widths > ", max_peakwidths, " : ", sum(peakwidths > max_peakwidths))
      combined.peaks <- combined.peaks[peakwidths >= min_peakwidths & peakwidths <= max_peakwidths]

      f <- list()
      for (i in seq_along(srtList)) {
        nm <- names(srtList)[i]
        srt <- srtList[[i]]
        f[[nm]] <- future(
          {
            counts <- FeatureMatrix(
              features = combined.peaks,
              fragments = Fragments(srt[["ATAC"]]),
              cells = colnames(srt[["ATAC"]]),
              verbose = FALSE
            )
            unified <- CreateChromatinAssay(
              counts = counts,
              genome = seqinfo(srt[["ATAC"]]),
              annotation = Annotation(srt[["ATAC"]]),
              fragments = Fragments(srt[["ATAC"]]),
              min.cells = 1
            )
            Key(unified) <- "atac_"
            srt[["ATAC"]] <- unified
            saveRDS(srt, file = paste0("./", nm, ".Seurat.rds"))
            srt
          },
          seed = TRUE
        )
      }
      srtList <- lapply(f, FUN = value)

      cat("Merge Seurat objects...\n")
      srtMerge <- merge(srtList[[1]], srtList[2:length(srtList)])
      saveRDS(srtMerge, file = paste0("./srtMerge.Seurat.rds"))
    }
  } else {
    saveRDS(srtList[[1]], file = paste0("./", names(srtList)[1], ".Seurat.rds"))
    file.copy(paste0("./", names(srtList)[1], ".Seurat.rds"), paste0("./srtMerge.Seurat.rds"))
  }
} else {
  cat("All tasks have been completed.")
}
