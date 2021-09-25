work_dir <- "~"
setwd(work_dir)

srt_list_Standard <- readRDS("/data/lab/HuangMingQian/scRNA-seq/ESC-PGC-GSCLC-new/analysis/bk/srt_list_fastMNN.rds")
if (!species %in% c("Human", "Mouse")) {
  stop("Auto-annotation only support the species: Human,Mouse.\nStop the process.")
}

Sys.setenv(RETICULATE_PYTHON = "~/Program/SystemTools/miniconda3/bin/python")

# Library -----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(invisible(lapply(
  c(
    "Seurat", "SingleCellExperiment", "scran", "reticulate", "stringr", "dplyr",
    "tidyr", "reshape2", "BiocParallel", "biomaRt", "rvest", "xml2",
    "scCATCH", "cellassign", "SingleR", "celaref", "rSuperCT", "scibet"
  ),
  require,
  character.only = TRUE
))))


# tool types: ---------------------------------------------------------------
# build_in: SingleR,rSuperCT
# need_marker: scCATCH,SCSA,cellassign,
# need_expression_matrix: SingleR,celaref


# parameters: -------------------------------------------------------------
srt <- srt_list_Standard[[1]]
# srt <- scCATCH::mouse_kidney_203_Seurat
species <- "Human"
ref_marker <- "CellMatch" # CellMatch,CellMarker,PanglaoDB
ref_rnaseq <- ""
tissues <- c(
  "Testis", "Germ", "Gonad", "Fetal gonad", "Seminal plasma", "Pluripotent stem cell", # CellMarker CellMatch PanglaoDB
  "Embryo", "Embryonic stem cell", "Embryoid body", "Primitive streak"
  # "Reproductive","Embryo","Other"
  # "All"
)
cancers <- NULL
threads <- 125
Ensembl_version <- 101

# scCATCH
scCATCH_cell_min_pct <- 0.25
scCATCH_logfc <- 0.25
scCATCH_pvalue <- 0.05

# cellassign

# SCSA
SCSA_dir <- "/home/zhanghao/Program/NGS/SingleCell/SCSA/"
SCSA_scrpit <- paste0(SCSA_dir, "/SCSA.py")
SCSA_db <- paste0(SCSA_dir, "/whole.db")

# SingleR
SingleR_method <- "cluster"
SingleR_ref <- "build-in"

# celaref

# SuperCT

# SciBet


# Load marker database  ---------------------------------------------------
ref <- get(ref_marker)
if (tissues[1] == "All") {
  tissues <- unique(ref[ref$speciesType == species & ref$cancerType == "Normal", "tissueType"])
}
if (!is.null(cancers)) {
  cell_type <- "Cancer cell"
  if (cancers[1] == "All") {
    cancers <- unique(ref[["cancerType"]])
  }
} else {
  cell_type <- "Normal cell"
}

# Annotation: scCATCH(marker-based) -----------------------------------------------------
if (!all(tissues %in% ref[["tissueType"]])) {
  stop(paste0(
    "Incorrect tissue type:", paste0(tissues[!tissues %in% ref[["tissueType"]]], collapse = ","),
    case_when(
      ref_marker == "CellMatch" ~ paste0("\nPlease refer to the tissue/organ types of CellMatch (", "https://github.com/ZJUFanLab/scCATCH/wiki", ")."),
      ref_marker == "CellMarker" ~ paste0("\nPlease refer to the tissue/organ types of CellMarker (", "http://biocc.hrbmu.edu.cn/CellMarker", ")."),
      ref_marker == "PanglaoDB" ~ paste0("\nPlease refer to the tissue/organ types of PanglaoDB (", "https://panglaodb.se/markers.html?cell_type", ")."),
      TRUE ~ paste0("Please refer to the tissue/organ types of ", ref_marker)
    )
  ))
}

if (ncol(srt) >= 10000 | nlevels(Idents(srt)) >= 15) {
  scCATCH_match_ref <- TRUE
} else {
  scCATCH_match_ref <- FALSE
}

clu_markers_srt <- findmarkergenes(
  object = srt,
  species = species,
  cluster = "All",
  ref = ref,
  match_ref = scCATCH_match_ref,
  cancer = cancers,
  tissue = tissues,
  cell_min_pct = scCATCH_cell_min_pct,
  logfc = scCATCH_logfc,
  pvalue = scCATCH_pvalue,
  BPPARAM = MulticoreParam()
)
clu_markers <- clu_markers_srt$clu_markers

clu_markers <- Markers_Wilcox
clu_markers[, "pct"] <- clu_markers[, "pct.1"]
clu_markers <- clu_markers[clu_markers$gene %in% ref$geneSymbol, ]
clu_ann <- scCATCH(
  object = clu_markers,
  species = species,
  ref = ref,
  cancer = cancers,
  tissue = tissues
)

srt[["scCATCH"]] <- plyr::mapvalues(
  x = as.character(Idents(srt)),
  from = clu_ann[["cluster"]],
  to = clu_ann[["cell_type"]]
) %>% as.factor()
DimPlot(srt,
  reduction = "StandardUMAP2d",
  group.by = c( "scCATCH"), label = TRUE, pt.size = 0.5
)

# Annotation: SCSA(marker-based) -----------------------------------------------------
if (nrow(srt@tools)) {

}

DefaultAssay(srt) <- "RNA"
Markers_MAST <- FindAllMarkers(
  object = srt, only.pos = T, min.pct = 0.25, logfc.threshold = log(2),
  test.use = "MAST"
)
Markers_Wilcox <- FindAllMarkers(
  object = srt, only.pos = T, min.pct = 0.25, logfc.threshold = log(2),
  test.use = "wilcox"
)
Markers_ROC <- FindAllMarkers(
  object = srt, only.pos = T, min.pct = 0.25, logfc.threshold = log(2),
  test.use = "roc"
)
srt@tools$FindAllMarkers <- setNames(
  object = list(Markers_MAST, Markers_ROC),
  nm = c("Markers_MAST", "Markers_ROC")
)

temp <- tempfile()
de <- srt@tools$DEtest_Standardclusters$AllMarkers_Wilcoxon
de$avg_logFC <- de$avg_log2FC
de$cluster <- de$group1
write.csv(de, file = temp, row.names = FALSE)
command <- paste("python3", SCSA_scrpit, "-d", SCSA_db, "-s seurat -i", temp, "-k", paste0('"', paste0(tissues, collapse = ","), '"'), "-g", species, "-E")
system(command)
unlink(temp)


# Annotation: SingleR(reference-based) -----------------------------------------------------
library(SingleR)
if (SingleR_ref == "build-in") {
  state <- 1
  while (state == 1) {
    state <- tryCatch(
      expr = {
        if (species == "Human") {
          hpca <- HumanPrimaryCellAtlasData()
          blueprint <- BlueprintEncodeData()
          dice <- DatabaseImmuneCellExpressionData()
          Monaco <- MonacoImmuneData()
          comom <- Reduce(intersect, c(
            lapply(list(hpca, blueprint, dice, Monaco), rownames)
          ))
          ref <- list(hpca[comom, ], blueprint[comom, ], dice[comom, ], Monaco[comom, ])
          labels <- list(hpca$label.main, blueprint$label.main, dice$label.main, Monaco$label.main)
        }
        if (species == "Mouse") {
          mrd <- MouseRNAseqData()
          igd <- ImmGenData()
          comom <- Reduce(intersect, c(
            lapply(list(mrd, igd), rownames)
          ))
          ref <- list(mrd[comom, ], igd[comom, ])
          labels <- list(mrd$label.main, igd$label.main)
        }
        state <- 0
      }
    )
  }

  # trained <- trainSingleR(ref = ref, labels = labels)
  # prediction <- classifySingleR(test = hESCs[common,], trained)

  prediction <- SingleR(
    test = GetAssayData(srt, slot = "data", assay = "RNA"),
    ref = ref,
    labels = labels,
    method = SingleR_method,
    clusters = Idents(srt),
    BPPARAM = MulticoreParam()
  )
  # plotScoreHeatmap(prediction)

  if (SingleR_method == "single") {
    srt[[paste0("SingleR_", SingleR_method)]] <- plyr::mapvalues(
      x = colnames(srt),
      from = rownames(prediction),
      to = prediction[["labels"]]
    ) %>% as.factor()
  }

  if (SingleR_method == "cluster") {
    srt[[paste0("SingleR_", SingleR_method)]] <- plyr::mapvalues(
      x = as.character(Idents(srt)),
      from = rownames(prediction),
      to = prediction[["labels"]]
    ) %>% as.factor()
  }

  DimPlot(srt,
    reduction = "umap",
    group.by = c("seurat_clusters", paste0("SingleR_", SingleR_method)), label = TRUE, pt.size = 0.5
  )
}


# Annotation: cellassign(marker-based) --------------------------------------------------
library(cellassign)
library(scran)
markerlist <- ref %>%
  filter(species == species, cellType == cell_type, tissueType %in% tissues) %>%
  group_by(shortname) %>%
  summarise(list = list(geneSymbol))
names(markerlist$list) <- markerlist$shortname
marker_mat <- marker_list_to_mat(markerlist$list)

if (nrow(marker_mat) > 2000) {
  stop(paste0("Number of markers(", nrow(marker_mat), ") is too large. Please reduce to less than 2000 markers for cellassign."))
} else {
  sce <- as.SingleCellExperiment(srt, assay = "RNA")
  sce <- computeSumFactors(sce, BPPARAM = MulticoreParam())
  sf <- sizeFactors(sce)
  keep <- intersect(rownames(marker_mat), rownames(sce))
  
  fit <- cellassign(
    exprs_obj = sce[keep, ],
    marker_gene_info = marker_mat[keep, ],
    s = sf,
    learning_rate = 1e-2,
    shrinkage = TRUE,
    threads = 2,
    verbose = TRUE
  )
  unique(fit$cell_type)
  srt[["cellassign"]] <- fit$cell_type
  DimPlot(srt,
          reduction = "umap",
          group.by = c("seurat_clusters", "cellassign"), label = TRUE, pt.size = 0.5
  )
}



# Annotation: Garnett(marker-based)  ----------------------------------------------------------------

# Annotation: celaref(reference-based) -----------------------------------------------------
# Paths to data files.
counts_filepath.query <- system.file("extdata", "sim_query_counts.tab", package = "celaref")
cell_info_filepath.query <- system.file("extdata", "sim_query_cell_info.tab", package = "celaref")
counts_filepath.ref <- system.file("extdata", "sim_ref_counts.tab", package = "celaref")
cell_info_filepath.ref <- system.file("extdata", "sim_ref_cell_info.tab", package = "celaref")

# Load data
toy_ref_se <- load_se_from_files(counts_file = counts_filepath.ref, cell_info_file = cell_info_filepath.ref)
toy_query_se <- load_se_from_files(counts_file = counts_filepath.query, cell_info_file = cell_info_filepath.query)

# Filter data
toy_ref_se <- trim_small_groups_and_low_expression_genes(toy_ref_se)
toy_query_se <- trim_small_groups_and_low_expression_genes(toy_query_se)

# Setup within-experiment differential expression
de_table.toy_ref <- contrast_each_group_to_the_rest(toy_ref_se, dataset_name = "ref")
de_table.toy_query <- contrast_each_group_to_the_rest(toy_query_se, dataset_name = "query")

# Plot and get group labels
make_ranking_violin_plot(de_table.test = de_table.toy_query, de_table.ref = de_table.toy_ref)
make_ref_similarity_names(de_table.toy_query, de_table.toy_ref)

# Annotation: CHETAH(reference-based)  ----------------------------------------------------------------

# Annotation: SciBet(reference-based) ------------------------------------------------------






