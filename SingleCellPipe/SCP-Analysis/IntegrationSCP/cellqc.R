#!/usr/bin/env Rscript

# Set parameters ---------------------------------------------------------
## parameters: global settings ---------------------------------------------
SCPwork_dir <- "/ssd/lab/HuangMingQian/scRNAseq/H0CellLine-project/NGSmodule_SCP_work/"
SCP_path <- "/home/zhanghao/Program/NGS/UniversalTools/NGSmodule/SingleCellPipe/"
Rscript_threads <- 120
Rscript_memory <- 500

## parameters: cell filtering ----------------------------------------------
db_method <- "scDblFinder" ## "scDblFinder,Scrublet,DoubletDetection,scds_cxds,scds_bcds,scds_hybrid" or NULL, comma-separated.
gene_threshold <- 1000 ## 1000. Cells that exceed this threshold will be kept.
UMI_threshold <- 3000 ## 3000. Cells that exceed this threshold will be kept.
mito_threshold <- 20 ## 20. Cells that exceed this threshold will be discarded.
ribo_threshold <- 50 ## 50. Cells that exceed this threshold will be discarded.
species <- "Homo_sapiens,Mus_musculus" ## "Homo_sapiens,Mus_musculus" or NULL, comma-separated.
species_gene_prefix <- "GRCh38-,mm10-" ## "GRCh38-,mm10-" or NULL, comma-separated.
species_percent <- 5 ## Cells that exceed this threshold will be discarded.
species_nCount <- 5000 ## Cells that exceed this threshold will be discarded.
exogenous_genes <- "GFP" ## Gene symbol or NULL, comma-separated.
features_inspect <- "POU5F1,SOX2,SOX17,SOX9,EOMES,PAX6,PRDM1,SOX1,OTX2,nCount_RNA,nFeature_RNA,percent.mito,percent.ribo" ## Gene symbol or other qc features, comma-separated.

# Library -----------------------------------------------------------------
reticulate::py_run_string("import matplotlib.pyplot as plt")
reticulate::import("scanpy")
suppressWarnings(suppressPackageStartupMessages(invisible(lapply(
  c(
    "Seurat", "SeuratDisk", "future", "dplyr", "ggplot2", "ggridges", "reshape2", "scuttle"
  ),
  require,
  character.only = TRUE
))))
source(paste0(SCP_path, "/SCP-workflow-funtcion.R"))
source(paste0(SCP_path, "/SCP-plot.R"))

main_dir <- paste0(SCPwork_dir, "/../NGSmodule_SCP_analysis/")
dir.create(main_dir, recursive = T, showWarnings = FALSE)
setwd(main_dir)

dir.create("CellQC/RNA", recursive = T, showWarnings = FALSE)
dir.create("CellQC/Velocity", recursive = T, showWarnings = FALSE)
dir.create("CellQC/Plot", recursive = T, showWarnings = FALSE)

db_method <- strsplit(db_method, split = ",") %>% unlist()
species <- strsplit(species, split = ",") %>% unlist()
species_gene_prefix <- strsplit(species_gene_prefix, split = ",") %>% unlist()
exogenous_genes <- strsplit(exogenous_genes, split = ",") %>% unlist()
features_inspect <- c(exogenous_genes, strsplit(features_inspect, split = ",") %>% unlist())

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
rawList <- list()
for (i in 1:length(samples)) {
  cat("++++++", as.character(samples[i]), "(Data Loading)", "++++++", "\n")
  rawList[[as.character(samples[i])]] <- LoadH5Seurat(paste0(main_dir, "/Prepare/RNA/", samples[i], ".h5Seurat"), verbose = FALSE)
}

## qc plot -----------------------------------------------------------------
meta_raw <- Reduce(function(x, y) rbind(x, y), lapply(rawList, function(x) x@meta.data))
meta_raw$CellName <- rownames(meta_raw)
qcvariable <- colnames(meta_raw)[colnames(meta_raw) != "CellName"]
df0 <- reshape2::melt(meta_raw[, qcvariable])
df0$orig.ident <- factor(df0$orig.ident, levels = samples)
nrow <- ceiling(sqrt(length(qcvariable) - 1))
ncol <- ceiling((length(qcvariable) - 1) / nrow)
p0 <- ggplot(df0, aes(x = variable, y = value)) +
  geom_violin(aes(fill = orig.ident), color = "black") +
  scale_fill_manual(name = "Samples", values = palette_zh(df0$orig.ident)) +
  facet_wrap(~variable, scales = "free", nrow = nrow) +
  labs(x = "", y = "") +
  theme_zh(
    aspect.ratio = 0.6,
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
ggsave(p0,
  filename = "CellQC/Plot/raw_QCvariable.vln.png",
  width = ncol * 5, height = nrow * 5
)

plist <- list()
df1 <- table(meta_raw$orig.ident) %>% reshape2::melt()
df1$Var1 <- factor(df1$Var1, levels = samples)
plist[["cellstat"]] <- ggplot(df1, aes(x = Var1, y = value)) +
  geom_bar(aes(fill = Var1), color = "black", stat = "identity") +
  geom_text(aes(label = value), vjust = -0.25) +
  scale_fill_manual(name = "Samples", values = palette_zh(df1$Var1)) +
  scale_y_continuous(limits = c(0, max(df1$value * 1.2))) +
  labs(title = "Cell Number", x = "", y = "") +
  theme_zh(
    aspect.ratio = 0.6,
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
if (length(species) == 2) {
  df2 <- meta_raw[, c("CellName", "orig.ident", paste0("nCount_RNA.", species))]
  df2$orig.ident <- factor(df2$orig.ident, levels = samples)
  p <- ggplot(df2, aes(
    x = .data[[paste0("nCount_RNA.", species)[1]]],
    y = .data[[paste0("nCount_RNA.", species)[2]]],
    fill = orig.ident
  )) +
    geom_point(color = "black", shape = 21, alpha = 0.5) +
    scale_fill_manual(name = "Samples", values = palette_zh(df2$orig.ident)) +
    guides(fill = guide_legend(override.aes = list(alpha = 1))) +
    labs(
      x = paste0("nCount_RNA.", species)[1],
      y = paste0("nCount_RNA.", species)[2]
    ) +
    theme_zh(
      aspect.ratio = 1,
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
  plist[[paste0(unlist(list(samples)), collapse = "-")]] <- p

  for (sp in species) {
    df3 <- meta_raw[, c("CellName", "orig.ident", paste0("percent.genome.", sp))]
    df3$orig.ident <- factor(df3$orig.ident, levels = samples)
    p <- ggplot(df3, aes(
      x = .data[[paste0("percent.genome.", sp)[1]]],
      y = orig.ident,
      fill = orig.ident
    )) +
      geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, gradient_lwd = 1, show.legend = T) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) +
      scale_y_discrete(limits = rev) +
      scale_fill_manual(name = "Group", values = alpha(palette_zh(df3$orig.ident), alpha = 0.8)) +
      labs(y = "") +
      theme_zh() +
      theme(
        aspect.ratio = 0.6, axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major = element_line()
      )
    plist[[sp]] <- p
  }
}
nrow <- ceiling(sqrt(length(plist)))
ncol <- ceiling(length(plist) / nrow)
p <- cowplot::plot_grid(plotlist = plist, nrow = nrow, ncol = ncol, align = "hv", axis = "tblr")
ggsave(p,
  filename = "CellQC/Plot/raw_cellstat.comb.png",
  width = ncol * 5 / PanelRatio(plist[[1]], width = 5, height = 5, units = "in"), height = nrow * 5
)

# Cell QC -----------------------------------
srt_list_QC <- list()
for (i in 1:length(samples)) {
  cat("++++++", samples[i], "(Cell QC)", "++++++", "\n")
  if (!file.exists(paste0("CellQC/RNA/", samples[i], ".qc.h5Seurat"))) {
    cat("Perform QC for the sample:", samples[i], "...\n")
    srt <- rawList[[samples[i]]]
    srt <- CellQC(
      srt = srt, db_method = db_method,
      UMI_threshold = UMI_threshold, gene_threshold = gene_threshold,
      mito_threshold = mito_threshold, ribo_threshold = ribo_threshold,
      species = species, species_gene_prefix = species_gene_prefix,
      species_percent = species_percent, species_nCount = species_nCount
    )
    cat("Save the Seurat object to h5Seurat...\n")
    SaveH5Seurat(
      object = srt,
      filename = paste0("CellQC/RNA/", samples[i], ".qc.h5Seurat"),
      overwrite = TRUE,
      verbose = FALSE
    )
  } else {
    cat("Load", samples[i], "qc data from the h5Seurat...\n")
    srt <- LoadH5Seurat(paste0("CellQC/RNA/", samples[i], ".qc.h5Seurat"), verbose = FALSE)
  }
  srt_list_QC[[samples[i]]] <- srt
}

srt_list_all <- c(srt_list_QC, "All_samples" = Reduce(function(x, y) merge(x, y), srt_list_QC))
for (sample in names(srt_list_all)) {
  cat("++++++", sample, "(QC with Standard_SCP)", "++++++", "\n")
  srt_tmp <- Standard_SCP(srt_list_all[[sample]])
  features <- c(
    features_inspect,
    paste0(features_inspect, rep(paste0(".", species), each = length(features_inspect)))
  )
  features <- c(
    features,
    unlist(sapply(paste0(species_gene_prefix, "-*", rep(features_inspect, each = length(species_gene_prefix)), "$"), function(pat) {
      rownames(srt_tmp)[grep(pattern = pat, x = rownames(srt_tmp))]
    }))
  )
  features <- c(
    "CellFilterng", "db_out", "qc_out", "umi_out", "gene_out", "mito_out", "ribo_out", "species_out",
    exogenous_genes, features
  )
  features <- features[features %in% c(colnames(srt_tmp@meta.data), rownames(srt_tmp))]
  p <- SummaryPlot(
    srt = srt_tmp, clusters = "Standardpcaclusters", groups = "orig.ident",
    reduction = "StandardpcaUMAP2d",
    features = features
  )
  ggsave(
    plot = p, filename = paste0("CellQC/Plot/", sample, ".StandardpcaUMAP2d.qc.png"),
    units = "mm", height = 200, width = 80 * max(c(2 + ceiling(length(features) / 2), length(unique(srt_tmp[["orig.ident", drop = T]])))),
    scale = 1.5, limitsize = FALSE
  )
}

## qc plot --------------------------------------------------------------------
meta_qc <- Reduce(function(x, y) rbind(x, y), lapply(srt_list_QC, function(x) x@meta.data))
meta_qc$CellName <- rownames(meta_qc)
df4 <- meta_qc[, c("CellName", "orig.ident", "db_out", "qc_out", "umi_out", "gene_out", "mito_out", "ribo_out", "species_out", "CellFilterng")]
df4$orig.ident <- factor(df4$orig.ident, levels = samples)
df4$CellFilterng <- factor(df4$CellFilterng, levels = c("Discarded", "Kept"))
p <- ggplot(df4, aes(x = orig.ident, fill = CellFilterng)) +
  geom_bar(color = "black") +
  scale_fill_manual(name = "", values = setNames(c("grey80", "steelblue"), c("Discarded", "Kept"))) +
  geom_text(aes(label = ..count.., group = CellFilterng, color = CellFilterng),
    stat = "count", position = position_stack(vjust = 0.5), color = "black", show.legend = F
  ) +
  scale_color_manual(values = setNames(c("white", "black"), c("FALSE", "TRUE"))) +
  labs(y = "Cell number", x = "") +
  theme_zh(aspect.ratio = 0.6, axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(plot = p, filename = "CellQC/Plot/qc_stat.bar.png", width = 7, height = 7)

filter_upset <- df4 %>%
  dplyr::select(CellName, orig.ident, db_out, qc_out, umi_out, gene_out, mito_out, ribo_out, species_out) %>%
  reshape2::melt(
    measure.vars = c("db_out", "qc_out", "umi_out", "gene_out", "mito_out", "ribo_out", "species_out"),
    variable.name = "filter_by",
    value.name = "filter_out"
  ) %>%
  dplyr::filter(filter_out == "Discarded") %>%
  group_by(CellName) %>%
  dplyr::summarize(
    orig.ident = orig.ident,
    Method_list = list(filter_by),
    Method_comb = paste(filter_by, collapse = ",")
  )

upset_list <- list()
library(ggupset)
for (i in samples) {
  filter_upset0 <- filter_upset %>%
    dplyr::filter(orig.ident == i)
  y_max <- max(table(pull(filter_upset0, "Method_comb")))
  p <- filter_upset0 %>%
    group_by(Method_comb) %>%
    mutate(count = n()) %>%
    ungroup() %>%
    filter(count %in% head(sort(unique(count), decreasing = T), 10)) %>%
    ggplot(aes(x = Method_list)) +
    geom_bar(aes(fill = ..count..), color = "black", width = 0.5) +
    geom_text(aes(label = ..count..), stat = "count", vjust = -0.5, hjust = 0, angle = 45) +
    labs(title = i, x = "", y = "Number") +
    scale_x_upset(sets = c("db_out", "qc_out", "umi_out", "gene_out", "mito_out", "ribo_out", "species_out")) +
    scale_y_continuous(limits = c(0, 1.3 * y_max)) +
    scale_fill_material(name = "Count", palette = "blue-grey") +
    theme_zh() +
    theme(
      aspect.ratio = 0.5,
      legend.position = "none"
    ) +
    theme_combmatrix(
      combmatrix.label.text = element_text(size = 12, color = "black"),
      combmatrix.label.extra_spacing = 6
    )
  upset_list[[i]] <- p
}
nrow <- ceiling(sqrt(length(upset_list)))
ncol <- ceiling(length(upset_list) / nrow)
p <- cowplot::plot_grid(plotlist = upset_list, nrow = nrow, ncol = ncol, align = "hv", axis = "tblr")
ggsave(p,
  filename = "CellQC/Plot/qc_stat.upset.png",
  width = ncol * 5, height = nrow * 5
)

# Cell Filtering -----------------------------------
main_gene <- rownames(srt_list_QC[[1]])[grep(pattern = paste0("^", species_gene_prefix[1]), x = rownames(srt_list_QC[[1]]))]
srt_list_filter <- list()
for (i in 1:length(samples)) {
  cat("++++++", samples[i], "(Cell Filtering)", "++++++", "\n")
  if (!file.exists(paste0("CellQC/RNA/", samples[i], ".filtered.h5Seurat"))) {
    cat("Perform the cell-filtering according to the QC result for the sample:", samples[i], "...\n")
    srt <- srt_list_QC[[samples[i]]]
    srt_filter <- subset(srt, cells = colnames(srt)[srt[["CellFilterng", drop = T]] == "Kept"], features = main_gene)
    srt_filter <- RenameFeatures(srt_filter, newnames = gsub(x = rownames(srt_filter), pattern = paste0("^", species_gene_prefix[1], "-*"), replacement = ""))
    srt_filter[["percent.mito"]] <- PercentageFeatureSet(object = srt_filter, pattern = paste0(c("^MT-", "^Mt-", "^mt-"), collapse = "|"))
    srt_filter[["percent.ribo"]] <- PercentageFeatureSet(object = srt_filter, pattern = paste0(c("^RP[SL]\\d+(\\w+)$", "^Rp[sl]\\d+(\\w+)$", "^rp[sl]\\d+(\\w+)$"), collapse = "|"))
    srt_filter[["log10_nFeature_RNA"]] <- log10(srt_filter[["nFeature_RNA", drop = T]])
    srt_filter[["log10_nCount_RNA"]] <- log10(srt_filter[["nCount_RNA", drop = T]])
    if (isTRUE(exogenous_genes %in% rownames(srt))) {
      for (g in exogenous_genes) {
        srt_filter[[g]] <- log1p(FetchData(object = srt, vars = g, slot = "counts", cells = colnames(srt_filter)))
      }
    }
    cat("Save the filtered RNA data to h5Seurat...\n")
    SaveH5Seurat(
      object = srt_filter,
      filename = paste0("CellQC/RNA/", samples[i], ".filtered.h5Seurat"),
      overwrite = TRUE,
      verbose = FALSE
    )
  } else {
    cat("Load", samples[i], "filtered RNA data from the h5Seurat....\n")
    srt_filter <- LoadH5Seurat(paste0("CellQC/RNA/", samples[i], ".filtered.h5Seurat"), verbose = FALSE)
  }
  srt_list_filter[[samples[i]]] <- srt_filter

  cat("Create the corresponding velocity object...\n")
  velocity <- LoadH5Seurat(paste0("Prepare/Velocity/", samples[i], ".h5Seurat"), verbose = FALSE)
  velocity_filter <- subset(velocity, cells = colnames(srt_filter), features = main_gene)
  velocity_filter <- RenameFeatures(velocity_filter, newnames = gsub(x = rownames(velocity_filter), pattern = paste0("^", species_gene_prefix[1], "-*"), replacement = ""))
  if (!(identical(rownames(srt_filter), rownames(velocity_filter)) &
    identical(colnames(srt_filter), colnames(velocity_filter)))) {
    for (assay in Seurat::Assays(velocity_filter)) {
      velocity_filter[[assay]] <- Seurat::CreateAssayObject(counts = velocity_filter[[assay]][rownames(srt_filter), colnames(srt_filter)])
      velocity_filter[[assay]]@key <- paste0(assay, "_")
    }
  }
  cat("Save the velocity data to h5Seurat...\n")
  SaveH5Seurat(
    object = velocity_filter,
    filename = paste0("CellQC/Velocity/", samples[i], ".filtered.velocity.h5Seurat"),
    overwrite = TRUE,
    verbose = FALSE
  )
}

## qc plot --------------------------------------------------------------------
meta_filter <- Reduce(function(x, y) rbind(x, y), lapply(srt_list_filter, function(x) x@meta.data))
meta_filter$CellName <- rownames(meta_filter)
df0 <- reshape2::melt(meta_filter[, qcvariable])
df0$orig.ident <- factor(df0$orig.ident, levels = samples)
nrow <- ceiling(sqrt(length(qcvariable) - 1))
ncol <- ceiling((length(qcvariable) - 1) / nrow)
p0 <- ggplot(df0, aes(x = variable, y = value)) +
  geom_violin(aes(fill = orig.ident), color = "black") +
  scale_fill_manual(name = "Samples", values = palette_zh(df0$orig.ident)) +
  facet_wrap(~variable, scales = "free", nrow = nrow) +
  labs(x = "", y = "") +
  theme_zh(
    aspect.ratio = 0.6,
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
ggsave(p0,
  filename = "CellQC/Plot/filter_QCvariable.vln.png",
  width = ncol * 5, height = nrow * 5
)

plist <- list()
df1 <- table(meta_filter$orig.ident) %>% reshape2::melt()
df1$Var1 <- factor(df1$Var1, levels = samples)
plist[["cellstat"]] <- ggplot(df1, aes(x = Var1, y = value)) +
  geom_bar(aes(fill = Var1), color = "black", stat = "identity") +
  geom_text(aes(label = value), vjust = -0.25) +
  scale_fill_manual(name = "Samples", values = palette_zh(df1$Var1)) +
  scale_y_continuous(limits = c(0, max(df1$value * 1.2))) +
  labs(title = "Cell Number", x = "", y = "") +
  theme_zh(
    aspect.ratio = 0.6,
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
if (length(species) == 2) {
  df2 <- meta_filter[, c("CellName", "orig.ident", paste0("nCount_RNA.", species))]
  df2$orig.ident <- factor(df2$orig.ident, levels = samples)
  p <- ggplot(df2, aes(
    x = .data[[paste0("nCount_RNA.", species)[1]]],
    y = .data[[paste0("nCount_RNA.", species)[2]]],
    fill = orig.ident
  )) +
    geom_point(color = "black", shape = 21, alpha = 0.5) +
    scale_fill_manual(name = "Samples", values = palette_zh(df2$orig.ident)) +
    guides(fill = guide_legend(override.aes = list(alpha = 1))) +
    labs(
      x = paste0("nCount_RNA.", species)[1],
      y = paste0("nCount_RNA.", species)[2]
    ) +
    theme_zh(
      aspect.ratio = 1,
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
  plist[[paste0(unlist(list(samples)), collapse = "-")]] <- p

  for (sp in species) {
    df3 <- meta_filter[, c("CellName", "orig.ident", paste0("percent.genome.", sp))]
    df3$orig.ident <- factor(df3$orig.ident, levels = samples)
    p <- ggplot(df3, aes(
      x = .data[[paste0("percent.genome.", sp)[1]]],
      y = orig.ident,
      fill = orig.ident
    )) +
      geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, gradient_lwd = 1, show.legend = T) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) +
      scale_y_discrete(limits = rev) +
      scale_fill_manual(name = "Group", values = alpha(palette_zh(df3$orig.ident), alpha = 0.8)) +
      labs(y = "") +
      theme_zh() +
      theme(
        aspect.ratio = 0.6, axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major = element_line()
      )
    plist[[sp]] <- p
  }
}
nrow <- ceiling(sqrt(length(plist)))
ncol <- ceiling(length(plist) / nrow)
p <- cowplot::plot_grid(plotlist = plist, nrow = nrow, ncol = ncol, align = "hv", axis = "tblr")
ggsave(p,
  filename = "CellQC/Plot/filter_cellstat.comb.png",
  width = ncol * 5 / PanelRatio(plist[[1]], width = 5, height = 5, units = "in"), height = nrow * 5
)

# The end -----------------------------------------------------------------
