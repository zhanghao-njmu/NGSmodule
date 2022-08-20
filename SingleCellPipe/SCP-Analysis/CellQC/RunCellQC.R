#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
analysis_dir <- as.character(args[1])
samples <- as.character(args[2])
simple_analysis <- as.logical(args[3])
db_method <- as.character(args[4])
outlier_cutoff <- as.character(args[5])
outlier_n <- as.numeric(args[6])
gene_threshold <- as.numeric(args[7])
UMI_threshold <- as.numeric(args[8])
mito_threshold <- as.numeric(args[9])
ribo_threshold <- as.numeric(args[10])
ribo_mito_ratio_min <- as.numeric(args[11])
ribo_mito_ratio_max <- as.numeric(args[12])
species <- as.character(args[13])
species_gene_prefix <- as.character(args[14])
species_percent <- as.numeric(args[15])
exogenous_genes <- as.character(args[16])
features_inspect <- as.character(args[17])


# # Set parameters ---------------------------------------------------------
# ## parameters: global settings ---------------------------------------------
# analysis_dir <- "/ssd/lab/ChenMengQi/scRNAseq/in_vitro_culture/analysis_zh/NGSmodule_SCP_analysis/"
# samples <- "" ## Leave blank or comma-separated names of samples to be analyzed. Leave blank means analyze all samples.
# simple_analysis <- TRUE ## Whether to do a simple analysis.
#
# ## parameters: cell filtering ----------------------------------------------
# db_method <- "scDblFinder" ## Doublet-calling methods used. Can be one of scDblFinder, Scrublet, DoubletDetection, scds_cxds, scds_bcds, scds_hybrid.
# outlier_cutoff <- "log10_nCount:lower:2.5,log10_nCount:higher:5,log10_nFeature:lower:2.5,log10_nFeature:higher:5,featurecount_dist:lower:2.5"
# outlier_n <- 1
# gene_threshold <- 1000. ## 1000. Minimum threshold for the cell gene count.
# UMI_threshold <- 3000 ## 3000. Minimum threshold for the cell UMI count.
# mito_threshold <- 20 ## 15. Maximum threshold for the count proportion of mitochondrial genes.
# ribo_threshold <- 20 ## 50. Maximum threshold for the count proportion of ribosomal genes.
# ribo_mito_ratio_threshold <- 1
# species <- "" ## Leave blank or comma-separated species names, e.g. "Homo_sapiens,Mus_musculus". The first is the species of interest.
# species_gene_prefix <- "" ## Leave blank or comma-separated prefixes, e.g. "GRCh38-,mm10-". The first is the species of interest.
# species_percent <- 95 ## Count proportion thresholds for species of interest.
# exogenous_genes <- "" ## Leave blank or or comma-separated gene symbol.
# features_inspect <- "nCount_RNA,nFeature_RNA,percent.mito,percent.ribo" ## Comma-separated gene symbol or other features.

# Library -----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(invisible(lapply(
  c(
    "SCP", "Seurat", "SeuratDisk", "dplyr", "reshape2",
    "ggplot2", "ggridges", "ggrepel", "ggupset"
  ),
  require,
  character.only = TRUE
))))
set.seed(11)
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
setwd(analysis_dir)
dir.create("CellQC/Plot", recursive = TRUE, showWarnings = FALSE)

db_method <- strsplit(db_method, split = ",") %>% unlist()
outlier_cutoff <- strsplit(outlier_cutoff, split = ",") %>% unlist()
species <- strsplit(species, split = ",") %>% unlist()
species_gene_prefix <- strsplit(species_gene_prefix, split = ",") %>% unlist()
exogenous_genes <- strsplit(exogenous_genes, split = ",") %>% unlist()
features_inspect <- c(exogenous_genes, strsplit(features_inspect, split = ",") %>% unlist())

if (length(species) == 0 || species == "") {
  species <- species_gene_prefix <- NULL
}
if (samples == "") {
  samples <- gsub(".h5Seurat", "", list.files(path = paste0(analysis_dir, "/Prepare"), pattern = ".h5Seurat"))
} else {
  samples <- strsplit(samples, split = ",") %>% unlist()
}

# Data Loading ------------------------------------------------
rawList <- list()
for (i in 1:length(samples)) {
  cat("++++++", as.character(samples[i]), "(Data Loading)", "++++++", "\n")
  rawList[[as.character(samples[i])]] <- LoadH5Seurat(paste0(analysis_dir, "/Prepare/", samples[i], ".h5Seurat"), verbose = FALSE)
}

## basic plot -----------------------------------------------------------------
meta_raw <- Reduce(function(x, y) rbind(x, y), lapply(rawList, function(x) x@meta.data))
meta_raw$CellName <- rownames(meta_raw)
qcvariable <- colnames(meta_raw)[colnames(meta_raw) != "CellName"]
df0 <- reshape2::melt(meta_raw[, qcvariable])
df0$orig.ident <- factor(df0$orig.ident, levels = samples)
ncol <- ceiling(sqrt(length(qcvariable) - 1))
df <- df0 %>%
  group_by(variable, orig.ident) %>%
  summarise(
    lower = quantile(value, 0.25),
    middle = quantile(value, 0.5),
    upper = quantile(value, 0.75),
    ymin = max(min(value), lower - 1.5 * (upper - lower)),
    ymax = min(max(value), upper + 1.5 * (upper - lower)),
    .groups = "keep"
  )
p <- ggplot(df, aes(x = variable, fill = orig.ident)) +
  geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
    stat = "identity", color = "black"
  ) +
  scale_fill_manual(name = "Samples", values = palette_scp(df0$orig.ident)) +
  facet_wrap(~variable, scales = "free", ncol = ncol) +
  labs(x = "", y = "") +
  theme_scp(
    aspect.ratio = 5 / length(samples),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90")
  )
p <- panel_fix(p, height = 2, width = max(2 / (5 / length(samples)), 3), save = "CellQC/Plot/raw_QCvariable.box.pdf")

df1 <- table(meta_raw$orig.ident) %>% reshape2::melt()
df1$Var1 <- factor(df1$Var1, levels = samples)
p <- ggplot(df1, aes(x = Var1, y = value)) +
  geom_bar(aes(fill = Var1), color = "black", stat = "identity") +
  geom_text_repel(aes(label = value),
    nudge_y = max(df1$value) * 0.1, min.segment.length = 1,
    point.size = NA, bg.color = "white", bg.r = 0.1, size = 4, color = "black"
  ) +
  scale_fill_manual(name = "Samples", values = palette_scp(df1$Var1)) +
  scale_y_continuous(limits = c(0, max(df1$value * 1.15))) +
  labs(title = "Cell Number", x = "", y = "") +
  theme_scp(
    aspect.ratio = 5 / length(samples),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.y = element_line(colour = "grey90")
  )
p <- panel_fix(p, height = 2, save = "CellQC/Plot/raw_cellnumber.bar.pdf")

# Cell QC -----------------------------------
srt_list_QC <- list()
for (i in 1:length(samples)) {
  cat("++++++", samples[i], "(Cell QC)", "++++++", "\n")
  if (!file.exists(paste0("CellQC/", samples[i], ".qc.h5Seurat"))) {
    cat("Perform QC for the sample:", samples[i], "...\n")
    srt <- rawList[[samples[i]]]
    srt <- RunCellQC(
      srt = srt, db_method = db_method,
      outlier_cutoff = outlier_cutoff, outlier_n = outlier_n,
      UMI_threshold = UMI_threshold, gene_threshold = gene_threshold,
      mito_threshold = mito_threshold, ribo_threshold = ribo_threshold,
      ribo_mito_ratio_range = c(ribo_mito_ratio_min, ribo_mito_ratio_max),
      species = species, species_gene_prefix = species_gene_prefix, species_percent = species_percent
    )
    cat("Save the Seurat object as h5Seurat...\n")
    SaveH5Seurat(
      object = srt,
      filename = paste0("CellQC/", samples[i], ".qc.h5Seurat"),
      overwrite = TRUE,
      verbose = FALSE
    )
  } else {
    cat("Load", samples[i], "qc data from the h5Seurat...\n")
    srt <- LoadH5Seurat(paste0("CellQC/", samples[i], ".qc.h5Seurat"), verbose = FALSE)
  }
  srt_list_QC[[samples[i]]] <- srt
}
if (!file.exists(paste0("CellQC/Merge.qc.h5Seurat"))) {
  srt_qc_merge <- Reduce(merge, srt_list_QC)
  for (qc in c("db_qc", "outlier_qc", "umi_qc", "gene_qc", "mito_qc", "ribo_qc", "ribo_mito_ratio_qc", "species_qc", "CellQC")) {
    srt_qc_merge[[qc]] <- factor(srt_qc_merge[[qc, drop = TRUE]], levels = c("Pass", "Fail"))
  }
  cat("Save all qc data as one merged h5Seurat...\n")
  SaveH5Seurat(
    object = srt_qc_merge,
    filename = "CellQC/Merge.qc.h5Seurat",
    overwrite = TRUE,
    verbose = FALSE
  )
}

## qc statistic plot --------------------------------------------------------------------
meta_qc <- Reduce(function(x, y) rbind(x, y), lapply(srt_list_QC, function(x) x@meta.data))
meta_qc$CellName <- rownames(meta_qc)
if (length(species) == 2) {
  df2 <- meta_qc[, c("CellName", "orig.ident", paste0("nCount_RNA.", species))]
  df2$orig.ident <- factor(df2$orig.ident, levels = samples)
  s1_q99 <- quantile(df2[, paste0("nCount_RNA.", species)[1]], 0.99)
  s2_q99 <- quantile(df2[, paste0("nCount_RNA.", species)[2]], 0.99)
  df2 <- df2[df2[[paste0("nCount_RNA.", species)[1]]] < s1_q99 & df2[[paste0("nCount_RNA.", species)[2]]] < s2_q99, ]
  p <- ggplot(df2, aes(
    x = .data[[paste0("nCount_RNA.", species)[1]]],
    y = .data[[paste0("nCount_RNA.", species)[2]]],
    fill = orig.ident
  )) +
    geom_point(color = "black", shape = 21) +
    scale_fill_manual(name = "Samples", values = palette_scp(df2$orig.ident)) +
    guides(fill = guide_legend(override.aes = list(alpha = 1, size = 4.5))) +
    labs(
      title = "Species-mixed cells",
      x = paste0("nCount_RNA.", species)[1],
      y = paste0("nCount_RNA.", species)[2]
    ) +
    theme_scp(
      aspect.ratio = 1,
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major = element_line(colour = "grey90")
    )
  p <- panel_fix(p, height = 3, save = "CellQC/Plot/raw_speciesmixed.point.pdf")

  for (sp in species) {
    df3 <- meta_qc[, c("CellName", "orig.ident", paste0("percent.genome.", sp))]
    df3$orig.ident <- factor(df3$orig.ident, levels = samples)
    p <- ggplot(df3, aes(
      x = .data[[paste0("percent.genome.", sp)[1]]],
      y = orig.ident,
      fill = orig.ident
    )) +
      geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, gradient_lwd = 1, show.legend = FALSE) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) +
      scale_y_discrete(limits = rev) +
      scale_fill_manual(name = "Group", values = alpha(palette_scp(df3$orig.ident), alpha = 0.8)) +
      labs(title = "Species RNA proportion distribution", y = "") +
      theme_scp(
        aspect.ratio = length(samples) / 10, axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major = element_line(colour = "grey90")
      )
    p <- panel_fix(p, width = 3, margin = 0.1, save = paste0("CellQC/Plot/raw_proportion.", sp, ".density.pdf"))
  }
}

df4 <- meta_qc[, c("CellName", "orig.ident", "db_qc", "outlier_qc", "umi_qc", "gene_qc", "mito_qc", "ribo_qc", "ribo_mito_ratio_qc", "species_qc", "CellQC")]
df4$orig.ident <- factor(df4$orig.ident, levels = samples)
df4$CellQC <- factor(df4$CellQC, levels = c("Fail", "Pass"))
p <- ggplot(df4, aes(x = orig.ident, fill = CellQC)) +
  geom_bar(color = "black") +
  scale_fill_manual(values = setNames(c("grey80", "steelblue"), c("Fail", "Pass"))) +
  geom_text_repel(aes(label = ..count.., group = CellQC, color = CellQC),
    stat = "count", position = position_stack(vjust = 0.5),
    point.size = NA, bg.color = "white", bg.r = 0.1, size = 4, color = "black", show.legend = FALSE
  ) +
  scale_color_manual(values = setNames(c("white", "black"), c("FALSE", "TRUE"))) +
  labs(title = "Quality Control summary", y = "Cell number", x = "") +
  theme_scp(
    aspect.ratio = 5 / length(samples),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_line(colour = "grey90")
  )
p <- panel_fix(p, height = 3, save = "CellQC/Plot/qc_stat.bar.pdf")

filter_upset <- df4 %>%
  dplyr::select(CellName, orig.ident, db_qc, outlier_qc, umi_qc, gene_qc, mito_qc, ribo_qc, ribo_mito_ratio_qc, species_qc) %>%
  reshape2::melt(
    measure.vars = c("db_qc", "outlier_qc", "umi_qc", "gene_qc", "mito_qc", "ribo_qc", "ribo_mito_ratio_qc", "species_qc"),
    variable.name = "filter_by",
    value.name = "filter_out"
  ) %>%
  dplyr::filter(filter_out == "Fail") %>%
  group_by(CellName) %>%
  summarise(
    orig.ident = orig.ident,
    Method_list = list(filter_by),
    Method_comb = paste(filter_by, collapse = ","),
    .groups = "keep"
  )

upset_list <- list()
for (i in samples) {
  filter_upset0 <- filter_upset %>%
    dplyr::filter(orig.ident == i)
  y_max <- max(table(pull(filter_upset0, "Method_comb")))
  p <- filter_upset0 %>%
    group_by(Method_comb) %>%
    mutate(count = n()) %>%
    ungroup() %>%
    filter(count %in% head(sort(unique(count), decreasing = TRUE), 10)) %>%
    ggplot(aes(x = Method_list)) +
    geom_bar(aes(fill = ..count..), color = "black", width = 0.5, show.legend = FALSE) +
    geom_text_repel(aes(label = ..count..),
      stat = "count", vjust = 0.5, hjust = 0.5, angle = 45,
      nudge_y = max(table(filter_upset0$Method_comb)) * 0.1, min.segment.length = 1,
      point.size = NA, bg.color = "white", bg.r = 0.1, size = 4, color = "black"
    ) +
    labs(title = i, x = "", y = "Number of cells with failed QC") +
    scale_x_upset(sets = c("db_qc", "outlier_qc", "umi_qc", "gene_qc", "mito_qc", "ribo_qc", "ribo_mito_ratio_qc", "species_qc")) +
    scale_y_continuous(limits = c(0, 1.2 * y_max)) +
    scale_fill_gradientn(name = "Count", colors = palette_scp(palette = "material-indigo")) +
    theme_scp(
      aspect.ratio = 0.6,
      axis.title.y = element_text(hjust = 1, angle = 90),
      panel.grid.major = element_line(colour = "grey90")
    ) +
    theme_combmatrix(
      combmatrix.label.text = element_text(size = 12, color = "black"),
      combmatrix.label.extra_spacing = 6
    )
  p <- panel_fix(p, height = 2, margin = 0.1)
  upset_list[[i]] <- p
}
units <- attr(p, "size")$units
width <- attr(p, "size")$width
height <- attr(p, "size")$height
nrow <- ceiling(sqrt(length(upset_list)))
ncol <- ceiling(length(upset_list) / nrow)
p_all <- cowplot::plot_grid(plotlist = upset_list, nrow = nrow, ncol = ncol, align = "hv", axis = "tblr")
ggsave(p_all,
  filename = "CellQC/Plot/qc_stat.upset.pdf", limitsize = FALSE,
  width = ncol * width, height = nrow * height, units = units
)

# Cell Filtering -----------------------------------
srt_list_filter <- list()
for (i in 1:length(samples)) {
  cat("++++++", samples[i], "(Cell Filtering)", "++++++", "\n")
  if (!file.exists(paste0("CellQC/", samples[i], ".filtered.h5Seurat"))) {
    cat("Perform the cell-filtering according to the QC result for the sample:", samples[i], "...\n")
    srt <- srt_list_QC[[samples[i]]]
    cell_pass <- colnames(srt)[srt[["CellQC", drop = TRUE]] == "Pass"]
    if (length(cell_pass) == 0) {
      warning("No cells passed quality control", immediate. = TRUE)
      next
    }
    main_gene <- rownames(srt)[grep(pattern = paste0("^", species_gene_prefix[1]), x = rownames(srt))]
    srt_filter <- srt[main_gene, cell_pass]
    srt_filter <- RenameFeatures(srt_filter, newnames = gsub(x = rownames(srt_filter), pattern = paste0("^", species_gene_prefix[1], "-*"), replacement = ""))
    srt_filter[["percent.mito"]] <- PercentageFeatureSet(object = srt_filter, pattern = paste0(c("^MT-", "^Mt-", "^mt-"), collapse = "|"))
    srt_filter[["percent.ribo"]] <- PercentageFeatureSet(object = srt_filter, pattern = paste0(c("^RP[SL]\\d+\\w{0,1}\\d*$", "^Rp[sl]\\d+\\w{0,1}\\d*$", "^rp[sl]\\d+\\w{0,1}\\d*$"), collapse = "|"))
    srt_filter[["log10_nFeature_RNA"]] <- log10(srt_filter[["nFeature_RNA", drop = TRUE]])
    srt_filter[["log10_nCount_RNA"]] <- log10(srt_filter[["nCount_RNA", drop = TRUE]])
    if (isTRUE(exogenous_genes %in% rownames(srt))) {
      for (g in exogenous_genes) {
        srt_filter[[g]] <- FetchData(object = srt, vars = g, slot = "counts", cells = colnames(srt_filter))
      }
    }
    cat("Save the filtered RNA data as h5Seurat...\n")
    SaveH5Seurat(
      object = srt_filter,
      filename = paste0("CellQC/", samples[i], ".filtered.h5Seurat"),
      overwrite = TRUE,
      verbose = FALSE
    )
  } else {
    cat("Load", samples[i], "filtered RNA data from the h5Seurat....\n")
    srt_filter <- LoadH5Seurat(paste0("CellQC/", samples[i], ".filtered.h5Seurat"), verbose = FALSE)
  }
  srt_list_filter[[samples[i]]] <- srt_filter
}
if (!file.exists(paste0("CellQC/Merge.filtered.h5Seurat"))) {
  srt_filter_merge <- Reduce(merge, srt_list_filter)
  for (qc in c("db_qc", "outlier_qc", "umi_qc", "gene_qc", "mito_qc", "ribo_qc", "ribo_mito_ratio_qc", "species_qc", "CellQC")) {
    srt_filter_merge[[qc]] <- factor(srt_filter_merge[[qc, drop = TRUE]], levels = c("Pass", "Fail"))
  }
  cat("Save all filtered data as one merged h5Seurat...\n")
  SaveH5Seurat(
    object = srt_filter_merge,
    filename = "CellQC/Merge.filtered.h5Seurat",
    overwrite = TRUE,
    verbose = FALSE
  )
}

## basic plot --------------------------------------------------------------------
meta_filter <- Reduce(function(x, y) rbind(x, y), lapply(srt_list_filter, function(x) x@meta.data))
meta_filter$CellName <- rownames(meta_filter)
df0 <- reshape2::melt(meta_filter[, qcvariable])
df0$orig.ident <- factor(df0$orig.ident, levels = samples)
ncol <- ceiling(sqrt(length(qcvariable) - 1))
df <- df0 %>%
  group_by(variable, orig.ident) %>%
  summarise(
    lower = quantile(value, 0.25),
    middle = quantile(value, 0.5),
    upper = quantile(value, 0.75),
    ymin = max(min(value), lower - 1.5 * (upper - lower)),
    ymax = min(max(value), upper + 1.5 * (upper - lower)),
    .groups = "keep"
  )
p <- ggplot(df, aes(x = variable, fill = orig.ident)) +
  geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
    stat = "identity", color = "black"
  ) +
  scale_fill_manual(name = "Samples", values = palette_scp(df0$orig.ident)) +
  facet_wrap(~variable, scales = "free", ncol = ncol) +
  labs(x = "", y = "") +
  theme_scp(
    aspect.ratio = 5 / length(samples),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90")
  )
p <- panel_fix(p, height = 2, width = max(2 / (5 / length(samples)), 3), save = "CellQC/Plot/filtered_QCvariable.box.pdf")

df1 <- table(meta_filter$orig.ident) %>% reshape2::melt()
df1$Var1 <- factor(df1$Var1, levels = samples)
p <- ggplot(df1, aes(x = Var1, y = value)) +
  geom_bar(aes(fill = Var1), color = "black", stat = "identity") +
  geom_text_repel(aes(label = value),
    nudge_y = max(df1$value) * 0.1, min.segment.length = 1,
    point.size = NA, bg.color = "white", bg.r = 0.1, size = 4, color = "black"
  ) +
  scale_fill_manual(name = "Samples", values = palette_scp(df1$Var1)) +
  scale_y_continuous(limits = c(0, max(df1$value * 1.15))) +
  labs(title = "Cell Number", x = "", y = "") +
  theme_scp(
    aspect.ratio = 5 / length(samples),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.y = element_line(colour = "grey90")
  )
p <- panel_fix(p, height = 2, save = "CellQC/Plot/filtered_cellnumber.bar.pdf")

# Simple Analysis --------------------------------------------------------
if (isTRUE(simple_analysis)) {
  if (file.exists("CellQC/Merge.qc.h5Seurat")) {
    srt_qc_merge <- LoadH5Seurat("CellQC/Merge.qc.h5Seurat", verbose = FALSE)
    srt_list_all <- c(srt_list_QC, "Merge" = srt_qc_merge)
  } else {
    srt_qc_merge <- Reduce(merge, srt_list_QC)
    for (qc in c("db_qc", "outlier_qc", "umi_qc", "gene_qc", "mito_qc", "ribo_qc", "ribo_mito_ratio_qc", "species_qc", "CellQC")) {
      srt_qc_merge[[qc]] <- factor(srt_qc_merge[[qc, drop = TRUE]], levels = c("Pass", "Fail"))
    }
    srt_list_all <- c(srt_list_QC, "Merge" = srt_qc_merge)
  }
  for (sample in names(srt_list_all)) {
    cat("++++++", sample, "(Simple Analysis)", "++++++", "\n")
    srt <- srt_list_all[[sample]]
    if (length(Reductions(srt)) == 0) {
      srt <- Standard_SCP(srt)
      cat("Save the Seurat object as h5Seurat...\n")
      SaveH5Seurat(
        object = srt,
        filename = paste0("CellQC/", sample, ".qc.tmp.h5Seurat"),
        overwrite = TRUE,
        verbose = FALSE
      )
      file.rename(
        from = paste0("CellQC/", sample, ".qc.tmp.h5Seurat"),
        to = paste0("CellQC/", sample, ".qc.h5Seurat")
      )
    } else {
      cat("Sample:", sample, "has been analyzed.\n")
    }
    features <- c(
      features_inspect,
      paste0(features_inspect, rep(paste0(".", species), each = length(features_inspect)))
    )
    if (length(species_gene_prefix) > 0) {
      features <- c(
        features,
        unlist(sapply(paste0(species_gene_prefix, "-*", rep(features_inspect, each = length(species_gene_prefix)), "$"), function(pat) {
          rownames(srt)[grep(pattern = pat, x = rownames(srt))]
        }))
      )
    }
    features <- c(
      "CellQC", "db_qc", "outlier_qc", "umi_qc", "gene_qc", "mito_qc", "ribo_qc", "ribo_mito_ratio_qc", "species_qc",
      exogenous_genes, features
    )
    features <- features[features %in% c(colnames(srt@meta.data), rownames(srt))]
    p <- SummaryPlot(
      srt = srt, group.by = c("orig.ident", "Standardclusters"), group.split.by = "orig.ident",
      features = features, reduction = "StandardUMAP2D"
    )
    ggsave(
      plot = p, filename = paste0("CellQC/Plot/analysis_", sample, ".umap.png"), limitsize = FALSE,
      units = attr(p, "size")$units, width = attr(p, "size")$width, height = attr(p, "size")$height
    )
  }
}

# The end -----------------------------------------------------------------
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}
