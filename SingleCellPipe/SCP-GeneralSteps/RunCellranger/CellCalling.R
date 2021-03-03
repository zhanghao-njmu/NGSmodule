#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
cellranger_dir <- as.character(args[1])
sample <- as.character(args[2])
threads <- as.character(args[3])
EmptyThreshold <- as.character(args[4])
CellLabel <- as.character(args[5])

### parameters: global settings ---------------------------------------------
# cellranger_dir <- "/ssd/lab/ZhangHao/test/SCP/NGSmodule_SCP_work/Testis-d50/Alignment-Cellranger/"
# sample <- "Testis-d50"
# threads <- 80
# EmptyThreshold <- "AUTO"
# CellLabel <- "GFP"
#
# cellranger_dir <- "/ssd/lab/HuangMingQian/scRNAseq/iPSC-ESC-iMeLC-PGCd6-CellLine-Testis/NGSmodule_SCP_work/ESC/Alignment-Cellranger/"
# sample <- "ESC"
# threads <- 80
# EmptyThreshold <- "AUTO"
# CellLabel <- "GFP"

# parameters: emptyDrops --------------------------------------------------
iters <- 1e+5
fdr <- 1e-3

# parameters: dropEst -----------------------------------------------------
dropest_quality <- 0.99


########################### Start the workflow ############################
library(SingleCellExperiment)
library(DropletUtils)
library(dropestr)
library(ggplot2)
library(ggrepel)
library(ggupset)
library(ggsci)
library(png)
library(grid)
library(gridExtra)
library(cowplot)
library(scales)
library(dplyr)
library(tidyr)
library(reshape2)

set.seed(11)
bpparam <- BiocParallel::MulticoreParam(workers = threads)
plotlist <- list()
color <- pal_material("blue-grey")(10)

# Method: Cellranger ------------------------------------------------------
path_raw <- paste0(cellranger_dir, "/", sample, "/outs/raw_feature_bc_matrix")
raw <- read10xCounts(path_raw)
colnames(raw) <- colData(raw)[["Barcode"]]
rownames(raw) <- rowData(raw)[["Symbol"]] %>% make.unique(sep = "@")
raw <- raw[, colSums(counts(raw)) > 0]

path_filtered <- paste0(cellranger_dir, "/", sample, "/outs/filtered_feature_bc_matrix")
filtered <- read10xCounts(path_filtered)
colnames(filtered) <- colData(filtered)[["Barcode"]]
rownames(filtered) <- rowData(filtered)[["Symbol"]] %>% make.unique(sep = "@")
filtered <- filtered[, colSums(counts(filtered)) > 0]

colData(raw)$Cellranger_v2 <- defaultDrops(counts(raw))
colData(raw)$Cellranger_v3 <- colData(raw)$Barcode %in% colData(filtered)$Barcode

if (EmptyThreshold != "AUTO") {
  EmptyThreshold <- as.numeric(EmptyThreshold)
} else {
  # EmptyThreshold <- max(sort(colSums(counts(raw)), decreasing = TRUE)[sum(colData(raw)$Cellranger_v2) * 10], 100)
  EmptyThreshold <- max(sort(colSums(counts(raw)), decreasing = TRUE)[sum(defaultDrops(counts(raw), lower.prop = 0.001))], 100)
}

if (CellLabel != "NULL" & !CellLabel %in% rownames(raw)) {
  warning(paste0("CellLabel:", CellLabel, " is not in the features for the sample ", sample, "\nCellLabel will set to 'NULL'."))
  CellLabel <- "NULL"
}
if (CellLabel == "NULL") {
  CellLabel <- "NULL"
  colData(raw)$CellLabel <- FALSE
} else {
  colData(raw)$CellLabel <- counts(raw)[CellLabel, ] > 0
  colData(raw)$CellLabel <- factor(colData(raw)$CellLabel, levels = c(TRUE, FALSE))
}

# Method: emptyDrops ------------------------------------------------------
bc_ranks <- barcodeRanks(counts(raw), lower = EmptyThreshold)
colData(raw)$nCount <- bc_ranks$total
colData(raw)$nCount_rank <- bc_ranks$rank
colData(raw)$nFeature <- colSums(counts(raw) > 0)
colData(raw)$nFeature_rank <- rank(-(colSums(counts(raw) > 0)))
colData(raw)$Fitted <- bc_ranks$fitted
colData(raw)$knee <- colData(raw)$nCount > metadata(bc_ranks)$knee
colData(raw)$inflection <- colData(raw)$nCount > metadata(bc_ranks)$inflection

vars <- c("nCount", "nCount_rank", "nFeature", "nFeature_rank")
vars_expand <- expand.grid(vars, vars,stringsAsFactors = FALSE)


df1 <- colData(raw) %>%
  as.data.frame() %>%
  group_by(Barcode) %>%
  summarise(
    Barcode = Barcode,
    nCount = nCount,
    nCount_rank = nCount_rank,
    nFeature = nFeature,
    nFeature_rank = nFeature_rank,
    nCount2 = nCount,
    nCount_rank2 = nCount_rank,
    nFeature2 = nFeature,
    nFeature_rank2 = nFeature_rank)
df2 <- reshape2::melt(df1,
  id.vars = c("Barcode", "nCount2", "nCount_rank2", "nFeature2", "nFeature_rank2"),
  variable.name = "var2", value.name = "value2"
)
colnames(df2) <- c("Barcode", "nCount", "nCount_rank", "nFeature", "nFeature_rank", "var2", "value2")
df3 <- reshape2::melt(df2,
  id.vars = c("Barcode", "var2", "value2"),
  variable.name = "var1", value.name = "value1"
)

p <- ggplot(df3, aes(x = value1, y = value2)) +
  geom_point(alpha = 0.5, shape = 16, color = "steelblue") +
  scale_y_continuous(
    trans = log10_trans(),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_x_continuous(
    trans = log10_trans(),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  annotation_logticks() +
  facet_wrap(var1 ~ var2,scales = "fixed") +
  labs(title = "nCount vs nFeature") +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.title = element_blank(),
    panel.grid.major = element_line()
  )
ggsave(p, filename = paste0(sample, ".Basic.png"), width = 11, height = 8)


if (CellLabel != "NULL") {
  p <- colData(raw) %>%
    as.data.frame() %>%
    ggplot(
      aes(x = nCount_rank, y = nCount)
    ) +
    geom_segment(aes(xend = nCount_rank * 2, yend = nCount * 2), color = "grey90", alpha = 0.01) +
    geom_point(
      aes(color = log10(nCount)),
      alpha = 0.5, shape = 16
    ) +
    scale_color_material(palette = "blue-grey", guide = FALSE) +
    geom_point(
      aes(x = nCount_rank * 2, y = nCount * 2, fill = CellLabel),
      color = "transparent",
      shape = 21
    ) +
    scale_fill_manual(values = setNames(c(alpha("red", 0.05), "transparent"), c(TRUE, FALSE))) +
    guides(fill = guide_legend(override.aes = list(fill = c("red", "transparent")))) +
    geom_line(aes(x = nCount_rank, y = Fitted), color = "red") +
    geom_hline(
      yintercept = metadata(bc_ranks)$knee,
      colour = "dodgerblue", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = metadata(bc_ranks)$inflection,
      colour = "forestgreen", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = EmptyThreshold,
      colour = "darkorchid", linetype = "dashed"
    ) +
    annotate("text",
      x = 1, y = metadata(bc_ranks)$knee,
      label = paste0(
        "Knee:",
        " UMI=", round(metadata(bc_ranks)$knee),
        " Rank=", sum(colData(raw)$nCount > metadata(bc_ranks)$knee)
      ),
      colour = "dodgerblue", size = 3, hjust = 0, vjust = -0.5
    ) +
    annotate("text",
      x = 1, y = metadata(bc_ranks)$inflection,
      label = paste0(
        "Inflection:",
        " UMI=", round(metadata(bc_ranks)$inflection),
        " Rank=", sum(colData(raw)$nCount > metadata(bc_ranks)$inflection)
      ),
      colour = "forestgreen", size = 3, hjust = 0, vjust = -0.5
    ) +
    annotate("text",
      x = 1, y = EmptyThreshold,
      label = paste0(
        "Empty threshold:",
        " UMI=", EmptyThreshold,
        " Rank=", sum(colData(raw)$nCount > EmptyThreshold)
      ),
      colour = "darkorchid", size = 3, hjust = 0, vjust = 1.5
    ) +
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks() +
    labs(title = "Barcode Rank Plot", x = "Barcode rank", y = "UMI counts") +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
} else {
  p <- colData(raw) %>% # [sample(1:nrow(colData(raw)),100000),]
    as.data.frame() %>%
    ggplot(
      aes(x = nCount_rank, y = nCount)
    ) +
    geom_point(
      aes(color = log10(nCount)),
      alpha = 0.5, shape = 16
    ) +
    scale_color_material(palette = "blue-grey", guide = FALSE) +
    geom_line(aes(x = nCount_rank, y = Fitted), color = "red") +
    geom_hline(
      yintercept = metadata(bc_ranks)$knee,
      colour = "dodgerblue", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = metadata(bc_ranks)$inflection,
      colour = "forestgreen", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = EmptyThreshold,
      colour = "darkorchid", linetype = "dashed"
    ) +
    annotate("text",
      x = 1, y = metadata(bc_ranks)$knee,
      label = paste0(
        "Knee:",
        " UMI=", round(metadata(bc_ranks)$knee),
        " Rank=", sum(colData(raw)$nCount > metadata(bc_ranks)$knee)
      ),
      colour = "dodgerblue", size = 3, hjust = 0, vjust = -0.5
    ) +
    annotate("text",
      x = 1, y = metadata(bc_ranks)$inflection,
      label = paste0(
        "Inflection:",
        " UMI=", round(metadata(bc_ranks)$inflection),
        " Rank=", sum(colData(raw)$nCount > metadata(bc_ranks)$inflection)
      ),
      colour = "forestgreen", size = 3, hjust = 0, vjust = -0.5
    ) +
    annotate("text",
      x = 1, y = EmptyThreshold,
      label = paste0(
        "Empty threshold:",
        " UMI=", EmptyThreshold,
        " Rank=", sum(colData(raw)$nCount > EmptyThreshold)
      ),
      colour = "darkorchid", size = 3, hjust = 0, vjust = 1.5
    ) +
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks() +
    labs(title = "Barcode Rank Plot", x = "Barcode rank", y = "UMI counts") +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
}
plotlist[["emptyDrops-nCount_rankFitted"]] <- p


mito_gene <- grep(x = rownames(raw), pattern = "(^MT-)|(^Mt-)|(^mt-)", perl = T)
ribo_gene <- grep(x = rownames(raw), pattern = "(^RP[SL]\\d+(\\w|)$)|(^Rp[sl]\\d+(\\w|)$)|(^rp[sl]\\d+(\\w|)$)", perl = T)
residual <- 1000
while (residual != 0) {
  cat("emptyDrops niters:", iters, "\n")
  emp_drops <- emptyDrops(counts(raw)[-c(mito_gene, ribo_gene), ],
    lower = EmptyThreshold, niters = iters,
    test.ambient = TRUE, retain = Inf, BPPARAM = bpparam
  )
  is_cell <- emp_drops$FDR < fdr
  is_cell[is.na(is_cell)] <- FALSE
  print(table(Limited = emp_drops$Limited, Significant = is_cell))
  residual <- table(Limited = emp_drops$Limited, Significant = is_cell)[2, 1]
  iters <- iters * 2
}

colData(raw)$EmpDropsLogProb <- emp_drops$LogProb
colData(raw)$EmpDropsPValue <- emp_drops$PValue
colData(raw)$EmpDropsLimited <- emp_drops$Limited
colData(raw)$EmpDropsFDR <- emp_drops$FDR
colData(raw)$EmptyDrops <- is_cell & colData(raw)$nCount > EmptyThreshold

p <- colData(raw) %>%
  as.data.frame() %>%
  dplyr::filter(
    nCount <= EmptyThreshold,
    nCount > 0
  ) %>%
  ggplot(aes(x = EmpDropsPValue)) +
  geom_histogram(aes(fill = ..count..), bins = 30) +
  scale_fill_material(name = "Count", palette = "blue-grey") +
  labs(
    x = "P-value", y = "UMI counts",
    title = "P-value distribution for empty droplets",
    subtitle = "the distribution of p-values should be approximately uniform"
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    legend.position = "none"
  )
plotlist[["emptyDrops-emptyPValue"]] <- p

if (CellLabel != "NULL") {
  p <- colData(raw) %>%
    as.data.frame() %>%
    arrange(desc(EmptyDrops)) %>%
    ggplot(aes(x = nCount, y = -EmpDropsLogProb)) +
    geom_point(aes(color = EmptyDrops),
      alpha = 0.5, shape = 16
    ) +
    scale_color_manual(
      name = "is cell",
      values = setNames(color[c(9, 2)], c(TRUE, FALSE))
    ) +
    geom_point(
      aes(x = nCount, y = -EmpDropsLogProb * 5, fill = CellLabel),
      color = "transparent",
      shape = 21
    ) +
    scale_fill_manual(values = setNames(c(alpha("red", 0.05), "transparent"), c(TRUE, FALSE))) +
    guides(fill = guide_legend(override.aes = list(fill = c("red", "transparent")))) +
    geom_vline(
      xintercept = metadata(bc_ranks)$knee,
      colour = "dodgerblue", linetype = "dashed"
    ) +
    geom_vline(
      xintercept = metadata(bc_ranks)$inflection,
      colour = "forestgreen", linetype = "dashed"
    ) +
    geom_vline(
      xintercept = EmptyThreshold,
      colour = "darkorchid", linetype = "dashed"
    ) +
    annotate("text",
      y = 1, x = metadata(bc_ranks)$knee,
      label = paste0(
        "Knee:",
        " UMI=", round(metadata(bc_ranks)$knee),
        " Rank=", sum(colData(raw)$nCount > metadata(bc_ranks)$knee)
      ),
      colour = "dodgerblue", size = 3, vjust = 1.5, hjust = 0, angle = 90
    ) +
    annotate("text",
      y = 1, x = metadata(bc_ranks)$inflection,
      label = paste0(
        "Inflection:",
        " UMI=", round(metadata(bc_ranks)$inflection),
        " Rank=", sum(colData(raw)$nCount > metadata(bc_ranks)$inflection)
      ),
      colour = "forestgreen", size = 3, vjust = 1.5, hjust = 0, angle = 90
    ) +
    annotate("text",
      y = 1, x = EmptyThreshold,
      label = paste0(
        "Empty threshold:",
        " UMI=", EmptyThreshold,
        " Rank=", sum(colData(raw)$nCount > EmptyThreshold)
      ),
      colour = "darkorchid", size = 3, vjust = -0.5, hjust = 0, angle = 90
    ) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks() +
    labs(
      x = "UMI counts", y = "-log10(Probability)",
      title = "The negative log-probability of being empty-dropets",
      subtitle = "droplets detected as cells should show up with \nlarge negative log-probabilities"
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
} else {
  p <- colData(raw) %>%
    as.data.frame() %>%
    arrange(desc(EmptyDrops)) %>%
    ggplot(aes(x = nCount, y = -EmpDropsLogProb)) +
    geom_point(aes(color = EmptyDrops),
      alpha = 0.5, shape = 16
    ) +
    scale_color_manual(
      name = "is cell",
      values = setNames(color[c(9, 2)], c(TRUE, FALSE))
    ) +
    geom_vline(
      xintercept = metadata(bc_ranks)$knee,
      colour = "dodgerblue", linetype = "dashed"
    ) +
    geom_vline(
      xintercept = metadata(bc_ranks)$inflection,
      colour = "forestgreen", linetype = "dashed"
    ) +
    geom_vline(
      xintercept = EmptyThreshold,
      colour = "darkorchid", linetype = "dashed"
    ) +
    annotate("text",
      y = 1, x = metadata(bc_ranks)$knee,
      label = paste0(
        "Knee:",
        " UMI=", round(metadata(bc_ranks)$knee),
        " Rank=", sum(colData(raw)$nCount > metadata(bc_ranks)$knee)
      ),
      colour = "dodgerblue", size = 3, vjust = 1.5, hjust = 0, angle = 90
    ) +
    annotate("text",
      y = 1, x = metadata(bc_ranks)$inflection,
      label = paste0(
        "Inflection:",
        " UMI=", round(metadata(bc_ranks)$inflection),
        " Rank=", sum(colData(raw)$nCount > metadata(bc_ranks)$inflection)
      ),
      colour = "forestgreen", size = 3, vjust = 1.5, hjust = 0, angle = 90
    ) +
    annotate("text",
      y = 1, x = EmptyThreshold,
      label = paste0(
        "Empty threshold:",
        " UMI=", EmptyThreshold,
        " Rank=", sum(colData(raw)$nCount > EmptyThreshold)
      ),
      colour = "darkorchid", size = 3, vjust = -0.5, hjust = 0, angle = 90
    ) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks() +
    labs(
      x = "UMI counts", y = "-log10(Probability)",
      title = "The negative log-probability of being empty-dropets",
      subtitle = "droplets detected as cells should show up with \nlarge negative log-probabilities"
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
}
plotlist[["emptyDrops-Probability_vs_counts"]] <- p

if (CellLabel != "NULL") {
  p <- colData(raw) %>%
    as.data.frame() %>%
    ggplot(
      aes(x = nCount_rank, y = nCount)
    ) +
    geom_segment(aes(xend = nCount_rank * 2, yend = nCount * 2), color = "grey90", alpha = 0.01) +
    geom_point(
      aes(color = EmptyDrops),
      alpha = 0.5, shape = 16
    ) +
    scale_color_manual(
      name = "is cell",
      values = setNames(color[c(9, 2)], c(TRUE, FALSE))
    ) +
    geom_point(
      aes(x = nCount_rank * 2, y = nCount * 2, fill = CellLabel),
      color = "transparent",
      shape = 21
    ) +
    scale_fill_manual(values = setNames(c(alpha("red", 0.05), "transparent"), c(TRUE, FALSE))) +
    guides(fill = guide_legend(override.aes = list(fill = c("red", "transparent")))) +
    geom_hline(
      yintercept = metadata(bc_ranks)$knee,
      colour = "dodgerblue", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = metadata(bc_ranks)$inflection,
      colour = "forestgreen", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = EmptyThreshold,
      colour = "darkorchid", linetype = "dashed"
    ) +
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks() +
    labs(title = "Barcode Rank Plot", subtitle = paste("cells:", sum(colData(raw)$EmptyDrops)), x = "Barcode rank", y = "UMI counts") +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
} else {
  p <- colData(raw) %>%
    as.data.frame() %>%
    ggplot(
      aes(x = nCount_rank, y = nCount)
    ) +
    geom_point(
      aes(color = EmptyDrops),
      alpha = 0.5, shape = 16
    ) +
    scale_color_manual(
      name = "is cell",
      values = setNames(color[c(9, 2)], c(TRUE, FALSE))
    ) +
    geom_hline(
      yintercept = metadata(bc_ranks)$knee,
      colour = "dodgerblue", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = metadata(bc_ranks)$inflection,
      colour = "forestgreen", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = EmptyThreshold,
      colour = "darkorchid", linetype = "dashed"
    ) +
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks() +
    labs(title = "Barcode Rank Plot", subtitle = paste("cells:", sum(colData(raw)$EmptyDrops)), x = "Barcode rank", y = "UMI counts") +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
}

plotlist[["emptyDrops-nCount_rank"]] <- p

p <- plot_grid(
  plotlist = plotlist[c(
    "emptyDrops-nCount_rankFitted",
    "emptyDrops-emptyPValue",
    "emptyDrops-Probability_vs_counts",
    "emptyDrops-nCount_rank"
  )],
  nrow = 2, ncol = 2, align = "hv", axis = "tblr"
)
title <- ggdraw() +
  draw_label(
    label = "Method: emptyDrops",
    fontface = "bold", x = 0, hjust = 0
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 7)
  )
p <- plot_grid(
  title, p,
  ncol = 1,
  rel_heights = c(0.05, 1)
)
ggsave(p, filename = paste0(sample, ".emptyDrops.png"), width = 11, height = 8)

# Method: dropEst ---------------------------------------------------------
path_dropEst <- paste0(cellranger_dir, "/", sample, "/dropEst/cell.counts.rds")
dropest <- readRDS(path_dropEst)
cell_num <- EstimateCellsNumber(dropest$aligned_umis_per_cell)
chr <- colnames(dropest$reads_per_chr_per_cells$Exon)
chr_mito <- chr[grep(x = chr, pattern = "(^chrM)|(^chrMT)|(^MT)|(^mt)|(^mito)|(^Mito)")]
if (length(chr_mito) != 0) {
  lq_cells_df <- PrepareLqCellsDataPipeline(data = dropest, mit.chromosome.name = chr_mito)
} else {
  lq_cells_df <- PrepareLqCellsDataPipeline(data = dropest)
}
scores <- ScoreQualityData(
  umi.counts = dropest$aligned_umis_per_cell,
  scoring.data = lq_cells_df, cell.number = cell_num
)
qualified <- intersect(
  x = names(sort(dropest$aligned_umis_per_cell, decreasing = TRUE)[1:cell_num$min]),
  y = names(scores[scores >= dropest_quality])
)
colData(raw)$Score <- scores[colData(raw)$Barcode]
colData(raw)$dropEst <- colData(raw)$Barcode %in% qualified & colData(raw)$nCount > EmptyThreshold

if (CellLabel != "NULL") {
  p <- colData(raw) %>%
    as.data.frame() %>%
    ggplot() +
    geom_point(aes(x = nCount_rank, y = Score, color = dropEst),
      alpha = 0.5, shape = 16
    ) +
    scale_color_manual(
      name = "is cell",
      values = setNames(color[c(9, 2)], c(TRUE, FALSE))
    ) +
    geom_point(
      aes(x = nCount_rank, y = Score - 0.1, fill = CellLabel),
      color = "transparent",
      shape = 21
    ) +
    scale_fill_manual(values = setNames(c(alpha("red", 0.05), "transparent"), c(TRUE, FALSE))) +
    guides(fill = guide_legend(override.aes = list(fill = c("red", "transparent")))) +
    geom_hline(yintercept = dropest_quality, linetype = "dashed") +
    annotate("text",
      x = 1, y = dropest_quality,
      label = paste("Score >=", dropest_quality),
      colour = "black", size = 3, hjust = 0, vjust = 1.5
    ) +
    geom_vline(
      xintercept = cell_num$min,
      colour = "red3", linetype = "dashed"
    ) +
    geom_vline(
      xintercept = cell_num$max,
      colour = "royalblue4", linetype = "dashed"
    ) +
    annotate("text",
      x = cell_num$min, y = 0,
      label = paste0(
        "High-Quality",
        " Rank=", cell_num$min
      ),
      colour = "red3", size = 3, hjust = 0, vjust = -0.5, angle = 90
    ) +
    annotate("text",
      x = cell_num$max, y = 0,
      label = paste0(
        "Low-Quality",
        " Rank=", cell_num$max
      ),
      colour = "royalblue4", size = 3, hjust = 0, vjust = 1.5, angle = 90
    ) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    labs(title = "Cell Score", x = "Barcode rank", y = "Quality score") +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
} else {
  p <- colData(raw) %>%
    as.data.frame() %>%
    ggplot() +
    geom_point(aes(x = nCount_rank, y = Score, color = dropEst),
      alpha = 0.5, shape = 16
    ) +
    scale_color_manual(
      name = "is cell",
      values = setNames(color[c(9, 2)], c(TRUE, FALSE))
    ) +
    geom_hline(yintercept = dropest_quality, linetype = "dashed") +
    annotate("text",
      x = 1, y = dropest_quality,
      label = paste("Score >=", dropest_quality),
      colour = "black", size = 3, hjust = 0, vjust = 1.5
    ) +
    geom_vline(
      xintercept = cell_num$min,
      colour = "red3", linetype = "dashed"
    ) +
    geom_vline(
      xintercept = cell_num$max,
      colour = "royalblue4", linetype = "dashed"
    ) +
    annotate("text",
      x = cell_num$min, y = 0,
      label = paste0(
        "High-Quality",
        " Rank=", cell_num$min
      ),
      colour = "red3", size = 3, hjust = 0, vjust = -0.5, angle = 90
    ) +
    annotate("text",
      x = cell_num$max, y = 0,
      label = paste0(
        "Low-Quality",
        " Rank=", cell_num$max
      ),
      colour = "royalblue4", size = 3, hjust = 0, vjust = 1.5, angle = 90
    ) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    labs(title = "Cell Score", x = "Barcode rank", y = "Quality score") +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
}
plotlist[["dropEst-CellScore"]] <- p

if (CellLabel != "NULL") {
  p <- colData(raw) %>%
    as.data.frame() %>%
    ggplot(
      aes(x = nCount_rank, y = nCount)
    ) +
    geom_segment(aes(xend = nCount_rank * 2, yend = nCount * 2), color = "grey90", alpha = 0.01) +
    geom_point(
      aes(color = dropEst),
      alpha = 0.5, shape = 16
    ) +
    scale_color_manual(
      name = "is cell",
      values = setNames(color[c(9, 2)], c(TRUE, FALSE))
    ) +
    geom_point(
      aes(x = nCount_rank * 2, y = nCount * 2, fill = CellLabel),
      color = "transparent",
      shape = 21
    ) +
    scale_fill_manual(values = setNames(c(alpha("red", 0.05), "transparent"), c(TRUE, FALSE))) +
    guides(fill = guide_legend(override.aes = list(fill = c("red", "transparent")))) +
    geom_vline(
      xintercept = cell_num$min,
      colour = "red3", linetype = "dashed"
    ) +
    geom_vline(
      xintercept = cell_num$max,
      colour = "royalblue4", linetype = "dashed"
    ) +
    annotate("text",
      x = cell_num$min, y = 1,
      label = paste0(
        "High-Quality:",
        " Rank=", cell_num$min
      ),
      colour = "red3", size = 3, hjust = 0, vjust = -0.5, angle = 90
    ) +
    annotate("text",
      x = cell_num$max, y = 1,
      label = paste0(
        "Low-Quality:",
        " Rank=", cell_num$max
      ),
      colour = "royalblue4", size = 3, hjust = 0, vjust = 1.5, angle = 90
    ) +
    geom_hline(
      yintercept = metadata(bc_ranks)$knee,
      colour = "dodgerblue", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = metadata(bc_ranks)$inflection,
      colour = "forestgreen", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = EmptyThreshold,
      colour = "darkorchid", linetype = "dashed"
    ) +
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks() +
    labs(title = "Barcode Rank Plot", subtitle = paste("cells:", sum(colData(raw)$dropEst)), x = "Barcode rank", y = "UMI counts") +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
} else {
  p <- colData(raw) %>%
    as.data.frame() %>%
    ggplot(
      aes(x = nCount_rank, y = nCount)
    ) +
    geom_point(
      aes(color = dropEst),
      alpha = 0.5, shape = 16
    ) +
    scale_color_manual(
      name = "is cell",
      values = setNames(color[c(9, 2)], c(TRUE, FALSE))
    ) +
    geom_vline(
      xintercept = cell_num$min,
      colour = "red3", linetype = "dashed"
    ) +
    geom_vline(
      xintercept = cell_num$max,
      colour = "royalblue4", linetype = "dashed"
    ) +
    annotate("text",
      x = cell_num$min, y = 1,
      label = paste0(
        "High-Quality:",
        " Rank=", cell_num$min
      ),
      colour = "red3", size = 3, hjust = 0, vjust = -0.5, angle = 90
    ) +
    annotate("text",
      x = cell_num$max, y = 1,
      label = paste0(
        "Low-Quality:",
        " Rank=", cell_num$max
      ),
      colour = "royalblue4", size = 3, hjust = 0, vjust = 1.5, angle = 90
    ) +
    geom_hline(
      yintercept = metadata(bc_ranks)$knee,
      colour = "dodgerblue", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = metadata(bc_ranks)$inflection,
      colour = "forestgreen", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = EmptyThreshold,
      colour = "darkorchid", linetype = "dashed"
    ) +
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks() +
    labs(title = "Barcode Rank Plot", subtitle = paste("cells:", sum(colData(raw)$dropEst)), x = "Barcode rank", y = "UMI counts") +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
}
plotlist[["dropEst-nCount_rank"]] <- p

p <- plot_grid(
  plotlist = plotlist[c(
    "dropEst-CellScore",
    "dropEst-nCount_rank"
  )],
  nrow = 1, align = "hv", axis = "tblr"
)
title <- ggdraw() +
  draw_label(
    label = "Method: dropEst",
    fontface = "bold", x = 0, hjust = 0
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 7)
  )
p <- plot_grid(
  title, p,
  ncol = 1,
  rel_heights = c(0.05, 1)
)
ggsave(p, filename = paste0(sample, ".dropEst.png"), width = 11, height = 8)

# Method: zUMIs ---------------------------------------------------------
library(inflection)
bccount <- colData(raw) %>%
  as.data.frame() %>%
  filter(nCount_rank < quantile(nCount_rank, 0.9)) %>%
  arrange(nCount_rank) %>%
  mutate(nCount_cumsum = cumsum(nCount))
if (length(unique(as.integer(quantile(bccount[, "nCount_rank"])))) != 5) {
  bccount[nrow(bccount), "nCount_rank"] <- bccount[nrow(bccount), "nCount_rank"] + 1
}
ntop <- floor(uik(x = bccount[, "nCount_rank"], y = bccount[, "nCount_cumsum"]))
cutoff_counts <- max(c(bccount[ntop, "nCount"], EmptyThreshold))
colData(raw)$zUMIs <- colData(raw)$nCount > cutoff_counts

p <- bccount %>%
  mutate(zUMIs = nCount > cutoff_counts) %>%
  arrange(nCount_rank) %>%
  mutate(nCount_rank_uniq = 1:n()) %>%
  ggplot(aes(x = nCount_rank_uniq, y = nCount)) +
  geom_ribbon(aes(ymin = 0, ymax = nCount, fill = zUMIs), alpha = 0.5) +
  geom_line() +
  geom_vline(xintercept = ntop, color = "red3") +
  scale_fill_manual(
    name = "is cell",
    values = setNames(color[c(9, 2)], c(TRUE, FALSE))
  ) +
  labs(title = "UMI counts vs barcode rank ", x = "Barcode rank (unique)", y = "UMI counts") +
  theme_classic() +
  theme(
    aspect.ratio = 1,
  )
plotlist[["zUMIs-UMIbar"]] <- p

if (CellLabel != "NULL") {
  p <- bccount %>%
    mutate(zUMIs = nCount > cutoff_counts) %>%
    ggplot(aes(x = nCount_rank, y = nCount_cumsum)) +
    geom_segment(aes(xend = nCount_rank, yend = nCount_cumsum + max(nCount_cumsum) * 0.1), color = "grey90", alpha = 0.01) +
    geom_point(aes(color = zUMIs),
      alpha = 0.5, shape = 16
    ) +
    scale_color_manual(
      name = "is cell",
      values = setNames(color[c(9, 2)], c(TRUE, FALSE))
    ) +
    geom_point(
      aes(x = nCount_rank, y = nCount_cumsum + max(nCount_cumsum) * 0.1, fill = CellLabel),
      color = "transparent",
      shape = 21
    ) +
    scale_fill_manual(values = setNames(c(alpha("red", 0.01), "transparent"), c(TRUE, FALSE))) +
    guides(fill = guide_legend(override.aes = list(fill = c("red", "transparent")))) +
    geom_vline(xintercept = ntop, color = "red3") +
    labs(title = "Cumulative distribution plot", x = "Barcode rank", y = "Cumulative UMI counts") +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
} else {
  p <- bccount %>%
    mutate(zUMIs = nCount > cutoff_counts) %>%
    ggplot(aes(x = nCount_rank, y = nCount_cumsum)) +
    geom_point(aes(color = zUMIs),
      alpha = 0.5, shape = 16
    ) +
    scale_color_manual(
      name = "is cell",
      values = setNames(color[c(9, 2)], c(TRUE, FALSE))
    ) +
    geom_vline(xintercept = ntop, color = "red3") +
    labs(title = "Cumulative distribution plot", x = "Barcode rank", y = "Cumulative UMI counts") +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
}

plotlist[["zUMIs-CumsumPlot"]] <- p

p <- colData(raw) %>%
  as.data.frame() %>%
  filter(nCount_rank < quantile(nCount_rank, 0.9)) %>%
  ggplot(aes(x = nCount)) +
  geom_histogram(aes(fill = ..count..), bins = 50, color = "black") +
  geom_vline(xintercept = cutoff_counts, color = "red3") +
  scale_fill_material(name = "Count", palette = "blue-grey") +
  scale_x_continuous(
    trans = log10_trans(),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  annotation_logticks(sides = "b") +
  labs(
    x = "UMI counts", y = "Density",
    title = "Density plot of the UMI counts",
    subtitle = "red line indicated the intact-cell threshold."
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    legend.position = "none"
  )
plotlist[["zUMIs-UMIdensity"]] <- p

if (CellLabel != "NULL") {
  p <- colData(raw) %>%
    as.data.frame() %>%
    ggplot(
      aes(x = nCount_rank, y = nCount)
    ) +
    geom_segment(aes(xend = nCount_rank * 2, yend = nCount * 2), color = "grey90", alpha = 0.01) +
    geom_point(
      aes(color = zUMIs),
      alpha = 0.5, shape = 16
    ) +
    scale_color_manual(
      name = "is cell",
      values = setNames(color[c(9, 2)], c(TRUE, FALSE))
    ) +
    geom_point(
      aes(x = nCount_rank * 2, y = nCount * 2, fill = CellLabel),
      color = "transparent",
      shape = 21
    ) +
    scale_fill_manual(values = setNames(c(alpha("red", 0.05), "transparent"), c(TRUE, FALSE))) +
    guides(fill = guide_legend(override.aes = list(fill = c("red", "transparent")))) +
    geom_vline(xintercept = ntop, color = "red3") +
    annotate("text",
      x = ntop, y = 1,
      label = paste0(
        "zUMIs threshold:",
        " UMI=", cutoff_counts,
        " Rank=", ntop
      ),
      colour = "red3", size = 3, vjust = -0.5, hjust = 0, angle = 90
    ) +
    geom_hline(
      yintercept = metadata(bc_ranks)$knee,
      colour = "dodgerblue", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = metadata(bc_ranks)$inflection,
      colour = "forestgreen", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = EmptyThreshold,
      colour = "darkorchid", linetype = "dashed"
    ) +
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks() +
    labs(title = "Barcode Rank Plot", subtitle = paste("cells:", sum(colData(raw)$zUMIs)), x = "Barcode rank", y = "UMI counts") +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
} else {
  p <- colData(raw) %>%
    as.data.frame() %>%
    ggplot(
      aes(x = nCount_rank, y = nCount)
    ) +
    geom_point(
      aes(color = zUMIs),
      alpha = 0.5, shape = 16
    ) +
    scale_color_manual(
      name = "is cell",
      values = setNames(color[c(9, 2)], c(TRUE, FALSE))
    ) +
    geom_vline(xintercept = ntop, color = "red3") +
    annotate("text",
      x = ntop, y = 1,
      label = paste0(
        "zUMIs threshold:",
        " UMI=", cutoff_counts,
        " Rank=", ntop
      ),
      colour = "red3", size = 3, vjust = -0.5, hjust = 0, angle = 90
    ) +
    geom_hline(
      yintercept = metadata(bc_ranks)$knee,
      colour = "dodgerblue", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = metadata(bc_ranks)$inflection,
      colour = "forestgreen", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = EmptyThreshold,
      colour = "darkorchid", linetype = "dashed"
    ) +
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks() +
    labs(title = "Barcode Rank Plot", subtitle = paste("cells:", sum(colData(raw)$zUMIs)), x = "Barcode rank", y = "UMI counts") +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
}
plotlist[["zUMIs-nCount_rank"]] <- p

p <- plot_grid(
  plotlist = plotlist[c(
    "zUMIs-UMIbar",
    "zUMIs-CumsumPlot",
    "zUMIs-UMIdensity",
    "zUMIs-nCount_rank"
  )],
  nrow = 2, align = "hv", axis = "tblr"
)
title <- ggdraw() +
  draw_label(
    label = "Method: zUMIs",
    fontface = "bold", x = 0, hjust = 0
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 7)
  )
p <- plot_grid(
  title, p,
  ncol = 1,
  rel_heights = c(0.05, 1)
)
ggsave(p, filename = paste0(sample, ".zUMIs.png"), width = 11, height = 8)

# Method: dropSplit ---------------------------------------------------------
library(dropSplit)
result <- dropSplit(counts = counts(raw), Redundancy_control = T)
meta_info <- result$meta_info
colData(raw)[rownames(meta_info), "RankMSE"] <- meta_info$RankMSE
colData(raw)[rownames(meta_info), "CellEntropy"] <- meta_info$CellEntropy
colData(raw)[rownames(meta_info), "CellRedundancy"] <- meta_info$CellRedundancy
colData(raw)[rownames(meta_info), "CellGini"] <- meta_info$CellGini
colData(raw)[rownames(meta_info), "dropSplitClass"] <- as.character(meta_info$dropSplitClass)
colData(raw)[rownames(meta_info), "dropSplitScore"] <- meta_info$dropSplitScore
colData(raw)[rownames(meta_info), "dropSplit"] <- ifelse(meta_info$dropSplitClass == "Cell", TRUE, FALSE)
cutoff_counts <- min(meta_info$nCount[meta_info$preDefinedClass == "Cell"])
Cell_rank <- max(meta_info$nCount_rank[meta_info$preDefinedClass == "Cell"])
Uncertain_rank <- max(meta_info$nCount_rank[meta_info$preDefinedClass == "Uncertain"])
Empty_rank <- max(meta_info$nCount_rank[meta_info$preDefinedClass == "Empty"])


if (CellLabel != "NULL") {
  p <- colData(raw) %>%
    as.data.frame() %>%
    ggplot(aes(x = nCount_rank, y = nFeature_rank)) +
    geom_point(
      aes(color = dropSplit),
      alpha = 0.5, shape = 16
    ) +
    scale_color_manual(
      name = "is cell",
      values = setNames(color[c(9, 2)], c(TRUE, FALSE))
    ) +
    geom_point(
      aes(x = nCount_rank, y = nFeature_rank * 10, fill = CellLabel),
      color = "transparent",
      shape = 21
    ) +
    scale_fill_manual(values = setNames(c(alpha("red", 0.05), "transparent"), c(TRUE, FALSE))) +
    guides(fill = guide_legend(override.aes = list(fill = c("red", "transparent")))) +
    geom_vline(
      xintercept = c(Cell_rank, Uncertain_rank, Empty_rank),
      color = c("red3", "forestgreen", "steelblue")
    ) +
    annotate("text",
      x = Cell_rank, y = 1,
      label = paste0(
        "preDefined-Cell threshold:",
        " Rank=", Cell_rank
      ),
      colour = "red3", size = 3, vjust = -0.5, hjust = 0, angle = 90
    ) +
    annotate("text",
      x = Uncertain_rank, y = 1,
      label = paste0(
        "preDefined-Uncertain threshold:",
        " Rank=", Uncertain_rank
      ),
      colour = "forestgreen", size = 3, vjust = -0.5, hjust = 0, angle = 90
    ) +
    annotate("text",
      x = Empty_rank, y = 1,
      label = paste0(
        "preDefined-Empty threshold:",
        " Rank=", Empty_rank
      ),
      colour = "steelblue", size = 3, vjust = -0.5, hjust = 0, angle = 90
    ) +
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks() +
    labs(title = "nCount-nFeature Rank Plot", subtitle = paste("cells:", sum(colData(raw)$dropSplit)), x = "nCount_rank", y = "nFeature_rank") +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
} else {
  p <- colData(raw) %>%
    as.data.frame() %>%
    ggplot(aes(x = nCount_rank, y = nFeature_rank)) +
    geom_point(
      aes(color = dropSplit),
      alpha = 0.5, shape = 16
    ) +
    scale_color_manual(
      name = "is cell",
      values = setNames(color[c(9, 2)], c(TRUE, FALSE))
    ) +
    geom_vline(
      xintercept = c(Cell_rank, Uncertain_rank, Empty_rank),
      color = c("red3", "forestgreen", "steelblue")
    ) +
    annotate("text",
      x = Cell_rank, y = 1,
      label = paste0(
        "preDefined-Cell threshold:",
        " Rank=", Cell_rank
      ),
      colour = "red3", size = 3, vjust = -0.5, hjust = 0, angle = 90
    ) +
    annotate("text",
      x = Uncertain_rank, y = 1,
      label = paste0(
        "preDefined-Uncertain threshold:",
        " Rank=", Uncertain_rank
      ),
      colour = "forestgreen", size = 3, vjust = -0.5, hjust = 0, angle = 90
    ) +
    annotate("text",
      x = Empty_rank, y = 1,
      label = paste0(
        "preDefined-Empty threshold:",
        " Rank=", Empty_rank
      ),
      colour = "steelblue", size = 3, vjust = -0.5, hjust = 0, angle = 90
    ) +
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks() +
    labs(title = "nCount-nFeature Rank Plot", subtitle = paste("cells:", sum(colData(raw)$dropSplit)), x = "nCount_rank", y = "nFeature_rank") +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
}
plotlist[["dropSplit-Rank"]] <- p

if (CellLabel != "NULL") {
  p <- colData(raw) %>%
    as.data.frame() %>%
    ggplot(aes(x = nCount_rank, y = RankMSE)) +
    geom_point(
      aes(color = dropSplit),
      alpha = 0.5, shape = 16
    ) +
    scale_color_manual(
      name = "is cell",
      values = setNames(color[c(9, 2)], c(TRUE, FALSE))
    ) +
    geom_point(
      aes(x = nCount_rank, y = RankMSE * 10, fill = CellLabel),
      color = "transparent",
      shape = 21
    ) +
    scale_fill_manual(values = setNames(c(alpha("red", 0.05), "transparent"), c(TRUE, FALSE))) +
    guides(fill = guide_legend(override.aes = list(fill = c("red", "transparent")))) +
    geom_vline(
      xintercept = c(Cell_rank, Uncertain_rank, Empty_rank),
      color = c("red3", "forestgreen", "steelblue")
    ) +
    annotate("text",
      x = Cell_rank, y = 1,
      label = paste0(
        "preDefined-Cell threshold:",
        " Rank=", Cell_rank
      ),
      colour = "red3", size = 3, vjust = -0.5, hjust = 0, angle = 90
    ) +
    annotate("text",
      x = Uncertain_rank, y = 1,
      label = paste0(
        "preDefined-Uncertain threshold:",
        " Rank=", Uncertain_rank
      ),
      colour = "forestgreen", size = 3, vjust = -0.5, hjust = 0, angle = 90
    ) +
    annotate("text",
      x = Empty_rank, y = 1,
      label = paste0(
        "preDefined-Empty threshold:",
        " Rank=", Empty_rank
      ),
      colour = "steelblue", size = 3, vjust = -0.5, hjust = 0, angle = 90
    ) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks() +
    labs(title = "Mean Squared Error of nCount/nFeature Rank", subtitle = paste("cells:", sum(colData(raw)$dropSplit)), x = "nCount_rank", y = "RankMSE") +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
} else {
  p <- colData(raw) %>%
    as.data.frame() %>%
    ggplot(aes(x = nCount_rank, y = RankMSE)) +
    geom_point(
      aes(color = dropSplit),
      alpha = 0.5, shape = 16
    ) +
    scale_color_manual(
      name = "is cell",
      values = setNames(color[c(9, 2)], c(TRUE, FALSE))
    ) +
    geom_vline(
      xintercept = c(Cell_rank, Uncertain_rank, Empty_rank),
      color = c("red3", "forestgreen", "steelblue")
    ) +
    annotate("text",
      x = Cell_rank, y = 1,
      label = paste0(
        "preDefined-Cell threshold:",
        " Rank=", Cell_rank
      ),
      colour = "red3", size = 3, vjust = -0.5, hjust = 0, angle = 90
    ) +
    annotate("text",
      x = Uncertain_rank, y = 1,
      label = paste0(
        "preDefined-Uncertain threshold:",
        " Rank=", Uncertain_rank
      ),
      colour = "forestgreen", size = 3, vjust = -0.5, hjust = 0, angle = 90
    ) +
    annotate("text",
      x = Empty_rank, y = 1,
      label = paste0(
        "preDefined-Empty threshold:",
        " Rank=", Empty_rank
      ),
      colour = "steelblue", size = 3, vjust = -0.5, hjust = 0, angle = 90
    ) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks() +
    labs(title = "Mean Squared Error of nCount/nFeature Rank within a window", subtitle = paste("cells:", sum(colData(raw)$dropSplit)), x = "nCount_rank", y = "RankMSE") +
    theme_classic() +
    theme(
      aspect.ratio = 1,
    )
}
plotlist[["dropSplit-RankMSE"]] <- p


p <- colData(raw) %>%
  as.data.frame() %>%
  subset(!is.na(CellEntropy)) %>%
  mutate(dropSplitClass = factor(dropSplitClass, levels = c("Cell", "Uncertain", "Empty", "Discarded"))) %>%
  ggplot(aes(x = nCount, y = CellEntropy)) +
  geom_point(
    aes(color = dropSplitClass),
    alpha = 0.5, shape = 16
  ) +
  scale_color_manual(
    name = "is cell",
    values = setNames(
      c("red3", "forestgreen", "steelblue", "grey80"),
      c("Cell", "Uncertain", "Empty", "Discarded")
    )
  ) +
  scale_x_continuous(
    trans = log10_trans(),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  annotation_logticks(sides = "b") +
  labs(title = "Cell Entropy Plot", subtitle = paste("cells:", sum(colData(raw)$dropSplit)), x = "nCount", y = "CellEntropy") +
  theme_classic() +
  theme(
    aspect.ratio = 1,
  )

plotlist[["dropSplit-CellEntropy"]] <- p


p <- colData(raw) %>%
  as.data.frame() %>%
  subset(!is.na(CellEntropy)) %>%
  ggplot(aes(x = nCount, y = CellEntropy)) +
  geom_point(
    aes(color = dropSplitScore),
    alpha = 0.5, shape = 16
  ) +
  scale_color_viridis_c(
    name = "dropSplitScore"
  ) +
  scale_x_continuous(
    trans = log10_trans(),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  annotation_logticks(sides = "b") +
  labs(title = "Cell Entropy Plot", subtitle = paste("cells:", sum(colData(raw)$dropSplit)), x = "nCount", y = "CellEntropy") +
  theme_classic() +
  theme(
    aspect.ratio = 1,
  )
plotlist[["dropSplit-CellEntropy2"]] <- p

p <- plot_grid(
  plotlist = plotlist[c(
    "dropSplit-Rank",
    "dropSplit-RankMSE",
    "dropSplit-CellEntropy",
    "dropSplit-CellEntropy2"
  )],
  nrow = 2, align = "hv", axis = "tblr"
)
title <- ggdraw() +
  draw_label(
    label = "Method: dropSplit",
    fontface = "bold", x = 0, hjust = 0
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 7)
  )
p <- plot_grid(
  title, p,
  ncol = 1,
  rel_heights = c(0.05, 1)
)
ggsave(p, filename = paste0(sample, ".dropSplit.png"), width = 11, height = 8)


# Summary -----------------------------------------------------------------
cell_rank <- colData(raw) %>%
  as.data.frame() %>%
  dplyr::select(
    Barcode, nCount_rank, nCount, Fitted, CellLabel,
    Cellranger_v2, Cellranger_v3, EmptyDrops, dropEst, zUMIs, dropSplit
  ) %>%
  arrange(nCount_rank) %>%
  reshape2::melt(
    measure.vars = c("Cellranger_v2", "Cellranger_v3", "EmptyDrops", "dropEst", "zUMIs", "dropSplit"),
    variable.name = "Method",
    value.name = "is_cell"
  ) %>%
  group_by(Method) %>%
  mutate(Method = paste0(Method, ": ", sum(is_cell)))
cell_rank[, "Method"] <- factor(pull(cell_rank, "Method"),
  levels = unique(pull(cell_rank, "Method"))
)

if (CellLabel != "NULL") {
  p <- ggplot(cell_rank, aes(x = nCount_rank, y = nCount)) +
    geom_segment(aes(xend = nCount_rank * 2, yend = nCount * 2), color = "grey90", alpha = 0.01) +
    geom_point(
      aes(color = is_cell),
      alpha = 0.5, shape = 16
    ) +
    scale_color_manual(
      values = setNames(color[c(9, 2)], c(TRUE, FALSE)),
      guide = FALSE
    ) +
    geom_point(
      aes(x = nCount_rank * 2, y = nCount * 2, fill = CellLabel),
      color = "transparent",
      shape = 21
    ) +
    scale_fill_manual(values = setNames(c(alpha("red", 0.05), "transparent"), c(TRUE, FALSE)), guide = FALSE) +
    geom_hline(
      yintercept = metadata(bc_ranks)$knee,
      colour = "dodgerblue", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = metadata(bc_ranks)$inflection,
      colour = "forestgreen", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = EmptyThreshold,
      colour = "darkorchid", linetype = "dashed"
    ) +
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks() +
    labs(title = "Barcode Rank Plot", x = "Barcode rank", y = "UMI counts") +
    facet_wrap(. ~ Method, nrow = 1) +
    theme_classic() +
    theme(
      aspect.ratio = 1
    )
} else {
  p <- ggplot(cell_rank, aes(x = nCount_rank, y = nCount)) +
    geom_point(
      aes(color = is_cell),
      alpha = 0.5, shape = 16
    ) +
    scale_color_manual(
      values = setNames(color[c(9, 2)], c(TRUE, FALSE)),
      guide = FALSE
    ) +
    geom_hline(
      yintercept = metadata(bc_ranks)$knee,
      colour = "dodgerblue", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = metadata(bc_ranks)$inflection,
      colour = "forestgreen", linetype = "dashed"
    ) +
    geom_hline(
      yintercept = EmptyThreshold,
      colour = "darkorchid", linetype = "dashed"
    ) +
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks() +
    labs(title = "Barcode Rank Plot", x = "Barcode rank", y = "UMI counts") +
    facet_wrap(. ~ Method, nrow = 1) +
    theme_classic() +
    theme(
      aspect.ratio = 1
    )
}
plotlist[["MethodCompare-nCount_rank"]] <- p

cell_count <- colData(raw) %>%
  as.data.frame() %>%
  dplyr::select(Barcode, Cellranger_v2, Cellranger_v3, EmptyDrops, dropEst, zUMIs, dropSplit) %>%
  reshape2::melt(
    measure.vars = c("Cellranger_v2", "Cellranger_v3", "EmptyDrops", "dropEst", "zUMIs", "dropSplit"),
    variable.name = "Method",
    value.name = "is_cell"
  ) %>%
  filter(is_cell == TRUE) %>%
  group_by(Method) %>%
  summarise(Method = Method, cellnum = n()) %>%
  distinct()

p <- ggplot(
  cell_count,
  aes(x = Method, y = cellnum, fill = cellnum, label = cellnum)
) +
  geom_bar(color = "black", stat = "identity") +
  geom_text(vjust = -0.5, hjust = 0, angle = 45) +
  scale_fill_material(name = "Count", palette = "blue-grey", guide = FALSE) +
  scale_y_continuous(limits = c(0, 1.2 * max(cell_count[, "cellnum"]))) +
  labs(title = "Cell-containing droplet identification", x = "", y = "Cell number") +
  theme_classic() +
  theme(
    aspect.ratio = 0.8,
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )
plotlist[["MethodCompare-Barplot"]] <- p

cell_upset <- colData(raw) %>%
  as.data.frame() %>%
  dplyr::select(Barcode, Cellranger_v2, Cellranger_v3, EmptyDrops, dropEst, zUMIs, dropSplit) %>%
  reshape2::melt(
    measure.vars = c("Cellranger_v2", "Cellranger_v3", "EmptyDrops", "dropEst", "zUMIs", "dropSplit"),
    variable.name = "Method",
    value.name = "is_cell"
  ) %>%
  dplyr::filter(is_cell == TRUE) %>%
  group_by(Barcode) %>%
  summarize(
    Method_list = list(Method),
    Method_comb = paste(Method, collapse = ","),
    Method_num = n()
  )
y_max <- max(table(pull(cell_upset, "Method_comb")))

p <- cell_upset %>%
  group_by(Method_comb) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  filter(count %in% head(sort(unique(count), decreasing = T), 10)) %>%
  ggplot(aes(x = Method_list)) +
  geom_bar(aes(fill = ..count..), color = "black", width = 0.5) +
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.5, hjust = 0, angle = 45) +
  labs(title = "Cell intersection among differnent methods", x = "", y = "Cell number") +
  scale_x_upset(sets = c("Cellranger_v2", "Cellranger_v3", "EmptyDrops", "dropEst", "zUMIs", "dropSplit")) +
  scale_y_continuous(limits = c(0, 1.2 * y_max)) +
  scale_fill_material(name = "Count", palette = "blue-grey") +
  theme_combmatrix(
    combmatrix.label.text = element_text(size = 10, color = "black"),
    combmatrix.label.extra_spacing = 6
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 0.6,
    legend.position = "none"
  )
plotlist[["MethodCompare-upsetR"]] <- p

p <- plot_grid(
  plotlist[["MethodCompare-nCount_rank"]],
  plot_grid(plotlist[["MethodCompare-Barplot"]],
    plotlist[["MethodCompare-upsetR"]],
    nrow = 1, rel_widths = c(0.5, 0.5)
  ),
  nrow = 2
)
title <- ggdraw() +
  draw_label(
    label = "Summary",
    fontface = "bold", x = 0, hjust = 0
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 7)
  )
p <- plot_grid(
  title, p,
  ncol = 1,
  rel_heights = c(0.05, 1)
)
ggsave(p, filename = paste0(sample, ".MethodCompare.png"), width = 11, height = 8)


# Output the report -------------------------------------------------------
saveRDS(raw, file = "raw.rds")
saveRDS(cell_upset, file = "cell_upset.rds")
write.table(x = cell_upset[["Barcode"]], file = "barcodes.tsv", row.names = FALSE, col.names = FALSE, quote = FALSE)

pdf(paste0(sample, ".CellCalling.pdf"), width = 11, height = 8)
invisible(lapply(paste0(sample, c(".Basic.png", ".emptyDrops.png", ".dropEst.png", ".zUMIs.png", ".dropSplit.png", ".MethodCompare.png")), function(x) {
  grid.arrange(rasterGrob(readPNG(x, native = FALSE),
    interpolate = FALSE
  ))
}))
invisible(dev.off())

if (file.exists("Rplots.pdf")) {
  invisible(file.remove("Rplots.pdf"))
}

cat("Cell Calling Done\n")
