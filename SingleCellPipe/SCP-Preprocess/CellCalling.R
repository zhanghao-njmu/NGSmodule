#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
cellranger_dir <- as.character(args[1])
sample <- as.character(args[2])
threads <- as.character(args[3])

# parameters: global settings ---------------------------------------------
# cellranger_dir <- "/data/lab/HuangMingQian/scRNA-seq/ESC-PGC-GSCLC-new/cellranger"
# sample <- "1_P31"
# threads <- 80

# parameters: emptyDrops --------------------------------------------------
empty_thresh <- 200
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


# Method: emptyDrops ------------------------------------------------------
bc_ranks <- barcodeRanks(counts(raw), lower = empty_thresh)
colData(raw)$UMIcounts <- bc_ranks$total
colData(raw)$BarcodeRank <- bc_ranks$rank
colData(raw)$Fitted <- bc_ranks$fitted
colData(raw)$knee <- colData(raw)$UMIcounts > metadata(bc_ranks)$knee
colData(raw)$inflection <- colData(raw)$UMIcounts > metadata(bc_ranks)$inflection
# colData(raw) %>%
#   as.data.frame() %>%
#   arrange(BarcodeRank) %>%
#   mutate(cumsum = cumsum(UMIcounts)) %>%
#   ggplot(aes(x = BarcodeRank, y = cumsum)) +
#   geom_point()

p <- colData(raw) %>%
  as.data.frame() %>%
  ggplot(
    aes(x = BarcodeRank, y = UMIcounts)
  ) +
  geom_point(
    aes(color = log10(UMIcounts)),
    alpha = 0.5, shape = 16
  ) +
  geom_line(aes(x = BarcodeRank, y = Fitted), color = "red") +
  geom_hline(
    yintercept = metadata(bc_ranks)$knee,
    colour = "dodgerblue", linetype = "dashed"
  ) +
  geom_hline(
    yintercept = metadata(bc_ranks)$inflection,
    colour = "forestgreen", linetype = "dashed"
  ) +
  geom_hline(
    yintercept = empty_thresh,
    colour = "darkorchid", linetype = "dashed"
  ) +
  annotate("text",
    x = 1, y = metadata(bc_ranks)$knee,
    label = paste0(
      "Knee:",
      " UMI=", round(metadata(bc_ranks)$knee),
      " Rank=", sum(colData(raw)$UMIcounts > metadata(bc_ranks)$knee)
    ),
    colour = "dodgerblue", size = 3, hjust = 0, vjust = -0.5
  ) +
  annotate("text",
    x = 1, y = metadata(bc_ranks)$inflection,
    label = paste0(
      "Inflection:",
      " UMI=", round(metadata(bc_ranks)$inflection),
      " Rank=", sum(colData(raw)$UMIcounts > metadata(bc_ranks)$inflection)
    ),
    colour = "forestgreen", size = 3, hjust = 0, vjust = -0.5
  ) +
  annotate("text",
    x = 1, y = empty_thresh,
    label = paste0(
      "Empty threshold:",
      " UMI=", empty_thresh,
      " Rank=", sum(colData(raw)$UMIcounts > empty_thresh)
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
  scale_color_material(palette = "blue-grey", guide = FALSE) +
  annotation_logticks() +
  labs(title = "Barcode Rank Plot", x = "Barcode rank", y = "UMI counts") +
  theme_classic() +
  theme(
    aspect.ratio = 1,
  )
plotlist[["emptyDrops-BarcodeRankFitted"]] <- p

cat("niters:", iters, "\n")
mito_gene <- grep(x = rownames(raw), pattern = "^(MT-)|(mt-)", perl = T)
ribo_gene <- grep(x = rownames(raw), pattern = "^(RP[SL]\\d+)|(rp[sl]\\d+)", perl = T)
keep <- !1:nrow(raw) %in% c(mito_gene, ribo_gene)
emp_drops <- emptyDrops(counts(raw)[keep, ],
  lower = empty_thresh, test.ambient = TRUE,
  niters = iters, BPPARAM = bpparam
)
is_cell <- emp_drops$FDR <= fdr
is_cell[is.na(is_cell)] <- FALSE
table(Limited = emp_drops$Limited, Significant = is_cell)

colData(raw)$EmpDropsLogProb <- emp_drops$LogProb
colData(raw)$EmpDropsPValue <- emp_drops$PValue
colData(raw)$EmpDropsLimited <- emp_drops$Limited
colData(raw)$EmpDropsFDR <- emp_drops$FDR
colData(raw)$EmptyDrops <- is_cell

p <- colData(raw) %>%
  as.data.frame() %>%
  dplyr::filter(
    UMIcounts <= empty_thresh,
    UMIcounts > 0
  ) %>%
  ggplot(aes(x = EmpDropsPValue)) +
  geom_histogram(aes(fill = ..count..)) +
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

p <- colData(raw) %>%
  as.data.frame() %>%
  ggplot(aes(x = UMIcounts, y = -EmpDropsLogProb)) +
  geom_point(aes(color = EmptyDrops),
    alpha = 0.5, shape = 16
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
    xintercept = empty_thresh,
    colour = "darkorchid", linetype = "dashed"
  ) +
  annotate("text",
    y = 1, x = metadata(bc_ranks)$knee,
    label = paste0(
      "Knee:",
      " UMI=", round(metadata(bc_ranks)$knee),
      " Rank=", sum(colData(raw)$UMIcounts > metadata(bc_ranks)$knee)
    ),
    colour = "dodgerblue", size = 3, vjust = 1.5, hjust = 0, angle = 90
  ) +
  annotate("text",
    y = 1, x = metadata(bc_ranks)$inflection,
    label = paste0(
      "Inflection:",
      " UMI=", round(metadata(bc_ranks)$inflection),
      " Rank=", sum(colData(raw)$UMIcounts > metadata(bc_ranks)$inflection)
    ),
    colour = "forestgreen", size = 3, vjust = 1.5, hjust = 0, angle = 90
  ) +
  annotate("text",
    y = 1, x = empty_thresh,
    label = paste0(
      "Empty threshold:",
      " UMI=", empty_thresh,
      " Rank=", sum(colData(raw)$UMIcounts > empty_thresh)
    ),
    colour = "darkorchid", size = 3, vjust = -0.5, hjust = 0, angle = 90
  ) +
  scale_color_manual(
    name = "is cell",
    values = setNames(color[c(9, 2)], c(TRUE, FALSE))
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
  labs(
    x = "UMI counts", y = "-log(Probability)",
    title = "The negative log-probability of being empty-dropets",
    subtitle = "droplets detected as cells should show up with \nlarge negative log-probabilities"
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
  )
plotlist[["emptyDrops-Probability_vs_counts"]] <- p

p <- colData(raw) %>%
  as.data.frame() %>%
  ggplot(
    aes(x = BarcodeRank, y = UMIcounts)
  ) +
  geom_point(
    aes(color = EmptyDrops),
    alpha = 0.5, shape = 16
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
    yintercept = empty_thresh,
    colour = "darkorchid", linetype = "dashed"
  ) +
  annotate("text",
    x = 1, y = metadata(bc_ranks)$knee,
    label = paste0(
      "Knee:",
      " UMI=", round(metadata(bc_ranks)$knee),
      " Rank=", sum(colData(raw)$UMIcounts > metadata(bc_ranks)$knee)
    ),
    colour = "dodgerblue", size = 3, hjust = 0, vjust = -0.5
  ) +
  annotate("text",
    x = 1, y = metadata(bc_ranks)$inflection,
    label = paste0(
      "Inflection:",
      " UMI=", round(metadata(bc_ranks)$inflection),
      " Rank=", sum(colData(raw)$UMIcounts > metadata(bc_ranks)$inflection)
    ),
    colour = "forestgreen", size = 3, hjust = 0, vjust = -0.5
  ) +
  annotate("text",
    x = 1, y = empty_thresh,
    label = paste0(
      "Empty threshold:",
      " UMI=", empty_thresh,
      " Rank=", sum(colData(raw)$UMIcounts > empty_thresh)
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
  scale_color_manual(
    name = "is cell",
    values = setNames(color[c(9, 2)], c(TRUE, FALSE))
  ) +
  annotation_logticks() +
  labs(title = "Barcode Rank Plot (emptyDrops method)", x = "Barcode rank", y = "UMI counts") +
  theme_classic() +
  theme(
    aspect.ratio = 1,
  )
plotlist[["emptyDrops-BarcodeRank"]] <- p

p <- plot_grid(
  plotlist = plotlist[c(
    "emptyDrops-BarcodeRankFitted",
    "emptyDrops-emptyPValue",
    "emptyDrops-Probability_vs_counts",
    "emptyDrops-BarcodeRank"
  )],
  nrow = 2, ncol = 2, align = "hv", axis = "tblr"
)
ggsave(p, filename = paste0(sample, ".emptyDrops.png"), width = 11, height = 8)

# Method: dropEst ---------------------------------------------------------
path_dropEst <- paste0(cellranger_dir, "/", sample, "/dropEst/cell.counts.rds")
dropest <- readRDS(path_dropEst)
cell_num <- EstimateCellsNumber(dropest$aligned_umis_per_cell)
chr <- colnames(dropest$reads_per_chr_per_cells$Exon)
chr_mito <- chr[grep(x = chr, pattern = "^(chrM)|(chrMT)|(MT)|(mt)|(mito)|(Mito)")]
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
colData(raw)$dropEst <- colData(raw)$Barcode %in% qualified

cell_score <- data.frame(colData(raw)[, c("Barcode", "BarcodeRank", "UMIcounts", "dropEst")], row.names = 1)
cell_score[, "Score"] <- scores[rownames(cell_score)]
p <- ggplot(cell_score) +
  geom_point(aes(x = BarcodeRank, y = Score, color = dropEst),
    alpha = 0.5, shape = 16
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
  scale_color_manual(
    name = "is cell",
    values = setNames(color[c(9, 2)], c(TRUE, FALSE))
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
plotlist[["dropEst-CellScore"]] <- p

p <- colData(raw) %>%
  as.data.frame() %>%
  ggplot(
    aes(x = BarcodeRank, y = UMIcounts)
  ) +
  geom_point(
    aes(color = dropEst),
    alpha = 0.5, shape = 16
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
    yintercept = empty_thresh,
    colour = "darkorchid", linetype = "dashed"
  ) +
  annotate("text",
    x = 1, y = metadata(bc_ranks)$knee,
    label = paste0(
      "Knee:",
      " UMI=", round(metadata(bc_ranks)$knee),
      " Rank=", sum(colData(raw)$UMIcounts > metadata(bc_ranks)$knee)
    ),
    colour = "dodgerblue", size = 3, hjust = 0, vjust = -0.5
  ) +
  annotate("text",
    x = 1, y = metadata(bc_ranks)$inflection,
    label = paste0(
      "Inflection:",
      " UMI=", round(metadata(bc_ranks)$inflection),
      " Rank=", sum(colData(raw)$UMIcounts > metadata(bc_ranks)$inflection)
    ),
    colour = "forestgreen", size = 3, hjust = 0, vjust = -0.5
  ) +
  annotate("text",
    x = 1, y = empty_thresh,
    label = paste0(
      "Empty threshold:",
      " UMI=", empty_thresh,
      " Rank=", sum(colData(raw)$UMIcounts > empty_thresh)
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
  scale_color_manual(
    name = "is cell",
    values = setNames(color[c(9, 2)], c(TRUE, FALSE))
  ) +
  annotation_logticks() +
  labs(title = "Barcode Rank Plot (dropEst method)", x = "Barcode rank", y = "UMI counts") +
  theme_classic() +
  theme(
    aspect.ratio = 1,
  )
plotlist[["dropEst-BarcodeRank"]] <- p

p <- plot_grid(
  plotlist = plotlist[c(
    "dropEst-CellScore",
    "dropEst-BarcodeRank"
  )],
  nrow = 1, align = "hv", axis = "tblr"
)
ggsave(p, filename = paste0(sample, ".dropEst.png"), width = 11, height = 8)


# Summary -----------------------------------------------------------------
cell_rank <- colData(raw) %>%
  as.data.frame() %>%
  dplyr::select(
    Barcode, BarcodeRank, UMIcounts, Fitted,
    Cellranger_v2, Cellranger_v3, EmptyDrops, dropEst, knee, inflection
  ) %>%
  arrange(BarcodeRank) %>%
  reshape2::melt(
    measure.vars = c("Cellranger_v2", "Cellranger_v3", "EmptyDrops", "dropEst", "knee", "inflection"),
    variable.name = "Method",
    value.name = "is_cell"
  ) %>%
  group_by(Method) %>%
  mutate(Method = paste0(Method, ": ", sum(is_cell)))
cell_rank[, "Method"] <- factor(pull(cell_rank, "Method"),
  levels = unique(pull(cell_rank, "Method"))
)
# n<- sample(c(1:nrow(cell_rank)),size=10000)

p <- ggplot(cell_rank, aes(x = BarcodeRank, y = UMIcounts)) +
  geom_point(
    aes(color = is_cell),
    alpha = 0.5, shape = 16
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
    yintercept = empty_thresh,
    colour = "darkorchid", linetype = "dashed"
  ) +
  scale_color_manual(
    values = setNames(color[c(9, 2)], c(TRUE, FALSE)),
    guide = FALSE
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
plotlist[["MethodCompare-BarcodeRank"]] <- p

cell_count <- colData(raw) %>%
  as.data.frame() %>%
  dplyr::select(Barcode, Cellranger_v2, Cellranger_v3, EmptyDrops, dropEst, knee, inflection) %>%
  reshape2::melt(
    measure.vars = c("Cellranger_v2", "Cellranger_v3", "EmptyDrops", "dropEst", "knee", "inflection"),
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
  dplyr::select(Barcode, Cellranger_v2, Cellranger_v3, EmptyDrops, dropEst, knee, inflection) %>%
  reshape2::melt(
    measure.vars = c("Cellranger_v2", "Cellranger_v3", "EmptyDrops", "dropEst", "knee", "inflection"),
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

p <- ggplot(cell_upset, aes(x = Method_list)) +
  geom_bar(aes(fill = ..count..), color = "black", width = 0.5) +
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.5, hjust = 0, angle = 45) +
  labs(title = "Cell intersection among differnent methods", x = "", y = "Cell number") +
  scale_x_upset(sets = c("Cellranger_v2", "Cellranger_v3", "EmptyDrops", "dropEst", "knee", "inflection")) +
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
  plotlist[["MethodCompare-BarcodeRank"]],
  plot_grid(plotlist[["MethodCompare-Barplot"]],
    plotlist[["MethodCompare-upsetR"]],
    nrow = 1, rel_widths = c(0.5, 0.5)
  ),
  nrow = 2
)
ggsave(p, filename = paste0(sample, ".MethodCompare.png"), width = 11, height = 8)



# Output the report -------------------------------------------------------
saveRDS(raw, file = "raw.rds")
saveRDS(cell_upset, file = "cell_upset.rds")

pdf(paste0(sample, ".CellCalling.pdf"), width = 11, height = 8)
invisible(lapply(paste0(sample, c(".emptyDrops.png", ".dropEst.png", ".MethodCompare.png")), function(x) {
  grid.arrange(rasterGrob(readPNG(x, native = FALSE),
    interpolate = FALSE
  ))
}))
invisible(dev.off())

if (file.exists("Rplots.pdf")) {
  invisible(file.remove("Rplots.pdf"))
}

cat("Successfully")
