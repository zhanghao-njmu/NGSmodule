#!/usr/bin/env Rscript
suppressWarnings(suppressPackageStartupMessages(invisible(lapply(
  c(
    "dplyr", "stringr", "ggplot2", "ggsci", "ggtree", "RColorBrewer", "cowplot",
    "aplot", "ggplotify", "edgeR", "sva", "limma", "patchwork", "ggrepel", "Rtsne",
    "plotly","plot3D","grid","ggforce"
  ),
  require,
  character.only = TRUE
))))
set.seed(11)

args <- commandArgs(trailingOnly = TRUE)
maindir <- args[1]
aligner <- args[2]
SampleInfoFile <- args[3]
script_path <- as.character(args[4])

# maindir <- "/data/database/SRR_collection/human/early_embyro/"
# aligner <- "hisat2"
# SampleInfoFile <- "/data/database/SRR_collection/human/early_embyro/temp_20200714173936.Sample_info.csv"
# script_path <- "/home/zhanghao/Program/NGS/UniversalTools/NGSmodule/Analysis/Quantification/BatchCorrection.R"

script_dir <- gsub(x = script_path, pattern = "BatchCorrection.R", replacement = "")
source(paste0(script_dir, "/BatchCorrection_function.R"))

count_file <- paste0(maindir, "/NGSmodule_analysis/Quantification/Quantification.", aligner, ".count.tab")

if (!file.exists(count_file)) {
  print(paste0("Can not find count file: ", count_file, "\nPlease run NGSmodule Quantification first!\n"))
  quit(status = 1)
}
if (!file.exists(SampleInfoFile)) {
  print(paste0("Can not find SampleInfoFile: ", SampleInfoFile, "\nPlease check the config parameters!\n"))
  quit(status = 1)
}

##### Load data #####
sample_info <- read.csv(SampleInfoFile, stringsAsFactors = F, header = T)
sample_info <- as.data.frame(sapply(sample_info, as.character))
sample_info <- sample_info[!is.na(sample_info[["Group"]]), ]
sample_info <- sample_info %>%
  group_by(SampleID) %>%
  mutate(RunID = paste0(RunID, collapse = ",")) %>%
  distinct() %>%
  as.data.frame()
rownames(sample_info) <- sample_info[, "SampleID"]
sample_info[, "Group"] <- factor(sample_info[, "Group"], levels = unique(sample_info[, "Group"]))
sample_info[, "BatchID"] <- factor(sample_info[, "BatchID"], levels = unique(sample_info[, "BatchID"]))

count_matrix_raw <- read.table(file = count_file, header = T, sep = "\t", row.names = 1, stringsAsFactors = F, check.names = F)
count_matrix <- count_matrix_raw[, which(str_detect(colnames(count_matrix_raw), pattern = paste0(".", aligner, ".count"))), drop = FALSE]
colnames(count_matrix) <- gsub(x = colnames(count_matrix), pattern = paste0(".", aligner, ".count"), replacement = "")
count_matrix <- count_matrix[, sample_info[["SampleID"]]]

logcpm <- cpm(count_matrix, log = TRUE, prior.count = 2)
logcpm_adj_raw <- logcpm

start <- which(str_detect(colnames(count_matrix_raw), ">>"))[1]
annotation_matrix <- count_matrix_raw[, start:ncol(count_matrix_raw)]

# Identify low-expressed genes --------------------------------------------
min_n <- max(min(table(sample_info[, "Group"])),2)
min_count <- 10
keep <- apply(count_matrix, 1, function(g) {
  sum(g >= min_count) >= min_n
})

# HVF calculation ---------------------------------------------------------
hvf <- data.frame(
  "gene" = rownames(count_matrix),
  "mean" = apply(count_matrix, 1, mean),
  "var" = apply(count_matrix, 1, var),
  stringsAsFactors = F
)
hvf <- hvf[hvf$var > 0, ]
hvf[, "log10_mean"] <- log10(hvf[, "mean"])
hvf[, "log10_var"] <- log10(hvf[, "var"])

fit <- loess(formula = log10_var ~ log10_mean, data = hvf, span = 0.3)
hvf[, "log10_var_exp"] <- fit$fitted
hvf[, "log10_var_residuals"] <- fit$residuals
hvf_gene <- hvf[keep,] %>% top_n(n = 1000,wt = log10_var_residuals) %>% pull("gene")

# qplot(x = log10_mean, y = log10_var, data = hvf) +
#   geom_line(aes(x = log10_mean, y = log10_var_exp), color = "red") +
#   geom_point(data = hvf[hvf_gene, ], aes(x = log10_mean, y = log10_var), color = "red") +
#   geom_abline(intercept = 0, slope = 1.5, color = "blue")

# select covariates for Combat and removeBatchEffect ----------------------
# covar <- c()
# for (col in c("Layout.PE.SE.")) {
#   col_batch <- unique(sample_info[, c("BatchID", col)])
#   if (min(table(col_batch[, col])) != 1 | length(table(col_batch[, col])) >= 3) {
#     covar <- c(covar, col)
#   }
# }
# if (length(covar) >= 1) {
#   covar_only_fmla <- as.formula(paste("~", paste(covar, collapse = "+")))
#   covar_only <- model.matrix(covar_only_fmla, data = sample_info)
#   covar_mod_fmla <- as.formula(paste("~", paste(covar, collapse = "+"), "+ Group"))
#   covar_mod <- model.matrix(covar_mod_fmla, data = sample_info)
# } else {
# covar_only <- NULL
# covar_mod <- model.matrix(~Group, data = sample_info)
# }

if (nlevels(sample_info[, "BatchID"]) >= 2) {
  # Construct design model -----------------------------------------------------
  mod1 <- model.matrix(~ 1 + as.factor(Group), data = sample_info)
  mod0 <- model.matrix(~1, data = sample_info)
  n.sv <- num.sv(dat = logcpm, mod = mod1, method = "be", seed = 11)

  # Apply raw ComBat funtion ------------------------------------------------
  cat("\n>>> Apply raw ComBat funtion\n")
  logcpm_adj_rawComBat <- ComBat(
    dat = logcpm,
    batch = as.factor(sample_info[["BatchID"]]),
    mod = mod1,
    prior.plots = FALSE,
    BPPARAM = bpparam("MulticoreParam")
  )

  # Apply ComBat-seq funtion ------------------------------------------------
  cat("\n>>> Apply ComBat-seq funtion\n")
  count_adj <- ComBat_seq_custom(
    counts = count_matrix,
    batch = as.factor(sample_info[["BatchID"]]),
    covar_mod = mod1
  )
  logcpm_adj_ComBatSeq <- cpm(count_adj, log = TRUE, prior.count = 2)

  # Apply SVA function ------------------------------------------------------
  cat("\n>>> Apply SVA funtion\n")
  svobj <- sva(dat = logcpm, mod = mod1, mod0 = mod0, n.sv = n.sv)
  logcpm_adj_SVA <- cleanY(y = logcpm, mod = mod1, svs = svobj$sv)

  # Apply removeBatchEffect function ----------------------------------------
  cat("\n>>> Apply removeBatchEffect funtion\n")
  logcpm_adj_removeBatchEffect <- removeBatchEffect(
    x = logcpm,
    batch = as.factor(sample_info[["BatchID"]]),
    design = mod1
  )
}

methods <- c("raw", "rawComBat", "ComBatSeq", "SVA", "removeBatchEffect")
pl <- list()
for (method in methods) {
  if (!exists(paste0("logcpm_adj_", method))) {
    next
  }
  cat("+++", method, "+++\n")
  logcpm_adj <- get(paste0("logcpm_adj_", method))
  write.table(
    x = cbind(data.frame(GeneID = rownames(logcpm_adj)), logcpm_adj),
    file = paste0(method, ".", aligner, ".log2CPM.tab"),
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
  )
  logcpm_adj <- logcpm_adj[hvf_gene, ]
  logcpm_adj_scale <- t(scale(t(logcpm_adj)))

  group_num <- length(unique(sample_info[, "Group"]))
  if (group_num <= 8) {
    group_color <- pal_nejm("default")(group_num)
  }
  if (group_num > 8 & group_num <= 20) {
    group_color <- pal_d3("category20")(group_num)
  }
  if (group_num > 20 & group_num <= 51) {
    group_color <- pal_igv("default")(group_num)
  }
  names(group_color) <- levels(sample_info[, "Group"])

  batch_num <- length(unique(sample_info[, "BatchID"]))
  if (batch_num <= 8) {
    batch_color <- pal_nejm("default")(batch_num)
  }
  if (batch_num > 8 & batch_num <= 20) {
    batch_color <- pal_d3("category20")(batch_num)
  }
  if (batch_num > 20 & batch_num <= 51) {
    batch_color <- pal_igv("default")(batch_num)
  }
  # batch_color <- colorRampPalette(brewer.pal(n = 11, name = "Spectral")[1:11])(length(unique(sample_info[["BatchID"]])))
  names(batch_color) <- levels(sample_info[["BatchID"]])

  plot_list <- list()
  ##### Dendrogram #####
  cat(">>> Hierarchical Clustering (colored by Groups)\n")
  hc <- hclust(dist(t(logcpm_adj_scale), method = "euclidean")) ## sqrt(sum((x_i - y_i)^2)). Large expressed genes -large vars.
  p1 <- ggtree(tr = hc) + layout_dendrogram()
  p2 <- ggplot(sample_info) +
    geom_tile(aes(x = SampleID, y = 1, fill = Group), color = ifelse(nrow(sample_info) <= 30, "black", "transparent")) +
    scale_fill_manual(values = group_color, name = "Group") +
    theme_void()
  p <- p2 %>% insert_top(p1, height = 8)
  plot_list[["Hierarchical_Clustering_colored_by_group"]] <- as.ggplot(aplotGrob(p))

  cat(">>> Hierarchical Clustering (colored by BatchID)\n")
  p2 <- ggplot(sample_info) +
    geom_tile(aes(x = SampleID, y = 1, fill = BatchID), color = ifelse(nrow(sample_info) <= 30, "black", "transparent")) +
    scale_fill_manual(values = batch_color, name = "Group") +
    theme_void()
  p <- p2 %>% insert_top(p1, height = 8)
  plot_list[["Hierarchical_Clustering_colored_by_batch"]] <- as.ggplot(aplotGrob(p))

  ##### PCA #####
  cat(">>> Principal Components Analysis\n")
  df_pca <- prcomp(t(logcpm_adj_scale), center = F, scale. = F)
  df_pca <- summary(df_pca)
  PoV <- round(df_pca$importance[2, ] * 100, 2)
  x <- df_pca$x[, "PC1"]
  y <- df_pca$x[, "PC2"]
  z <- df_pca$x[, "PC3"]

  df <- data.frame(
    x = x, y = y, z = z, sample = names(x),
    SampleID = sample_info[names(x), "SampleID"],
    Group = sample_info[names(x), "Group"],
    Batch = sample_info[names(x), "BatchID"],
    stringsAsFactors = FALSE
  )

  p <- ggplot(data = df, aes(x = x, y = y, fill = Group)) +
    geom_point(aes(SampleID = SampleID, Group = Group, Batch = Batch), shape = 21, alpha = 0.8, size = 3) +
    geom_rug(aes(color = Group), show.legend = FALSE) +
    # geom_mark_ellipse(aes(fill = Group, color = Group),show.legend = FALSE) +
    labs(title = "Principal Components Analysis", x = paste0("PC1(", PoV[1], "%)"), y = paste0("PC2(", PoV[2], "%)")) +
    scale_fill_manual(values = group_color) +
    scale_color_manual(values = group_color) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      panel.grid.major = element_line()
    )
  if (nrow(sample_info) < 20) {
    p <- p + geom_label_repel(
      aes(label = sample),
      size = 2.5, color = "white",
      min.segment.length = 0, segment.color = "black", segment.alpha = 0.8
    )
  }
  plot_list[["PCA2d_colored_by_group"]] <- p

  grob <- as.grob(~ {
    scatter3D(df[["x"]], df[["y"]], df[["z"]],
      type = "h",
      colvar = as.integer(df[["Group"]]),
      bg = alpha(group_color[df[["Group"]]], alpha = 0.8),
      col = group_color,
      colkey = list(at = c(1, 2, 3), plot = FALSE),
      xlab = paste0("PC1(", PoV[1], "%)"), ylab = paste0("PC2(", PoV[2], "%)"), zlab = paste0("PC3(", PoV[3], "%)"),
      pch = 21, cex = 1.2, cex.axis = 0.8,
      ticktype = "simple",
      bty = "b2", theta = 25, phi = 15, r = sqrt(10), d = 3,
      box = TRUE
    )
    # legend("topright", legend = names(group_color), pch = 16, col = group_color, cex = 0.9, bg = "white", inset = c(-0.18, 0))
  })
  p <- as.ggplot(grob, scale = 1.2) + theme(aspect.ratio = 1)
  plot_list[["PCA3d_colored_by_group"]] <- p


  p <- ggplot(data = df, aes(x = x, y = y, fill = Batch)) +
    geom_point(aes(SampleID = SampleID, Group = Group, Batch = Batch), shape = 21, alpha = 0.8, size = 3) +
    geom_rug(aes(color = Batch), show.legend = FALSE) +
    # geom_mark_ellipse(aes(fill = Batch, color = Batch), show.legend = FALSE) +
    labs(title = "Principal Components Analysis", x = paste0("PC1(", PoV[1], "%)"), y = paste0("PC2(", PoV[2], "%)")) +
    scale_fill_manual(values = batch_color) +
    scale_color_manual(values = batch_color) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      panel.grid.major = element_line()
    )
  if (nrow(sample_info) < 20) {
    p <- p + geom_label_repel(
      aes(label = sample),
      size = 2.5, color = "white",
      min.segment.length = 0, segment.color = "black", segment.alpha = 0.8
    )
  }
  plot_list[["PCA2d_colored_by_batch"]] <- p

  grob <- as.grob(~ {
    scatter3D(df[["x"]], df[["y"]], df[["z"]],
      type = "h",
      colvar = as.integer(df[["Batch"]]),
      bg = alpha(batch_color[df[["Batch"]]], alpha = 0.8),
      col = batch_color,
      colkey = list(at = c(1, 2, 3), plot = FALSE),
      xlab = paste0("PC1(", PoV[1], "%)"), ylab = paste0("PC2(", PoV[2], "%)"), zlab = paste0("PC3(", PoV[3], "%)"),
      pch = 21, cex = 1.2, cex.axis = 0.8,
      ticktype = "simple",
      bty = "b2", theta = 25, phi = 15, r = sqrt(10), d = 3,
      box = TRUE
    )
    # legend("topright", legend = names(batch_color), pch = 16, col = batch_color, cex = 0.9, bg = "white", inset = c(-0.18, 0))
  })
  p <- as.ggplot(grob, scale = 1.2) + theme(aspect.ratio = 1)
  plot_list[["PCA3d_colored_by_batch"]] <- p

  PC1_top <- sort(df_pca$rotation[, "PC1"], decreasing = T)[c(1:100, (nrow(df_pca$rotation) - 99):nrow(df_pca$rotation))]
  PC2_top <- sort(df_pca$rotation[, "PC2"], decreasing = T)[c(1:100, (nrow(df_pca$rotation) - 99):nrow(df_pca$rotation))]
  PC1_sd <- df_pca$sdev[1]
  PC2_sd <- df_pca$sdev[2]


  title <- ggdraw() +
    draw_label(
      label = paste("Batch-correction method:", method),
      fontface = "bold", x = 0, hjust = 0
    ) +
    theme(
      plot.margin = margin(0, 0, 0, 7)
    )
  res <- plot_grid(
    title, plot_grid(plotlist = plot_list, nrow = 2, byrow = FALSE),
    ncol = 1,
    rel_heights = c(0.05, 1)
  )
  pl[[method]] <- res
}


pdf("BatchCorrected.pdf", width = 14, height = 7)
invisible(lapply(pl, print))
invisible(dev.off())




##### check whether the unwanted file exists and remove it #####
if (file.exists("Rplots.pdf")) {
  invisible(file.remove("Rplots.pdf"))
}
