#!/usr/bin/env Rscript
suppressWarnings(suppressPackageStartupMessages(invisible(lapply(
  c(
    "limma", "edgeR", "data.table", "gplots", "stringr", "ComplexHeatmap",
    "ggsci", "ggpubr", "RColorBrewer", "circlize", "ggrepel", "GGally",
    "factoextra", "nord", "aplot", "ggtree"
  ),
  require,
  character.only = TRUE
))))

args <- commandArgs(trailingOnly = TRUE)
maindir <- args[1]
aligner <- args[2]
SampleInfoFile <- args[3]

######### example 1 ############
# maindir <- "/data/lab/ZhangHao/tmp/"
# aligner <- "hisat2"
# SampleInfoFile <- "/data/lab/ZhangHao/tmp/temp_20200413205117.Sample_info.csv"
################################

######### example 2 ############
# setwd("/data/lab/HeXi/Rpl39l/RNC-seq/NGSmodule_analysis/Quantification/postQuantificationQC/")
# maindir <- "/data/lab/HeXi/Rpl39l/RNC-seq/"
# aligner <- "hisat2"
# SampleInfoFile <- "/data/lab/HeXi/Rpl39l/RNC-seq/Sample_info.csv"
################################

######### example 3 ############
# setwd("/data/database/SRR_collection/human/early_embyro/NGSmodule_analysis/Quantification/postQuantificationQC/")
# maindir <- "/data/database/SRR_collection/human/early_embyro"
# aligner <- "hisat2"
# SampleInfoFile <- "/data/database/SRR_collection/human/early_embyro/temp_20200714173936.Sample_info.csv"
################################

##### source the function #####
script_dir <- gsub(x = script_path,pattern = "postQuantificationQC.R",replacement = "")
source(paste0(script_dir,"/postQuantificationQC_function.R"))

count_file <- paste0(maindir, "/NGSmodule_analysis/Quantification/Quantification.", aligner, ".count.tab")

if (!file.exists(count_file)) {
  print(paste0("Can not find count file: ", count_file, "\nPlease run NGSmodule Quantification first!\n"))
  quit(status = 1)
}
if (!file.exists(SampleInfoFile)) {
  print(paste0("Can not find SampleInfoFile: ", SampleInfoFile, "\nPlease check the config parameters!\n"))
  quit(status = 1)
}

sample_info <- read.csv(SampleInfoFile, stringsAsFactors = F, header = T)
sample_info <- unique(sample_info[, colnames(sample_info) != "RunID"])
rownames(sample_info) <- sample_info[, "SampleID"]

count_matrix_raw <- read.table(file = count_file, header = T, sep = "\t", row.names = 1, stringsAsFactors = F, check.names = F)
count_matrix <- count_matrix_raw[, which(str_detect(colnames(count_matrix_raw), pattern = paste0(".", aligner, ".count"))), drop = FALSE]
colnames(count_matrix) <- gsub(x = colnames(count_matrix), pattern = paste0(".", aligner, ".count"), replacement = "")
count_matrix <- count_matrix[, sample_info[["SampleID"]]]
count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]

logcpm <- cpm(count_matrix, log = TRUE, prior.count = 2)
logcpm_scale <- t(scale(t(logcpm)))

anno_start <- which(str_detect(colnames(count_matrix_raw), ">>"))[1]
anno_matrix <- count_matrix_raw[, anno_start:ncol(count_matrix_raw)]

QCpath <- getwd()

##### main NGSmodule QC #####
setwd(QCpath)
group_num <- length(unique(sample_info[, "Group"]))
if (group_num <= 5) {
  col_color <- nord("mountain_forms")[1:group_num]
}
if (group_num > 5 & group_num <= 8) {
  col_color <- pal_nejm("default")(group_num)
}
if (group_num > 8 & group_num <= 20) {
  col_color <- pal_d3("category20")(group_num)
}
if (group_num > 20 & group_num <= 51) {
  col_color <- pal_igv("default")(group_num)
}
names(col_color) <- unique(sample_info[, "Group"])
batch_color <- colorRampPalette(brewer.pal(n = 11, name = "Spectral")[1:11])(length(unique(sample_info[["BatchID"]])))
names(batch_color) <- unique(sample_info[["BatchID"]])

plot_list <- list()

##### Dendrogram #####
hc <- hclust(dist(t(logcpm_scale), method = "euclidean"))

cat(">>> Hierarchical Clustering (colored by Groups)\n")
p1 <- ggtree(tr = hc) + layout_dendrogram()
p2 <- ggplot(sample_info) +
  geom_tile(aes(x = SampleID, y = 1, fill = Group), color = ifelse(nrow(sample_info) <= 30, "black", "transparent")) +
  scale_fill_manual(values = col_color, name = "Group") +
  theme_void()
p <- p2 %>% insert_top(p1, height = 8)
plot_list[["Hierarchical_Clustering_colored_by_group"]] <- p


cat(">>> Hierarchical Clustering (colored by BatchID)\n")
p2 <- ggplot(sample_info) +
  geom_tile(aes(x = SampleID, y = 1, fill = BatchID), color = ifelse(nrow(sample_info) <= 30, "black", "transparent")) +
  scale_fill_manual(values = batch_color, name = "Group") +
  theme_void()
p <- p2 %>% insert_top(p1, height = 8)
plot_list[["newpage"]] <- ggplot()
plot_list[["Hierarchical_Clustering_colored_by_batch"]] <- p

##### correlation heatmap #####
cat(">>> Correlation Heatmap\n")
logcpm_cor <- round(cor(logcpm, method = "spearman"), 2)
q1 <- quantile(logcpm_cor, 0.99)
q2 <- quantile(logcpm_cor, 0.01)
color_palette <- colorRampPalette(brewer.pal(11, "RdBu"))(100)

df_bottom_annotation <- HeatmapAnnotation(
  "Groups" = anno_simple(
    x = sample_info[colnames(logcpm), "Group"],
    col = col_color[sample_info[colnames(logcpm), "Group"]],
    gp = gpar(col = ifelse(nrow(sample_info) <= 30, "black", "transparent"))
  ),
  "Batches" = anno_simple(
    x = sample_info[colnames(logcpm), "BatchID"],
    col = batch_color,
    gp = gpar(col = ifelse(nrow(sample_info) <= 30, "black", "transparent"))
  ),
  border = TRUE,
  show_legend = FALSE,
  show_annotation_name = TRUE,
  annotation_name_side = "left"
)
df_right_annotation <- HeatmapAnnotation(
  "Groups" = anno_simple(
    x = sample_info[colnames(logcpm), "Group"],
    col = col_color[sample_info[colnames(logcpm), "Group"]],
    gp = gpar(col = ifelse(nrow(sample_info) <= 30, "black", "transparent"))
  ),
  "Batches" = anno_simple(
    x = sample_info[colnames(logcpm), "BatchID"],
    col = batch_color,
    gp = gpar(col = ifelse(nrow(sample_info) <= 30, "black", "transparent"))
  ),
  border = TRUE,
  show_legend = FALSE,
  show_annotation_name = TRUE,
  annotation_name_side = "top",
  which = "row"
)

r_max <- 0.408 / ncol(logcpm_cor)
if (nrow(sample_info) <= 30) {
  ht <- Heatmap(logcpm_cor,
    col = colorRamp2(seq(q1, q2, length = 100), color_palette),
    rect_gp = gpar(type = "none"),
    cell_fun = function(j, i, x, y, w, h, fill) {
      p <- logcpm_cor[i, j]
      perc <- (p - min(logcpm_cor)) / (max(logcpm_cor) - min(logcpm_cor))
      grid.circle(x, y,
        r = r_max / 2 * (1 + perc),
        gp = gpar(col = "black", lwd = 1, fill = fill)
      )
      grid.rect(x, y,
        width = w, height = h,
        gp = gpar(col = "grey", lwd = 1, fill = "transparent")
      )
    },
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 10),
    row_names_gp = gpar(fontsize = 10),
    border = T,
    bottom_annotation = df_bottom_annotation,
    right_annotation = df_right_annotation,
    heatmap_legend_param = list(
      title = "Spearman's correlation coefficient",
      title_gp = gpar(fontsize = 12, fontfamily = "sans"),
      title_position = "lefttop",
      grid_height = unit(4, "mm"),
      grid_width = unit(4, "mm"),
      border = "black",
      labels_gp = gpar(fontsize = 8),
      legend_direction = "horizontal",
      legend_width = unit(4, "cm")
    )
  )
} else {
  ht <- Heatmap(logcpm_cor,
    col = colorRamp2(seq(q1, q2, length = 100), color_palette),
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_names_gp = gpar(fontsize = 10),
    row_names_gp = gpar(fontsize = 10),
    border = T,
    bottom_annotation = df_bottom_annotation,
    right_annotation = df_right_annotation,
    heatmap_legend_param = list(
      title = "Spearman's correlation coefficient",
      title_gp = gpar(fontsize = 12, fontfamily = "sans"),
      title_position = "lefttop",
      grid_height = unit(4, "mm"),
      grid_width = unit(4, "mm"),
      border = "black",
      labels_gp = gpar(fontsize = 8),
      legend_direction = "horizontal",
      legend_width = unit(4, "cm")
    )
  )
}
grob <- grid.grabExpr(ComplexHeatmap::draw(ht, heatmap_legend_side = "top"))
p <- as_ggplot(grob) + theme(aspect.ratio = 1)
plot_list[["Correlation_Heatmap"]] <- p

##### correlation paired scatter  #####
if (nrow(sample_info) <= 6) {
  cat(">>> Correlation Paired Scatter\n")
  cor_uniq <- sort(unique(c(round(cor(logcpm, method = "spearman"), 2))))
  color_cor <- colorRampPalette(brewer.pal(n = 9, name = "YlOrBr")[2:9])(length(cor_uniq))
  names(color_cor) <- cor_uniq

  p <- ggpairs(as.data.frame(logcpm),
    title = "Multivariate correlation analysis",
    xlab = paste0("log2CPM"),
    ylab = paste0("log2CPM"),
    upper = list(continuous = custom_cor),
    lower = list(continuous = custom_point),
    diag = list(continuous = custom_density)
  ) +
    theme(
      aspect.ratio = 1,
      strip.background = element_rect(fill = "transparent", color = "black"),
      strip.text = element_text(size = 12)
    )
  plot_list[["Correlation_PairedScatter"]] <- p
}

##### PCA #####
cat(">>> Principal Components Analysis\n")
df_pca <- prcomp(t(logcpm_adj_scale), center = F, scale. = F)
df_pca <- summary(df_pca)
PoV <- round(df_pca$importance[2, ] * 100, 2)
x <- df_pca$x[, "PC1"]
y <- df_pca$x[, "PC2"]
m <- max(c(abs(x), abs(y)))

df <- data.frame(
  x = x, y = y, sample = names(x),
  Group = sample_info[names(x), "Group"],
  Batch = sample_info[names(x), "BatchID"],
  stringsAsFactors = FALSE
)

p <- ggplot(data = df, aes(x = x, y = y, Group = Group, fill = Group, label = sample)) +
  geom_point(shape = 21, alpha = 0.8) +
  geom_rug(aes(color = Group), show.legend = FALSE) +
  labs(title = "Principal Components Analysis", x = paste0("PC1(", PoV[1], "%)"), y = paste0("PC2(", PoV[2], "%)")) +
  scale_fill_manual(values = col_color) +
  scale_color_manual(values = col_color) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  xlim(-m, m) +
  ylim(-m, m) +
  theme_pubr(border = T, legend = "right") +
  theme(
    aspect.ratio = 1,
    panel.grid.major = element_line()
  )
if (nrow(sample_info) <= 30) {
  p <- p + geom_label_repel(
    size = 2.5, color = "white",
    min.segment.length = 0, segment.color = "black", segment.alpha = 0.8
  )
}
plot_list[["PCA_colored_by_group"]] <- p

p <- ggplot(data = df, aes(x = x, y = y, Batch = Batch, fill = Batch, label = sample)) +
  geom_point(shape = 21, alpha = 0.8) +
  geom_rug(aes(color = Batch), show.legend = FALSE) +
  labs(title = "Principal Components Analysis", x = paste0("PC1(", PoV[1], "%)"), y = paste0("PC2(", PoV[2], "%)")) +
  scale_fill_manual(values = col_color) +
  scale_color_manual(values = col_color) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  xlim(-m, m) +
  ylim(-m, m) +
  theme_pubr(border = T, legend = "right") +
  theme(
    aspect.ratio = 1,
    panel.grid.major = element_line()
  )
if (nrow(sample_info) <= 30) {
  p <- p + geom_label_repel(
    size = 2.5, color = "white",
    min.segment.length = 0, segment.color = "black", segment.alpha = 0.8
  )
}
plot_list[["PCA_colored_by_batch"]] <- p



##### PCA top features #####
cat(">>> PCA top features\n")
PC1_top <- sort(df_pca$rotation[, "PC1"], decreasing = T)[c(1:100, (nrow(df_pca$rotation) - 99):nrow(df_pca$rotation))]
PC2_top <- sort(df_pca$rotation[, "PC2"], decreasing = T)[c(1:100, (nrow(df_pca$rotation) - 99):nrow(df_pca$rotation))]

mat_PC1 <- logcpm_scale[names(PC1_top), ]
mat_PC2 <- logcpm_scale[names(PC2_top), ]

color_palette <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(100))

df_bottom_annotation <- HeatmapAnnotation(
  foo1 = anno_simple(
    x = sample_info[colnames(logcpm), "Group"],
    col = col_color[sample_info[colnames(logcpm), "Group"]],
    gp = gpar(col = ifelse(nrow(sample_info) <= 30, "black", "transparent"))
  ),
  border = TRUE,
  show_legend = F,
  show_annotation_name = F
)

ht_legend <- list(
  title = paste0("Z-score"),
  title_gp = gpar(fontsize = 10, fontfamily = "sans"),
  title_position = "topleft",
  grid_height = unit(3, "mm"),
  grid_width = unit(3, "mm"),
  border = "black",
  labels_gp = gpar(fontsize = 10),
  legend_direction = "horizontal",
  legend_width = unit(3, "cm")
)

ht1 <- Heatmap(mat_PC1,
  column_title = "PC1 top features",
  col = colorRamp2(seq(-max(abs(mat_PC1)), max(abs(mat_PC1)), length = 100), color_palette),
  show_row_names = FALSE,
  show_column_names = ifelse(nrow(sample_info) >= 30, FALSE, TRUE),
  cluster_rows = F,
  column_title_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 8),
  border = T,
  bottom_annotation = df_bottom_annotation,
  heatmap_legend_param = ht_legend
)

ht2 <- Heatmap(mat_PC2,
  column_title = "PC2 top features",
  col = colorRamp2(seq(-max(abs(mat_PC2)), max(abs(mat_PC2)), length = 100), color_palette),
  show_row_names = FALSE,
  show_column_names = ifelse(nrow(sample_info) >= 30, FALSE, TRUE),
  cluster_rows = F,
  column_title_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 8),
  border = T,
  bottom_annotation = df_bottom_annotation,
  heatmap_legend_param = ht_legend
)


grob <- grid.grabExpr(ComplexHeatmap::draw(ht1 + ht2, heatmap_legend_side = "top"))
p <- as_ggplot(grob)

plot_list[["PCA_Heatmap"]] <- p


##### output report #####
pdf("./postQuantificationQC.pdf", width = 8, height = 5)
invisible(lapply(plot_list, print))
invisible(dev.off())


##### check whether the unwanted file exists and remove it #####
if (file.exists("Rplots.pdf")) {
  invisible(file.remove("Rplots.pdf"))
}
