work_dir <- "/data/lab/LiLaiHua/scRNA-seq/Gonadal_ridge/analysis_1002/"
setwd(work_dir)

srt_list_Standard <- readRDS("srt_list_Standard.rds")
srt_list_SCT <- readRDS("srt_list_SCT.rds")
srt_list_fastMNN <- readRDS("srt_list_fastMNN.rds")
srt_list_Harmony <- readRDS("srt_list_Harmony.rds")

#### choose srt object #####

for (srt_name in c("srt_list_Standard","srt_list_SCT","srt_list_fastMNN","srt_list_Harmony")) {
  srt_use <- get(srt_name)[[1]]
  DefaultAssay(srt_use) <- "RNA"
  srt_use$batch <- srt_use$orig.ident
  project_name <- srt_use@project.name
  p <- summary_plot(
    srt = srt_use, return_list = F,features = "DDX4", color_by = "seurat_clusters", reduction = "umap", split_by = "batch", palette = "nejm",
    do_save = T, file_save = paste0(srt_name, ".summary.png")
  )
}


##### test features #####
if (!dir.exists("result/general")) {
  dir.create("result/general", recursive = T)
}
srt_use <- srt_list_Standard[[1]]
DefaultAssay(srt_use) <- "RNA"
srt_use$batch <- factor(srt_use$batch, levels = c("5d", "7d"))

p0 <- DimPlot(object = srt_use, group.by = "Phase", reduction = "umap", label = TRUE, combine = F) %>%
  custom_Dim_plot(srt = srt_use, group_by = "Phase", title = "All cells", combine = T, palette = "npg")
ggsave(plot = p0, filename = "result/general/Cellcycle_umap.png", width = 5, height = 5)
p123 <- lapply(SplitObject(srt_use, split.by = "batch"), function(srt) {
  title <- paste0(srt$batch %>% unique(), " (", ncol(srt), ")")
  p <- DimPlot(object = srt, group.by = "Phase", reduction = "umap", label = TRUE, combine = F)
  p <- custom_Dim_plot(plotlist = p, srt = srt, group_by = "Phase", title = title, combine = T, palette = "npg")
  return(p)
})
p123 <- cowplot::plot_grid(plotlist = p123, nrow = 1)
ggsave(plot = p123, filename = "result/general/Cellcycle_umap_byBatch.png", width = 15, height = 5)
p4 <- stat_barplot(srt = srt_use, x = "batch", color_by = "Phase", position = "fill", palette = "npg")
p5 <- stat_barplot(srt = srt_use, x = "seurat_clusters", color_by = "Phase", position = "fill", palette = "npg")
p <- ggarrange(p4, p5, ncol = 1)
ggsave(plot = p, filename = "result/general/Cellcycle_statbar_byBatch.png", width = 5, height = 10)

p0 <- FeaturePlot(
  object = srt_use, features = "DDX4", reduction = "umap",
  coord.fixed = TRUE, order = T, combine = F
) %>%
  custom_Feature_plot(srt = srt_use, ncol = 1, labels = "DDX4", gradient_use = c("gold", "red3"), na_value = "grey")
ggsave(plot = p0, filename = "result/general/DDX4_umap.png", width = 5, height = 5)
p1 <- FeaturePlot(
  object = srt_use, features = "DDX4", reduction = "umap", split.by = "batch",
  coord.fixed = TRUE, order = T, combine = F
) %>%
  custom_Feature_plot(srt = srt_use, features = "DDX4", right_label = "batch", ncol = 3, gradient_use = c("gold", "red3"), na_value = "grey")
ggsave(plot = p1, filename = "result/general/DDX4_umap_byBatch.png", width = 15, height = 5)
p2 <- FeaturePlot(
  object = srt_use, features = "DDX4", reduction = "umap", split.by = "ident",
  coord.fixed = TRUE, order = T, combine = F
) %>%
  custom_Feature_plot(srt = srt_use, features = "DDX4", right_label = "ident", ncol = 3, gradient_use = c("gold", "red3"), na_value = "grey")
ggsave(plot = p2, filename = "result/general/ZBTB16_umap_byCluster.png", width = 3 * 5, height = 4 * 5)
p3 <- VlnPlot(object = srt_use, features = c("DDX4"), pt.size = 0.2, combine = F) %>%
  custom_Vln_plot(alpha = 0.1, combine = T, ncol = 1)
ggsave(plot = p3, filename = "result/general/DDX4_vln.png", width = 1 * 5, height = 1 * 5)

p4 <- VlnPlot(object = srt_use, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2, combine = F) %>%
  custom_Vln_plot(alpha = 0.1, combine = T, ncol = 3)
ggsave(plot = p4, filename = "result/general/QC_vln_byCluster.png", width = 3 * 5, height = 1 * 5)

p5 <- VlnPlot(object = srt_use, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by = "batch", ncol = 3, pt.size = 0.2, combine = F) %>%
  custom_Vln_plot(alpha = 0.1, combine = T, ncol = 3)
ggsave(plot = p5, filename = "result/general/QC_vln_byBatch.png", width = 3 * 5, height = 1 * 5)
##### plot known cell marker #####
DefaultAssay(srt_use) <- "RNA"
marker <- read.table("/data/lab/HuangMingQian/scRNA-seq/ESC-PGC-GSCLC-new/analysis/cellmarker3.txt", sep = "\t", header = T, stringsAsFactors = FALSE, row.names = 1)

marker <- data.frame(cell="cell",marker = c("AR","AMH","SOX9","GATA4","GDNF","ZBTB16","UTF1","NANOS3","KIT","STRA8","SPO11","MLH1","RAD51","SYCP3","PRM1"))
rownames(marker) <- marker$marker
genelist <- rownames(marker)
# genelist<- rownames(marker)[which(marker$cell=="SSC")]
# genelist <- c("EPCAM","ITGA6")
# genelist <- c("POU5F1", "NANOG", "SOX2", "TFAP2C", "SOX17", "PRDM1", "TBXT", "EOMES", "GATA6", "GATA4", "CDX2", "KRT7", "GATA3")
genelist <- genelist[which(genelist %in% rownames(srt_use))]

p0 <- FeaturePlot(
  object = srt_use, features = genelist, reduction = "umap",
  coord.fixed = TRUE, order = T, combine = F
) %>%
  custom_Feature_plot(
    srt = srt_use, features = genelist, labels = marker[genelist, "cell"],
    gradient_use = c("gold", "red3"), na_value = "grey", combine = F
  )

p1 <- FeaturePlot(
  object = srt_use, features = genelist, reduction = "umap", split.by = "batch",
  coord.fixed = TRUE, order = T, combine = F
) %>%
  custom_Feature_plot(
    srt = srt_use, ncol = length(unique(srt_use$batch)), features = genelist, labels = marker[genelist, "cell"], right_label = "batch",
    gradient_use = c("gold", "red3"), na_value = "grey", combine = FALSE
  )
p <- plot_grid(plotlist = c(p0, p1), ncol = 3, byrow = FALSE)
ggsave(plot = p, filename = "result/general/ClassicalMarker_FGC_batch.png", width = 150 * 4, height = 120 * ceiling(length(genelist) / 3) * 2, units = "mm", limitsize = F)


##### qiaojie data
setwd("/data/database/HECA/FGC/analysis")
library(Seurat)
library(future)
library(intrinsicDimension)
library(dplyr)
library(stringr)
library(ggsci)
library(ggplot2)
library(cowplot)

plan(multiprocess, workers = 80, gc = TRUE) # stop with the command 'future:::ClusterRegistry("stop")'
plan()
set.seed(11)
source("/data/lab/LiLaiHua/scRNA-seq/Gonadal_ridge/analysis_zh/scRNA-SeuratWorkflow-function.R")
source("/home/zhanghao/Documents/pipeline/Single_cell/customize_Seurat_FeaturePlot.R")
FGC <- read.table(file = "/data/database/HECA/FGC/FGC_umi_counts.tsv", row.names = 1, sep = "\t", header = T, stringsAsFactors = F)
meta <- openxlsx::read.xlsx("/data/database/HECA/FGC/meta.xlsx", sheet = 3, rowNames = T)
FGC <- FGC[, rownames(meta)]

srt <- CreateSeuratObject(counts = FGC, project = "FGC", meta.data = meta)
srt[["batch"]] <- srt[["orig.ident"]] <- Idents(object = srt) <- factor("FGC")
srt[["Cluster"]] <- factor(srt[["Cluster", drop = TRUE]], levels = unique(srt[["Cluster", drop = TRUE]]))
srt[["CellType"]] <- factor(srt[["CellType", drop = TRUE]], levels = unique(srt[["CellType", drop = TRUE]]))
srt[["week"]] <- str_extract(string = colnames(srt), pattern = "(?<=_)\\d+(?=W_)") %>% as.numeric()
srt[["week"]] <- factor(srt[["week", drop = TRUE]], levels = sort(unique(srt[["week", drop = TRUE]])))
srt[["sex"]] <- factor(str_extract(string = colnames(srt), pattern = "^[MF](?=_)"), levels = c("M", "F"))
srt[["percent.mt"]] <- PercentageFeatureSet(object = srt, pattern = "^MT-")

srt <- Standard_SCP(
  sc = srt, nHVF = 4000, maxPC = 100, resolution = 1,
  cc_S_genes = Seurat::cc.genes.updated.2019$s.genes, cc_G2M_genes = Seurat::cc.genes.updated.2019$g2m.genes,
  exogenous_genes = "", assay = "RNA"
)
DefaultAssay(srt)
p <- DimPlot(
  object = srt, group.by = c("seurat_clusters", "Cluster", "CellType", "week", "sex"), cols = palette_uniform(group = 17, "npg"),
  label = T, reduction = "umap", repel = T, combine = FALSE
)
p <- lapply(p, function(x) x + theme(aspect.ratio = 1))
p <- plot_grid(plotlist = p, nrow = 2, align = "hv", axis = "tblr")
ggsave(plot = p, filename = "DimPlot.png", width = 150 * 4, height = 150 * 2, units = "mm", limitsize = F)

FeaturePlot(srt, features = c("DDX4")) %>% ggplotly()

marker1 <- read.table("/data/lab/HuangMingQian/scRNA-seq/ESC-PGC-GSCLC-new/analysis/cellmarker_FGC.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
marker2 <- read.table("/data/lab/HuangMingQian/scRNA-seq/ESC-PGC-GSCLC-new/analysis/cellmarker3.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
marker3 <- read.table("/data/lab/HuangMingQian/scRNA-seq/ESC-PGC-GSCLC-new/analysis/cellmarker_3layer.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
markers <- unique(rbind(marker1, marker2, marker3))
markers <- markers[markers$marker %in% rownames(srt), ]
p0 <- FeaturePlot(
  object = srt, features = markers$marker, reduction = "umap",
  coord.fixed = TRUE, order = T, combine = F
) %>%
  custom_Feature_plot(
    srt = srt, ncol = 7, features = markers$marker, labels = markers$cell,
    gradient_use = c("lightgrey", "red3"), na_value = "grey", combine = T
  )
ggsave(plot = p0, filename = "markers.png", width = 100 * 7, height = 100 * ceiling(nrow(markers) / 7), units = "mm", limitsize = F)



####################################################




topn <- selected_feature %>%
  group_by(cluster) %>%
  top_n(10, power)
p1 <- FeaturePlot(
  object = srt_use, features = topn$gene, reduction = "umap", split.by = "batch",
  coord.fixed = TRUE, order = F, combine = F
) %>%
  custom_Feature_plot(
    srt = srt_use, nrow = 6, features = topn$gene, labels = topn$cluster, right_label = "batch",
    gradient_use = c("gold", "red3"), na_value = "grey", combine = T
  )
ggsave(plot = p1, filename = "result/general/FGCMarker_batch.png", height = 90 * 6, width = 100 * ceiling(length(topn$gene) / 3) * 2.5, units = "mm", limitsize = F)



##### find DE genes #####
DefaultAssay(srt_use) <- "RNA"
Markers.lr <- FindAllMarkers(
  object = srt_use, only.pos = T, min.pct = 0.25, test.use = "LR", latent.vars = c("batch", "nCount_RNA"),
  logfc.threshold = 0.25
)
Markers.roc <- FindAllMarkers(
  object = srt_use, only.pos = T, min.pct = 0.25, test.use = "roc", ## default power>0.4
  logfc.threshold = 0.25
)
# Markers.lr.filter <- Markers.lr[which(Markers.lr$p_val_adj<0.05&Markers.lr$avg_logFC>=log(2)),]
# Markers.roc.filter <- Markers.roc[which(Markers.roc$power>=0.5),]
Markers.lr.roc.filter <- merge(
  x = Markers.lr, by.x = c("cluster", "gene"), all.x = F,
  y = Markers.roc[, c("myAUC", "power", "cluster", "gene")], by.y = c("cluster", "gene"), all.y = F
)
Markers.lr.roc.filter.pos <- Markers.lr.roc.filter[which(Markers.lr.roc.filter$avg_logFC > 0), ]
openxlsx::write.xlsx(Markers.lr.roc.filter.pos, file = "Markers.lr.roc.filter.pos.xlsx")
# saveRDS(Markers.lr,"Markers.lr.rds")
# saveRDS(Markers.roc,"Markers.roc.rds")
Markers.lr <- readRDS("Markers.lr.rds")
Markers.roc <- readRDS("Markers.roc.rds")

## plot top marker genes
filter <- table(Markers.lr.roc.filter.pos$gene)
Markers.lr.roc.filter.unique <- Markers.lr.roc.filter.pos[which(Markers.lr.roc.filter.pos$gene %in% names(filter[which(filter == 1)])), ]
top <- Markers.lr.roc.filter.unique %>%
  group_by(cluster) %>%
  top_n(4, avg_logFC)
top <- top[with(top, order(cluster)), ]

p0 <- DotPlot(
  object = srt_use, features = top$gene,
  col.min = -2, col.max = 2, dot.min = 0, dot.scale = 6,
  scale.by = "radius", scale.min = NA, scale.max = NA,
  x.lab.rot = FALSE
) +
  coord_flip() + scale_colour_distiller(palette = "OrRd", direction = 1)
ggsave(
  plot = p0, filename = "result/general/findMarkers_heatmap.png",
  units = "mm", width = 150, height = 110, scale = 1.5, limitsize = FALSE
)

p1 <- FeaturePlot(
  object = srt_use, features = top$gene, reduction = "umap",
  coord.fixed = TRUE, order = T, combine = F
) %>%
  custom_Feature_plot(
    srt = srt_use, ncol = 4, features = top$gene, labels = paste("Cluster:", top$cluster),
    gradient_use = c("gold", "red3"), na_value = "grey", combine = T
  )
ggsave(plot = p1, filename = "result/general/findMarkers_umap.png", width = 150 * 4, height = 100 * ceiling(length(top$gene) / 4), units = "mm", limitsize = F)

p2 <- VlnPlot(object = srt_use, features = top$gene, pt.size = 0.2, combine = F) %>%
  custom_Vln_plot(alpha = 0.1, label = paste("Cluster:", top$cluster), combine = T, ncol = 3)
ggsave(plot = p2, filename = "result/general/findMarkers_VlnPlot.png", width = 200 * 3, height = 100 * ceiling(length(top$gene) / 4), units = "mm", limitsize = F)

##### "EPCAM","ITGA6" analysis #####
srt_use[["DoubleMarker_GeomMean"]] <- FetchData(srt_use, vars = c("EPCAM", "ITGA6")) %>%
  apply(c(1, 2), expm1) %>%
  as.data.frame(stringsAsFactors = F) %>%
  apply(1, function(x) log1p(exp(mean(log(x)))))
# gghistogram(data = data.frame(x=srt_use[["DoubleMarker_GeomMean",drop=T]],stringsAsFactors = F),x="x")

p0 <- FeaturePlot(
  object = srt_use, features = c("EPCAM", "ITGA6", "DoubleMarker_GeomMean"), reduction = "umap",
  coord.fixed = TRUE, order = T, combine = F
) %>%
  custom_Feature_plot(
    srt = srt_use, ncol = 3, features = c("EPCAM", "ITGA6", "DoubleMarker_GeomMean"),
    gradient_use = c("gold", "red3"), na_value = "grey", combine = T
  )
ggsave(plot = p0, filename = "result/double_marker_test/DoubleMarker_umap.png", width = 3 * 5, height = 1 * 5)

p1 <- FeaturePlot(
  object = srt_use, features = c("EPCAM", "ITGA6", "DoubleMarker_GeomMean"), reduction = "umap", split.by = "batch",
  coord.fixed = TRUE, order = T, combine = F
) %>%
  custom_Feature_plot(
    srt = srt_use, ncol = 3, features = c("EPCAM", "ITGA6", "DoubleMarker_GeomMean"), right_label = "batch",
    gradient_use = c("gold", "red3"), na_value = "grey", combine = T
  )
ggsave(plot = p1, filename = "result/double_marker_test/DoubleMarker_umap_byBatch.png", width = 3 * 5, height = 3 * 5)

p2 <- FeaturePlot(
  object = srt_use, features = c("EPCAM", "ITGA6", "DoubleMarker_GeomMean"), reduction = "umap", split.by = "ident",
  coord.fixed = TRUE, order = T, combine = F
) %>%
  custom_Feature_plot(
    srt = srt_use, ncol = 3, features = c("EPCAM", "ITGA6", "DoubleMarker_GeomMean"), right_label = "ident",
    gradient_use = c("gold", "red3"), na_value = "grey", combine = T
  )
ggsave(plot = p2, filename = "result/double_marker_test/DoubleMarker_umap_byCluster.png", width = 3 * 5, height = 7 * 5, limitsize = F)

p3 <- VlnPlot(object = srt_use, features = c("EPCAM", "ITGA6", "DoubleMarker_GeomMean"), pt.size = 0.2, combine = F) %>%
  custom_Vln_plot(alpha = 0, combine = T, ncol = 3, AddBox = T, palette = "nejm")
ggsave(plot = p3, filename = "result/double_marker_test/DoubleMarker_vln.png", width = 3 * 5, height = 1 * 5)

p4 <- FeaturePlot(
  object = srt_use, features = c("EPCAM", "ITGA6", "DoubleMarker_GeomMean", "ZBTB16"), reduction = "umap", split.by = "batch",
  cells = WhichCells(srt_use, batch = "Pgc6d", expression = DoubleMarker_GeomMean > 0),
  coord.fixed = TRUE, order = T, combine = F
) %>%
  custom_Feature_plot(
    srt = srt_use, ncol = 4, features = c("EPCAM", "ITGA6", "DoubleMarker_GeomMean", "ZBTB16"), right_label = "batch",
    gradient_use = c("gold", "red3"), na_value = "grey", combine = T
  )
ggsave(plot = p4, filename = "result/double_marker_test/DoubleMarker_umap_byBatch_pos_ZBTB16.png", width = 4 * 5, height = 3 * 5, limitsize = F)

p5 <- FeaturePlot(
  object = srt_use, features = c("EPCAM", "ITGA6", "DoubleMarker_GeomMean", "ZBTB16"), reduction = "umap", split.by = "batch",
  cells = WhichCells(srt_use, batch = "Pgc6d", expression = DoubleMarker_GeomMean == 0),
  coord.fixed = TRUE, order = T, combine = F
) %>%
  custom_Feature_plot(
    srt = srt_use, ncol = 4, features = c("EPCAM", "ITGA6", "DoubleMarker_GeomMean", "ZBTB16"), right_label = "batch",
    gradient_use = c("gold", "red3"), na_value = "grey", combine = T
  )
ggsave(plot = p5, filename = "result/double_marker_test/DoubleMarker_umap_byBatch_neg_ZBTB16.png", width = 4 * 5, height = 3 * 5, limitsize = F)

srt_sub <- subset(srt_use, subset = DoubleMarker_GeomMean == 0)
srt_sub <- SetIdent(object = srt_sub, value = srt_sub$batch)
p6 <- VlnPlot(object = srt_sub, features = c("EPCAM", "ITGA6", "DoubleMarker_GeomMean", "ZBTB16"), pt.size = 0.2, combine = F) %>%
  custom_Vln_plot(alpha = 0, combine = T, ncol = 2, vln_alpha = 0.1)
ggsave(plot = p6, filename = "result/double_marker_test/DoubleMarker_vln_byBatch_neg_ZBTB16.png", width = 2 * 5, height = 2 * 4, limitsize = F)

srt_sub <- subset(srt_use, subset = DoubleMarker_GeomMean > 0)
srt_sub <- SetIdent(object = srt_sub, value = srt_sub$batch)
p7 <- VlnPlot(object = srt_sub, features = c("EPCAM", "ITGA6", "DoubleMarker_GeomMean", "ZBTB16"), pt.size = 0.2, combine = F) %>%
  custom_Vln_plot(alpha = 0, combine = T, ncol = 2, vln_alpha = 0.1)
ggsave(plot = p7, filename = "result/double_marker_test/DoubleMarker_vln_byBatch_pos_ZBTB16.png", width = 2 * 5, height = 2 * 4, limitsize = F)


meta <- c()
for (batchname in unique(srt_use[["batch", drop = T]])) {
  ZBTB16_pos <- WhichCells(srt_use, expression = ZBTB16 > 0 & batch == batchname)
  DM_pos <- WhichCells(srt_use, expression = DoubleMarker_GeomMean > 0 & batch == batchname)
  venn_list <- list("DM_pos" = DM_pos, "ZBTB16_pos" = ZBTB16_pos)
  overlap <- calculate.overlap(venn_list)
  # venn.plot <- venn.diagram(x= venn_list,filename = paste0("result/double_marker_test/DoubleMarker_ZBTB16_venn_",batchname,".png"),imagetype="png",height = 2000,width = 2000,
  #                           scaled=T,force.unique=T,print.mode="raw",
  #                           # hyper.test=T,total.population=length(WhichCells(srt_use,expression=batch==batchname)),
  #                           lwd=2,lty=c(1,2),fill= brewer.pal(9, "Set1")[1:2],
  #                           cex=1.5,fontface = "bold",fontfamily = "serif",
  #                           cat.pos=c(280,80),cat.dist=c(0.2,0.3),cat.cex=1.5,cat.fontface="bold",cat.fontfamily="serif",cat.default.pos="outer",
  #                           margin = 0.3,
  #                           ext.text = TRUE,ext.pos = c(0,-180),ext.dist = c(-0.05,0.1),ext.length = c(1,1))
  venn.output <- venn(data = venn_list, show.plot = FALSE)
  intersections <- attributes(venn.output)$intersections
  intersections.df <- stack(intersections)
  intersections.vec <- as.character(intersections.df$ind)
  names(intersections.vec) <- intersections.df$values
  meta <- c(meta, intersections.vec)
}
srt_use <- AddMetaData(object = srt_use, metadata = meta, col.name = "DM_ZBTB16_intersections")
srt_use$DM_ZBTB16_intersections[which(is.na(srt_use$DM_ZBTB16_intersections))] <- "Others"
srt_use$DM_ZBTB16_intersections <- factor(x = srt_use$DM_ZBTB16_intersections, levels = c("DM_pos", "ZBTB16_pos", "DM_pos:ZBTB16_pos", "Others"))
p <- VlnPlot(
  object = srt_use, features = c("DoubleMarker_GeomMean", "ZBTB16"), pt.size = 0.1, combine = F,
  group.by = "DM_ZBTB16_intersections", split.by = "batch"
) %>%
  custom_Vln_plot(alpha = 0.1, combine = T, ncol = 2, vln_alpha = 0.1)
ggsave(plot = p, filename = "result/double_marker_test/DoubleMarker_vln_intersection.png", width = 2 * 5, height = 1 * 5)

srt_use <- SetIdent(object = srt_use, value = srt_use$batch)
p <- VlnPlot(object = srt_use, features = c("DoubleMarker_GeomMean", "ZBTB16", "EPCAM", "ITGA6"), pt.size = 0.2, combine = F) %>%
  custom_Vln_plot(alpha = 0, combine = T, ncol = 2, vln_alpha = 0.1, AddBox = T, palette = "nejm")
srt_use <- SetIdent(object = srt_use, value = srt_use$seurat_clusters)
ggsave(plot = p, filename = "result/double_marker_test/DoubleMarker_vln_byBatch.png", width = 2 * 5, height = 2 * 5)

p <- stat_barplot(srt_use, color_by = "DM_ZBTB16_intersections", position = "fill", palette = "nejm")
ggsave(plot = p, filename = "result/double_marker_test/DoubleMarker_bar_intersection.png", width = 8, height = 6)

saveRDS(object = srt_use, file = paste0(srt_integrated@project.name, ".regress_umi.rds"))

##### scmap #####
ref_sce <- SingleCellExperiment(
  assays = list(normcounts = as.matrix(QiaoJie_FGC)),
  rowData = data.frame(feature_symbol = rownames(QiaoJie_FGC), stringsAsFactors = F),
  colData = QiaoJie_FGC_meta[colnames(QiaoJie_FGC), "Cell_type", drop = FALSE]
)
logcounts(ref_sce) <- log1p(normcounts(ref_sce))
ref_sce <- selectFeatures(ref_sce, n_features = 500, suppress_plot = F)
query_sce <- SingleCellExperiment(
  assays = list(logcounts = as.matrix(srt_use@assays$integrated@data)),
  rowData = data.frame(feature_symbol = rownames(srt_use), stringsAsFactors = F),
  colData = data.frame(
    Cell_type = colnames(srt_use) %>% gsub(pattern = "_.*", replacement = "", perl = T),
    row.names = colnames(srt_use), stringsAsFactors = F
  )
)

### scmap-cluster ###
ref_sce <- indexCluster(ref_sce, cluster_col = "Cell_type")
# heatmap(as.matrix(metadata(ref_sce)$scmap_cluster_index))
nrow(metadata(ref_sce)$scmap_cluster_index)
scmapCluster_results <- scmapCluster(projection = query_sce, index_list = list(mFGC = metadata(ref_sce)$scmap_cluster_index), threshold = 0.65)
table(scmapCluster_results$scmap_cluster_labs)

srt_use[["FGC_assigned"]] <- scmapCluster_results$scmap_cluster_labs
srt_use[["FGC_assigned_sim"]] <- scmapCluster_results$scmap_cluster_siml
p <- stat_barplot(srt_use, x = "seurat_clusters", color_by = "FGC_assigned", position = "fill", palette = "npg")
ggsave(plot = p, filename = "result/FGC_compare/FGC_assigned_statbar_byCluster.png", width = 7, height = 5)
p <- stat_barplot(srt_use, x = "batch", color_by = "FGC_assigned", position = "fill", palette = "npg")
ggsave(plot = p, filename = "result/FGC_compare/FGC_assigned_statbar_byBatch.png", width = 7, height = 5)

p0 <- DimPlot(object = srt_use, group.by = "FGC_assigned", reduction = "umap", label = T, combine = F) %>%
  custom_Dim_plot(srt = srt_use, group_by = "FGC_assigned", title = "All cells", combine = T, palette = "npg")
ggsave(plot = p0, filename = "result/FGC_compare/FGC_assigned_umap.png", width = 18, height = 5)
p123 <- lapply(SplitObject(srt_use, split.by = "batch"), function(srt) {
  title <- paste0(srt$batch %>% unique(), " (", ncol(srt), ")")
  p <- DimPlot(object = srt, group.by = "FGC_assigned", reduction = "umap", label = F, combine = F)
  p <- custom_Dim_plot(plotlist = p, srt = srt, group_by = "FGC_assigned", title = title, combine = T, palette = "npg")
  return(p)
})
p123 <- cowplot::plot_grid(plotlist = p123, nrow = 1)
ggsave(plot = p123, filename = "result/FGC_compare/FGC_assigned_umap_byBatch.png", width = 18, height = 5)

index <- metadata(ref_sce)$scmap_cluster_index
ref.matrix <- metadata(ref_sce)$scmap_cluster_index
tmp <- setFeatures(query_sce, rownames(index))
query.matrix <- tmp[rowData(tmp)$scmap_features, ] %>%
  logcounts() %>%
  as.matrix()
query.matrix <- query.matrix[order(rownames(query.matrix)), ]

res_cos <- proxy::simil(t(index), t(query.matrix), method = "cosine")
res_cos <- matrix(res_cos, ncol = nrow(t(index)), byrow = TRUE)
max_inds1 <- max.col(res_cos)
maxs1 <- rowMaxs(res_cos)

res_pearson <- cor(index, query.matrix, method = "pearson") %>% t()
max_inds2 <- max.col(res_pearson)
maxs2 <- rowMaxs(res_pearson)

res_spearman <- cor(index, query.matrix, method = "spearman") %>% t()
max_inds3 <- max.col(res_spearman)
maxs3 <- rowMaxs(res_spearman)

srt_use[["mFGC#1_cosine_similarity"]] <- res_cos[, 1]
srt_use[["mFGC#2_cosine_similarity"]] <- res_cos[, 2]
srt_use[["mFGC#3_cosine_similarity"]] <- res_cos[, 3]
srt_use[["mSoma#1_cosine_similarity"]] <- res_cos[, 4]
srt_use[["mSoma#2_cosine_similarity"]] <- res_cos[, 5]
srt_use[["mSoma#3_cosine_similarity"]] <- res_cos[, 6]
srt_use[["mSoma#4_cosine_similarity"]] <- res_cos[, 7]
p <- FeaturePlot(
  object = srt_use, features = c("mFGC#1_cosine_similarity", "mFGC#2_cosine_similarity", "mFGC#3_cosine_similarity"), reduction = "umap", # split.by = "batch",
  coord.fixed = TRUE, order = T, combine = F
) %>%
  custom_Feature_plot(
    srt = srt_use, ncol = 3, features = c("mFGC#1_cosine", "mFGC#2_cosine", "mFGC#3_cosine"), # right_label = "batch",
    gradient_use = c("black", "red3"), min_exp_cutoff = 0, max_exp = 1, na_value = "grey", combine = T
  )

for (cell_type in c("mFGC#1", "mFGC#2", "mFGC#3", "mSoma#1", "mSoma#2", "mSoma#3", "mSoma#4")) {
  p <- FeaturePlot(
    object = srt_use, features = paste0(cell_type, "_cosine_similarity"), reduction = "umap", # split.by = "batch",
    coord.fixed = TRUE, order = T, combine = F
  ) %>%
    custom_Feature_plot(
      srt = srt_use, ncol = 1, features = paste0(cell_type, "_cosine_similarity"), # right_label = "batch",
      gradient_use = c("black", "grey70", "red3"), gradient_values = c(0, 0.65 / 0.75, 1), min_exp_cutoff = 0, max_exp = 0.75,
      na_value = "black", legend_title = "Cosine similarity", combine = T
    )
  ggsave(plot = p, filename = paste0("result/FGC_compare/", cell_type, "_similarity_dim.png"), width = 6, height = 5)
  p <- VlnPlot(object = srt_use, features = paste0(cell_type, "_cosine_similarity"), pt.size = 0.1, combine = F) %>%
    custom_Vln_plot(alpha = 0, combine = T, ncol = 1, vln_alpha = 0.1, y_cutoff = 0.65, pass = "h", AddBox = T, palette = "nejm")
  ggsave(plot = p, filename = paste0("result/FGC_compare/", cell_type, "_similarity_vln.png"), width = 12, height = 7)
}

## similarity heatmap for each cluster
c1 <- aggregate(x = srt_use$`mFGC#1_cosine_similarity`, by = list(srt_use$seurat_clusters), mean)
c2 <- aggregate(x = srt_use$`mFGC#2_cosine_similarity`, by = list(srt_use$seurat_clusters), mean)
c3 <- aggregate(x = srt_use$`mFGC#3_cosine_similarity`, by = list(srt_use$seurat_clusters), mean)
c4 <- aggregate(x = srt_use$`mSoma#1_cosine_similarity`, by = list(srt_use$seurat_clusters), mean)
c5 <- aggregate(x = srt_use$`mSoma#2_cosine_similarity`, by = list(srt_use$seurat_clusters), mean)
c6 <- aggregate(x = srt_use$`mSoma#3_cosine_similarity`, by = list(srt_use$seurat_clusters), mean)
c7 <- aggregate(x = srt_use$`mSoma#4_cosine_similarity`, by = list(srt_use$seurat_clusters), mean)
sim_matrix <- Reduce(function(x, y) merge(x, y, by = 1, all = TRUE), list(c1, c2, c3, c4, c5, c6, c7))
colnames(sim_matrix) <- c("Cluster", "mFGC#1", "mFGC#2", "mFGC#3", "mSoma#1", "mSoma#2", "mSoma#3", "mSoma#4")
rownames(sim_matrix) <- sim_matrix[, "Cluster"]
sim_matrix[, "Cluster"] <- NULL
Heatmap(sim_matrix, cluster_columns = F)

# saveRDS(object = srt_use,file = paste0(srt_integrated@project.name,".regress_umi.rds"))


##### singleR ####
singler <- CreateSinglerObject(
  counts = srt_use@assays$RNA@data, annot = NULL, project.name = "Pgc6d-PgcE-PgcL", min.genes = 0,
  technology = "10X", species = "Human", citation = "",
  ref.list = list(), normalize.gene.length = F, variable.genes = "de",
  fine.tune = T, do.signatures = T, clusters = srt_use@active.ident, do.main.types = T,
  reduce.file.size = T, numCores = 30
)
singler$seurat <- srt_use
singler$meta.data$orig.ident <- srt_use@meta.data$seurat_clusters
singler$meta.data$xy <- srt_use@reductions$tsne@cell.embeddings
singler$meta.data$clusters <- srt_use@active.ident
singler.new <- convertSingleR2Browser(singler, use.singler.cluster.annot = F)
saveRDS(singler.new, file = paste0(singler.new@project.name, ".singler.rds"))

names(singler$singler) <- c("a1", "a2")
png(filename = "test.png", width = 2000, height = 3000)
corrplot(singler$singler$a1$SingleR.clusters.main$scores, is.corr = FALSE, method = "square")
dev.off()


### velocyto ###
Pgc6d_loom <- ReadVelocity(file = "../PGC6d/velocyto/PGC6d.loom")
PgcE_loom <- ReadVelocity(file = "../E1/velocyto/E1.loom")
PgcL_loom <- ReadVelocity(file = "../F1/velocyto/F1.loom") %>% as.Seurat()
DefaultAssay(PgcE_loom)

colSums(Pgc6d_loom$spliced)
