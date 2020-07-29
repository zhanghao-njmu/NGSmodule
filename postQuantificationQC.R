#!/usr/bin/env Rscript
suppressWarnings(suppressPackageStartupMessages(invisible(lapply(
  c("limma","edgeR","data.table","gplots","stringr","ComplexHeatmap",
    "ggsci","ggpubr","RColorBrewer","circlize","ggrepel","GGally",
    "factoextra","nord"),
  require, character.only = TRUE))))

args <- commandArgs( trailingOnly = TRUE )
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


count_file <- paste0(maindir,"/NGSmodule_analysis/Quantification/Quantification.",aligner,".count.tab")

if (!file.exists(count_file)) {
  print(paste0("Can not find count file: ",count_file,"\nPlease run NGSmodule Quantification first!\n"))
  quit(status=1)
}
if (!file.exists(SampleInfoFile)) {
  print(paste0("Can not find SampleInfoFile: ",SampleInfoFile,"\nPlease check the config parameters!\n"))
  quit(status=1)
}

sample_info<- read.csv(SampleInfoFile,stringsAsFactors = F,header = T)
rownames(sample_info) <- sample_info[,2]
sample_info <- sample_info[order(sample_info[,3]),]

count_matrix_raw <- read.table(file = count_file, header = T,sep = "\t",row.names = 1,stringsAsFactors = F,check.names=F)
count_matrix <- count_matrix_raw[,which(str_detect(colnames(count_matrix_raw),pattern = paste0(".",aligner,".count"))),drop=FALSE]
colnames(count_matrix) <- gsub(x=colnames(count_matrix),pattern = paste0(".",aligner,".count"),replacement = "")
count_matrix <- count_matrix[,sample_info[,2]]
count_matrix <- count_matrix[rowSums(count_matrix)>0,]

logcpm_file <- paste0(maindir,"/NGSmodule_analysis/Quantification/Quantification.",aligner,".log2CPM.tab")
if (!file.exists(logcpm_file)) {
  cat(paste0("Can not find logcpm_file: ",logcpm_file,"\nWill do transformation on the count_matrix into log2CPM!\n"))
  logcpm <- cpm(count_matrix,log = TRUE,prior.count = 2)
}else{
  logcpm <- read.table(file = logcpm_file, header = T,sep = "\t",row.names = 1,stringsAsFactors = F,check.names=F)
}
logcpm_scale <- t(scale(t(logcpm)))

start<- which(str_detect(colnames(count_matrix_raw),">>"))[1]
annotation_matrix <- count_matrix_raw[,start:ncol(count_matrix_raw)]

QCpath<- getwd()

##### limma and edgeR: MDS plot and heatmap ##### 
wd_path <- paste0(QCpath,"/Other")
dir.create(path = wd_path, showWarnings = FALSE)
setwd(wd_path)

df <- DGEList(counts=count_matrix)
df <- calcNormFactors(df)
if (ncol(df)>=3) {
  pdf('edgeR_MDS_plot.pdf')
  MDSdata <- plotMDS(df)
  invisible(dev.off())
  write.table(MDSdata$distance.matrix, 'edgeR_MDS_distance_matrix.txt', quote=FALSE,append = FALSE)
  MDSxy = MDSdata$cmdscale.out
  colnames(MDSxy) = c(paste(MDSdata$axislabel, '1'), paste(MDSdata$axislabel, '2'))
  line <- "#plot_type: 'scatter'"
  write(line,'edgeR_MDS_Aplot_coordinates_mqc.csv',append = FALSE)
  suppressWarnings(write.table(MDSxy, 'edgeR_MDS_Aplot_coordinates_mqc.csv', sep = ",",
                               quote=FALSE, row.names=T, col.names=NA, append=TRUE))
}

# Calculate the Pearsons correlation between samples
# Plot a heatmap of correlations
pdf('edgeR_log2CPM_sample_correlation_heatmap.pdf')
hmap <- heatmap.2(as.matrix(cor(logcpm, method="pearson")),
                  key.title="Pearson's Correlation", trace="none",
                  dendrogram="row", margin=c(9, 9))
invisible(dev.off())

# Write correlation values to file
line <-"#plot_type: 'heatmap'"
write(line,'edgeR_log2CPM_sample_correlation_mqc.csv',append = FALSE)
suppressWarnings(write.table(hmap$carpet, 'edgeR_log2CPM_sample_correlation_mqc.csv',sep = ",",
                             quote=FALSE, row.names=T, col.names=NA, append=TRUE))

# Plot the heatmap dendrogram
pdf('edgeR_log2CPM_sample_distances_dendrogram.pdf')
hmap <- heatmap.2(as.matrix(dist(t(logcpm))))
plot(hmap$rowDendrogram, main="Sample Pearson's Correlation Clustering")
invisible(dev.off())


##### main NGSmodule QC #####
setwd(QCpath)
group_num <- length(unique(sample_info[,3]))
if(group_num<=5){
  col_color <- nord("mountain_forms")[1:group_num]
} 
if(group_num>5 & group_num<=8){
  col_color <- pal_nejm("default")(group_num)
}
if(group_num>8 & group_num<=16){
  col_color <- pal_simpsons("springfield")(group_num)
}
if(group_num>16 & group_num<=20){
  col_color <- pal_d3("category20")(group_num)
}
if(group_num>20 & group_num<=51){
  col_color <- pal_igv("default")(group_num)
}

names(col_color)<- unique(sample_info[,3])

plot_list <- list()

##### Dendrogram #####
cat(">>> Hierarchical Clustering\n")
hc<- hclust(dist(t(logcpm_scale),method = "euclidean")) ## sqrt(sum((x_i - y_i)^2)). Large expressed genes -large vars. 
p <- fviz_dend(hc, rect = FALSE, cex = 0.6,horiz=T,palette="simpsons",
               label_cols=col_color[sample_info[colnames(logcpm),3]][hc$order])+
  labs(y="")
plot_list[["Hierarchical_Clustering"]] <- p
# ggsave(p0,filename = "Hierarchical_Clustering.pdf",width = max(ncol(quant_matrix)*0.3,3),height = max(ncol(quant_matrix)*0.3,3))

##### correlation heatmap #####
cat(">>> Correlation Heatmap\n")
logcpm_cor<- round(cor(logcpm,method = "pearson"),2) ### cannot use a row-global scale for paired correlation
q1 <- quantile(logcpm_cor,0.99)
q2 <- quantile(logcpm_cor,0.01)
color_palette <- colorRampPalette(brewer.pal(11,"RdBu"))(100)

df_bottom_annotation <- HeatmapAnnotation(foo1 = anno_simple(x = sample_info[colnames(logcpm),3],
                                                             col = col_color[sample_info[colnames(logcpm),3]],
                                                             gp=gpar(col="black")),
                                          border =TRUE,
                                          show_legend = F,
                                          show_annotation_name=F)
df_right_annotation <- HeatmapAnnotation(foo1 = anno_simple(x = sample_info[colnames(logcpm),3],
                                                            col = col_color[sample_info[colnames(logcpm),3]],
                                                            gp=gpar(col="black")),
                                         border =TRUE,
                                         show_legend = F,
                                         show_annotation_name=F,
                                         which = "row")

r_max <- 0.408/ncol(logcpm_cor)
ht<-Heatmap(logcpm_cor, 
            col = colorRamp2(seq(q1, q2, length = 100), color_palette),
            rect_gp = gpar(type = "none"),
            cell_fun = function(j, i, x, y, w, h, fill) {
              p = logcpm_cor[i, j]
              perc= (p-min(logcpm_cor))/(max(logcpm_cor)-min(logcpm_cor))
              grid.circle(x, y, r = r_max/2*(1+perc),
                          gp=gpar(col="black", lwd=1, fill=fill))
              grid.rect(x, y, width = w, height = h,
                        gp=gpar(col="grey", lwd=1, fill="transparent"))
            },
            column_names_gp = gpar(fontsize=10),
            row_names_gp = gpar(fontsize=10),
            border = T,
            bottom_annotation = df_bottom_annotation,
            right_annotation = df_right_annotation,
            heatmap_legend_param=list(title="Pearson's correlation coefficient",
                                      title_gp = gpar(fontsize = 10, fontfamily="sans"),
                                      title_position="lefttop",
                                      grid_height = unit(4, "mm"),
                                      grid_width = unit(4, "mm"),
                                      border="black",
                                      labels_gp = gpar(fontsize = 8),
                                      legend_direction = "horizontal",
                                      legend_width = unit(4, "cm"))
)

grob <- grid.grabExpr(ComplexHeatmap::draw(ht,heatmap_legend_side = "top"))
p <- as_ggplot(grob)+theme(aspect.ratio = 1)
plot_list[["Correlation_Heatmap"]] <- p
# ggsave(p1,filename = "Correlation_Heatmap.pdf",width = max(ncol(quant_matrix)*0.3,3),height = max(ncol(quant_matrix)*0.3,3))

##### correlation paired scatter  #####
cat(">>> Correlation Paired Scatter\n")
cor_uniq<- sort(unique(c(round(cor(logcpm,method = "pearson"),2))))
color_cor<- colorRampPalette(brewer.pal(n = 9, name = "YlOrBr")[2:9])(length(cor_uniq))
names(color_cor) <- cor_uniq

custom_point <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping)+
    stat_density2d(geom="tile", aes(fill=..density..^0.25, alpha=1), contour=FALSE,show.legend = F)+
    stat_density2d(geom="tile", aes(fill=..density..^0.25,alpha=ifelse(..density..^0.25<0.4,0,1)), contour=FALSE,show.legend = F)+
    geom_point(shape=16,color="#E18727FF",size=1,alpha=0.4)+
    geom_rug(alpha=0.4)+
    geom_smooth(method = "lm", formula = "y ~ x",color = "red",alpha=0.5)+
    geom_abline(intercept = 0,slope = 1,size=1,alpha=0.5,color="black",linetype=2)+
    scale_fill_gradientn(colours = colorRampPalette(c("white", "#BC3C29FF"))(256))+
    # scale_x_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x)))+
    # scale_y_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x)))+
    theme_pubr(border=T)+
    theme(axis.text = element_text(size=10))
}

custom_density <- function(data,mapping,...){
  ggplot(data = data, mapping = mapping)+
    geom_histogram(color="black",fill="#E18727FF",bins = 30)+
    # scale_x_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x)))+
    theme_pubr(border=T)+
    theme(axis.text = element_text(size=10))
}

custom_cor <- function(data,mapping,...){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  cor<- round(cor(x = x,y = y,method = "pearson"),2)
  color<- color_cor[as.character(cor)]
  ggplot()+labs(x = NULL, y = NULL)+
    geom_text(aes(x=0.5,y=0.5,
                  label=paste0("Corr:\n",cor)),
              size=5,fontface=2,color="white")+
    theme_pubr(border=T)+
    theme(panel.background = element_rect(fill= color))
}

p <- ggpairs(as.data.frame(logcpm), 
             title = "Multivariate correlation analysis",
             xlab = paste0("log2CPM"),
             ylab = paste0("log2CPM"),
             upper = list(continuous = custom_cor),
             lower = list(continuous = custom_point),
             diag = list(continuous = custom_density))+
  theme(aspect.ratio = 1,
        strip.background = element_rect(fill="transparent",color = "black"),
        strip.text = element_text(size=12))
# plot_list[["Correlation_PairedScatter"]] <- p

# ggsave(p,filename = "Correlation_PairedScatter.png",width = ncol(logcpm)*1.5,height = ncol(logcpm)*1.5)


##### PCA #####
cat(">>> Principal Components Analysis\n")
df_pca <- prcomp(t(logcpm_scale),center = F,scale. = F)
df_pca <- summary(df_pca)
PoV <- round(df_pca$importance[2,]*100,2)
x <- df_pca$x[,"PC1"]
y <- df_pca$x[,"PC2"]
m<- max(c(abs(x),abs(y)))

df <- data.frame(x=x,y=y,sample = names(x),group = sample_info[names(x),3],stringsAsFactors = FALSE)

p<- ggplot(data = df,aes(x=x,y=y,group=group,fill=group,label=sample))+
  geom_point(shape=21)+
  geom_rug(aes(color=group))+
  geom_label_repel(size=2.5,color = "white",
                   min.segment.length = 0,segment.color = "black",segment.alpha = 0.8)+
  labs(title = "Principal Components Analysis",x=paste0("PC1(",PoV[1],"%)"),y=paste0("PC2(",PoV[2],"%)"))+
  scale_fill_manual(values = col_color)+
  scale_color_manual(values = col_color)+
  xlim(-m,m)+ylim(-m,m)+
  theme_pubr(border = T,legend = "none")+
  theme(aspect.ratio = 1,
        panel.grid.major = element_line())

plot_list[["PCA"]] <- p
# ggsave(p3,filename = "PCA_scatter.pdf",width = 5,height = 5)



##### PCA top features #####
cat(">>> PCA top features\n")
PC1_top<- sort(df_pca$rotation[,"PC1"],decreasing = T)[c(1:100,(nrow(df_pca$rotation)-99):nrow(df_pca$rotation))]
PC2_top<- sort(df_pca$rotation[,"PC2"],decreasing = T)[c(1:100,(nrow(df_pca$rotation)-99):nrow(df_pca$rotation))]

mat_PC1<- logcpm_scale[names(PC1_top),]
mat_PC2<- logcpm_scale[names(PC2_top),]

color_palette <- rev(colorRampPalette(brewer.pal(11,"Spectral"))(100))

df_bottom_annotation <- HeatmapAnnotation(foo1 = anno_simple(x = sample_info[colnames(logcpm),3],
                                                             col = col_color[sample_info[colnames(logcpm),3]],
                                                             gp=gpar(col="black")),
                                          border =TRUE,
                                          show_legend = F,
                                          show_annotation_name=F)

ht_legend <- list(title=paste0("Z-score"),
                  title_gp = gpar(fontsize = 10, fontfamily="sans"),
                  title_position="topleft",
                  grid_height = unit(3, "mm"),
                  grid_width = unit(3, "mm"),
                  border="black",
                  labels_gp = gpar(fontsize = 10),
                  legend_direction = "horizontal",
                  legend_width = unit(3, "cm"))

ht1<-Heatmap(mat_PC1, column_title="PC1 top features",
             col = colorRamp2(seq(-max(abs(mat_PC1)), max(abs(mat_PC1)), length = 100), color_palette),
             show_row_names = F,
             cluster_rows = F,
             column_title_gp = gpar(fontsize=10),
             column_names_gp = gpar(fontsize=8),
             border = T,
             bottom_annotation = df_bottom_annotation,
             heatmap_legend_param=ht_legend)

ht2<-Heatmap(mat_PC2,column_title="PC2 top features",
             col = colorRamp2(seq(-max(abs(mat_PC2)), max(abs(mat_PC2)), length = 100), color_palette),
             show_row_names = F,
             cluster_rows = F,
             column_title_gp = gpar(fontsize=10),
             column_names_gp = gpar(fontsize=8),
             border = T,
             bottom_annotation = df_bottom_annotation,
             heatmap_legend_param=ht_legend)


grob <- grid.grabExpr(ComplexHeatmap::draw(ht1+ht2,heatmap_legend_side = "top"))
p <- as_ggplot(grob)+theme(aspect.ratio = 28/(ncol(logcpm)+12))

plot_list[["PCA_Heatmap"]] <- p

# ggsave(p,filename = "PCA_Heatmap.pdf",width = max(ncol(quant_matrix_log2_scale)*0.3,3),height = 7)


##### enrichment analysis for important genes in PCA (ENSEMBL ID only) #####
# if (length(grep(pattern = "^ENS",x = names(PC1_top200),perl = TRUE,useBytes = ))==200) {
#   PC1_top200_id <- names(PC1_top200)
#   PC2_top200_id <- names(PC2_top200)
#   
#   
#   
#   
#   
# }

##### output report #####
pdf("./postQuantificationQC.pdf",width = 8,height = 5)
invisible(lapply(plot_list, print))
invisible(dev.off())


##### check whether the unwanted file exists and remove it ##### 
if (file.exists("Rplots.pdf")) {
  invisible(file.remove("Rplots.pdf"))
}
