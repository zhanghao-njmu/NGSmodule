#!/usr/bin/env Rscript
suppressWarnings(suppressPackageStartupMessages(invisible(lapply(
  c("BiocParallel","edgeR","DESeq2","dplyr","stringr","scales","RColorBrewer",
    "ggpubr","ggsci","ggforce","reshape2","VennDiagram","gridExtra",
    "gplots","openxlsx","ggalluvial","ggfittext","ComplexHeatmap","circlize","nord","ggupset"), require, character.only = TRUE))))
register(MulticoreParam())

args <- commandArgs( trailingOnly = TRUE )
maindir <- args[1]
aligner <- args[2]
SampleInfoFile <- args[3]
comparison <- args[4]
max_padj <- as.numeric(args[5])
min_fc <- as.numeric(args[6])
min_count <- as.numeric(args[7])
DGEs_multi_compare <- as.numeric(args[8])
script_path <- as.character(args[9])

######## example 1 ############
# setwd("/data/lab/ZhangHao/tmp/NGSmodule_analysis/DifferentialExpression/")
# maindir <- "/data/lab/ZhangHao/tmp/"
# aligner <- "hisat2"
# SampleInfoFile <- "/data/lab/ZhangHao/tmp/temp_20200413205117.Sample_info.csv"
# comparison <- "Hom,WT"
# max_padj <- 0.05
# min_fc <- 2
# min_count <- 10
# DGEs_multi_compare <- 1
# script_path <- "/data/lab/ZhangHao/NGSmodule/DifferentialExpression.R"
##############################

######### example 2 ############
# setwd("/data/lab/HeXi/Rpl39l/RNC-seq/NGSmodule_analysis/DifferentialExpression/")
# maindir <- "/data/lab/HeXi/Rpl39l/RNC-seq/"
# aligner <- "hisat2"
# SampleInfoFile <- "/data/lab/HeXi/Rpl39l/RNC-seq/Sample_info.csv"
# comparison <- "Hom-Input,WT-Input;Hom-80S,WT-80S;Hom-SPS,WT-SPS;Hom-LPS,WT-LPS"
# max_padj <- 0.05
# min_fc <- 2
# min_count <- 10
# DGEs_multi_compare <- 1
# script_path <- "/data/lab/ZhangHao/NGSmodule/DifferentialExpression.R"
##############################

##### source plot function #####

script_dir <- gsub(x = script_path,pattern = "DifferentialExpression.R",replacement = "")
source(paste0(script_dir,"/DifferentialExpression_function.R"))

##### Load data ######

count_file <- paste0(maindir,"/NGSmodule_analysis/Quantification/Quantification.",aligner,".count.tab")

if (!file.exists(count_file)) {
  print(paste0("Can not find count file: ",count_file,"\nPlease run NGSmodule Quantification first!\n"))
  quit(status=1)
}
if (!file.exists(SampleInfoFile)) {
  print(paste0("Can not find SampleInfoFile: ",SampleInfoFile,"\nPlease check the config parameters!\n"))
  quit(status=1)
}

counts_matrix_raw <- read.table(file = count_file, header = T,sep = "\t",row.names = 1,stringsAsFactors = F,check.names=F)
counts_matrix <- counts_matrix_raw[,which(str_detect(colnames(counts_matrix_raw),pattern = paste0(".",aligner,".count"))),drop=FALSE]
colnames(counts_matrix) <- gsub(x=colnames(counts_matrix),pattern = paste0(".",aligner,".count"),replacement = "")

start<- which(str_detect(colnames(counts_matrix_raw),">>Annotation"))[1]
annotation_matrix <- counts_matrix_raw[,start:ncol(counts_matrix_raw)]

sample_info<- read.csv(SampleInfoFile,stringsAsFactors = FALSE)
sample_info <- sample_info[order(sample_info[,3]),]



comparison_list <- lapply(as.list(str_split(string = comparison,pattern = ";")[[1]]), function(x){str_split(x,",")[[1]]})
comparison_list_name <- gsub(x = str_split(string = comparison,pattern = ";")[[1]],pattern = ",",replacement = "_vs_")
names(comparison_list) <- comparison_list_name


###### START #####

res_list <- list()
res_plot <- list()
for (compare in comparison_list) {
  label <- paste0(compare,collapse = "/")
  if (!all(compare%in%sample_info[,3])) {
    print(paste0("Comparison is wrong : ",label,"\nPlease check the Group colmn of SampleInfoFile and comparison parameters in the config file!\n"))
    quit(status=1)
  }
  
  cat("+++++++++++++++","Comparing",label,"+++++++++++++++\n\n")
  sample1 <- sample_info[which(sample_info[,3]==compare[1]),2]
  sample2 <- sample_info[which(sample_info[,3]==compare[2]),2]
  compare_counts <- counts_matrix[,c(sample1,sample2)]
  compare_groups <- factor(c(rep("Group1",length(sample1)),rep("Group2",length(sample2))))
  
  lib_size <- colSums(compare_counts)/1000000
  general_p1<- count_bar(lib_size = lib_size,label = label)
  general_p2<- var_mean_plot(count_matrix = compare_counts,group = compare_groups,label = label)
  
  res_plot[[paste0("general_p1-",label)]] <- general_p1
  res_plot[[paste0("general_p2-",label)]] <- general_p2
  
  wb <- createWorkbook(title="NGSmodule DifferentialExpression analysis")
  if (min(length(sample1),length(sample2))==1) {
    
    ##### edgeR #####
    genelist <- DGEList(counts=compare_counts,group=relevel(compare_groups,ref = "Group2"))
    keep <- filterByExpr(genelist,min.count = min_count)
    cat("*** edgeR Number of genes tested: ",sum(keep),"\n",sep = "")
    genelist_filter <- genelist[keep,, keep.lib.sizes=FALSE]
    genelist_norm <- calcNormFactors(genelist_filter,method = "TMM") ### cannot used for global differential expression
    edgeR_log_norm <- as.data.frame(cpm(genelist_norm, log=TRUE))
    edgeR_log_norm <- cbind(rownames(edgeR_log_norm),edgeR_log_norm)
    colnames(edgeR_log_norm)[1] <- "RowID"
    genelist_bcv <- genelist_norm
    bcv <- 0.2
    res <- exactTest(genelist_bcv, dispersion=bcv^2)
    res <- topTags(res, n = nrow(res$table))$table
    res[["DifferentialExpression"]] <- ifelse(res$FDR<max_padj & abs(res$logFC)>log2(min_fc),
                                              ifelse(res$logFC>0,"Up-regulated","Down-regulated"),"Non-significant")
    edgeR_up <- rownames(res)[which(res[["DifferentialExpression"]]=="Up-regulated")]
    edgeR_down <- rownames(res)[which(res[["DifferentialExpression"]]=="Down-regulated")]
    edgeR_nonsig <- rownames(res)[which(res[["DifferentialExpression"]]=="Non-significant")]
    edgeR_drop <- setdiff(x = rownames(compare_counts),y = c(edgeR_up,edgeR_down,edgeR_nonsig))
    cat("*** edgeR Number of genes tested/droped: ",sum(keep),"/",(nrow(compare_counts)-sum(keep)),"\n",sep = "")
    cat(">>> edgeR All DGEs: ",length(c(edgeR_up,edgeR_down)),"  Up-/Down-regulated: ",length(edgeR_up),"/",length(edgeR_down),"\n",
        ">>> edgeR Non-significant: ",length(edgeR_nonsig),"\n\n",
        sep = "")
    res <- merge(x=res,y=annotation_matrix,by="row.names",all.x=T)
    colnames(res)[which(colnames(res)=="Row.names")] <- "RowID"
    rownames(res) <- res[,"RowID"]
    res_list[[paste0("edgeR_res@",paste0(compare,collapse = "_vs_"))]] <- res[with(res,order(FDR)),]
    res_list[[paste0("edgeR_LogNorm@",paste0(compare,collapse = "_vs_"))]] <- edgeR_log_norm
    
    addWorksheet(wb, sheetName = "edgeR_LogNorm")
    writeDataTable(wb = wb,
                   x = edgeR_log_norm,
                   sheet="edgeR_LogNorm",
                   withFilter=F,
                   tableStyle = "TableStyleLight8")
    
    sheetname <- paste0("edgeR@U",length(edgeR_up),"-D",length(edgeR_down),"-N",length(edgeR_nonsig))
    addWorksheet(wb, sheetName = sheetname)
    writeDataTable(wb = wb,
                   x = res[with(res,order(FDR)),],
                   sheet=sheetname,
                   withFilter=F,
                   tableStyle = "TableStyleLight8")
    saveWorkbook(wb, file = paste0("./DGEs_data/",paste0(compare,collapse = "_vs_"),".xlsx"), overwrite = TRUE)
    
    edgeR_p1 <- MDplot(count_matrix = compare_counts,
                       group = compare_groups,
                       up = edgeR_up,down = edgeR_down,nonsig = edgeR_nonsig,drop = edgeR_drop,
                       label = paste("edgeR:",label))
    edgeR_p2 <- pvalue_distibution(pvalue = res$PValue,padj = res$FDR,label = paste("edgeR:",label))
    edgeR_p3 <- volcano_plot(count_matrix = compare_counts,
                             res = res,col_padj = "FDR",col_log2fc = "logFC",col_DEgroup = "DifferentialExpression",
                             DEgroup_level = c("Up-regulated","Down-regulated","Non-significant"),
                             max_padj = max_padj,min_fc = min_fc,label = paste("edgeR:",label))
    if (length(c(edgeR_up,edgeR_down)>0)) {
      row_split <- setNames(object = c(rep(paste0("Up-regulated (",length(edgeR_up),")"),length(edgeR_up)),
                                       rep(paste0("Down-regulated (",length(edgeR_down),")"),length(edgeR_down))),
                            nm = c(edgeR_up,edgeR_down))
      row_split <- factor(row_split,levels = c(paste0("Up-regulated (",length(edgeR_up),")"),
                                               paste0("Down-regulated (",length(edgeR_down),")")))
      column_split <- setNames(object = rep(compare,table(compare_groups)),nm = colnames(edgeR_log_norm[,-1]))
      column_split <- factor(column_split,levels = unique(column_split))
      edgeR_p4 <- heatmap(LogNormCounts = edgeR_log_norm[,-1],
                          genes=c(edgeR_up,edgeR_down),
                          res=res,res_log2fc_col="logFC",res_padj_col="FDR",
                          row_split=row_split,column_split=column_split,
                          label = paste("edgeR:",label))
    }else{
      edgeR_p4 <- NULL
    }
    edgeR_p5 <- pie_chart(res = res,annotation_matrix = annotation_matrix,label = paste("edgeR:",label))
    edgeR_p6 <- gene_test_plot(test = c(edgeR_up,edgeR_down,edgeR_nonsig),drop_low = edgeR_drop,
                               annotation_matrix = annotation_matrix,label = paste("edgeR:",label))
    
    for(i in 1:6){
      var_name <- paste0("edgeR_p",i)
      if (!is.null(get(var_name))) {
        res_plot[[paste0(var_name,"-",label)]] <- get(var_name)
      }else{
        next
      }
    }
    
    edgeR_minCPM <- min(rowMeans(cpm(compare_counts))[keep])
    general_p3<- cpm_gene_curve(count_matrix = compare_counts,cpm_seq = seq(0.1,50,0.1),label = label,
                                minCPM = edgeR_minCPM,minCPM_label = "edgeR")
    res_plot[[paste0("general_p3-",label)]] <- general_p3
    
  }else{
    
    ##### edgeR #####
    genelist <- DGEList(counts=compare_counts,group=relevel(compare_groups,ref = "Group2"))
    keep <- filterByExpr(genelist,min.count = min_count)
    genelist_filter <- genelist[keep,, keep.lib.sizes=FALSE]
    genelist_norm <- calcNormFactors(genelist_filter,method = "TMM") ### cannot used for global differential expression
    edgeR_log_norm <- as.data.frame(cpm(genelist_norm, log=TRUE))
    edgeR_log_norm <- cbind(rownames(edgeR_log_norm),edgeR_log_norm)
    colnames(edgeR_log_norm)[1] <- "RowID"
    design <- model.matrix(~0 + factor(compare_groups))
    colnames(design) <- levels(factor(compare_groups))
    genelist_Disp <- estimateDisp(genelist_norm, design, robust = TRUE)
    fit <- glmQLFit(genelist_Disp, design,robust = TRUE)
    res <- glmQLFTest(fit, contrast = makeContrasts("Group1-Group2", levels = design))
    res <- topTags(res, n = nrow(res$table))$table
    res[["DifferentialExpression"]] <- ifelse(res$FDR<max_padj & abs(res$logFC)>log2(min_fc),
                                              ifelse(res$logFC>0,"Up-regulated","Down-regulated"),"Non-significant")
    edgeR_up <- rownames(res)[which(res[["DifferentialExpression"]]=="Up-regulated")]
    edgeR_down <- rownames(res)[which(res[["DifferentialExpression"]]=="Down-regulated")]
    edgeR_nonsig <- rownames(res)[which(res[["DifferentialExpression"]]=="Non-significant")]
    edgeR_drop <- setdiff(x = rownames(compare_counts),y = c(edgeR_up,edgeR_down,edgeR_nonsig))
    cat("*** edgeR Number of genes tested/droped: ",sum(keep),"/",(nrow(compare_counts)-sum(keep)),"\n",sep = "")
    cat(">>> edgeR All DGEs: ",length(c(edgeR_up,edgeR_down)),"  Up-/Down-regulated: ",length(edgeR_up),"/",length(edgeR_down),"\n",
        ">>> edgeR Non-significant: ",length(edgeR_nonsig),"\n\n",
        sep = "")
    res <- merge(x=res,y=annotation_matrix,by="row.names",all.x=T)
    colnames(res)[which(colnames(res)=="Row.names")] <- "RowID"
    rownames(res) <- res[,"RowID"]
    res_list[[paste0("edgeR_res@",paste0(compare,collapse = "_vs_"))]] <- res[with(res,order(FDR)),]
    res_list[[paste0("edgeR_LogNorm@",paste0(compare,collapse = "_vs_"))]] <- edgeR_log_norm
    
    addWorksheet(wb, sheetName = "edgeR_LogNorm")
    writeDataTable(wb = wb,
                   x = edgeR_log_norm,
                   sheet="edgeR_LogNorm",
                   withFilter=F,
                   tableStyle = "TableStyleLight8")
    
    sheetname <- paste0("edgeR@U",length(edgeR_up),"-D",length(edgeR_down),"-N",length(edgeR_nonsig))
    addWorksheet(wb, sheetName = sheetname)
    writeDataTable(wb = wb,
                   x = res[with(res,order(FDR)),],
                   sheet=sheetname,
                   withFilter=F,
                   tableStyle = "TableStyleLight8")
    
    edgeR_p1 <- MDplot(count_matrix = compare_counts,
                       group = compare_groups,
                       up = edgeR_up,down = edgeR_down,nonsig = edgeR_nonsig,drop = edgeR_drop,
                       label = paste("edgeR:",label))
    edgeR_p2 <- pvalue_distibution(pvalue = res$PValue,padj = res$FDR,label = paste("edgeR:",label))
    edgeR_p3 <- volcano_plot(count_matrix = compare_counts,
                             res = res,col_padj = "FDR",col_log2fc = "logFC",col_DEgroup = "DifferentialExpression",
                             DEgroup_level = c("Up-regulated","Down-regulated","Non-significant"),
                             max_padj = max_padj,min_fc = min_fc,label = paste("edgeR:",label))
    if (length(c(edgeR_up,edgeR_down)>0)) {
      row_split <- setNames(object = c(rep(paste0("Up-regulated (",length(edgeR_up),")"),length(edgeR_up)),
                                       rep(paste0("Down-regulated (",length(edgeR_down),")"),length(edgeR_down))),
                            nm = c(edgeR_up,edgeR_down))
      row_split <- factor(row_split,levels = c(paste0("Up-regulated (",length(edgeR_up),")"),
                                               paste0("Down-regulated (",length(edgeR_down),")")))
      column_split <- setNames(object = rep(compare,table(compare_groups)),nm = colnames(edgeR_log_norm[,-1]))
      column_split <- factor(column_split,levels = unique(column_split))
      edgeR_p4 <- heatmap(LogNormCounts = edgeR_log_norm[,-1],
                          genes=c(edgeR_up,edgeR_down),
                          res=res,res_log2fc_col="logFC",res_padj_col="FDR",
                          row_split=row_split,column_split=column_split,
                          label = paste("edgeR:",label))
    }else{
      edgeR_p4 <- NULL
    }
    
    edgeR_p5 <- pie_chart(res = res,annotation_matrix = annotation_matrix,label = paste("edgeR:",label))
    edgeR_p6 <- gene_test_plot(test = c(edgeR_up,edgeR_down,edgeR_nonsig),drop_low = edgeR_drop,
                               annotation_matrix = annotation_matrix,label = paste("edgeR:",label))
    
    for(i in 1:6){
      var_name <- paste0("edgeR_p",i)
      if (!is.null(get(var_name))) {
        res_plot[[paste0(var_name,"-",label)]] <- get(var_name)
      }else{
        next
      }
    }
    
    ##### DESeq2 #####
    colData <- data.frame(row.names=colnames(compare_counts),condition=compare_groups)
    dds <- DESeqDataSetFromMatrix(countData = compare_counts,colData = colData, design = ~ condition)
    dds <- dds[keep,]
    dds <- suppressMessages(DESeq(dds, parallel = TRUE))
    DESeq2_log_norm <- as.data.frame(cpm(counts(dds,normalized = T), log=TRUE))
    DESeq2_log_norm <- cbind(rownames(DESeq2_log_norm),DESeq2_log_norm)
    colnames(DESeq2_log_norm)[1] <- "RowID"
    dds_res <- results(dds,alpha = 0.05,contrast=c("condition", "Group1", "Group2"))  ### Group1/Group2
    res <- as.data.frame(dds_res)
    res <- na.omit(res)
    res[["DifferentialExpression"]] <- ifelse(res$padj<max_padj & abs(res$log2FoldChange)>log2(min_fc),
                                              ifelse(res$log2FoldChange>0,"Up-regulated","Down-regulated"),"Non-significant")
    DESeq2_up <- rownames(res)[which(res[["DifferentialExpression"]]=="Up-regulated")]
    DESeq2_down <- rownames(res)[which(res[["DifferentialExpression"]]=="Down-regulated")]
    DESeq2_nonsig <- rownames(res)[which(res[["DifferentialExpression"]]=="Non-significant")]
    DESeq2_drop <- setdiff(x = rownames(compare_counts),y = c(DESeq2_up,DESeq2_down,DESeq2_nonsig))
    DESeq2_drop_outlier <- rownames(dds_res)[is.na(dds_res[,"pvalue"])&is.na(dds_res[,"padj"])]
    DESeq2_drop_low <- DESeq2_drop[which(!DESeq2_drop%in%DESeq2_drop_outlier)]
    cat("*** DESeq2 Number of genes tested/droped: ",sum(length(c(DESeq2_up,DESeq2_down,DESeq2_nonsig))),"/",(nrow(compare_counts)-sum(length(c(DESeq2_up,DESeq2_down,DESeq2_nonsig)))),"\n",sep = "")
    cat(">>> DESeq2 All DGEs: ",length(c(DESeq2_up,DESeq2_down)),"  Up-/Down-regulated: ",length(DESeq2_up),"/",length(DESeq2_down),"\n",
        ">>> DESeq2 Non-significant: ",length(DESeq2_nonsig),"\n\n",
        sep = "")
    res <- merge(x=res,y=annotation_matrix,by="row.names",all.x=T)
    colnames(res)[which(colnames(res)=="Row.names")] <- "RowID"
    rownames(res) <- res[,"RowID"]
    res_list[[paste0("DESeq2_res@",paste0(compare,collapse = "_vs_"))]] <- res[with(res,order(padj)),]
    res_list[[paste0("DESeq2_LogNorm@",paste0(compare,collapse = "_vs_"))]] <- DESeq2_log_norm
    
    addWorksheet(wb, sheetName = "DESeq2_LogNorm")
    writeDataTable(wb = wb,
                   x = DESeq2_log_norm,
                   sheet="DESeq2_LogNorm",
                   withFilter=F,
                   tableStyle = "TableStyleLight8")
    
    sheetname <- paste0("DESeq2@U",length(DESeq2_up),"-D",length(DESeq2_down),"-N",length(DESeq2_nonsig))
    addWorksheet(wb, sheetName = sheetname)
    writeDataTable(wb = wb,
                   x = res[with(res,order(padj)),],
                   sheet=sheetname,
                   withFilter=F,
                   tableStyle = "TableStyleLight8")
    
    DESeq2_p1 <- MDplot(count_matrix = compare_counts,
                        group = compare_groups,
                        up = DESeq2_up,down = DESeq2_down,nonsig = DESeq2_nonsig,drop = DESeq2_drop,
                        label = paste("DESeq2:",label))
    DESeq2_p2 <- pvalue_distibution(pvalue = dds_res$pvalue,padj = dds_res$padj,label = paste("DESeq2:",label))
    DESeq2_p3 <- filterNumRej_plot(dds_res = dds_res,max_padj = max_padj,count_matrix = compare_counts,label = paste("DESeq2:",label))
    DESeq2_p4 <- volcano_plot(count_matrix = compare_counts,
                              res = res,col_padj = "padj",col_log2fc = "log2FoldChange",col_DEgroup = "DifferentialExpression",
                              DEgroup_level = c("Up-regulated","Down-regulated","Non-significant"),
                              max_padj = max_padj,min_fc = min_fc,label = paste("DESeq2:",label))
    if (length(c(DESeq2_up,DESeq2_down)>0)) {
      row_split <- setNames(object = c(rep(paste0("Up-regulated (",length(DESeq2_up),")"),length(DESeq2_up)),
                                       rep(paste0("Down-regulated (",length(DESeq2_down),")"),length(DESeq2_down))),
                            nm = c(DESeq2_up,DESeq2_down))
      row_split <- factor(row_split,levels = c(paste0("Up-regulated (",length(DESeq2_up),")"),
                                               paste0("Down-regulated (",length(DESeq2_down),")")))
      column_split <- setNames(object = rep(compare,table(compare_groups)),nm = colnames(DESeq2_log_norm[,-1]))
      column_split <- factor(column_split,levels = unique(column_split))
      DESeq2_p5 <- heatmap(LogNormCounts = DESeq2_log_norm[,-1],
                           genes=c(DESeq2_up,DESeq2_down),
                           res=res,res_log2fc_col="log2FoldChange",res_padj_col="padj",
                           row_split=row_split,column_split=column_split,
                           label = paste("DESeq2:",label))
    }else{
      DESeq2_p5 <- NULL
    }
    DESeq2_p6 <- pie_chart(res = res,annotation_matrix = annotation_matrix,label = paste("DESeq2:",label))
    DESeq2_p7 <- gene_test_plot(test = c(DESeq2_up,DESeq2_down,DESeq2_nonsig),
                                drop_low = DESeq2_drop_low,drop_outlier = DESeq2_drop_outlier,
                                annotation_matrix = annotation_matrix,label = paste("DESeq2:",label))
    
    for(i in 1:7){
      var_name <- paste0("DESeq2_p",i)
      if (!is.null(get(var_name))) {
        res_plot[[paste0(var_name,"-",label)]] <- get(var_name)
      }else{
        next
      }
    }
    
    ##### scaling factor #####
    edgeR_sf <- setNames(object = genelist_norm$samples$norm.factors,nm = rownames(genelist_norm$samples))
    DESeq2_sf <- dds$sizeFactor
    general_p3 <- scale_factor_compare(edgeR_sf = edgeR_sf,DESeq2_sf = DESeq2_sf,lib_size=lib_size)
    res_plot[[paste0("general_p3-",label)]] <- general_p3
    
    ##### cpm-gene curve #####
    edgeR_minCPM <- min(rowMeans(cpm(compare_counts))[c(edgeR_up,edgeR_down,edgeR_nonsig)])
    DESeq2_minCPM <- min(rowMeans(cpm(compare_counts))[c(DESeq2_up,DESeq2_down,DESeq2_nonsig)])
    
    general_p4<- cpm_gene_curve(count_matrix = compare_counts,cpm_seq = seq(0.1,50,0.1),label = label,
                                minCPM = c(edgeR_minCPM,DESeq2_minCPM),minCPM_label = factor(c("edgeR","DESeq2"),levels = c("edgeR","DESeq2")))
    res_plot[[paste0("general_p4-",label)]] <- general_p4
    
    ##### Compare DGEs result #####
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
    venn_list = list("DESeq2"=c(DESeq2_up,DESeq2_down),"edgeR"=c(edgeR_up,edgeR_down))
    venn_p1 <- venn.diagram(x= venn_list,filename = NULL,
                            height = 1500,width = 1800, resolution = 500,
                            category.names =paste0(c("DESeq2","edgeR"),"\n(",lapply(venn_list, length) %>% unlist(),")"),
                            scaled=T,force.unique=T,print.mode="raw",
                            lwd=2,lty=1,
                            col= pal_nejm(palette = "default")(8)[1:2],
                            fill= alpha(pal_nejm(palette = "default")(8)[1:2],0.3),
                            cex=0.8,
                            fontface = "bold",
                            fontfamily = "sans",
                            cat.cex = 1,
                            cat.default.pos = "outer",
                            cat.dist = c(0.07, 0.07),
                            cat.fontfamily = "sans",
                            cat.col = pal_nejm(palette = "default")(8)[1:2],
                            margin=0.3
    )
    venn_p1 <- grobTree(venn_p1)
    venn_p1 <- addGrob(gTree = venn_p1,child = textGrob(label = paste0("DGEs of ",label), x = unit(0.5, "npc"), y = unit(0.75, "npc")))
    venn_p1 <- as_ggplot(venn_p1)+theme(aspect.ratio = 1)
    
    venn_list = list("DESeq2-Up"=DESeq2_up,"DESeq2-Down"=DESeq2_down,"edgeR-Up"=edgeR_up,"edgeR-Down"=edgeR_down)
    venn_p2 <- venn.diagram(x= venn_list,filename = NULL,
                            height = 1500,width = 1800, resolution = 500,
                            category.names =paste0(c("DESeq2-Up","DESeq2-Down","edgeR-Up","edgeR-Down"),"\n(",lapply(venn_list, length) %>% unlist(),")"),
                            scaled=T,force.unique=T,print.mode="raw",
                            lwd=2,lty=1,
                            col= pal_nejm(palette = "default")(8)[1:4],
                            fill= alpha(pal_nejm(palette = "default")(8)[1:4],0.3),
                            cex=0.8,
                            fontface = "bold",
                            fontfamily = "sans",
                            cat.cex = 1,
                            cat.default.pos = "outer",
                            cat.dist = c(0.24,0.24,0.14,0.14),
                            cat.fontfamily = "sans",
                            cat.col = pal_nejm(palette = "default")(8)[1:4],
                            margin=0.3
    )
    venn_p2 <- grobTree(venn_p2)
    venn_p2 <- addGrob(gTree = venn_p2,child = textGrob(label = paste0("DGEs of ",label), x = unit(0.5, "npc"), y = unit(0.83, "npc")))
    venn_p2 <- as_ggplot(venn_p2)+theme(aspect.ratio = 1)
    
    res_plot[[paste0("venn_p1-",label)]] <- venn_p1
    res_plot[[paste0("venn_p2-",label)]] <- venn_p2
    
    venn_output <- venn(data=venn_list,show.plot = F)
    intersections<- stack(attributes(venn_output)$intersections)
    colnames(intersections) <- c("RowID","Class")
    intersections[,"Algorithm"] <- gsub(x = intersections[,"Class"],pattern = "(-Up)|(-Down)",replacement = "")
    res_list[[paste0("Compare@",paste0(compare,collapse = "_vs_"))]] <- intersections
    
    sheetname <- "Compare_Algorithm"
    addWorksheet(wb, sheetName = sheetname)
    writeDataTable(wb = wb,
                   x = intersections,
                   sheet=sheetname,
                   withFilter=F,
                   tableStyle = "TableStyleLight8")
    
    DGEs_stat_1<- as.data.frame(table(intersections[,"Algorithm"]))
    colnames(DGEs_stat_1) <- c("Algorithm","Algorithm_Count")
    DGEs_stat_2<- as.data.frame(table(intersections[,"Class"]))
    colnames(DGEs_stat_2) <- c("Class","Class_Count")
    DGEs_stat_1[(nrow(DGEs_stat_1)+1):nrow(DGEs_stat_2),] <- NA
    DGEs_stat<- cbind(DGEs_stat_1,DGEs_stat_2)
    
    addWorksheet(wb, sheetName = "Stat")
    writeDataTable(wb = wb,
                   x = DGEs_stat,
                   sheet="Stat",
                   withFilter=F,
                   tableStyle = "TableStyleLight8")
    saveWorkbook(wb, file = paste0("./DGEs_data/DGEs@",paste0(compare,collapse = "_vs_"),".xlsx"), overwrite = TRUE)
    
  }
  res_plot_sectect <- res_plot[which(str_detect(string = names(res_plot),pattern = label))]
  
  # res_plot_Grob<- suppressWarnings(sapply(res_plot_sectect, ggplotGrob))
  # class(res_plot_Grob) <- c("arrangelist", class(res_plot_Grob))
  # suppressMessages(suppressWarnings(
  #   ggsave(filename = paste0("./DGEs_plot/DGEs_plot-",paste0(compare,collapse = "_vs_"),".pdf"), plot = res_plot_Grob,
  #          height = 5,width = 8)
  # ))
  
  pdf(paste0("./DGEs_plot/DGEs_plot-",paste0(compare,collapse = "_vs_"),".pdf"),width = 8,height = 5)
  invisible(lapply(res_plot_sectect, print))
  invisible(dev.off())
}
saveRDS(object = res_list,file = "./DGEs_rds/res_list.rds")
saveRDS(object = res_plot,file = "./DGEs_rds/res_plot.rds")




##### Compare DGEs among different comparison #####
if (length(comparison_list)>=2 & DGEs_multi_compare==1) {
  cat("\n+++++++++++++++ Comparing DGEs ... +++++++++++++++\n\n")
  compare_plot <- list()
  compare_list <- list()
  
  # UpsetR ------------------------------------------------------------------
  pal <- pal_material("blue-grey",n = 10,reverse = T)(10)
  
  ## edgeR
  res_edgeR <- res_list[str_detect(names(res_list),"edgeR_res")]
  res_edgeR<- lapply(names(res_edgeR), function(x){
    df <- res_edgeR[[x]]
    df[,"compare"] <- gsub(x = x,pattern = "edgeR_res@",replacement = "")
    df[,"compare_DE"] <- paste(df[,"compare"],df[,"DifferentialExpression"],sep = ":")
    return(df)
  })
  res_edgeR <- bind_rows(res_edgeR)
  res_edgeR <- res_edgeR[res_edgeR[,"DifferentialExpression"]!="Non-significant",]
  p1 <- res_edgeR %>% group_by(RowID) %>% 
    mutate(upset = list(compare)) %>%
    ggplot(aes(x = upset)) +
    geom_bar(aes(fill=..count..),color="black",width = 0.5) +
    geom_text(aes(label=..count..),stat="count",vjust=-0.5)+
    labs(title="Upset plot for edgeR DGEs",x="",y="Count")+
    scale_x_upset()+
    scale_y_continuous(expand = expansion(0.1))+
    scale_fill_material(name="Count",palette = "blue-grey")+
    theme(axis.line = element_line(),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          panel.background=element_rect(fill = "white",colour = NA),
          panel.grid.major = element_line(colour = "grey92"))+
    theme_combmatrix(combmatrix.label.text = element_text(size=10,color="black"),
                     combmatrix.label.extra_spacing=6)

  p2 <- res_edgeR %>% group_by(RowID) %>% 
    mutate(upset = list(compare_DE)) %>%
    ggplot(aes(x = upset)) +
    geom_bar(aes(fill=..count..),color="black",width = 0.5) +
    geom_text(aes(label=..count..),stat="count",vjust=-0.5)+
    labs(title="Upset plot for edgeR DGEs",x="",y="Count")+
    scale_x_upset()+
    scale_y_continuous(expand = expansion(0.1))+
    scale_fill_material(name="Count",palette = "blue-grey")+
    theme(axis.line = element_line(),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          panel.background=element_rect(fill = "white",colour = NA),
          panel.grid.major = element_line(colour = "grey92"))+
    theme_combmatrix(combmatrix.label.text = element_text(size=10,color="black"),
                     combmatrix.label.extra_spacing=6)
  
  compare_plot[[paste0("Compare_edgeR@upset1")]] <- p1
  compare_plot[[paste0("Compare_edgeR@upset1")]] <- p2
  
   ## DESeq2
  if (sum(str_detect(string = names(res_list),"DESeq2"))>0) {
    res_DESeq2 <- res_list[str_detect(names(res_list),"DESeq2_res")]
    res_DESeq2<- lapply(names(res_DESeq2), function(x){
      df <- res_DESeq2[[x]]
      df[,"compare"] <- gsub(x = x,pattern = "DESeq2_res@",replacement = "")
      df[,"compare_DE"] <- paste(df[,"compare"],df[,"DifferentialExpression"],sep = ":")
      return(df)
    })
    res_DESeq2 <- bind_rows(res_DESeq2)
    res_DESeq2 <- res_DESeq2[res_DESeq2[,"DifferentialExpression"]!="Non-significant",]
    p1 <- res_DESeq2 %>% group_by(RowID) %>% 
      mutate(upset = list(compare)) %>%
      ggplot(aes(x = upset)) +
      geom_bar(aes(fill=..count..),color="black",width = 0.5) +
      geom_text(aes(label=..count..),stat="count",vjust=-0.5)+
      labs(title="Upset plot for DESeq2 DGEs",x="",y="Count")+
      scale_x_upset()+
      scale_y_continuous(expand = expansion(0.1))+
      scale_fill_material(name="Count",palette = "blue-grey")+
      theme(axis.line = element_line(),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            panel.background=element_rect(fill = "white",colour = NA),
            panel.grid.major = element_line(colour = "grey92"))+
      theme_combmatrix(combmatrix.label.text = element_text(size=10,color="black"),
                       combmatrix.label.extra_spacing=6)
    
    p2 <- res_DESeq2 %>% group_by(RowID) %>% 
      mutate(upset = list(compare_DE)) %>%
      ggplot(aes(x = upset)) +
      geom_bar(aes(fill=..count..),color="black",width = 0.5) +
      geom_text(aes(label=..count..),stat="count",vjust=-0.5)+
      labs(title="Upset plot for DESeq2 DGEs",x="",y="Count")+
      scale_x_upset()+
      scale_y_continuous(expand = expansion(0.1))+
      scale_fill_material(name="Count",palette = "blue-grey")+
      theme(axis.line = element_line(),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            panel.background=element_rect(fill = "white",colour = NA),
            panel.grid.major = element_line(colour = "grey92"))+
      theme_combmatrix(combmatrix.label.text = element_text(size=10,color="black"),
                       combmatrix.label.extra_spacing=6)
    
    compare_plot[[paste0("Compare_DESeq2@upset1")]] <- p1
    compare_plot[[paste0("Compare_DESeq2@upset1")]] <- p2
    
    }
  
  combn<- combn(names(comparison_list), 2,simplify = F)
  for (n in combn) {
    cat(">>> ",paste(n,collapse = " & "),"\n")
    wb <- createWorkbook()
    ##### edgeR #####
    dge_res <- paste0("edgeR_res@",n)
    
    df1 <- res_list[[dge_res[1]]]
    start<- which(str_detect(colnames(df1),">>Annotation"))[1]
    df1 <- df1[,1:(start-1)]
    df1 <- cbind(paste0(">>",n[1]),df1)
    colnames(df1)[1] <- paste0(">>",n[1])
    
    df2 <- res_list[[dge_res[2]]]
    start<- which(str_detect(colnames(df2),">>Annotation"))[1]
    df2 <- df2[,1:(start-1)]
    df2 <- cbind(paste0(">>",n[2]),df2)
    colnames(df2)[1] <- paste0(">>",n[2])
    
    df<- merge(x=df1,y=df2,by="RowID",all=TRUE)
    df[is.na(df$DifferentialExpression.x),"DifferentialExpression.x"] <- "Droped"
    df[is.na(df$DifferentialExpression.y),"DifferentialExpression.y"] <- "Droped"
    df[,"DifferentialExpression.x"] <- factor(df[,"DifferentialExpression.x"] ,levels = c("Down-regulated","Non-significant","Up-regulated","Droped"))
    df[,"DifferentialExpression.y"] <- factor(df[,"DifferentialExpression.y"] ,levels = c("Down-regulated","Non-significant","Up-regulated","Droped"))
    df <- merge(x=df,y=annotation_matrix,by.x="RowID",by.y="row.names",all.x=T)
    compare_list[[paste0("Compare_edgeR@",n,collapse = "&")]] <- df
    
    addWorksheet(wb, sheetName = "edgeR")
    writeDataTable(wb = wb,
                   x = df,
                   sheet="edgeR",
                   withFilter=F,
                   tableStyle = "TableStyleLight8")
    
    label <- gsub(x = n,pattern = "_vs_",replacement = "/")
    df[,"label"] <- paste0(label[1],":",df[,c("DifferentialExpression.x")],"\n",label[2],":",df[,"DifferentialExpression.y"])
    
    plots<- dges_compare(data = df,col_log2fc=c("logFC.x","logFC.y"),col_label="label",comparsion=label,min_fc=min_fc,tool="edgeR")
    compare_plot[[paste0("Compare_edgeR_logfc@",n,collapse = "&")]] <- plots[[1]]
    compare_plot[[paste0("Compare_edgeR_bubble@",n,collapse = "&")]] <- plots[[2]]
    
    ##### DESeq2 #####
    if (sum(str_detect(string = names(res_list),"DESeq2_res"))>0) {
      dge_res <- paste0("DESeq2_res@",n)
      
      df1 <- res_list[[dge_res[1]]]
      start<- which(str_detect(colnames(df1),">>Annotation"))[1]
      df1 <- df1[,1:(start-1)]
      df1 <- cbind(paste0(">>",n[1]),df1)
      colnames(df1)[1] <- paste0(">>",n[1])
      
      df2 <- res_list[[dge_res[2]]]
      start<- which(str_detect(colnames(df2),">>Annotation"))[1]
      df2 <- df2[,1:(start-1)]
      df2 <- cbind(paste0(">>",n[2]),df2)
      colnames(df2)[1] <- paste0(">>",n[2])
      
      df<- merge(x=df1,y=df2,by="RowID",all=TRUE)
      df[is.na(df$DifferentialExpression.x),"DifferentialExpression.x"] <- "Droped"
      df[is.na(df$DifferentialExpression.y),"DifferentialExpression.y"] <- "Droped"
      df[,"DifferentialExpression.x"] <- factor(df[,"DifferentialExpression.x"] ,levels = c("Down-regulated","Non-significant","Up-regulated","Droped"))
      df[,"DifferentialExpression.y"] <- factor(df[,"DifferentialExpression.y"] ,levels = c("Down-regulated","Non-significant","Up-regulated","Droped"))
      df <- merge(x=df,y=annotation_matrix,by.x="RowID",by.y="row.names",all.x=T)
      compare_list[[paste0("Compare_DESeq2@",n,collapse = "&")]] <- df
      
      addWorksheet(wb, sheetName = "DESeq2")
      writeDataTable(wb = wb,
                     x = df,
                     sheet="DESeq2",
                     withFilter=F,
                     tableStyle = "TableStyleLight8")
      
      label <- gsub(x = n,pattern = "_vs_",replacement = "/")
      df[,"label"] <- paste0(label[1],":",df[,c("DifferentialExpression.x")],"\n",label[2],":",df[,"DifferentialExpression.y"])
      
      plots<- dges_compare(data = df,col_log2fc=c("log2FoldChange.x","log2FoldChange.y"),col_label="label",comparsion=label,min_fc=min_fc,tool="DESeq2")
      compare_plot[[paste0("Compare_DESeq2_logfc@",n,collapse = "&")]] <- plots[[1]]
      compare_plot[[paste0("Compare_DESeq2_bubble@",n,collapse = "&")]] <- plots[[2]]
      
    }
    saveWorkbook(wb, file = paste0("./DGEs_data/Compare@",paste0(n,collapse = "&"),".xlsx"), overwrite = TRUE)
  }
  
  ##### edgeR alluvial plot
  p <- alluvial_plot(res_list = res_list,annotation_matrix = annotation_matrix,tool = "edgeR",stratum_col = "DifferentialExpression")
  compare_plot[[paste0("Compare_edgeR_alluvial@",n,collapse = "&")]] <- p
  
  ##### DEseq2 alluvial plot
  if (sum(str_detect(string = names(res_list),"DESeq2"))>0) {
    p <- alluvial_plot(res_list = res_list,annotation_matrix = annotation_matrix,tool = "DESeq2",stratum_col = "DifferentialExpression")
    compare_plot[[paste0("Compare_DESeq2_alluvial@",n,collapse = "&")]] <- p
  }
  
  # compare_plot_Grob<- suppressWarnings(sapply(compare_plot, ggplotGrob))
  # class(compare_plot_Grob) <- c("arrangelist", class(compare_plot_Grob))
  # suppressMessages(suppressWarnings(
  #   ggsave(filename = "./DGEs_plot/DGEs_multi_compare.pdf", plot = compare_plot_Grob,height = 5,width = 8)
  # ))
  pdf("./DGEs_plot/DGEs_multi_compare.pdf",width = 8,height = 5)
  invisible(lapply(compare_plot, print))
  invisible(dev.off())
  
  saveRDS(object = compare_list,file = "./DGEs_rds/compare_list.rds")
  saveRDS(object = compare_plot,file = "./DGEs_rds/compare_plot.rds")
}



##### check whether the unwanted file exists and remove it ##### 
if (file.exists("Rplots.pdf")) {
  invisible(file.remove("Rplots.pdf"))
}


