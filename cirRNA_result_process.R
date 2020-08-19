rm(list = ls())
maindir="F:/Lab/SHZ/circ/circRNA_result/"
setwd(maindir)
# AllFiles <- list.files(path = "./rawdata/",recursive=T)
# AllSamples <- list.dirs(path = "./rawdata/",recursive=F,full.names = F)
# dir.create("./circRNA_result",showWarnings = FALSE)

##### Load library #####
library(stringr)
library(openxlsx)
library(dplyr)
library(VennDiagram)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(gplots)
library(ggplotify)
library(cowplot)
library(data.table)
library(tidyr)
library(GenomicFeatures)
library(genomation)
library(VariantAnnotation)
library(refGenome)
library(edgeR)
library(ggrepel)
library(AnnotationHub)

make.unique.2 = function(x, sep='.'){
  ave(x, x, FUN=function(a){if(length(a) > 1){paste(a, 1:length(a), sep=sep)} else {a}})
}

##### circ table #####
sampleinfo <- read.xlsx("./sampleMeta.xlsx")
rownames(sampleinfo) <- sampleinfo[["Sample.Name"]]
sampleinfo[,"newname"] <- make.unique.2(sampleinfo[["Group.Name"]])
normal <- sampleinfo[which(sampleinfo[["Group.Name"]]=="Normal"),"newname"]
degeneration <- sampleinfo[which(sampleinfo[["Group.Name"]]=="Degeneration"),"newname"]


##### CIRI2(Standard pipeline) #####
CiriStandard_flies <- AllFiles %>% grep(x = .,pattern = ".*CIRI_trim/.*.ciri$",perl = T,value=T)
CiriStandard_sample <- gsub(x = CiriStandard_flies,pattern = "\\/.*",replacement = "",perl = T)
CiriStandard_df<- lapply(1:length(CiriStandard_flies), function(x){
  file <- paste0("./rawdata/",CiriStandard_flies[x])
  read.table(file = file,header = T,sep = "\t",stringsAsFactors = F,quote = "",fill = T,comment.char = "")
})
names(CiriStandard_df) <- CiriStandard_sample
### CIRI criterion and preprocessing
CiriStandard_df <- lapply(names(CiriStandard_df), function(name){
  x <- CiriStandard_df[[name]]
  x <- x[which(x[["X.junction_reads"]]>=2),]
  x <- x[,c("chr","circRNA_start","circRNA_end","strand","X.junction_reads")]
  colnames(x) <- c("chrom","start","end","strand","n_reads")
  x[,c("start","end","n_reads")] <- lapply(x[,c("start","end","n_reads")], as.numeric)
  x[,"CircID"] <- paste(x[["chrom"]],x[["start"]],x[["end"]],x[["strand"]],sep = ":")
  x[,"Tool"] <- "CIRI2"
  x[,"Sample"] <- name
  return(x)
})
CiriStandard_df <- rbindlist(CiriStandard_df)


##### CircExplorer2 #####
CircExplorer2_flies <- AllFiles %>% grep(x = .,pattern = ".*CIRCexplorer2/annotate/circularRNA_known.txt$",perl = T,value=T)
CircExplorer2_sample <- gsub(x = CircExplorer2_flies,pattern = "\\/.*",replacement = "",perl = T)
CircExplorer2_df<- lapply(1:length(CircExplorer2_flies), function(x){
  file <- paste0("./rawdata/",CircExplorer2_flies[x])
  read.table(file = file,header = F,col.names = c("chrom","start","end","CircularRNA/JunctionReads","score","strand","thickStart","thickEnd","itemRgb",
                                                  "exonCount","exonSizes","exonOffsets","readNumber","circType","geneName","isoformName","index","flankIntron"),
             sep = "\t",stringsAsFactors = F,quote = "",fill = T,comment.char = "")
})
names(CircExplorer2_df) <- CircExplorer2_sample
### CircExplorer2 criterion and preprocessing
CircExplorer2_df <- lapply(names(CircExplorer2_df), function(name){
  x <- CircExplorer2_df[[name]]
  x <- x[which(x[["readNumber"]]>=2),]
  x <- x[,c("chrom","start","end","strand","readNumber")]
  colnames(x) <- c("chrom","start","end","strand","n_reads")
  x[,c("start","end","n_reads")] <- lapply(x[,c("start","end","n_reads")], as.numeric)
  x[["start"]] <- x[["start"]]+1
  x[,"CircID"] <- paste(x[["chrom"]],x[["start"]],x[["end"]],x[["strand"]],sep = ":")
  x[,"Tool"] <- "CircExplorer2"
  x[,"Sample"] <- name
  return(x)
})
CircExplorer2_df <- rbindlist(CircExplorer2_df)


##### find_circ #####
FindCirc_flies <- AllFiles %>% grep(x = .,pattern = ".*find_circ/circ_candidates.bed$",perl = T,value=T)
FindCirc_sample <- gsub(x = FindCirc_flies,pattern = "\\/.*",replacement = "",perl = T)
FindCirc_df<- lapply(1:length(FindCirc_flies), function(x){
  file <- paste0("./rawdata/",FindCirc_flies[x])
  read.table(file = file,header = F,col.names = c("chrom","start","end","name","n_reads","strand","n_uniq","uniq_bridges","best_qual_left","best_qual_right",
                                                  "tissues","tiss_counts","edits","anchor_overlap","breakpoints","signal","strandmatch","category"),
             sep = "\t",stringsAsFactors = F,quote = "",fill = T,comment.char = "")
})
names(FindCirc_df) <- FindCirc_sample
### find_circ criterion and preprocessing
FindCirc_df <- lapply(names(FindCirc_df), function(name){
  x <- FindCirc_df[[name]]
  x <- x[which(x[["n_reads"]]>=2),]
  x <- x[,c("chrom","start","end","strand","n_reads")]
  x[,c("start","end","n_reads")] <- lapply(x[,c("start","end","n_reads")], as.numeric)
  x[["start"]] <- x[["start"]]+1
  x[,"CircID"] <- paste(x[["chrom"]],x[["start"]],x[["end"]],x[["strand"]],sep = ":")
  x[,"Tool"] <- "find_circ"
  x[,"Sample"] <- name
  return(x)
})
FindCirc_df <- rbindlist(FindCirc_df)


##### KNIFE #####
KNIFE_flies <- AllFiles %>% grep(x = .,pattern = ".*KNIFE/output/circReads/combinedReports/.*circJuncProbs.txt$",perl = T,value=T)
KNIFE_sample <- gsub(x = KNIFE_flies,pattern = "\\/.*",replacement = "",perl = T)
KNIFE_df<- lapply(1:length(KNIFE_flies), function(x){
  file <- paste0("./rawdata/",KNIFE_flies[x])
  read.table(file = file,header = T,
             sep = "\t",stringsAsFactors = F,quote = "",fill = T,comment.char = "")
})
names(KNIFE_df) <- KNIFE_sample
### KNIFE criterion and preprocessing
KNIFE_df <- lapply(names(KNIFE_df), function(name){
  x <- KNIFE_df[[name]]
  x <- x[which(x[["orig_count"]]>=2 & x[["orig_posterior"]]>=0.9),]
  b<- lapply(x[["junction"]], function(x){str_split(string = x,pattern = "\\||:")%>% unlist()})%>%unlist()%>%matrix(ncol = 7,byrow = T) %>%as.data.frame(stringsAsFactors=F)
  x[,c("chrom","gene_start","start","gene_end","end","circ_type","strand")] <- b
  x <- x[,c("chrom","start","end","strand","orig_count")]
  x<-  x %>% mutate(start2=if_else(as.numeric(start)>as.numeric(end),end,start),
                                       end2=if_else(as.numeric(start)>as.numeric(end),start,end))
  x <- x[,c("chrom","start2","end2","strand","orig_count")]
  colnames(x) <- c("chrom","start","end","strand","n_reads")
  x[,c("start","end","n_reads")] <- lapply(x[,c("start","end","n_reads")], as.numeric)
  x[,"CircID"] <- paste(x[["chrom"]],x[["start"]],x[["end"]],x[["strand"]],sep = ":")
  x[,"Tool"] <- "KNIFE"
  x[,"Sample"] <- name
  return(x)
})
KNIFE_df <- rbindlist(KNIFE_df)


##### CircMaker #####
CircMaker_flies <- AllFiles %>% grep(x = .,pattern = ".*CircMaker/Detection_Result/Brief_sum.txt$",perl = T,value=T)
CircMaker_sample <- gsub(x = CircMaker_flies,pattern = "\\/.*",replacement = "",perl = T)
CircMaker_df<- lapply(1:length(CircMaker_flies), function(x){
  file <- paste0("./rawdata/",CircMaker_flies[x])
  read.table(file = file,header = F,col.names = c("chrom","start","end","n_reads","strand","circ_type"),
             sep = " ",stringsAsFactors = F,quote = "",fill = T,comment.char = "")
})
names(CircMaker_df) <- CircMaker_sample
### CircMaker criterion and preprocessing
CircMaker_df <- lapply(names(CircMaker_df), function(name){
  x <- CircMaker_df[[name]]
  x <- x[which(x[["n_reads"]]>=2),]
  x <- x[,c("chrom","start","end","strand","n_reads")]
  x<-  x %>% mutate(start2=if_else(as.numeric(start)>as.numeric(end),end,start),
                    end2=if_else(as.numeric(start)>as.numeric(end),start,end))
  x <- x[,c("chrom","start2","end2","strand","n_reads")]
  colnames(x) <- c("chrom","start","end","strand","n_reads")
  x[,c("start","end","n_reads")] <- lapply(x[,c("start","end","n_reads")], as.numeric)
  x[,"CircID"] <- paste(x[["chrom"]],x[["start"]],x[["end"]],x[["strand"]],sep = ":")
  x[,"Tool"] <- "CircMaker"
  x[,"Sample"] <- name
  return(x)
})
CircMaker_df <- rbindlist(CircMaker_df)


##### all data #####
Circ_df<- rbind(CiriStandard_df,CircExplorer2_df,FindCirc_df,KNIFE_df,CircMaker_df)%>%as.data.frame(stringsAsFactors=FALSE)
Circ_df[,c("Sample","Group")] <- sampleinfo[Circ_df[["Sample"]],c("newname","Group.Name")]
#saveRDS(Circ_df,"./Circ_df.rds")

Circ_df <- readRDS("./Circ_df.rds")
##### Venn plot for the results #####
Circ_list <- list()
ggList <- list()

for(sample in sampleinfo[["newname"]]){
  venn_list = list("CIRI2"= Circ_df[which(Circ_df[["Sample"]]==sample & Circ_df[["Tool"]]=="CIRI2"),"CircID"],
                   "CircExplorer2"= Circ_df[which(Circ_df[["Sample"]]==sample & Circ_df[["Tool"]]=="CircExplorer2"),"CircID"],
                   "find_circ"= Circ_df[which(Circ_df[["Sample"]]==sample & Circ_df[["Tool"]]=="find_circ"),"CircID"],
                   "KNIFE"= Circ_df[which(Circ_df[["Sample"]]==sample & Circ_df[["Tool"]]=="KNIFE"),"CircID"],
                   "CircMaker"= Circ_df[which(Circ_df[["Sample"]]==sample & Circ_df[["Tool"]]=="CircMaker"),"CircID"])
  overlap <- calculate.overlap(venn_list)
  venn.plot <- venn.diagram(x= venn_list,filename = NULL,
                            height = 1500,width = 1800, resolution = 500, 
                            category.names = lapply(venn_list, length) %>% unlist() %>% paste0(names(.),"\n(",.,")"),
                            scaled=T,force.unique=T,print.mode="raw",
                            lwd=2,lty=1,
                            col= pal_nejm(palette = "default")(8)[1:5],
                            fill= alpha(pal_nejm(palette = "default")(8)[1:5],0.3),
                            cex=0.5,
                            fontface = "bold",
                            fontfamily = "sans",
                            cat.cex = 0.8,
                            cat.default.pos = "outer",
                            cat.dist = c(0.26, 0.35, 0.26,0.29,0.30),
                            cat.fontfamily = "sans",
                            cat.col = pal_nejm(palette = "default")(8)[1:5],
                            margin=0.3
                            )
  venn.plot <- grobTree(venn.plot)
  venn.plot <- addGrob(gTree = venn.plot,child = textGrob(label = paste("Sample:",sample), x = unit(0.5, "npc"), y = unit(0.9, "npc")))
  venn.plot <- as.ggplot(venn.plot)+theme(aspect.ratio = 1)
  ggList[[sample]] <- venn.plot

  venn.output <- venn(data=venn_list,show.plot = F)
  intersections<- attributes(venn.output)$intersections
  intersections.df<- stack(intersections)
  intersections.df <- cbind(intersections.df,sample,stringsAsFactors=FALSE)
  colnames(intersections.df) <- c("CircID","Tool","Sample")
  intersections.df[,"nTools"] <- str_count(intersections.df[["Tool"]],":")+1
  intersections.df <- intersections.df %>% mutate(Tool = strsplit(as.character(Tool), ":")) %>% unnest(Tool)
  intersections.df <- merge(x = intersections.df,by.x = c("CircID","Tool","Sample"),
                            y=Circ_df[,c("CircID","Tool","Sample","Group","n_reads")],by.y = c("CircID","Tool","Sample")) %>% as.data.frame(stringsAsFactors=FALSE)
  intersections.df <- aggregate(cbind(Tool,n_reads)~.,data=intersections.df,FUN = function(x){paste0(x,collapse = ":")})
  Circ_list[[sample]] <- intersections.df
}

##### save circ venn plot #####
p<- plot_grid(plotlist = ggList)
num_plots <- length(ggList)
cols <- ceiling(sqrt(num_plots))
rows <- ceiling(num_plots/cols)
ggsave(filename = "./circ_perSample.pdf",plot = p,width = cols*3.5,height = rows*3.5)

##### make Circ_table #####
Circ_table <- rbindlist(Circ_list)

Circ_table[["Sample"]] <- factor(Circ_table[["Sample"]],levels = c(normal,degeneration))
Circ_table[["Group"]] <- factor(Circ_table[["Group"]],levels = c("Normal","Degeneration"))
Circ_table<- Circ_table %>%group_by(CircID,Sample)%>%
  mutate(meanReads=n_reads %>% str_split(":") %>% unlist() %>% as.numeric() %>% mean() %>% round(digits = 2),
         width=CircID %>% str_split(":") %>% unlist() %>% as.numeric() %>% (function(x)(x[3]-x[2])))
Circ_table <- Circ_table %>% group_by(CircID) %>%
  mutate(Sample_combine=paste(as.character(Sample)[order(match(as.character(Sample),c(normal,degeneration)))],collapse = ":"),
         Group_combine= paste(as.character(Group)[order(match(as.character(Group),c("Normal","Degeneration")))],collapse = ":"),
         nNormal=sum(Group=="Normal"),
         nDegeneration=sum(Group=="Degeneration"),
         nTools_combine=paste(nTools[order(match(as.character(Sample),c(normal,degeneration)))],collapse = ":"),
         meanReads_combine=paste(meanReads[order(match(as.character(Sample),c(normal,degeneration)))],collapse = ":")
  ) %>%                                                         
  as.data.frame(stringsAsFactors=FALSE)


write.xlsx(x = Circ_table,file =  "./Circ_table.xlsx",
  asTable = T,colWidths="auto",tableStyle = "TableStyleLight8",firstRow=T)
saveRDS(Circ_table,file = "./Circ_table.rds")

##### circRNA filtering #####
Circ_table <- readRDS("Circ_table.rds")

p1 <- ggplot(data = Circ_table,aes(x=nTools))+
  geom_bar(aes(fill=..count..),color="black",alpha=0.7)+
  geom_line(aes(y=..count..),stat = "count",size=1)+
  geom_point(aes(y=..count..),stat = "count",shape=21,fill="black",color="white",size=4)+
  geom_text(aes(label=..count..),vjust=1.8,stat='count',color="white",size=3.5)+
  labs(title = "CircRNA detected filtering by nTools",y="CircRNA number")+
  scale_fill_gradient2(name="Count")+
  theme_classic2()+
  theme(aspect.ratio=0.8,
        panel.border=element_rect(color = "black",size = 1,fill = "transparent"),
        panel.grid.major=element_line(size=0.5,colour = "grey80",linetype = 2),
        line=element_line(size = 1),
        axis.line = element_blank(),
        axis.text=element_text(colour = "black",size = 10),
        axis.title=element_text(colour = "black",size = 12),
        legend.title=element_text(colour = "black",size = 10),
        legend.position = "right")

br <- sort(unique(c(0:6,log2(3))))
br_label <- br
br_label[which(br_label==(log2(3)))] <- "log2(3)"
col_x <- ifelse(br==(log2(3)),"red","black")

p2 <- ggplot(data = Circ_table,aes(x=log2(meanReads),fill = ..x..))+
  geom_histogram(color="black",alpha=0.7)+
  geom_vline(xintercept = log2(3),color="red3",linetype=2)+
  labs(title = "CircRNA detected filtering by nTools",y="CircRNA number")+
  scale_fill_gradient2(name="log2(meanReads)")+
  scale_x_continuous(limits = c(0,6),breaks = br,labels = br_label)+
  theme_classic2()+
  theme(aspect.ratio=0.8,
        panel.border=element_rect(color = "black",size = 1,fill = "transparent"),
        panel.grid.major=element_line(size=0.5,colour = "grey80",linetype = 2),
        line=element_line(size = 1),
        axis.line = element_blank(),
        axis.text=element_text(colour = "black",size = 10),
        axis.title=element_text(colour = "black",size = 12),
        axis.text.x = element_text(colour = col_x,angle = 45,hjust = 1),
        legend.title=element_text(colour = "black",size = 10),
        legend.position = "right")


p <- cowplot::plot_grid(plotlist = list(p1,p2))
ggsave(plot = p,filename = "./Circ_ntools_reads.bar.pdf",width = 10,height = 4)

Circ_table <- Circ_table %>% group_by(CircID) %>% filter(max(nTools)>=2 & max(meanReads)>=3) %>% as.data.frame(stringsAsFactors = F)

##### custom analysis #####
p1 <- ggplot(data = Circ_table,aes(x=Sample,fill=Group))+
  geom_bar(color="black")+
  geom_text(aes(label=..count..),vjust=1.2,stat='count',color="white",size=3)+
  labs(title = paste0("Number of circRNAs in the all samples(",length(unique(Circ_table[["CircID"]])),")"),y="Count by sample")+
  scale_fill_nejm()+
  theme_classic2()+
  theme(aspect.ratio=0.8,
        panel.border=element_rect(color = "black",size = 1,fill = "transparent"),
        panel.grid.major=element_line(size=0.5,colour = "grey80",linetype = 2),
        line=element_line(size = 1),
        axis.line = element_blank(),
        axis.text=element_text(colour = "black",size = 10),
        axis.title=element_text(colour = "black",size = 12),
        axis.text.x = element_text(colour = "black",size = 10,angle = 45,hjust = 1),
        legend.title=element_text(colour = "black",size = 10),
        legend.position = c(0.8,0.85),
        legend.background = element_rect(fill=alpha(colour = "grey90",alpha = 0.3),color = "black"))

dat<- Circ_table%>% group_by(Sample) %>% summarise(Count=n(),Group=unique(Group))
p2 <- ggplot(data = dat,aes(x=Group,y=Count,fill=Group))+
  geom_boxplot(width=0.5,size=0.8,show.legend = F)+
  stat_summary(fun.y=median, colour="white", geom="text", show.legend = FALSE, 
               vjust=-0.5, aes( label=round(..y.., digits=1)))+
  labs(title = "Number of circRNAs in the two groups",y="Count by group")+
  stat_compare_means(method = "t.test",label.x.npc=0.26,color="red3",size=5)+
  scale_fill_nejm()+
  theme_classic2()+
  theme(aspect.ratio=0.8,
        panel.border=element_rect(color = "black",size = 1,fill = "transparent"),
        panel.grid.major=element_line(size=0.5,colour = "grey80",linetype = 2),
        line=element_line(size = 1),
        axis.line = element_blank(),
        axis.text=element_text(colour = "black",size = 10),
        axis.title=element_text(colour = "black",size = 12),
        axis.text.x = element_text(colour = "black",size = 10,angle = 45,hjust = 1),
        legend.title=element_text(colour = "black",size = 10),
        legend.position = c(0.8,0.85),
        legend.background = element_rect(fill=alpha(colour = "grey90",alpha = 0.3),color = "black"))

p3 <- ggplot(data = Circ_table,aes(x=Sample,y=log2(meanReads),fill=Group))+
  geom_boxplot()+
  stat_summary(fun.y=median, colour="white", geom="text", show.legend = FALSE, 
               vjust=-0.5, aes( label=round(..y.., digits=1)))+
  labs(title = "Value of mean supported reads per circRNA across the samples",y="log2(MeanReads)")+
  scale_fill_nejm()+
  theme_classic2()+
  theme(aspect.ratio=0.8,
        panel.border=element_rect(color = "black",size = 1,fill = "transparent"),
        panel.grid.major=element_line(size=0.5,colour = "grey80",linetype = 2),
        line=element_line(size = 1),
        axis.line = element_blank(),
        axis.text=element_text(colour = "black",size = 10),
        axis.title=element_text(colour = "black",size = 12),
        axis.text.x = element_text(colour = "black",size = 10,angle = 45,hjust = 1),
        legend.title=element_text(colour = "black",size = 10),
        legend.position = c(0.8,0.85),
        legend.background = element_rect(fill=alpha(colour = "grey90",alpha = 0.3),color = "black"))
p4 <- ggplot(data = Circ_table,aes(x=Group,y=log2(meanReads),fill=Group))+
  geom_boxplot(width=0.5,size=0.8,show.legend = F)+
  stat_summary(fun.y=median, colour="white", geom="text", show.legend = FALSE, 
               vjust=-0.5, aes( label=round(..y.., digits=1)))+
  labs(title = "Value of mean supported reads samples in the two groups",y="log2(MeanReads)")+
  stat_compare_means(method = "t.test",label.x.npc=0.26,color="red3",size=5)+
  scale_fill_nejm()+
  theme_classic2()+
  theme(aspect.ratio=0.8,
        panel.border=element_rect(color = "black",size = 1,fill = "transparent"),
        panel.grid.major=element_line(size=0.5,colour = "grey80",linetype = 2),
        line=element_line(size = 1),
        axis.line = element_blank(),
        axis.text=element_text(colour = "black",size = 10),
        axis.title=element_text(colour = "black",size = 12),
        axis.text.x = element_text(colour = "black",size = 10,angle = 45,hjust = 1),
        legend.title=element_text(colour = "black",size = 10),
        legend.position = c(0.8,0.85),
        legend.background = element_rect(fill=alpha(colour = "grey90",alpha = 0.3),color = "black"))

p <- cowplot::plot_grid(plotlist = list(p1,p2,p3,p4),nrow = 2)
ggsave(p,filename = "./Circ_stat.barbox.pdf",width = 12,height = 10)


##### DE test ##### 
Circ_table_melt <- dcast(Circ_table,CircID+width+Sample_combine+Group_combine+nNormal+nDegeneration+nTools_combine+meanReads_combine~Sample,value.var="meanReads")
Circ_table_melt[is.na(Circ_table_melt)] <- 0
sampleinfo[["Group.Name"]] <- factor(x = sampleinfo[["Group.Name"]],levels = c("Normal","Degeneration"))

dge <- edgeR::DGEList(counts = Circ_table_melt[,sampleinfo[["newname"]]], group = sampleinfo[["Group.Name"]])
dge <- edgeR::calcNormFactors(dge, method = "none",na.rm = TRUE)
dge <- edgeR::estimateCommonDisp(dge)
dge <- edgeR::estimateTagwiseDisp(dge)
statistics <- edgeR::topTags(exactTest(dge),n = nrow(dge$counts),adjust.method = "BH", sort.by = "none")
edgerRes <- dplyr::bind_cols(data.frame(Circ_table_melt$CircID,statistics$table, dge$counts)) %>% 
  dplyr::select(-.data$logCPM) %>% 
  dplyr::rename(log2FC = .data$logFC, 
                pvalue = .data$PValue, 
                p.adjust = .data$FDR, 
                CircID = .data$Circ_table_melt.CircID) %>% 
  dplyr::mutate(CircID = as.character(.data$CircID))
edgerRes[,"Mean(Normal)"] <- apply(edgerRes[,normal], 1,mean)
edgerRes[,"Mean(Degeneration)"] <- apply(edgerRes[,degeneration], 1,mean)

##### paramaters
project_name = "Differentially expressed circRNA"
group_num <- 2
group1 <- normal ### control group
group2 <- degeneration ### experiment group
fc_cutoff <- 4
paired_group <- FALSE
adjust_method <- "BH"
stat_method <- "p.adjust" ### p-value or p.adjust
stat_cutoff <- 0.01
feature_label <- ""
feature_column <- ""
theme_use <- NULL
main_title=""
facet_x_title="Degeneration vs Normal" #Pachytene spermatocyte / Elongated spermatids
facet_y_title=""

diff_method <- function(dat_row) {
  if (dat_row["log2FC"] < log2(1 / fc_cutoff) & dat_row[stat_method] < stat_cutoff) {
    "Down-regulated"
  } else if (dat_row["log2FC"] > log2(fc_cutoff) &dat_row[stat_method] < stat_cutoff) {
    "Up-regulated"
  } else{
    "No significance"
  }
}
discrete_level <- c("Up-regulated","Down-regulated","No significance")

edgerRes[,"MaxMean"] <- apply(edgerRes[,c("Mean(Normal)","Mean(Degeneration)")], 1,max)
edgerRes[,"Diff_group"] <- apply(edgerRes[,c("log2FC",stat_method)], 1, diff_method)
edgerRes[,"Diff_group"] <- factor(edgerRes[,"Diff_group"] ,levels = discrete_level)
(dat_table<- table(edgerRes[,"Diff_group"]))
saveRDS(object = edgerRes,"./edgerRes.rds")

edgerRes <- readRDS("./edgerRes.rds")
dat <- edgerRes

maxFC <-
  edgerRes %>% filter(.data[["log2FC"]] != Inf &
                   .data[["log2FC"]] != (-Inf)) %>% 
  pull("log2FC") %>% 
  abs() %>% 
  max() %>% 
  ceiling()

theme_use<- pal_nejm(palette = "default")(8)
summary_label <- c(paste0("Up(",dat_table[1],")"),
                   paste0("Down(",dat_table[2],")"),
                   paste0("No significance(",dat_table[3],")"))
dat[["summary_label"]] <- NA
dat[["summary_label"]][1] <- summary_label

dat[["facet_x"]] <- facet_x_title
dat[["facet_y"]] <- facet_y_title
dat[["feature_label"]] <- NA
dat[["feature_label"]][which(dat[[feature_column]]%in%feature_label)] <- 
  dat[[feature_column]][which(dat[[feature_column]]%in%feature_label)]

br <- sort(unique(c(0:5,-log10(stat_cutoff))))
br_label <- br
br_label[which(br_label==(-log10(stat_cutoff)))] <- paste0(stat_method,"=",stat_cutoff)
col_y <- ifelse(br==(-log10(stat_cutoff)),"red","black")

p1<- ggplot(data = dat,aes(x=dat[,"log2FC"],y=-log10(dat[,stat_method]),group=dat[,"Diff_group"]))+
  geom_point(aes(fill=dat[,"Diff_group"],alpha=dat[,"Diff_group"],size=log2(dat[,"MaxMean"])),shape=21,position = "jitter")+
  geom_vline(xintercept = log2(fc_cutoff),linetype="longdash",color=theme_use[1],size=1,show.legend = F)+
  geom_vline(xintercept = -log2(fc_cutoff),linetype="longdash",color=theme_use[2],size=1,show.legend = F)+
  geom_hline(yintercept = -log10(stat_cutoff),linetype="longdash",color="black",size=1,show.legend = F)+
  geom_rug(aes(color=dat[,"Diff_group"]),show.legend = F)+
  geom_label_repel(aes(label=feature_label),min.segment.length = 0)+
  scale_x_continuous(breaks = seq(-maxFC,maxFC,1),limits = c(-maxFC,maxFC))+
  scale_y_continuous(breaks = br,limits = c(0,5),labels = br_label)+
  scale_alpha_manual(breaks=discrete_level,values=c(1,1,0.5),guide=F)+
  scale_color_manual(breaks=discrete_level,values = c(theme_use[1],theme_use[2],"grey80"))+
  scale_fill_manual(breaks=discrete_level,values = c(theme_use[1],theme_use[2],"grey80"),labels=summary_label)+
  scale_size_continuous(range = c(0.3,4))+
  labs(x="Log2(Fold Change)",y=paste0("-Log10(",stat_method,")"),title = main_title)+
  facet_grid(.~facet_x)+
  guides(fill=guide_legend(title = "Groups",override.aes = list(size = 3),order = 1),
         size=guide_legend(title = "Log2(Max Expression)",override.aes = list(color="black",shape=16)))+
  theme_classic2()+
  theme(aspect.ratio=0.9,
        panel.border=element_rect(color = "black",size = 1,fill = "transparent"),
        panel.grid.major=element_line(size=0.5,colour = "grey80",linetype = 2),
        line=element_line(size = 1),
        axis.line = element_blank(),
        axis.text=element_text(colour = "black",size = 10),
        axis.title=element_text(colour = "black",size = 12),
        axis.text.y = element_text(colour = col_y),
        legend.title=element_text(colour = "black",size = 10))
ggsave(p1,filename = paste0("./Volcano_1.pdf"),width = 7,height =7,limitsize = F)



##### filter circRNA by group #####
venn_list = list("Normal_highLevel"= unique(Circ_table[which(Circ_table[["nNormal"]]>=5),"CircID"]),
                 "Degeneration_highLevel"=  unique(Circ_table[which(Circ_table[["nDegeneration"]]>=5),"CircID"]),
                 "Normal_lowLevel"= unique(Circ_table[which(Circ_table[["nNormal"]]>=1 & Circ_table[["nNormal"]]<5),"CircID"]),
                 "Degeneration_lowLevel"=  unique(Circ_table[which(Circ_table[["nDegeneration"]]>=1 & Circ_table[["nDegeneration"]]<5),"CircID"])
)
overlap <- calculate.overlap(venn_list)
venn.plot <- venn.diagram(x= venn_list,filename = NULL,
                          height = 1500,width = 1800, resolution = 500, 
                          category.names = lapply(venn_list, length) %>% unlist() %>% paste0(names(.),"\n(",.,")"),
                          scaled=T,force.unique=T,print.mode="raw",
                          lwd=2,lty=1,
                          col= pal_nejm(palette = "default")(8)[1:4],
                          fill= alpha(pal_nejm(palette = "default")(8)[1:4],0.3),
                          cex=0.8,
                          fontface = "bold",
                          fontfamily = "sans",
                          cat.cex = 0.8,
                          cat.default.pos = "outer",
                          cat.dist = c(0.3, 0.2,0.2,0.3),
                          cat.fontfamily = "sans",
                          cat.col = pal_nejm(palette = "default")(8)[1:4],
                          margin=0.2
)
venn.plot <- grobTree(venn.plot)
venn.plot <- addGrob(gTree = venn.plot,child = textGrob(label = paste0("CircRNA Number (",length(unique(unlist(venn_list))),")"), x = unit(0.5, "npc"), y = unit(0.9, "npc")))
venn.plot <- as.ggplot(venn.plot)+theme(aspect.ratio = 1)
ggsave(plot = venn.plot,filename = "./Normal_vs_Degeneration.pdf",width = 5,height = 3)

venn.output <- venn(data=venn_list,show.plot = F)
intersections<- attributes(venn.output)$intersections
intersections.df<- stack(intersections)
colnames(intersections.df) <- c("CircID","Intersections")
intersections.df <- merge(x=intersections.df,by.x = "CircID",y = unique(Circ_table[,c("CircID","width","Sample_combine","Group_combine","nNormal","nDegeneration","nTools_combine","meanReads_combine")]),by.y = "CircID",all.x = T)


##### opt1. make annotation for the circ #####
TxDb<- makeTxDbFromGFF(file = "./Homo_sapiens.GRCh38.94.gtf",format = "gtf")
gtf = ensemblGenome()
read.gtf(gtf, filename="./Homo_sapiens.GRCh38.94.gtf",useBasedir=FALSE)
genes = gtf@ev$genes[ ,c("gene_id","gene_name","start","end","strand","gene_biotype")]


dat <- matrix(unlist(str_split(unique(intersections.df[["CircID"]]),pattern = ":")),ncol = 4,byrow = T) %>% as.data.frame(stringsAsFactors=F)
colnames(dat) <- c("Chr","Start","End","Strand")
dat_gr <- makeGRangesFromDataFrame(dat)
act<- isActiveSeq(TxDb)
act[which(!names(act)%in%unique(dat[["Chr"]]))] <- FALSE
isActiveSeq(TxDb) <- act
TxDb_df <- as.data.frame(transcripts(TxDb))

dat_loc_gr <- locateVariants(dat_gr, TxDb, AllVariants(promoter = PromoterVariants(),
                                                       intergenic = IntergenicVariants(upstream = 1e+09L, downstream = 1e+09L)))

dat_loc_df <- as.data.frame(dat_loc_gr,stringsAsFactors=F)
dat_loc_df[,"CircID"] <- apply(dat_loc_df, 1, function(x){
  paste(as.character(x[c(1:3,5)]),collapse = ":")
})
dat_loc_df <- merge(x=dat_loc_df,by.x="TXID",y=TxDb_df[,c("tx_id","tx_name")],by.y="tx_id",all.x=T)
dat_loc_df<- merge(x = dat_loc_df,by.x = "GENEID", y = genes, by.y = "gene_id",all.x = TRUE)
dat_loc_df <- unique(dat_loc_df[,c("CircID","width","LOCATION","tx_name","GENEID","gene_name","gene_biotype")])
dat_loc_df<- dat_loc_df %>% group_by(CircID) %>%
  mutate(tx_name=paste(unique(tx_name),collapse = ":"),
         LOCATION=paste(unique(LOCATION),collapse = ":"),
         GENEID=paste(unique(GENEID),collapse = ":"),
         gene_name=paste(unique(gene_name),collapse = ":"),
         gene_biotype=paste(unique(gene_biotype),collapse = ":"))%>%
  unique()

nrow(dat_loc_df)
length(dat_loc_df[["CircID"]])

intersections.df <- merge(x = intersections.df,by.x="CircID",y = dat_loc_df,by.y="CircID",all.x=T)
write.xlsx(x = intersections.df,file =  "./Circ_intersection.xlsx",
           asTable = T,colWidths="auto",tableStyle = "TableStyleLight8",firstRow=T)


##### opt2. make annotation for the circ #####                                        
dat <- matrix(unlist(str_split(unique(intersections.df[["CircID"]]),pattern = ":")),ncol = 4,byrow = T) %>% as.data.frame(stringsAsFactors=F)
colnames(dat) <- c("Chr","Start","End","Strand")
dat_gr <- makeGRangesFromDataFrame(dat)

genes_gr = gffToGRanges("./Homo_sapiens.GRCh38.94.gtf")
names(genes_gr) <- 
# what regions overlap what genes?
overlapGenes <- findOverlaps(dat_gr, genes_gr)

# Return any genes with an overlap.
# Convert the resulting "Hits" object to a data frame
# and use it as an index
overlapGenes.df <- as.data.frame(overlapGenes)
a <- as.data.frame(dat_gr[overlapGenes.df$queryHits])
b <- as.data.frame(genes_gr[overlapGenes.df$subjectHits])[,c("gene_id","gene_name")]
res <- cbind(a,b) %>% unique()

# extract the regions that hit genes
regionsWithHits <- as.data.frame(dat_gr[overlapGenes.df$queryHits])
# add the names of the genes
regionsWithHits$genes <- genes_gr[overlapGenes.df$subjectHits]$gene_name

# Some segments might not overlap any gene.
# Find the distance to the nearest gene
distanceToNearest(dat_gr, genes_gr)


## 

res[,"CircID"] <- apply(res, 1, function(x){
  paste(trimws(as.character(x[c(1:3,5)])),collapse = ":")
})

length(which(intersections.df[["CircID"]] %in% res[,"CircID"]))

intersections.df <- merge(x = intersections.df,by.x="CircID",y=res[,c("CircID","gene_id","gene_name")],by.y="CircID",all.x=T)
write.xlsx(x = intersections.df,file =  "./Circ_intersection.xlsx",
           asTable = T,colWidths="auto",tableStyle = "TableStyleLight8",firstRow=T)







