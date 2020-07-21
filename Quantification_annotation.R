#!/usr/bin/env Rscript

library(refGenome)
library(data.table)


args <- commandArgs( trailingOnly = TRUE )
work_dir <- args[1]
gtfFile <- args[2]
aligner <- args[3]
species <- args[4]
database <- args[5]

######### example 1 ############
# work_dir <- "/data/lab/HeXi/Rpl39l/RNC-seq/NGSpipe_work/"
# gtfFile <- "/data/database/iGenomes/mouse/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf"
# aligner <- "hisat2"
# species <- "mouse"
# database <- "Ensembl"
##############################

gtf = ensemblGenome()
read.gtf(gtf, filename=gtfFile,useBasedir=FALSE)
genetable = getGeneTable(gtf)
genetable <- cbind(">>Annotation.gtf",genetable)
colnames(genetable)[1] <- ">>Annotation.gtf"

if (species!="") {
  library(AnnotationDbi)
  org <- switch(species,"human"="org.Hs.eg.db","mouse"="org.Mm.eg.db","rhesus"="org.Mmu.eg.db","fly"="org.Dm.eg.db")
  library(org,character.only = T)
  idtype <- switch(database,"Ensembl"="ENSEMBL","NCBI"="ENTREZID","UCSC"="SYMBOL")
  keys <- keys(get(org),keytype = idtype)
  columns_select <- c("SYMBOL","ALIAS","GENENAME","ENTREZID","ENSEMBL")
  bioc_anno <- AnnotationDbi::select(get(org), keys=keys, keytype = idtype,columns =columns_select)
  bioc_anno <- aggregate(x = bioc_anno,by=list(bioc_anno[[idtype]]), function(x){paste0(unique(x),collapse = ";")})
  bioc_anno <- bioc_anno[,-1]
  bioc_anno <- cbind(">>Annotation.org_eg_db",bioc_anno)
  colnames(bioc_anno)[1] <- ">>Annotation.org_eg_db"
}


for (type in c("count","rpkm","fpkm","tpm")) {
  files<- list.files(work_dir,recursive = T,full.names = T) %>% grep(x = .,pattern = paste0("Quantification/.*",aligner,".",type,"$"),perl = T,value=T) %>%sort()
  if (length(files)!=0) {
    df<- lapply(1:length(files), function(x){
      read.table(file = files[x],header = T,sep = "\t",stringsAsFactors = F,quote = "",fill = T,comment.char = "",check.names = F)
    })
    input<- Reduce(function(x, y) merge(x, y,by=1, all=TRUE), df)
    input_gtf<- merge(x = input,by.x = "GeneID", y = genetable, by.y = "gene_id",all.x = TRUE)
    if (species!="") {
      output<- merge(x = input_gtf,by.x = "GeneID", y = bioc_anno, by.y = idtype,all.x = TRUE)
    }else{
      output <- input_gtf
    }
    write.table(x = output,file=paste("Quantification",".",aligner,".",type,".tab",sep=""),sep = "\t",row.names = F)
  }else{
    cat("Warning: No .",type," file in ",work_dir,"\n")
    next
  }
}
