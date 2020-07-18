#!/usr/bin/env Rscript
args <- commandArgs( trailingOnly = TRUE )
maindir <- args[1]
version <- args[2]
  
#maindir="/data/database/SortmeRNA/"
#version=98

orgs <- c("Homo_sapiens","Mus_musculus","Macaca_fascicularis","Macaca_mulatta","Drosophila_melanogaster")

#fastq_screen_dir <- "/data/database/FastQ_Screen/FastQ_Screen_Genomes/"

##### Main process #####
setwd(maindir)
version_rmbk <- gsub(pattern = " ",replacement = "",x = version)
library(biomaRt)
library(Biostrings)
# listMarts()
# listEnsembl()
# ensembl_datasets <- listDatasets(useMart("ensembl",version = version))

for (org in orgs) {
  ensenbl_dataset <- switch(org,
                            "Homo_sapiens" = "hsapiens_gene_ensembl",
                            "Mus_musculus" = "mmusculus_gene_ensembl",
                            "Macaca_fascicularis" = "mfascicularis_gene_ensembl",
                            "Macaca_mulatta" = "mmulatta_gene_ensembl",
                            "Drosophila_melanogaster" = "dmelanogaster_gene_ensembl"
  )
  archives<- listEnsemblArchives()
  url <- archives[which(archives$version==version),"url"]
  
  mart <- useMart(biomart = "ensembl",
                  dataset = ensenbl_dataset,
                  host    = url)

  # attributes<- listAttributes(mart)
  # filters <- listFilters(mart)
  
  result <- getBM(mart = mart,
                  filters = "biotype",
                  values = c("rRNA","Mt_rRNA","Mt_tRNA"),
                  attributes = c("gene_biotype",
                                 "external_gene_name",
                                 "ensembl_gene_id",
                                 "ensembl_transcript_id",
                                 "transcript_exon_intron"))
  result[["Header"]] <- paste(result[["gene_biotype"]],
                              result[["external_gene_name"]],
                              result[["ensembl_gene_id"]],
                              result[["ensembl_transcript_id"]],
                              sep = ";"
  )
  result_list <- split.data.frame(x = result ,result[["gene_biotype"]])
  
  lapply(names(result_list), function(x){
    res<- result_list[[x]]
    d<- DNAStringSet(res[["transcript_exon_intron"]])
    names(d) <- res[["Header"]]
    prefix <- paste(x,org,version_rmbk,sep = ".")
    writeXStringSet(d,filepath = paste0(prefix,".fa"),format = "fasta")

    
    ##### fastq_screen
    #    BTdir <- paste0(fastq_screen_dir,"/",x,"_",org)
    #    if (!dir.exists(BTdir)) {
    #      dir.create(BTdir)
    #    }
    #    command2 <- paste0("bowtie2-build --threads 10 ",prefix,".fa ",BTdir,"/",x,"_",org)
    #    system(command2)
  })
}










