#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
maindir <- args[1]
version <- args[2]

maindir <- "/reference/SortmeRNA/"
version <- 101

species <- c("Homo_sapiens", "Mus_musculus", "Macaca_fascicularis", "Macaca_mulatta", "Drosophila_melanogaster")
sequence_type<-c("rRNA", "Mt_rRNA", "Mt_tRNA")
# fastq_screen_dir <- "/data/database/FastQ_Screen/FastQ_Screen_Genomes/"

##### Main process #####
setwd(maindir)
version_rmbk <- gsub(pattern = " ", replacement = "", x = version)
library(biomaRt)
library(Biostrings)

for (s in species) {
  species_split <- unlist(strsplit(s, split = "_"))
  species_dataset <- paste0(tolower(substring(species_split[1], 1, 1)), species_split[2], "_gene_ensembl")
  archives <- listEnsemblArchives()
  url <- archives[which(archives$version == version), "url"]

  mart <- useMart(
    biomart = "ensembl",
    dataset = species_dataset,
    host = url
  )
  attributes <- listAttributes(mart)
  filters <- listFilters(mart)
  filtervalues <- listFilterValues(mart, filter = "biotype")

  result <- getBM(
    mart = mart,
    filters = "biotype",
    values = sequence_type,
    attributes = c(
      "gene_biotype",
      "external_gene_name",
      "ensembl_gene_id",
      "ensembl_transcript_id",
      "transcript_exon_intron"
    )
  )
  result[["Header"]] <- paste(result[["gene_biotype"]],
    result[["external_gene_name"]],
    result[["ensembl_gene_id"]],
    result[["ensembl_transcript_id"]],
    sep = ";"
  )
  result_list <- split.data.frame(x = result, result[["gene_biotype"]])

  invisible(lapply(names(result_list), function(x) {
    res <- result_list[[x]]
    d <- DNAStringSet(res[["transcript_exon_intron"]])
    names(d) <- res[["Header"]]
    prefix <- paste(x, s, version_rmbk, sep = ".")
    writeXStringSet(d, filepath = paste0(prefix, ".fa"), format = "fasta")


    ##### fastq_screen
    #    BTdir <- paste0(fastq_screen_dir,"/",x,"_",s)
    #    if (!dir.exists(BTdir)) {
    #      dir.create(BTdir)
    #    }
    #    command2 <- paste0("bowtie2-build --threads 10 ",prefix,".fa ",BTdir,"/",x,"_",s)
    #    system(command2)
  }))
}
