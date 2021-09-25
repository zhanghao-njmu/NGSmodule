# prepare reference  ------------------------------------------------------
marker_dir <- "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/markers_data/"
dir.create(marker_dir, recursive = T, showWarnings = FALSE)

SCP_path <- "/home/zhanghao/Program/NGS/UniversalTools/NGSmodule/SingleCellPipe/"
source(paste0(SCP_path, "/SCP-workflow-funtcion.R"))
source(paste0(SCP_path, "/SCP-GeneralSteps/RunSeurat/RunSeurat-function.R"))
source(paste0(SCP_path, "/SCP-plot.R"))


# CellMatch ------------------------------------------------------
temp <- tempfile()
state <- 1
while (state != 0) {
  state <- tryCatch(expr = {
    download.file("https://raw.githubusercontent.com/ZJUFanLab/scCATCH/master/R/sysdata.rda", temp, method = "wget")
    stat <- 0
  }, error = function(error) {
    message(error)
    return(1)
  })
}
load(file = temp)
CellMatch <- CellMatch[!CellMatch$geneSymbol %in% c("", NA), ]
CellMatch <- CellMatch %>%
  mutate(third = NA, second = NA, first = NA, shortname = cellName) %>%
  as.data.frame()
saveRDS(CellMatch, file = paste0(marker_dir,"/CellMatch.rds"))
unlink(temp)


# CellMarker ------------------------------------------------------
temp <- tempfile()
state <- 1
while (state != 0) {
  state <- tryCatch(expr = {
    download.file("http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/all_cell_markers.txt", temp, method = "wget")
    stat <- 0
  }, error = function(error) {
    message(error)
    return(1)
  })
}
CellMarker <- read.table(file = temp, sep = "\t", header = T, fill = TRUE) %>%
  mutate(
    geneSymbol = strsplit(gsub(x = as.character(geneSymbol), pattern = "\\[|\\]| ", replacement = "", perl = TRUE), ","),
    third = NA, second = NA, first = NA, shortname = cellName
  ) %>%
  unnest(c(geneSymbol)) %>%
  as.data.frame()
CellMarker <- CellMarker[!CellMarker$geneSymbol %in% c("", NA), ]
saveRDS(CellMarker,file = paste0(marker_dir,"/CellMarker.rds"))
unlink(temp)


# PanglaoDB ------------------------------------------------------
temp <- tempfile()
state <- 1
while (state != 0) {
  state <- tryCatch(expr = {
    download.file("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz", temp, method = "libcurl")
    stat <- 0
  }, error = function(error) {
    message(error)
    return(1)
  })
}
PanglaoDB <- read.table(file = temp, sep = "\t", header = T, fill = TRUE) %>%
  mutate(
    species = strsplit(species, " "),
    cellType = "Normal cell", cancerType = "Normal", PMID = NA, third = NA, second = NA, first = NA, shortname = cell.type
  ) %>%
  unnest(c(species)) %>%
  as.data.frame()
colnames(PanglaoDB) <- plyr::mapvalues(
  x = colnames(PanglaoDB),
  from = c("species", "official.gene.symbol", "cell.type", "organ"),
  to = c("speciesType", "geneSymbol", "cellName", "tissueType")
)
PanglaoDB[, "speciesType"] <- plyr::mapvalues(
  x = PanglaoDB[, "speciesType"],
  from = c("Hs", "Mm"),
  to = c("Human", "Mouse")
)
PanglaoDB[is.na(PanglaoDB[, "tissueType"]), "tissueType"] <- "Other"

MmGenes <- PanglaoDB[PanglaoDB$speciesType == "Mouse","geneSymbol"]
out <- GeneConvert(
  geneID = MmGenes,
  geneID_from_IDtype = c("ensembl_symbol", "entrez_symbol"),
  geneID_to_IDtype = "ensembl_symbol",
  species_from = "Homo_sapiens",
  species_to = "Mus_musculus"
)
geneID_res <- unique(out$geneID_res[out$geneID_res[, "to_geneID"] != "", c("from_geneID", "to_geneID")])
geneID_res <- geneID_res %>%
  group_by(from_geneID) %>%
  summarise(
    from_geneID = from_geneID,
    to_geneID = paste0(to_geneID[to_geneID != ""], collapse = ";")
  ) %>%
  distinct() %>%
  as.data.frame()
rownames(geneID_res) <- make.unique(geneID_res$from_geneID)
PanglaoDB[PanglaoDB$speciesType == "Mouse", "geneSymbol"] <- geneID_res[PanglaoDB[PanglaoDB$speciesType == "Mouse", "geneSymbol"], "to_geneID"]
PanglaoDB <- PanglaoDB %>%
  mutate(geneSymbol = strsplit(geneSymbol, ";")) %>%
  unnest(c(geneSymbol)) %>%
  distinct() %>%
  as.data.frame()
PanglaoDB <- PanglaoDB[!PanglaoDB$geneSymbol %in% c("", NA), ]
saveRDS(PanglaoDB,file = paste0(marker_dir,"/PanglaoDB.rds"))
unlink(temp)



# tmp ---------------------------------------------------------------------
Ensembl_version <- 102
library(biomaRt)
archives <- listEnsemblArchives(https = FALSE)
url <- archives[which(archives$version == Ensembl_version), "url"]
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = url)
MmGenes <- getBM(
  mart = mart,
  attributes = c("external_gene_name", "mmusculus_homolog_associated_gene_name"),
  filters = c("external_gene_name"),
  values = list(PanglaoDB %>% filter(speciesType == "Mouse") %>% pull("geneSymbol"))
)
MmGenes <- rbind(MmGenes, data.frame(
  external_gene_name = PanglaoDB[, "geneSymbol"],
  mmusculus_homolog_associated_gene_name = ""
))
MmGenes <- MmGenes %>%
  group_by(external_gene_name) %>%
  summarise(
    external_gene_name = external_gene_name,
    mmusculus_homolog_associated_gene_name = paste0(mmusculus_homolog_associated_gene_name[mmusculus_homolog_associated_gene_name != ""], collapse = ";")
  ) %>%
  distinct() %>%
  as.data.frame()
rownames(MmGenes) <- MmGenes[, "external_gene_name"]
PanglaoDB[PanglaoDB$speciesType == "Mouse", "geneSymbol"] <- MmGenes[PanglaoDB[PanglaoDB$speciesType == "Mouse", "geneSymbol"], "mmusculus_homolog_associated_gene_name"]
PanglaoDB <- PanglaoDB %>%
  mutate(geneSymbol = strsplit(geneSymbol, ";")) %>%
  unnest(c(geneSymbol)) %>%
  distinct() %>%
  as.data.frame()
PanglaoDB <- PanglaoDB[!PanglaoDB$geneSymbol %in% c("", NA), ]
saveRDS(PanglaoDB,file = "PanglaoDB.rds")
unlink(temp)
