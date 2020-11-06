# prepare reference  ------------------------------------------------------
# CellMatch
temp <- tempfile()
state <- 1
while (state != 0) {
  state <- tryCatch(expr = {
    download.file("https://raw.githubusercontent.com/ZJUFanLab/scCATCH/master/R/sysdata.rda", temp, method = "auto")
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
# saveRDS(CellMatch,file = "CellMatch.rds")
unlink(temp)


# CellMarker
temp <- tempfile()
state <- 1
while (state != 0) {
  state <- tryCatch(expr = {
    download.file("http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/all_cell_markers.txt", temp, method = "auto")
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
# saveRDS(CellMarker,file = "CellMarker.rds")
unlink(temp)


# PanglaoDB
temp <- tempfile()
state <- 1
while (state != 0) {
  state <- tryCatch(expr = {
    download.file("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz", temp, method = "auto")
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
# saveRDS(PanglaoDB,file = "PanglaoDB.rds")
unlink(temp)
