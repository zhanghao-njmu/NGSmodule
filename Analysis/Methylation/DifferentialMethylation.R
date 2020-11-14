#!/usr/bin/env Rscript

library(methylKit)
library(limma)
library(bsseq)
library(dplyr)
library(stringr)
library(data.table)

NGSmodule_work <- "/data/huangmingqian/BSseq/NGSmodule_work"
SampleInfoFile <- "/data/huangmingqian/BSseq/temp_20201029014752.Sample_info.csv"
contrasts <- c("ESC,iMeLC;iMeLC,PGCLCd6;PGCLCd6,CellLine1;PGCLCd6,CellLine2;")
build <- "GRCh38"
threads <- 100
tiles_win <- 1000

contrasts_groups <- contrasts %>%
  str_split("[;,]") %>%
  unlist() %>%
  .[. != ""] %>%
  unique()

samples <- list.dirs(path = NGSmodule_work, full.names = FALSE, recursive = FALSE)
sample_info <- read.csv(SampleInfoFile, stringsAsFactors = F, header = T)
sample_info <- as.data.frame(sapply(sample_info, as.character))
sample_info <- sample_info %>%
  filter(SampleID %in% samples & Group %in% contrasts_groups) %>%
  group_by(SampleID) %>%
  mutate(RunID = paste0(RunID, collapse = ",")) %>%
  distinct() %>%
  as.data.frame()
rownames(sample_info) <- sample_info[, "SampleID"]
sample_info[, "Group"] <- factor(sample_info[, "Group"], unique(sample_info[, "Group"]))
sample_info <- sample_info[order(sample_info[, "Group"]), ]
sample_info <- sample_info[!is.na(sample_info[["Group"]]), ]


# Load data ---------------------------------------------------------------
# cov <- readBismarkCoverage(
#   location = paste0(NGSmodule_work, "/", sample_info[["SampleID"]], "/bismark_bowtie2/bismark_methylation_extractor/", sample_info[["SampleID"]], "_bismark_bt2_pe.deduplicated.bismark.cov.gz"),
#   sample.id = sample_info[["SampleID"]], assembly = build, treatment = sample_info[["Group"]], context = "CpG", min.cov = 3
# )
CpG <- readBismarkCytosineReport(
  location = paste0(NGSmodule_work, "/", sample_info[["SampleID"]], "/bismark_bowtie2/bismark_methylation_extractor/", sample_info[["SampleID"]], "_bismark_bt2_pe.deduplicated.CpG_report.txt.gz"),
  sample.id = sample_info[["SampleID"]], assembly = build, treatment = sample_info[["Group"]], context = "CpG", min.cov = 3
)

# base-resolution ---------------------------------------------------------
CpG_filtered <- filterByCoverage(CpG,
  lo.count = 3, lo.perc = NULL,
  hi.count = NULL, hi.perc = 99.9
)
CpG_unite <- unite(CpG_filtered, destrand = FALSE,mc.cores = threads)
CpG_unite_meth <- percMethylation(CpG_unite, rowids = TRUE)
CpG_unite_meth_logistic <- log2((CpG_unite_meth + 1) / (100 - CpG_unite_meth + 1))

# region-resolution -------------------------------------------------------
tiles <- tileMethylCounts(CpG, win.size = tiles_win, step.size = tiles_win, cov.bases = ceiling(tiles_win/10), mc.cores = threads)
tiles_filtered <- filterByCoverage(tiles,
  lo.count = 10, lo.perc = NULL,
  hi.count = NULL, hi.perc = 99.9
)
tiles_unite <- unite(tiles_filtered, destrand = FALSE,mc.cores = threads)
tiles_unite_meth <- percMethylation(tiles_unite, rowids = TRUE)
tiles_unite_meth_logistic <- log2((tiles_unite_meth + 1) / (100 - tiles_unite_meth + 1))








design <- model.matrix(~ 0 + Group, data = sample_info)
colnames(design) <- gsub(x = colnames(design), pattern = "Group", replacement = "")
limma_contrasts <- contrasts %>%
  str_split(";") %>%
  unlist() %>%
  .[. != ""] %>%
  gsub(x = ., pattern = ",", replacement = "-") %>%
  makeContrasts(contrasts = ., levels = design)

fit <- lmFit(CpG_unite_meth_logistic, design = design)
fit2 <- contrasts.fit(fit, limma_contrasts)
fit2 <- eBayes(fit2)

Diff <- lapply(setNames(colnames(limma_contrasts), colnames(limma_contrasts)), function(x) {
  topTable(fit2, adjust.method = "BH", coef = x, sort.by = "P", number = Inf)
})










# getMethylationStats(meth_CpG[[1]], plot = T, both.strands = FALSE)
# getCoverageStats(meth_CpG[[1]], plot = TRUE, both.strands = FALSE)


getCorrelation(meth, plot = TRUE)
clusterSamples(meth, dist = "correlation", method = "ward", plot = TRUE)
PCASamples(meth)
as <- assocComp(mBase = meth, data.frame(batch_id = sample_info[["BatchID"]]))
newObj <- removeComp(meth, comp = 1)

tiles <- tileMethylCounts(meth_CpG, win.size = 1000, step.size = 1000, cov.bases = 10)
myDiff <- calculateDiffMeth(meth)
res_limma <- limma.meth(methylBase.obj = meth)
res_DSS <- run.DSS(methylBase.obj = meth, difference = 5)






############ test ###################

# list the files

# load packages required for analysis
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")
list.files(dataDirectory, recursive = TRUE)
targets <- read.metharray.sheet(dataDirectory, pattern = "SampleSheet.csv")
targets
rgSet <- read.metharray.exp(targets = targets)
targets$ID <- paste(targets$Sample_Group, targets$Sample_Name, sep = ".")
sampleNames(rgSet) <- targets$ID
rgSet
