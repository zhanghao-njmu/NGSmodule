#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
rfile <- as.character(args[1])
gfile <- as.character(args[2])
mfile <- as.character(args[3])
ploid <- as.numeric(args[4])
sample <- as.character(args[5])

# rfile <- "/ssd/lab/CuiYiQiang/WGS/mus_1N/NGSmodule_work/SC_112/bowtie2/CNV/HMMcopy/SC_112.bowtie2.w1000000.wig"
# gfile <- "/archive/reference/iGenomes/Mus_musculus/UCSC/mm10/Sequence/GemIndex/Mappability/150mer/windows/win1000000/genome.win1000000.gc.wig"
# mfile <- "/archive/reference/iGenomes/Mus_musculus/UCSC/mm10/Sequence/GemIndex/Mappability/150mer/windows/win1000000/genome.win1000000.150mer.gem.wig"
# ploid <- 2
# sample <- "SC_112.bowtie2.HMMcopy"

library(HMMcopy)
library(DNAcopy)
library(dplyr)
library(stringr)
library(scales)
library(ggplot2)
library(ggsci)
library(aplot)

uncorrected <- wigsToRangedData(readfile = rfile, gcfile = gfile, mapfile = mfile)
corrected <- correctReadcount(uncorrected)

data(coriell)
ob <- coriell$Coriell.05296
coriell$Chromosome,coriell$Position,




rawratios <- read.table(input, header = TRUE)
ratios <- rawratios[order(rawratios$chrmarker), ]

cols <- names(ratios)
for (i in 3)
{
  # eval
  ratiostmp <- CNA(, ratios$chr, ratios$pos, data.type = "logratio", sampleid = cols[i], presorted = TRUE)
  smooth.ratiostmp <- smooth.CNA(ratiostmp)
  # 	seg.smooth.ratiostmp<-segment(smooth.ratiostmp, verbose=1);
  seg.smooth.ratiostmp <- segment(smooth.ratiostmp, verbose = 1, alpha = 0.001, undo.splits = "sdundo", undo.SD = 1)
  # change alpha for different sensitivity

  png(file = paste(cols[i], "png", sep = "."), width = 1500, height = 480, bg = "white")
  plot(seg.smooth.ratiostmp, ylim = c(-3, 3))
  dev.off()
  write.table(segments.summary(seg.smooth.ratiostmp), quote = FALSE, sep = "\t", file = paste(cols[i], "seg", sep = "."))
}



