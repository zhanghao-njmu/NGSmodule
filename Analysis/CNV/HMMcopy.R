#!/usr/bin/env Rscript

# evalue <- 0.995
# rfile <- "/ssd/llh_haploid/NGSmodule_work/m-1_FDSW202470982-1r/bwa/CNV/HMMcopy/m-1_FDSW202470982-1r.bwa.1M_input.wig"
# gfile <- "/data/database/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/GenmapIndex/windows/1000000/genome_main.w1000000.gc.wig"
# mfile <- "/data/database/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/GenmapIndex/windows/1000000/genome_main.w1000000.130mer.genmap.wig"
# ploid <- 2
# samp <- "test"

library(HMMcopy)
library(scales)

args <- commandArgs(TRUE)
evalue <- as.numeric(args[1])
rfile <- as.character(args[2])
gfile <- as.character(args[3])
mfile <- as.character(args[4])
ploid <- as.numeric(args[5])
samp <- as.character(args[6])


x1_untrimmed_reads <- wigsToRangedData(readfile = rfile, gcfile = gfile, mapfile = mfile)
corrected <- correctReadcount(x1_untrimmed_reads)
param <- HMMsegment(corrected, getparam = TRUE)
param$e <- evalue

segments <- HMMsegment(corrected, param = param)


lengths <- aggregate(corrected$end, by = list(chr = corrected$chr), max)
chr_order <- c(paste0("chr", 1:25), "chrX", "chrY")
lengths <- lengths[order(match(as.character(lengths$chr), chr_order)), ]
lengths$offset <- c(0, cumsum(as.numeric(lengths$x))[-nrow(lengths)])
corrected$state <- segments$state
corrected$sample <- samp
write.csv(x = corrected, file = "HMMcopy.corrected.csv", row.names = F)

pdf(paste(samp, ".pdf", sep = ""), height = 3.2, width = 12)
data <- merge(corrected, lengths[, c(1, 3)], by.x = "chr", by.y = "chr") 
dat2 <- data.frame(x_position = data$start + data$offset, y_position = data$copy)
dat2 <- dat2[with(dat2, order(x_position)), ]
x_position <- dat2$x_position
y_position <- dat2$y_position
a <- data.frame()
n <- 0
chr <- 1
for (i in 1:length(x_position)) {
  if (chr <= 20) {
    if (lengths$offset[chr] <= x_position[i] & x_position[i] <= lengths$offset[chr + 1]) {
      a[i, c("x_position")] <- x_position[i]
      a[i, c("chr")] <- chr
    } else {
      chr <- chr + 1
      a[i, c("x_position")] <- x_position[i]
      a[i, c("chr")] <- chr
    }
  } else {
    a[i, c("x_position")] <- x_position[i]
    a[i, c("chr")] <- chr
  }
  if (chr %% 2 == 0) {
    a[i, c("color")] <- "#4494C1"
  } else {
    a[i, c("color")] <- "#E15F47"
  }
}
range <- quantile(data$copy, na.rm = TRUE, prob = c(0.01, 0.99)) # Optional range, tweak this as desired


segs <- merge(segments$segs, lengths[, c(1, 3)], all.x = TRUE)
segs$x_start <- segs$start + segs$offset
segs$x_end <- segs$end + segs$offset
segs$y_median <- 2^segs$median * ploid # 2^segs$median*2 assumed to be diploid
segs <- segs[with(segs, order(x_start)), ]
segs$pre_median <- c(segs$y_median[1], segs$y_median[1:(length(segs$y_median) - 1)])

# 2^y_position*2 assumed to be diploid
plot(x = x_position, y = 2^y_position * ploid, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = alpha(a$color, 0.7), pch = 20, ylim = c(0, 8), cex = 0.4)
title(main = samp)
title(xlab = list("Chromosome"), ylab = list("CopyNumber"), mgp = c(2, 1, 1))
axis(2, at = 1:8, las = 2)
abline(v = lengths$offset, lty = 6, col = "#8F8F8F", lwd = 1.2) # Separate chromosomes
abline(h = c(0:8), lty = 3, col = "#BABABA")

lengths$median_offset <- lengths$offset + (lengths$x / 2)
text(lengths$median_offset, -1, labels = gsub("chr", "", lengths$chr), xpd = NA) # Label chromosomes

# Inserting the segmentation lines if desired...
segments(x0 = segs$x_start, y0 = segs$y_median, x1 = segs$x_end, col = "#363636", lwd = 1.5) # Note, segments is a function... not the output, unfortunately naming...
segments(x0 = segs$x_start, y0 = segs$pre_median, y1 = segs$y_median, col = "#363636", lwd = 1.5) # Note, segments is a function... not the output, unfortunately naming...

dev.off()

write.table(segments$segs, file = "segments.txt", sep = "\t")

# n=matrix(c(as.character(corrected$space),corrected$copy),ncol=2)
# write.table(n, file="data.txt", sep="\t",row.names=FALSE, col.names=FALSE)
