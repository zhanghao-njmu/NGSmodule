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
corrected[, "chr"] <- gsub(x = corrected[, chr], pattern = "chr", replacement = "", perl = T)
chr_order <- intersect(c(1:30, "X", "Y"),corrected[["chr"]])
corrected <- dplyr::filter(corrected, chr %in% chr_order)
corrected[, "chr"] <- factor(
  x = corrected[, chr],
  levels = chr_order
)

CNA.object <- CNA(log2(corrected$reads/median(corrected$reads)),
  chrom = corrected$chr, maploc = corrected$start,
  data.type = "logratio", sampleid = sample,presorted = TRUE
)
smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha = 0.001, undo.splits = "sdundo", undo.SD = 3, verbose = 0)

corrected <- as.data.frame(corrected)
corrected[, "sample"] <- sample
corrected[, "DNAcopy"] <- segment.smoothed.CNA.object$data[[sample]]
corrected[, "CN_predict"] <- 2^corrected[, "DNAcopy"] * ploid
corrected[corrected[, "CN_predict"] > 4 * ploid & !is.na(corrected[, "CN_predict"]), "CN_predict"] <- 4 * ploid

color <- pal_material("blue-grey")(10)[c(3, 9, 1, 4)]
chr_info <- corrected %>%
  group_by(chr) %>%
  dplyr::summarise(
    chr = unique(chr),
    length = max(end),
    reads = sum(reads)
  )
chr_info[, "offset"] <- c(0, cumsum(as.numeric(chr_info$length))[-nrow(chr_info)])
chr_info[, "chr_cum_median"] <- chr_info[, "offset"] + chr_info[, "length"] / 2
chr_info[, "chr_cum_end"] <- chr_info[, "offset"] + chr_info[, "length"]
chr_info[, "chr_color"] <- rep(c(color[1], color[2]), times = ceiling(nrow(chr_info) / 2))[1:nrow(chr_info)]

corrected <- merge(corrected, chr_info, by.x = "chr", by.y = "chr")
corrected[, "cum_pos"] <- corrected[, "start"] + corrected[, "offset"]
corrected <- arrange(corrected, cum_pos)


colnames(segment.smoothed.CNA.object$output) <- c("ID", "chr", "start", "end", "num.mark", "median")
segs <- merge(segment.smoothed.CNA.object$output, chr_info, all.x = TRUE, by = "chr")
segs[, "cum_start"] <- segs[, "start"] + segs[, "offset"]
segs[, "cum_end"] <- segs[, "end"] + segs[, "offset"]
segs[, "CN_predict"] <- 2^segs[, "median"] * ploid # 2 ^ segs[,"median"] * 2 assumed to be diploid
segs[segs[, "CN_predict"] > 4 * ploid & !is.na(segs[, "CN_predict"]), "CN_predict"] <- 4 * ploid
segs <- arrange(segs, cum_start)
segs[, "previous_CN_predict"] <- lag(segs[, "CN_predict"], default = segs[, "CN_predict"][1])


# plot the cnv -----------------------------------------------------------
p1 <- ggplot() +
  # geom_rect(
  #   data = chr_info,
  #   aes(
  #     xmin = offset, xmax = chr_cum_end,
  #     ymin = 0, ymax = 8, fill = chr_color
  #   ),
  #   color = NA, alpha = 0.5
  # ) +
  geom_vline(xintercept = pull(chr_info, "offset"), linetype = 1, color = "grey80", size = 0.5) +
  geom_point(
    data = corrected,
    aes(x = cum_pos, y = CN_predict, color = chr_color),
    shape = 20, alpha = 1
  ) +
  geom_segment(
    data = segs,
    aes(x = cum_start, y = CN_predict, xend = cum_end, yend = CN_predict), size = 1
  ) +
  geom_segment(
    data = segs,
    aes(x = cum_start, y = previous_CN_predict, xend = cum_start, yend = CN_predict), size = 1
  ) +
  scale_color_identity() +
  scale_fill_manual(values = setNames(color[c(3, 4)], color[c(1, 2)]), guide = FALSE) +
  scale_y_continuous(breaks = seq(0, 8, 1),limits = c(0,4 * ploid)) +
  scale_x_continuous(breaks = pull(chr_info, "chr_cum_median"), labels = pull(chr_info, "chr")) +
  labs(title = sample, y = "Copy Number\n(assumed to be diploid)") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(colour = "grey80", linetype = 2),
    panel.border = element_rect(fill = "transparent", color = "black", size = 1),
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

p2 <- ggplot(corrected) +
  geom_histogram(
    aes(y = CN_predict),
    bins = 50, color = "black", fill = color[4]
  ) +
  scale_y_continuous(breaks = seq(0, 8, 1),limits = c(0,4 * ploid)) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_line(colour = "grey80", linetype = 2)
  )

p3 <- ggplot() +
  geom_vline(xintercept = pull(chr_info, "offset"), linetype = 1, color = "grey80", size = 0.5) +
  geom_point(
    data = subset(corrected, reads.x < mean(reads.x)*4),
    aes(x = cum_pos, y = reads.x, color = chr_color),
    shape = 20, alpha = 1
  ) +
  scale_color_identity() +
  scale_fill_manual(values = setNames(color[c(3, 4)], color[c(1, 2)]), guide = FALSE) +
  scale_y_continuous(limits = c(0,4*mean(corrected[["reads.x"]]))) +
  scale_x_continuous(breaks = pull(chr_info, "chr_cum_median"), labels = pull(chr_info, "chr")) +
  labs(title = sample, y = "Number of reads per bin") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(colour = "grey80", linetype = 2),
    panel.border = element_rect(fill = "transparent", color = "black", size = 1),
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

p4 <- ggplot(data = subset(corrected, reads.x < mean(reads.x)*4)) +
  geom_histogram(
    aes(y = reads.x),
    bins = 50, color = "black", fill = color[4]
  ) +
  scale_y_continuous(limits = c(0,4*mean(corrected[["reads.x"]]))) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_line(colour = "grey80", linetype = 2)
  )

plot <- aplot::plot_list(gglist = list(p1, p2, p3, p4), nrow = 2, ncol = 2, widths = rep(c(0.9, 0.1), 2))

##### output report #####
pdf(paste0(sample, ".plot.pdf"), width = nrow(chr_info) / 2 * 1.1 + 1, height = 3)
invisible(print(plot))
invisible(dev.off())

saveRDS(corrected, file = paste0(sample, ".corrected.rds"))
saveRDS(segs, file = paste0(sample, ".segs.rds"))
saveRDS(chr_info, file = paste0(sample, ".chr_info.rds"))
saveRDS(list(p1, p2, p3, p4), file = paste0(sample, ".plotlist.rds"))


##### check whether the unwanted file exists and remove it #####
if (file.exists("Rplots.pdf")) {
  invisible(file.remove("Rplots.pdf"))
}


