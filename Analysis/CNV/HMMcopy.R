#!/usr/bin/env Rscript


# rfile <- "/data/lab/LiLaiHua/10-examples/NGSmodule_work/SCG_C1/bwa/CNV/HMMcopy/SCG_C1.bwa.w1000000.wig"
# gfile <- "/data/database/iGenomes/Macaca_fascicularis/UCSC/Macaca_fascicularis_5.0_Add_rhesus_chY/Sequence/GemIndex/windows/1000000/genome.w1000000.gc.wig"
# mfile <- "/data/database/iGenomes/Macaca_fascicularis/UCSC/Macaca_fascicularis_5.0_Add_rhesus_chY/Sequence/GemIndex/windows/1000000/genome.w1000000.150mer.gem.wig"
# ploid <- 2
# sample <- "test"

args <- commandArgs(TRUE)
rfile <- as.character(args[1])
gfile <- as.character(args[2])
mfile <- as.character(args[3])
ploid <- as.numeric(args[4])
sample <- as.character(args[5])

library(HMMcopy)
library(dplyr)
library(scales)
library(ggplot2)
library(ggsci)
library(aplot)

uncorrected <- wigsToRangedData(readfile = rfile, gcfile = gfile, mapfile = mfile)
corrected <- correctReadcount(uncorrected)
corrected[, "chr"] <- gsub(x = corrected[, chr], pattern = "chr", replacement = "", perl = T)
corrected <- dplyr::filter(corrected,!str_detect(chr,pattern = "(M)|(MT)|(Mt)|(mt)"))
chr_order <- c(1:30, "X", "Y")
chr_uniq <- unique(corrected[, chr])
chr_levels <- c(
  chr_order[chr_order %in% chr_uniq],
  chr_uniq[!chr_uniq %in% chr_order]
)
corrected[, "chr"] <- factor(
  x = corrected[, chr],
  levels = chr_levels
)

segments <- HMMsegment(corrected, maxiter = 1e4)
corrected <- as.data.frame(corrected)
corrected[, "sample"] <- sample
corrected[, "CN_predict"] <- 2^corrected[, "copy"] * ploid
corrected[, "state"] <- segments$state

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

segs <- merge(segments[["segs"]], chr_info, all.x = TRUE, by = "chr")
segs[, "cum_start"] <- segs[, "start"] + segs[, "offset"]
segs[, "cum_end"] <- segs[, "end"] + segs[, "offset"]
segs[, "CN_predict"] <- 2^segs[, "median"] * ploid # 2 ^ segs[,"median"] * 2 assumed to be diploid
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
    data = subset(corrected, CN_predict < 8),
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
  scale_y_continuous(breaks = seq(0, 8, 1)) +
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

p2 <- ggplot(data = subset(corrected, CN_predict < 8)) +
  geom_histogram(
    aes(y = CN_predict),
    bins = 50, color = "black", fill = color[4]
  ) +
  scale_y_continuous(breaks = seq(0, 8, 1)) +
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
    data = subset(corrected, CN_predict < 8 & reads.x < quantile(reads.x, 0.999)),
    aes(x = cum_pos, y = reads.x, color = chr_color),
    shape = 20, alpha = 1
  ) +
  scale_color_identity() +
  scale_fill_manual(values = setNames(color[c(3, 4)], color[c(1, 2)]), guide = FALSE) +
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

p4 <- ggplot(data = subset(corrected, CN_predict < 8 & reads.x < quantile(reads.x, 0.999))) +
  geom_histogram(
    aes(y = reads.x),
    bins = 50, color = "black", fill = color[4]
  ) +
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
