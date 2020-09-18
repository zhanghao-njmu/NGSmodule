#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
vcfFile <- as.character(args[1])
CNV_prefix <- as.character(args[2])
sample <- as.character(args[3])

library(vcfR)
library(dplyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggsci)
library(ggforce)
library(cowplot)
library(png)
library(grid)
library(gridExtra)

# vcfFile <- system.file("extdata", "pinf_sc50.vcf.gz",
#                         package = "pinfsc50")
# vcf <- read.vcfR(vcfFile, verbose = F)
# sample_index <- "P17777us22" # 2
# sample_index <- "P13626" # 3

vcf <- read.vcfR(vcfFile, verbose = F)
sample_index <- 1
plotlist <- list()
color <- pal_material("blue-grey")(10)[c(3, 9, 1, 4)]

chr_info <- readRDS(file = paste0(CNV_prefix, ".chr_info.rds"))
p1 <- readRDS(file = paste0(CNV_prefix, ".p1.rds"))
p2 <- readRDS(file = paste0(CNV_prefix, ".p2.rds"))

##### Allele balance #####
chrom <- vcf@fix[, "CHROM"]
pos <- as.numeric(vcf@fix[, "POS"])
qual <- as.numeric(vcf@fix[, "QUAL"])
mq <- extract.info(vcf, element = "MQ", as.numeric = TRUE)
dp <- extract.gt(vcf, element = "DP", as.numeric = TRUE)[, sample_index, drop = F]
gq <- extract.gt(vcf, element = "GQ", as.numeric = TRUE)[, sample_index, drop = F]
ad <- extract.gt(vcf, element = "AD")[, sample_index, drop = F]
allele1 <- masplit(ad, sort = 1, record = 1)
allele2 <- masplit(ad, sort = 1, record = 2)
ad1 <- allele1 / (allele1 + allele2)
ad2 <- allele2 / (allele1 + allele2)
gt <- extract.gt(vcf, element = "GT")[, sample_index, drop = F]
hets <- is.het(gt)

alleles <- cbind(
  gsub(x = chrom, pattern = "chr", replacement = "", perl = T),
  pos,
  gt[, 1],
  hets[, 1],
  allele1,
  allele2,
  ad1,
  ad2,
  dp[, 1],
  qual,
  mq
) %>%
  as.data.frame() %>%
  na.omit()
colnames(alleles) <- c(
  "chr", "Position", "GeneType", "is_het",
  "Allele1", "Allele2", "Major_allele", "Minor_allele", "DP", "QUAL", "MQ"
)
chr_order <- c(1:30, "X", "Y", "MT")
chr_uniq <- unique(pull(chr_info, "chr"))
chr_levels <- c(
  chr_order[chr_order %in% chr_uniq],
  chr_uniq[!chr_uniq %in% chr_order]
)
alleles[, "chr"] <- factor(
  x = alleles[, "chr"],
  levels = chr_levels
)
alleles[, "is_het"] <- factor(
  x = alleles[, "is_het"],
  levels = c(FALSE, TRUE)
)
numeric_col <- c("Position", "Allele1", "Allele2", "Major_allele", "Minor_allele", "DP", "QUAL", "MQ")
alleles[, numeric_col] <- sapply(alleles[, numeric_col], as.numeric)


genetype_summary <- function(GeneType_table) {
  x <- GeneType_table
  is_het <- c(is.het(as.matrix(names(x))))
  ratio <- NULL
  if (all(is_het)) {
    het <- sum(x[is_het])
    hom <- 0
    ratio <- Inf
  }
  if (all(!is_het)) {
    het <- 0
    hom <- sum(x[!is_het])
    ratio <- (-Inf)
  }
  if (is.null(ratio)) {
    het <- sum(x[is_het])
    hom <- sum(x[!is_het])
    ratio <- as.character(round(het / hom, digits = 3))
  }
  x <- sort(x, decreasing = T)[1:min(4, length(x))]
  n <- names(x)
  x <- x[!is.na(x)]
  label <- paste0(paste0("Het/Hom=", ratio, "\n"), paste0(" ", n, ": ", x, collapse = "\n"))

  return(setNames(c(het, hom, as.numeric(ratio), label), c("het", "hom", "ratio", "label")))
}

# Allele Depth distribution -----------------------------------------------
AD_q15 <- max(
  quantile(alleles[["Allele1"]], probs = 0.15),
  quantile(alleles[["Allele2"]], probs = 0.15)
)
AD_q95 <- max(
  quantile(alleles[["Allele1"]], probs = 0.95),
  quantile(alleles[["Allele2"]], probs = 0.95)
)
ad_df <- melt(alleles, measure.vars = c("Allele1", "Allele2"), variable.name = "Allele")
p_AD <- ggplot() +
  geom_histogram(
    data = subset(ad_df, value > 0 & value < AD_q95),
    aes(x = value, fill = Allele), color = "black", alpha = 0.8
  ) +
  geom_vline(xintercept = AD_q15, color = "red3", linetype = 2, size = 1) +
  annotate("text", x = AD_q15, y = Inf, vjust = 1, hjust = 0, label = paste("q15:", AD_q15), color = "red3") +
  scale_fill_manual(name = "Allele", values = color[c(1, 2)]) +
  labs(title = paste("Allele depth(AD):", sample), x = "", y = "Count") +
  theme_classic() +
  theme(
    aspect.ratio = 0.8,
    panel.grid.major = element_line(colour = "grey80")
  )

# Depth distribution ------------------------------------------------------
DP_q15 <- quantile(alleles[, "DP"], probs = 0.15)
DP_q95 <- quantile(alleles[, "DP"], probs = 0.95)
p_DP <- ggplot() +
  geom_histogram(
    data = subset(alleles, DP > 0 & DP < DP_q95),
    aes(x = DP, fill = is_het), color = "black", alpha = 0.8
  ) +
  geom_vline(xintercept = DP_q15, color = "red3", linetype = 2, size = 1) +
  annotate("text", x = DP_q15, y = Inf, vjust = 1, hjust = 0, label = paste("q15:", DP_q15), color = "red3") +
  scale_fill_manual(name = "Heterozygous", values = color[c(1, 2)]) +
  labs(title = paste("Depth(DP):", sample), x = "", y = "Count") +
  theme_classic() +
  theme(
    aspect.ratio = 0.8,
    panel.grid.major = element_line(colour = "grey80")
  )

# QUAL distribution -------------------------------------------------------
QUAL_q15 <- quantile(alleles[, "QUAL"], probs = 0.15)
QUAL_q95 <- quantile(alleles[, "QUAL"], probs = 0.95)
p_QUAL <- ggplot() +
  geom_histogram(
    data = subset(alleles, QUAL > 0 & QUAL < QUAL_q95),
    aes(x = QUAL, fill = is_het), color = "black", alpha = 0.8
  ) +
  geom_vline(xintercept = QUAL_q15, color = "red3", linetype = 2, size = 1) +
  annotate("text", x = QUAL_q15, y = Inf, vjust = 1, hjust = 0, label = paste("q15:", QUAL_q15), color = "red3") +
  scale_fill_manual(name = "Heterozygous", values = color[c(1, 2)]) +
  labs(title = paste("Quality(QUAL):", sample), x = "", y = "Count") +
  theme_classic() +
  theme(
    aspect.ratio = 0.8,
    panel.grid.major = element_line(colour = "grey80")
  )

# Mapping quality distribution --------------------------------------------
MQ_q15 <- quantile(alleles[, "MQ"], probs = 0.15)
p_MQ <- ggplot() +
  geom_histogram(
    data = alleles,
    aes(x = MQ, fill = is_het), color = "black", alpha = 0.8
  ) +
  geom_vline(xintercept = MQ_q15, color = "red3", linetype = 2, size = 1) +
  annotate("text",
    x = MQ_q15, y = Inf, vjust = 1, hjust = 1,
    label = paste("q15:", MQ_q15), color = "red3"
  ) +
  scale_fill_manual(name = "Heterozygous", values = color[c(1, 2)]) +
  labs(title = paste("Mapping quality:", sample), y = "Count") +
  theme_classic() +
  theme(
    aspect.ratio = 0.8,
    panel.grid.major = element_line(colour = "grey80")
  )

plotlist[["Alleles_QC"]] <- plot_grid(p_AD, p_DP, p_QUAL, p_MQ, nrow = 2)

# subset alleles ----------------------------------------------------------
alleles <- subset(alleles, (Allele1 >= AD_q15 | Allele2 >= AD_q15) &
  DP > DP_q15 & QUAL > QUAL_q15 & MQ > MQ_q15)

# all alleles -------------------------------------------------------------
p_Chr <- ggplot(alleles) +
  geom_bar(aes(x = chr, fill = is_het)) +
  scale_fill_manual(name = "Heterozygous", values = color[c(1, 2)], drop = F) +
  scale_x_discrete(drop = FALSE) +
  labs(title = paste("Alleles on each chromosome:", sample), y = "Count") +
  theme_classic() +
  theme(
    aspect.ratio = 0.8,
    axis.title.x = element_blank(),
    strip.background = element_rect(fill = "grey85")
  )

alleles_melt <- melt(alleles, measure.vars = c("Major_allele", "Minor_allele"), variable.name = "Type")
p0 <- ggplot(data = alleles_melt) +
  geom_histogram(aes(x = value, y = ..density.., fill = Type), color = "black", bins = 50, alpha = 0.8) +
  geom_density(aes(x = value), fill = "grey80", color = "black", alpha = 0.5) +
  scale_x_continuous(
    breaks = c(0, 1 / 4, 1 / 3, 1 / 2, 2 / 3, 3 / 4, 1),
    labels = c("0", "1/4", "1/3", "1/2", "1/3", "1/4", "1")
  ) +
  scale_fill_manual(
    name = "Type",
    values = setNames(
      color[c(2, 1)],
      c("Major_allele", "Minor_allele")
    )
  ) +
  labs(x = "Minor/Major allele frequency", y = "Density") +
  theme_classic() +
  theme(
    aspect.ratio = 0.8
  )

p_all <- p0 + annotate("text",
  label = genetype_summary(table(alleles[, "GeneType"]))["label"],
  x = 0.1, y = Inf,
  hjust = 0, vjust = 1.1
) +
  labs(title = paste("All alleles:", sample))
plotlist[["All_alleles_frequency"]] <- plot_grid(p_all, p_Chr, nrow = 1)

gt_chr <- alleles %>%
  group_by(chr) %>%
  select("GeneType") %>%
  table() %>%
  as.data.frame.matrix()
gt_label_chr <- lapply(1:nrow(gt_chr), function(i) {
  genetype_summary(gt_chr[i, ])
}) %>%
  bind_rows()
gt_label_chr[, "het"] <- as.numeric(pull(gt_label_chr, "het"))
gt_label_chr[, "hom"] <- as.numeric(pull(gt_label_chr, "hom"))
gt_label_chr[, "ratio"] <- as.numeric(pull(gt_label_chr, "ratio"))
gt_label_chr[, "chr"] <- rownames(gt_chr)
gt_label_chr[, "chr"] <- factor(
  x = pull(gt_label_chr, "chr"),
  levels = chr_levels
)
gt_label_chr <- gt_label_chr %>% mutate(color = case_when(
  chr %in% c("X") ~ "red3",
  chr %in% c("Y") ~ "royalblue3",
  TRUE ~ "black"
))
gt_label_chr <- gt_label_chr %>% mutate(fill = case_when(
  chr %in% c("X") ~ "red3",
  chr %in% c("Y") ~ "royalblue3",
  TRUE ~ "#90A4ADFF"
))
gt_label_chr <- gt_label_chr %>%
  mutate(
    log10ratio = case_when(
      log10(ratio) < (-3) ~ -3,
      log10(ratio) > 3 ~ 3,
      TRUE ~ log10(ratio)
    ),
    hjust = case_when(
      log10(ratio) <= 0 ~ 0,
      log10(ratio) > 0 ~ 1,
      is.na(log10(ratio)) ~ 0
    )
  )
gt_label_chr <- merge(x = gt_label_chr, y = chr_info, by = "chr", all.x = TRUE)
gt_label_chr <- gt_label_chr %>%
  mutate(
    het_norm = het / reads,
    hom_norm = hom / reads
  )
gt_label_chr_melt <- melt(gt_label_chr, measure.vars = c("hom_norm", "het_norm"), variable.name = "Type")


p <- p0 + facet_wrap(. ~ chr, nrow = 4, ) +
  geom_text(
    data = gt_label_chr, aes(label = label),
    x = 0.1, y = Inf, size = 3,
    hjust = 0, vjust = 1.1
  ) +
  labs(title = paste("All alleles:", sample))
plotlist[["All_alleles_frequency_bychr"]] <- p

p_ratio1 <- ggplot(gt_label_chr, aes(x = log10ratio, y = reorder(chr, desc(chr)), fill = fill, hjust = hjust)) +
  geom_col(color = "black") +
  geom_vline(xintercept = 0, color = "black", size = 1) +
  geom_vline(xintercept = c(-1, 1), color = c("royalblue4", "red3"), size = 1, linetype = 2) +
  geom_text(aes(x = ifelse(log10ratio < 0, 0.15, -0.15), label = as.character(round(log10ratio, digits = 3)), color = color), size = 5) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_continuous(breaks = seq(-3, 3, 1), limits = c(-3, 3)) +
  labs(title = paste("Het/Hom ratio:", sample), x = "log10(Het/Hom)", y = "Chromosome") +
  theme_classic() +
  theme(
    aspect.ratio = 0.8,
    panel.grid.major.y = element_line(colour = "grey80")
  )
p_ratio2 <- ggplot(data = gt_label_chr_melt, aes(x = chr, y = value, fill = Type)) +
  geom_col() +
  scale_fill_manual(
    name = "Heterozygous",
    values = setNames(color[c(1, 2)], c("hom_norm", "het_norm")),
    labels = c("FALSE", "TRUE")
  ) +
  labs(title = paste("Het/Hom ratio:", sample), y = "Average number of alleles per reads") +
  theme_classic() +
  theme(
    aspect.ratio = 0.8,
    axis.title.x = element_blank(),
    panel.grid.major = element_line(colour = "grey80")
  )
plotlist[["HetHom_ratio"]] <- plot_grid(p_ratio1, p_ratio2, nrow = 1)

# het alleles -------------------------------------------------------------
het_alleles <- subset(alleles, Allele1 > AD_q15 & Allele2 > AD_q15 & is_het == TRUE)

p_Chr <- ggplot(het_alleles) +
  geom_bar(aes(x = chr, fill = is_het)) +
  scale_fill_manual(name = "Heterozygous", values = color[c(1, 2)], drop = F) +
  scale_x_discrete(drop = FALSE) +
  labs(title = paste("Alleles on each chromosome:", sample), y = "Count") +
  theme_classic() +
  theme(
    aspect.ratio = 0.8,
    axis.title.x = element_blank(),
    strip.background = element_rect(fill = "grey85")
  )

het_alleles_melt <- melt(het_alleles, measure.vars = c("Major_allele", "Minor_allele"), variable.name = "Type")
p0 <- ggplot(data = het_alleles_melt) +
  geom_histogram(aes(x = value, y = ..density.., fill = Type),
    color = "black", bins = 50, alpha = 0.8
  ) +
  geom_density(aes(x = value), fill = "grey80", color = "black", alpha = 0.5) +
  scale_x_continuous(
    breaks = c(0, 1 / 4, 1 / 3, 1 / 2, 2 / 3, 3 / 4, 1),
    labels = c("0", "1/4", "1/3", "1/2", "1/3", "1/4", "1")
  ) +
  scale_fill_manual(
    name = "Type",
    values = setNames(
      color[c(2, 1)],
      c("Major_allele", "Minor_allele")
    )
  ) +
  labs(x = "Minor/Major allele frequency", y = "Density") +
  theme_classic() +
  theme(
    aspect.ratio = 0.8
  )

p_het <- p0 + annotate("text",
  label = genetype_summary(table(het_alleles[, "GeneType"]))["label"],
  x = 0.1, y = Inf,
  hjust = 0, vjust = 1.1
) +
  labs(title = paste("AD(q15)-filtered heterozygous:", sample), subtitle = paste("Allele1 >", AD_q15, "&", "Allele2 >", AD_q15))
plotlist[["Het_alleles_frequency"]] <- plot_grid(p_het, p_Chr, nrow = 1)

gt_chr <- het_alleles %>%
  group_by(chr) %>%
  select("GeneType") %>%
  table() %>%
  as.data.frame.matrix()
gt_label_chr <- lapply(1:nrow(gt_chr), function(i) {
  genetype_summary(gt_chr[i, ])
}) %>%
  bind_rows()
gt_label_chr[, "chr"] <- rownames(gt_chr)
gt_label_chr[, "chr"] <- factor(
  x = pull(gt_label_chr, "chr"),
  levels = chr_levels
)

p <- p0 + facet_wrap(. ~ chr, nrow = 4) +
  geom_text(
    data = gt_label_chr, aes(label = label),
    x = 0.1, y = Inf, size = 3,
    hjust = 0, vjust = 1.1
  ) +
  labs(title = paste("AD(q15)-filtered heterozygous:", sample), subtitle = paste("Allele1 >", AD_q15, "&", "Allele2 >", AD_q15))
plotlist[["Het_alleles_frequency_bychr"]] <- p


# Intergrated analysis ----------------------------------------------------
all_df <- merge(x = alleles, y = chr_info, by = "chr", all.x = T)
all_df[, "cum_pos"] <- all_df[, "Position"] + all_df[, "offset"]
all_df <- subset(all_df, cum_pos < max(chr_info[, "chr_cum_end"]))
all_df <- arrange(all_df, cum_pos)

het_df <- merge(x = het_alleles, y = chr_info, by = "chr", all.x = T)
het_df[, "cum_pos"] <- het_df[, "Position"] + het_df[, "offset"]
het_df <- subset(het_df, cum_pos < max(chr_info[, "chr_cum_end"]))
het_df <- arrange(het_df, cum_pos)


p3 <- ggplot(all_df) +
  geom_vline(xintercept = pull(chr_info, "offset"), linetype = 1, color = "grey80", size = 0.5) +
  geom_point(
    aes(x = cum_pos, y = Minor_allele, color = chr_color),
    shape = 20, alpha = 1, size = 0.1
  ) +
  geom_point(
    aes(x = cum_pos, y = Major_allele, color = chr_color),
    shape = 20, alpha = 1, size = 0.1
  ) +
  scale_color_identity() +
  scale_fill_manual(values = setNames(color[c(3, 4)], color[c(1, 2)]), guide = FALSE) +
  scale_x_continuous(breaks = pull(chr_info, "chr_cum_median"), labels = pull(chr_info, "chr")) +
  scale_y_continuous(
    breaks = c(0, 1 / 4, 1 / 3, 1 / 2, 2 / 3, 3 / 4, 1),
    labels = c("0", "1/4", "1/3", "1/2", "1/3", "1/4", "1")
  ) +
  labs(y = "Allele frequency\n(all alleles)") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(colour = "grey80", linetype = 2),
    panel.border = element_rect(fill = "transparent", color = "black", size = 1),
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

p4 <- ggplot(melt(all_df, measure.vars = c("Major_allele", "Minor_allele"), variable.name = "Type")) +
  geom_histogram(
    aes(x = ..density.., y = value),
    bins = 50, color = "black", fill = color[4]
  ) +
  geom_density(aes(y = value), fill = "grey80", color = "black", alpha = 0.5) +
  scale_y_continuous(
    breaks = c(0, 1 / 4, 1 / 3, 1 / 2, 2 / 3, 3 / 4, 1),
    labels = c("0", "1/4", "1/3", "1/2", "1/3", "1/4", "1")
  ) +
  labs(x = "Density") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_line(colour = "grey80", linetype = 2)
  )

p5 <- ggplot(het_df) +
  geom_vline(xintercept = pull(chr_info, "offset"), linetype = 1, color = "grey80", size = 0.5) +
  geom_point(
    aes(x = cum_pos, y = Minor_allele, color = chr_color),
    shape = 20, alpha = 1, size = 0.1
  ) +
  geom_point(
    aes(x = cum_pos, y = Major_allele, color = chr_color),
    shape = 20, alpha = 1, size = 0.1
  ) +
  scale_color_identity() +
  scale_fill_manual(values = setNames(color[c(3, 4)], color[c(1, 2)]), guide = FALSE) +
  scale_x_continuous(breaks = pull(chr_info, "chr_cum_median"), labels = pull(chr_info, "chr")) +
  scale_y_continuous(
    breaks = c(0, 1 / 4, 1 / 3, 1 / 2, 2 / 3, 3 / 4, 1),
    labels = c("0", "1/4", "1/3", "1/2", "1/3", "1/4", "1")
  ) +
  labs(y = "Allele frequency\n(AD-filtered heterozygous)", x = "Chromosome") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(colour = "grey80", linetype = 2),
    panel.border = element_rect(fill = "transparent", color = "black", size = 1),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

p6 <- ggplot(melt(het_df, measure.vars = c("Major_allele", "Minor_allele"), variable.name = "Type")) +
  geom_histogram(
    aes(x = ..density.., y = value),
    bins = 50, color = "black", fill = color[4]
  ) +
  geom_density(aes(y = value), fill = "grey80", color = "black", alpha = 0.5) +
  scale_y_continuous(
    breaks = c(0, 1 / 4, 1 / 3, 1 / 2, 2 / 3, 3 / 4, 1),
    labels = c("0", "1/4", "1/3", "1/2", "1/3", "1/4", "1")
  ) +
  labs(x = "Density") +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_line(colour = "grey80", linetype = 2)
  )

plot <- aplot::plot_list(list(p1, p2, p3, p4, p5, p6),
  nrow = 3, ncol = 2,
  widths = rep(c(0.9, 0.1), 3)
)
ggsave(plot, filename = paste0(sample, ".plot.png"), width = nrow(chr_info) / 2, height = 2 * 3)


##### output report #####
pdf(paste0(sample, ".SNV2ploidy.pdf"), width = 11, height = 8)
invisible(lapply(plotlist, print))
thePlot <- rasterGrob(readPNG(paste0(sample, ".plot.png"), native = FALSE),
  interpolate = FALSE
)
grid.arrange(thePlot)
invisible(dev.off())

##### check whether the unwanted file exists and remove it #####
if (file.exists("Rplots.pdf")) {
  invisible(file.remove("Rplots.pdf"))
}








#
# i <- 2
# winsize <- windows*10
# bin_width <- 0.02
#
# Peaks1 <- freq_peak(myMat = as.matrix(het_df[, "Major_allele"]),pos = het_df[,"cum_pos"], winsize = 1e6, bin_width = bin_width)
# Peaks2 <- freq_peak(myMat = as.matrix(all_df[, "Minor_allele"]),pos = all_df[,"cum_pos"], winsize = 1e6, bin_width = bin_width)

# layout(matrix(1:2, nrow = 1), widths = c(4, 1))
# par(mar = c(5, 4, 4, 0))
#
# mySample <- "test"
# plot(chr1_het_alleles$Position, chr1_het_alleles[, "Major_allele"],
#   ylim = c(0, 1), type = "n", yaxt = "n",
#   main = mySample, xlab = "POS", ylab = "Allele balance"
# )
# axis(
#   side = 2, at = c(0, 0.25, 0.333, 0.5, 0.666, 0.75, 1),
#   labels = c(0, "1/4", "1/3", "1/2", "2/3", "3/4", 1), las = 1
# )
# abline(h = c(0.25, 0.333, 0.5, 0.666, 0.75), col = 8)
# points(chr1_het_alleles$Position,chr1_het_alleles[, "Major_allele"], pch = 20, col = "#A6CEE344")
# points(chr1_het_alleles$Position,  chr1_het_alleles[, "Minor_allele"], pch = 20, col = "#1F78B444")
# segments(
#   x0 = myPeaks1$wins[, "START_pos"], y0 = myPeaks1$peaks[, 1],
#   x1 = myPeaks1$wins[, "END_pos"], lwd = 3
# )
# segments(
#   x0 = myPeaks1$wins[, "START_pos"], y0 = myPeaks2$peaks[, 1],
#   x1 = myPeaks1$wins[, "END_pos"], lwd = 3
# )
# bp1 <- hist(chr1_het_alleles[, "Major_allele"], breaks = seq(0, 1, by = bin_width), plot = FALSE)
# # Load libraries
# library(vcfR)
# library(pinfsc50)
#
# # Determine file locations
# vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz",
#                         package = "pinfsc50")
#
# # Read data into memory
# vcf <- read.vcfR(vcf_file, verbose = FALSE)
#
# ad <- extract.gt(vcf, element = 'AD')
#
# allele1 <- masplit(ad, record = 1)
# allele2 <- masplit(ad, record = 2)
#
# ad1 <- allele1 / (allele1 + allele2)
# ad2 <- allele2 / (allele1 + allele2)
#
# dp <- allele1
# # sums <- apply(dp, MARGIN=2, quantile, probs=c(0.15, 0.95), na.rm=TRUE)
# sums <- apply(dp, MARGIN = 2, quantile, probs = c(0.1, 0.9), na.rm = TRUE)
#
# par(mfrow = c(4, 3))
# par(mar = c(2, 2, 1, 1))
# par(oma = c(1, 1, 0, 0))
#
# for (i in 1:12) {
#   hist(allele1[, i], breaks = seq(0, 1e3, by = 1), xlim = c(0, 100), col = 8, main = "", xlab = "", ylab = "")
#   title(main = colnames(allele1)[i])
#   abline(v = sums[, i], col = 2)
# }
# title(xlab = "Depth", line = 0, outer = TRUE, font = 2)
# title(ylab = "Count", line = 0, outer = TRUE, font = 2)
