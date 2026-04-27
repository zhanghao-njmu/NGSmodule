#!/usr/bin/env Rscript
# Self-contained project-level QC for a counts matrix.
#
# Usage: Rscript qc.R <counts_matrix.tab> <sample_info.csv> <out_dir>
#
# Produces in <out_dir>/:
#   library_sizes.tab           per-sample total counts
#   detection_rates.tab         per-sample fraction of genes with count > 0
#   sample_correlation.tab      Spearman correlation between samples (matrix)
#   library_sizes.pdf           bar plot of library sizes
#   sample_correlation.pdf      heatmap of correlation matrix
#   pca.tab                     PCA coordinates (PC1/PC2) per sample
#   pca.pdf                     scatterplot of PC1 vs PC2 with Group colour
#
# Uses only base R + grDevices — no Bioconductor dependencies. Heavy
# packages (limma/edgeR/ComplexHeatmap) from the legacy postQuantificationQC.R
# can be added by users in their own R script if needed.

suppressPackageStartupMessages({
  library(stats)
  library(grDevices)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: qc.R <counts.matrix.tab> <sample_info.csv> <out_dir>")
}
counts_path <- args[1]
sample_info_path <- args[2]
out_dir <- args[3]

if (!file.exists(counts_path)) {
  stop(sprintf("Counts matrix not found: %s", counts_path))
}
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Read counts. The MergeCounts pipeline emits:
#   Geneid Chr Start End Strand Length sample1 sample2 ...
counts_full <- read.table(counts_path, header = TRUE, sep = "\t",
                          check.names = FALSE, stringsAsFactors = FALSE)
meta_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
sample_cols <- setdiff(colnames(counts_full), meta_cols)
if (length(sample_cols) == 0) {
  stop("No sample columns found in counts matrix.")
}
counts <- as.matrix(counts_full[, sample_cols, drop = FALSE])
rownames(counts) <- counts_full$Geneid
storage.mode(counts) <- "double"

# Optional sample_info (RunID, SampleID, Layout, Group, Batch ...)
group_map <- setNames(rep("default", length(sample_cols)), sample_cols)
if (file.exists(sample_info_path)) {
  si <- read.csv(sample_info_path, stringsAsFactors = FALSE)
  if (all(c("SampleID", "Group") %in% colnames(si))) {
    si <- si[!duplicated(si$SampleID), c("SampleID", "Group")]
    matches <- match(sample_cols, si$SampleID)
    keep <- !is.na(matches)
    if (any(keep)) group_map[sample_cols[keep]] <- si$Group[matches[keep]]
  }
}

# 1. Library sizes
lib <- data.frame(
  SampleID = sample_cols,
  Group = group_map[sample_cols],
  TotalCount = colSums(counts),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write.table(lib, file.path(out_dir, "library_sizes.tab"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 2. Detection rate (fraction of genes with count > 0)
det <- data.frame(
  SampleID = sample_cols,
  GenesDetected = colSums(counts > 0),
  TotalGenes = nrow(counts),
  DetectionRate = colSums(counts > 0) / nrow(counts),
  stringsAsFactors = FALSE
)
write.table(det, file.path(out_dir, "detection_rates.tab"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 3. Sample-sample correlation (Spearman, on log1p)
if (length(sample_cols) >= 2) {
  log_counts <- log1p(counts)
  cor_mat <- suppressWarnings(cor(log_counts, method = "spearman"))
  cor_df <- data.frame(SampleID = rownames(cor_mat), cor_mat,
                       check.names = FALSE)
  write.table(cor_df, file.path(out_dir, "sample_correlation.tab"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# 4. PCA (only on genes with non-zero variance; works for ≥2 samples)
if (length(sample_cols) >= 2) {
  log_counts <- log1p(counts)
  v <- apply(log_counts, 1, var)
  log_filt <- log_counts[v > 0, , drop = FALSE]
  if (nrow(log_filt) > 1) {
    pca <- prcomp(t(log_filt), scale. = TRUE)
    pcs <- pca$x[, 1:min(2, ncol(pca$x)), drop = FALSE]
    pca_df <- data.frame(
      SampleID = sample_cols,
      Group = group_map[sample_cols],
      pcs,
      stringsAsFactors = FALSE
    )
    write.table(pca_df, file.path(out_dir, "pca.tab"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

# Plots — base R only, so any R install works.
plot_lib_sizes <- function() {
  pdf(file.path(out_dir, "library_sizes.pdf"), width = max(6, length(sample_cols) * 0.5), height = 5)
  par(mar = c(8, 5, 3, 2))
  groups <- factor(lib$Group)
  cols <- rainbow(length(levels(groups)))[as.integer(groups)]
  barplot(lib$TotalCount, names.arg = lib$SampleID, las = 2,
          col = cols, ylab = "Total counts", main = "Per-sample library size")
  dev.off()
}
try(plot_lib_sizes(), silent = TRUE)

if (exists("cor_mat") && is.matrix(cor_mat)) {
  try({
    pdf(file.path(out_dir, "sample_correlation.pdf"),
        width = max(6, ncol(cor_mat) * 0.6), height = max(5, nrow(cor_mat) * 0.6))
    par(mar = c(8, 8, 2, 2))
    image(seq_len(ncol(cor_mat)), seq_len(nrow(cor_mat)), t(cor_mat[nrow(cor_mat):1, ]),
          axes = FALSE, xlab = "", ylab = "", col = hcl.colors(50, "Blues 3", rev = TRUE),
          main = "Spearman correlation (log1p counts)")
    axis(1, at = seq_len(ncol(cor_mat)), labels = colnames(cor_mat), las = 2)
    axis(2, at = seq_len(nrow(cor_mat)), labels = rev(rownames(cor_mat)), las = 2)
    dev.off()
  }, silent = TRUE)
}

if (exists("pca_df")) {
  try({
    pdf(file.path(out_dir, "pca.pdf"), width = 6, height = 5)
    par(mar = c(5, 5, 3, 2))
    groups <- factor(pca_df$Group)
    cols <- rainbow(length(levels(groups)))[as.integer(groups)]
    plot(pca_df$PC1, if ("PC2" %in% colnames(pca_df)) pca_df$PC2 else rep(0, nrow(pca_df)),
         pch = 19, col = cols, xlab = "PC1",
         ylab = if ("PC2" %in% colnames(pca_df)) "PC2" else "",
         main = "PCA (log1p counts, scaled)")
    text(pca_df$PC1, if ("PC2" %in% colnames(pca_df)) pca_df$PC2 else rep(0, nrow(pca_df)),
         pca_df$SampleID, cex = 0.7, pos = 3)
    if (length(levels(groups)) > 1) {
      legend("topright", legend = levels(groups), fill = rainbow(length(levels(groups))), bty = "n")
    }
    dev.off()
  }, silent = TRUE)
}

cat(sprintf("postQuantificationQC: wrote %d output files in %s\n",
            length(list.files(out_dir)), out_dir))
