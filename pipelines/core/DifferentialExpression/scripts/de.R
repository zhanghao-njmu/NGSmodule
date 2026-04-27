#!/usr/bin/env Rscript
# Differential expression with edgeR (quasi-likelihood F-tests).
#
# Usage:
#   Rscript de.R <counts.matrix.tab> <sample_info.csv> <out_dir> \
#                <max_padj> <min_log2fc> <min_count> <group_a> <group_b>
#
# Produces in <out_dir>/:
#   de_results_<group_a>_vs_<group_b>.tab    full table (all genes)
#   de_significant_<a>_vs_<b>.tab            filtered (padj < max_padj, |logFC| >= min_log2fc)
#   ma_<a>_vs_<b>.pdf                        MA plot
#   volcano_<a>_vs_<b>.pdf                   volcano plot
#   summary.json                             counts of up/down/ns + thresholds

suppressPackageStartupMessages({
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("edgeR is required. Install with: BiocManager::install('edgeR')")
  }
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop("Usage: de.R <counts.matrix.tab> <sample_info.csv> <out_dir> <max_padj> <min_log2fc> <min_count> <group_a> <group_b>")
}
counts_path     <- args[1]
sample_info_path <- args[2]
out_dir         <- args[3]
max_padj        <- as.numeric(args[4])
min_log2fc      <- as.numeric(args[5])
min_count       <- as.numeric(args[6])
group_a         <- args[7]
group_b         <- args[8]

if (!file.exists(counts_path))     stop(sprintf("counts matrix not found: %s", counts_path))
if (!file.exists(sample_info_path)) stop(sprintf("SampleInfo not found: %s", sample_info_path))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Read counts (MergeCounts format: Geneid Chr Start End Strand Length sample1 sample2 ...)
cm <- read.table(counts_path, header = TRUE, sep = "\t",
                 check.names = FALSE, stringsAsFactors = FALSE)
meta_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
sample_cols <- setdiff(colnames(cm), meta_cols)
counts <- as.matrix(cm[, sample_cols, drop = FALSE])
rownames(counts) <- cm$Geneid
storage.mode(counts) <- "integer"

# Read sample info
si <- read.csv(sample_info_path, stringsAsFactors = FALSE)
si <- si[!duplicated(si$SampleID), c("SampleID", "Group"), drop = FALSE]
rownames(si) <- si$SampleID

# Restrict to samples actually present in the counts matrix
si <- si[si$SampleID %in% sample_cols, , drop = FALSE]
counts <- counts[, si$SampleID, drop = FALSE]

# Restrict to the two groups being compared
keep <- si$Group %in% c(group_a, group_b)
si <- si[keep, , drop = FALSE]
counts <- counts[, si$SampleID, drop = FALSE]

if (nrow(si) < 2) {
  stop(sprintf("Need ≥2 samples in groups [%s, %s]; got %d", group_a, group_b, nrow(si)))
}

# Build edgeR DGEList; sample groups as factor with $group_a as the reference.
group <- factor(si$Group, levels = c(group_a, group_b))
y <- DGEList(counts = counts, group = group)

# Drop low-count genes (CPM-based filter equivalent to filterByExpr but explicit)
keep_genes <- rowSums(counts >= min_count) >= min(table(group))
y <- y[keep_genes, , keep.lib.sizes = FALSE]

# TMM normalisation + dispersion estimation
y <- calcNormFactors(y)
design <- model.matrix(~ group)
y <- estimateDisp(y, design)

# QL F-test
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)  # group_b vs group_a
res <- topTags(qlf, n = Inf)$table

# Order by FDR ascending; write full and filtered tables
res$Geneid <- rownames(res)
out_full <- file.path(out_dir, sprintf("de_results_%s_vs_%s.tab", group_b, group_a))
write.table(res[, c("Geneid", setdiff(colnames(res), "Geneid"))],
            file = out_full, sep = "\t", quote = FALSE, row.names = FALSE)

sig <- res[res$FDR < max_padj & abs(res$logFC) >= min_log2fc, , drop = FALSE]
out_sig <- file.path(out_dir, sprintf("de_significant_%s_vs_%s.tab", group_b, group_a))
write.table(sig[, c("Geneid", setdiff(colnames(sig), "Geneid"))],
            file = out_sig, sep = "\t", quote = FALSE, row.names = FALSE)

n_up <- sum(sig$logFC >  0)
n_dn <- sum(sig$logFC <  0)

# MA plot
try({
  pdf(file.path(out_dir, sprintf("ma_%s_vs_%s.pdf", group_b, group_a)), width = 7, height = 5)
  par(mar = c(5, 5, 3, 2))
  is_sig <- res$FDR < max_padj & abs(res$logFC) >= min_log2fc
  plot(res$logCPM, res$logFC,
       col = ifelse(is_sig, ifelse(res$logFC > 0, "red", "blue"), "grey60"),
       pch = 16, cex = 0.5,
       xlab = "log2 CPM", ylab = sprintf("log2 FC (%s / %s)", group_b, group_a),
       main = sprintf("MA plot — %s vs %s (n_up=%d, n_down=%d)", group_b, group_a, n_up, n_dn))
  abline(h = c(-min_log2fc, 0, min_log2fc), col = "black", lty = c(2, 1, 2))
  legend("topright", legend = c("up", "down", "ns"), pch = 16,
         col = c("red", "blue", "grey60"), bty = "n")
  dev.off()
}, silent = TRUE)

# Volcano
try({
  pdf(file.path(out_dir, sprintf("volcano_%s_vs_%s.pdf", group_b, group_a)), width = 7, height = 5)
  par(mar = c(5, 5, 3, 2))
  is_sig <- res$FDR < max_padj & abs(res$logFC) >= min_log2fc
  plot(res$logFC, -log10(res$FDR + 1e-300),
       col = ifelse(is_sig, ifelse(res$logFC > 0, "red", "blue"), "grey60"),
       pch = 16, cex = 0.5,
       xlab = sprintf("log2 FC (%s / %s)", group_b, group_a), ylab = "-log10 FDR",
       main = sprintf("Volcano — %s vs %s", group_b, group_a))
  abline(v = c(-min_log2fc, min_log2fc), col = "black", lty = 2)
  abline(h = -log10(max_padj),           col = "black", lty = 2)
  legend("topright", legend = c("up", "down", "ns"), pch = 16,
         col = c("red", "blue", "grey60"), bty = "n")
  dev.off()
}, silent = TRUE)

# Machine-readable summary
summary <- list(
  comparison = sprintf("%s_vs_%s", group_b, group_a),
  reference  = group_a,
  test_group = group_b,
  thresholds = list(max_padj = max_padj, min_log2fc = min_log2fc, min_count = min_count),
  n_genes_tested = nrow(res),
  n_significant  = nrow(sig),
  n_up           = n_up,
  n_down         = n_dn
)
# Write JSON without dependencies
to_json <- function(x, indent = 0) {
  pad <- paste(rep("  ", indent), collapse = "")
  if (is.list(x)) {
    items <- mapply(function(k, v) sprintf('%s"%s": %s', paste0(pad, "  "), k, to_json(v, indent + 1)),
                    names(x), x, SIMPLIFY = TRUE)
    sprintf("{\n%s\n%s}", paste(items, collapse = ",\n"), pad)
  } else if (is.numeric(x) && length(x) == 1) {
    if (is.na(x)) "null" else as.character(x)
  } else {
    sprintf('"%s"', x)
  }
}
writeLines(to_json(summary), con = file.path(out_dir, "summary.json"))

cat(sprintf("DifferentialExpression: %d up, %d down (FDR<%g, |log2FC|>=%g)\n",
            n_up, n_dn, max_padj, min_log2fc))
