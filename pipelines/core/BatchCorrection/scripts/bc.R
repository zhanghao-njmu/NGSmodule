#!/usr/bin/env Rscript
# Batch correction with sva::ComBat_seq (preferred for raw counts) or
# fallback to limma::removeBatchEffect on log-CPM when sva is missing.
#
# Usage: Rscript bc.R <counts.matrix.tab> <sample_info.csv> <out_dir> <method>
#   method ∈ {combat_seq, limma}
#
# Produces in <out_dir>/:
#   counts.matrix.corrected.tab    same shape as input + corrected sample cols
#   pca_before.pdf                 PCA of input matrix coloured by Batch+Group
#   pca_after.pdf                  PCA of corrected matrix coloured by Batch+Group
#   summary.json                   method, n_genes, n_batches, group_levels

suppressPackageStartupMessages(library(stats))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: bc.R <counts.matrix.tab> <sample_info.csv> <out_dir> <method>")
}
counts_path <- args[1]
sample_info_path <- args[2]
out_dir <- args[3]
method <- args[4]

if (!file.exists(counts_path))      stop(sprintf("counts not found: %s", counts_path))
if (!file.exists(sample_info_path)) stop(sprintf("SampleInfo not found: %s", sample_info_path))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cm <- read.table(counts_path, header = TRUE, sep = "\t",
                 check.names = FALSE, stringsAsFactors = FALSE)
meta_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
sample_cols <- setdiff(colnames(cm), meta_cols)
counts <- as.matrix(cm[, sample_cols, drop = FALSE])
rownames(counts) <- cm$Geneid
storage.mode(counts) <- "integer"

si <- read.csv(sample_info_path, stringsAsFactors = FALSE)
si <- si[!duplicated(si$SampleID), , drop = FALSE]
rownames(si) <- si$SampleID
if (!"Batch" %in% colnames(si)) {
  stop("SampleInfoFile is missing a 'Batch' column — required for BatchCorrection")
}
si <- si[si$SampleID %in% sample_cols, , drop = FALSE]
counts <- counts[, si$SampleID, drop = FALSE]

batch <- factor(si$Batch)
group <- factor(if ("Group" %in% colnames(si)) si$Group else rep("default", nrow(si)))

if (length(levels(batch)) < 2) {
  stop(sprintf("Need ≥2 batches; got %d (%s)",
               length(levels(batch)), paste(levels(batch), collapse = ",")))
}

# --- PCA helper --------------------------------------------------------
pca_plot <- function(mat, file_path, title) {
  v <- apply(mat, 1, var)
  filt <- mat[v > 0, , drop = FALSE]
  if (nrow(filt) < 2 || ncol(filt) < 2) return(invisible(NULL))
  pc <- prcomp(t(filt), scale. = TRUE)
  px <- pc$x[, 1, drop = TRUE]
  py <- if (ncol(pc$x) >= 2) pc$x[, 2] else rep(0, length(px))
  pdf(file_path, width = 7, height = 5); on.exit(dev.off(), add = TRUE)
  par(mar = c(5, 5, 3, 8), xpd = TRUE)
  cols <- rainbow(length(levels(batch)))[as.integer(batch)]
  pchs <- (15 + as.integer(group)) %% 25 + 1
  plot(px, py, col = cols, pch = pchs, cex = 1.6,
       xlab = "PC1", ylab = "PC2", main = title)
  text(px, py, colnames(mat), cex = 0.6, pos = 3)
  legend("topright", inset = c(-0.25, 0), legend = c(
    paste0("Batch ", levels(batch)),
    paste0("Group ", levels(group))
  ), col = c(rainbow(length(levels(batch))), rep("black", length(levels(group)))),
     pch = c(rep(19, length(levels(batch))), (15 + seq_along(levels(group))) %% 25 + 1),
     bty = "n", cex = 0.7)
}

pca_plot(log1p(counts), file.path(out_dir, "pca_before.pdf"),
         "PCA before correction (log1p)")

# --- Run chosen method ------------------------------------------------
corrected <- counts
method <- tolower(method)

if (method == "combat_seq") {
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("method=combat_seq requires the 'sva' package; install via BiocManager::install('sva')")
  }
  message("Running sva::ComBat_seq on raw counts")
  sv <- sva::ComBat_seq(counts = counts, batch = batch,
                        group = if (length(levels(group)) >= 2) group else NULL)
  corrected <- as.matrix(sv)
  storage.mode(corrected) <- "integer"
} else if (method == "limma") {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("method=limma requires the 'limma' package; install via BiocManager::install('limma')")
  }
  message("Running limma::removeBatchEffect on log1p counts (fallback)")
  log_counts <- log1p(counts)
  if (length(levels(group)) >= 2) {
    design <- model.matrix(~ group)
    log_corrected <- limma::removeBatchEffect(log_counts, batch = batch, design = design)
  } else {
    log_corrected <- limma::removeBatchEffect(log_counts, batch = batch)
  }
  # Round back to integer counts (clamped at 0)
  corrected <- pmax(0L, round(expm1(log_corrected)))
  storage.mode(corrected) <- "integer"
} else {
  stop(sprintf("Unknown method '%s'; valid: combat_seq, limma", method))
}

pca_plot(log1p(corrected), file.path(out_dir, "pca_after.pdf"),
         "PCA after correction (log1p)")

# --- Write corrected matrix in MergeCounts schema --------------------
out <- cm
for (s in sample_cols) out[[s]] <- corrected[, s]
write.table(out, file.path(out_dir, "counts.matrix.corrected.tab"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Summary JSON -----------------------------------------------------
summary <- list(
  method = method,
  n_genes = nrow(counts),
  n_samples = ncol(counts),
  n_batches = length(levels(batch)),
  batches = paste(levels(batch), collapse = ","),
  group_levels = paste(levels(group), collapse = ",")
)
to_json <- function(x) {
  if (is.list(x)) {
    items <- sprintf('  "%s": %s', names(x), sapply(x, to_json))
    paste0("{\n", paste(items, collapse = ",\n"), "\n}")
  } else if (is.numeric(x) && length(x) == 1) {
    if (is.na(x)) "null" else as.character(x)
  } else {
    sprintf('"%s"', x)
  }
}
writeLines(to_json(summary), con = file.path(out_dir, "summary.json"))

cat(sprintf("BatchCorrection: method=%s, %d genes × %d samples, %d batches\n",
            method, nrow(corrected), ncol(corrected), length(levels(batch))))
