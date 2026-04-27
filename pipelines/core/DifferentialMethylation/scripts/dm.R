#!/usr/bin/env Rscript
# Differential methylation: pairwise CpG-level test between two groups.
#
# Usage:
#   Rscript dm.R <cov_glob> <sample_info.csv> <out_dir> \
#                <group_a> <group_b> <min_cov> <max_qval> <min_meth_diff>
#
# Inputs:
#   cov_glob       — glob pattern that resolves to per-sample bismark.cov.gz
#                     files (one per sample); the script discovers them and
#                     pairs against sample_info[SampleID].
#   sample_info    — CSV with at least SampleID + Group columns.
#   out_dir        — output directory.
#   group_a, group_b — Group labels to compare (matched against Group column).
#   min_cov        — drop CpG positions with <min_cov reads in any sample.
#   max_qval       — FDR threshold for "significant".
#   min_meth_diff  — minimum |% methylation difference| for "significant".
#
# Outputs:
#   dm_results_<b>_vs_<a>.tab     all CpG positions tested
#   dm_significant_<b>_vs_<a>.tab filtered (qval<max_qval & |diff|>=min_meth_diff)
#   summary.json                  thresholds + counts of hyper/hypo-methylated CpGs
#
# Engine: pure base R (Welch's t-test on logit-transformed methylation
# rates per CpG). Users can swap in methylKit via --dm_script.

suppressPackageStartupMessages({
  library(stats)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop("Usage: dm.R <cov_glob> <sample_info.csv> <out_dir> <group_a> <group_b> <min_cov> <max_qval> <min_meth_diff>")
}
cov_glob       <- args[1]
sample_info_pt <- args[2]
out_dir        <- args[3]
group_a        <- args[4]
group_b        <- args[5]
min_cov        <- as.integer(args[6])
max_qval       <- as.numeric(args[7])
min_meth_diff  <- as.numeric(args[8])

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Resolve coverage files via Sys.glob
cov_files <- Sys.glob(cov_glob)
if (length(cov_files) == 0) {
  stop(sprintf("no .bismark.cov.gz files match: %s", cov_glob))
}

# Read sample info
si <- read.csv(sample_info_pt, stringsAsFactors = FALSE)
si <- si[!duplicated(si$SampleID), c("SampleID", "Group"), drop = FALSE]
si <- si[si$Group %in% c(group_a, group_b), , drop = FALSE]
if (nrow(si) < 2) {
  stop(sprintf("Need ≥2 samples in groups [%s, %s]; got %d", group_a, group_b, nrow(si)))
}

# Pair coverage files to samples (by basename containing SampleID)
read_cov <- function(path) {
  # Format: chr<TAB>start<TAB>end<TAB>%meth<TAB>numCs<TAB>numTs
  df <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                   col.names = c("chr", "start", "end", "pct_meth", "numCs", "numTs"))
  df$pos_id <- paste(df$chr, df$start, sep = ":")
  df$cov <- df$numCs + df$numTs
  df <- df[df$cov >= min_cov, , drop = FALSE]
  df
}

per_sample <- list()
for (i in seq_len(nrow(si))) {
  sid <- si$SampleID[i]
  match_idx <- grep(sid, basename(cov_files), fixed = TRUE)
  if (length(match_idx) == 0) {
    stop(sprintf("no cov file for sample %s in: %s", sid, paste(cov_files, collapse = ", ")))
  }
  cov_path <- cov_files[match_idx[1]]
  message(sprintf("[%s] %s", sid, cov_path))
  per_sample[[sid]] <- read_cov(cov_path)
}

# Compute the intersection of CpG positions covered in every sample
all_pos <- Reduce(intersect, lapply(per_sample, function(d) d$pos_id))
message(sprintf("CpG positions covered in all %d samples: %d", length(per_sample), length(all_pos)))

if (length(all_pos) == 0) {
  stop("No CpG positions covered in all samples after filtering.")
}

# Build a pos × sample matrix of % methylation
meth_mat <- matrix(NA_real_, nrow = length(all_pos), ncol = nrow(si),
                   dimnames = list(all_pos, si$SampleID))
for (sid in si$SampleID) {
  d <- per_sample[[sid]]
  d <- d[d$pos_id %in% all_pos, , drop = FALSE]
  rownames(d) <- d$pos_id
  meth_mat[all_pos, sid] <- d[all_pos, "pct_meth"]
}

# Logit transform with Laplace smoothing
logit <- function(p, eps = 0.5) log2((p + eps) / (100 - p + eps))
logit_mat <- logit(meth_mat)

# Per-CpG Welch's t-test
group_vec <- si$Group
idx_a <- which(group_vec == group_a)
idx_b <- which(group_vec == group_b)
if (length(idx_a) < 2 || length(idx_b) < 2) {
  message(sprintf("WARNING: group %s has %d samples, group %s has %d. T-test needs ≥2 per group.",
                  group_a, length(idx_a), group_b, length(idx_b)))
}

n_pos <- nrow(logit_mat)
pvals <- rep(NA_real_, n_pos)
mean_a <- rep(NA_real_, n_pos)
mean_b <- rep(NA_real_, n_pos)
diff_pct <- rep(NA_real_, n_pos)

for (i in seq_len(n_pos)) {
  a <- logit_mat[i, idx_a]
  b <- logit_mat[i, idx_b]
  mean_a[i] <- mean(meth_mat[i, idx_a], na.rm = TRUE)
  mean_b[i] <- mean(meth_mat[i, idx_b], na.rm = TRUE)
  diff_pct[i] <- mean_b[i] - mean_a[i]
  if (length(idx_a) >= 2 && length(idx_b) >= 2 && var(a, na.rm = TRUE) + var(b, na.rm = TRUE) > 0) {
    pvals[i] <- tryCatch(
      stats::t.test(a, b, var.equal = FALSE)$p.value,
      error = function(e) NA_real_
    )
  }
}

qvals <- p.adjust(pvals, method = "BH")

# Build the table
parts <- strsplit(all_pos, ":", fixed = TRUE)
chrs   <- vapply(parts, `[`, character(1), 1)
starts <- as.integer(vapply(parts, `[`, character(1), 2))

res <- data.frame(
  chr = chrs, start = starts,
  mean_meth_a = mean_a, mean_meth_b = mean_b,
  meth_diff = diff_pct,
  pvalue = pvals, qvalue = qvals,
  stringsAsFactors = FALSE
)
colnames(res)[3:4] <- c(paste0("mean_meth_", group_a), paste0("mean_meth_", group_b))
res <- res[order(res$qvalue, na.last = TRUE), , drop = FALSE]

out_full <- file.path(out_dir, sprintf("dm_results_%s_vs_%s.tab", group_b, group_a))
write.table(res, file = out_full, sep = "\t", quote = FALSE, row.names = FALSE)

sig <- res[!is.na(res$qvalue) & res$qvalue < max_qval &
           abs(res$meth_diff) >= min_meth_diff, , drop = FALSE]
out_sig <- file.path(out_dir, sprintf("dm_significant_%s_vs_%s.tab", group_b, group_a))
write.table(sig, file = out_sig, sep = "\t", quote = FALSE, row.names = FALSE)

n_hyper <- sum(sig$meth_diff > 0)
n_hypo  <- sum(sig$meth_diff < 0)

# JSON summary (no jsonlite dependency)
summary <- list(
  comparison = sprintf("%s_vs_%s", group_b, group_a),
  reference  = group_a, test_group = group_b,
  thresholds = list(min_cov = min_cov, max_qval = max_qval, min_meth_diff = min_meth_diff),
  n_positions_tested = nrow(res),
  n_significant      = nrow(sig),
  n_hyper            = n_hyper,
  n_hypo             = n_hypo
)
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

cat(sprintf("DifferentialMethylation: %d hyper, %d hypo (qval<%g, |meth_diff|>=%g%%)\n",
            n_hyper, n_hypo, max_qval, min_meth_diff))
