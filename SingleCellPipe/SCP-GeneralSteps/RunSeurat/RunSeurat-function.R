db_scDblFinder <- function(srt, ...) {
  if (class(srt) != "Seurat") {
    stop("'srt' is not a Seurat object.",
      call. = FALSE
    )
  }
  db_rate <- ncol(srt) / 1000 * 0.01
  require(scDblFinder)
  sce <- as.SingleCellExperiment(srt, assay = "RNA")
  sce <- scDblFinder(sce, dbr = db_rate)
  srt[["db.scDblFinder_score"]] <- sce[["scDblFinder.score"]]
  srt[["db.scDblFinder_class"]] <- sce[["scDblFinder.class"]]
  db_out <- colnames(srt)[srt[["scDblFinder_class"]] == "doublet"]
  return(list(srt = srt, db_out = db_out))
}

db_scds <- function(srt, method = "hybrid", ...) {
  if (class(srt) != "Seurat") {
    stop("'srt' is not a Seurat object.",
      call. = FALSE
    )
  }
  db_rate <- ncol(srt) / 1000 * 0.01
  require(scds)
  sce <- as.SingleCellExperiment(srt, assay = "RNA")
  sce <- cxds_bcds_hybrid(sce)
  srt[["db.cxds_score"]] <- sce[["cxds_score"]]
  srt[["db.bcds_score"]] <- sce[["bcds_score"]]
  srt[["db.hybrid_score"]] <- sce[["hybrid_score"]]
  ntop <- ceiling(db_rate * ncol(sce))
  db_out <- names(sort(srt[[paste0(method, "_score"), drop = T]], decreasing = T)[1:ntop])
  srt[[paste0("db.",method, "_class")]] <- "singlet"
  srt[[paste0("db.",method, "_class")]][db_out, ] <- "doublet"
  return(list(srt = srt, db_out = db_out))
}

db_Scrublet <- function(srt, ...) {
  if (class(srt) != "Seurat") {
    stop("'srt' is not a Seurat object.",
      call. = FALSE
    )
  }
  db_rate <- ncol(srt) / 1000 * 0.01
  require(reticulate)
  scr <- reticulate::import("scrublet")
  raw_counts <- r_to_py(t(as.matrix(GetAssayData(object = srt, assay = "RNA", slot = "counts"))))
  scrub <- scr$Scrublet(raw_counts, expected_doublet_rate = db_rate)
  res <- scrub$scrub_doublets()
  doublet_scores <- res[[1]]
  predicted_doublets <- res[[2]]

  srt[["db.Scrublet_score"]] <- doublet_scores
  srt[["db.Scrublet_class"]] <- sapply(predicted_doublets, function(i) {
    switch(as.character(i),
      "FALSE" = "singlet", "TRUE" = "doublet"
    )
  })
  db_out <- colnames(srt)[srt[["db.Scrublet_class"]] == "doublet"]
  return(list(srt = srt, db_out = db_out))
}

db_DoubletDetection <- function(srt, ...) {
  if (class(srt) != "Seurat") {
    stop("'srt' is not a Seurat object.",
      call. = FALSE
    )
  }
  db_rate <- ncol(srt) / 1000 * 0.01
  require(reticulate)
  raw_counts <- r_to_py(t(as.matrix(GetAssayData(object = srt, assay = "RNA", slot = "counts"))))
  doubletdetection <- reticulate::import("doubletdetection")
  clf <- doubletdetection$BoostClassifier()
  labels <- clf$fit(raw_counts)$predict()
  scores <- clf$doublet_score()

  srt[["db.DoubletDetection_score"]] <- scores
  srt[["db.DoubletDetection_class"]] <- sapply(labels, function(i) {
    switch(as.character(i),
      "0" = "singlet", "1" = "doublet"
    )
  })
  db_out <- colnames(srt)[srt[["db.DoubletDetection_class"]] == "doublet"]
  return(list(srt = srt, db_out = db_out))
}

DoubletIdentification <- function(srt, db_method = "scDblFinder") {
  if (class(srt) != "Seurat") {
    stop("'srt' is not a Seurat object.",
      call. = FALSE
    )
  }
  if (db_method %in% c("scDblFinder", "Scrublet", "DoubletDetection", "scds_cxds", "scds_bcds", "scds_hybrid")) {
    methods <- unlist(strsplit(db_method, "_"))
    method1 <- methods[1]
    method2 <- methods[2]
    args1 <- c(mget(names(formals()), sys.frame(sys.nframe())), method = method2)
    args2 <- c(as.list(match.call()), method = method2)
    for (n in names(args2)) {
      args1[[n]] <- args2[[n]]
    }
    tryCatch(expr = {
      res <- base::do.call(
        what = paste0("db_", method1),
        args = args1
      )
    }, error = function(e) {
      message(e)
      res <- NA
    })
    return(res)
  } else {
    stop(paste(db_method, "is not a suppoted doublet identification method!"),
      call. = FALSE
    )
  }
}


# aplotGrob function from aplot
aplotGrob <- function(x) {
  mp <- x$plotlist[[1]]
  if (length(x$plotlist) == 1) {
    return(ggplotGrob(mp))
  }

  for (i in x$layout[, x$main_col]) {
    if (is.na(i)) next
    if (i == 1) next
    x$plotlist[[i]] <- suppressMessages(x$plotlist[[i]] + xlim2(mp))
  }
  for (i in x$layout[x$main_row, ]) {
    if (is.na(i)) next
    if (i == 1) next
    x$plotlist[[i]] <- suppressMessages(x$plotlist[[i]] + ylim2(mp))
  }

  idx <- as.vector(x$layout)
  idx[is.na(idx)] <- x$n + 1
  x$plotlist[[x$n + 1]] <- ggplot() +
    theme_void() # plot_spacer()
  plotlist <- x$plotlist[idx]

  pp <- plotlist[[1]] + theme_no_margin()
  for (i in 2:length(plotlist)) {
    pp <- pp + (plotlist[[i]] + theme_no_margin())
  }

  res <- pp + plot_layout(
    byrow = F, ncol = ncol(x$layout),
    widths = x$width,
    heights = x$height,
    guides = "collect"
  )
  patchworkGrob(res)
}
