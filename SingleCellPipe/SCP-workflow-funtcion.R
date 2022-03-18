db_scDblFinder <- function(srt, db_rate = ncol(srt) / 1000 * 0.01, ...) {
  if (class(srt) != "Seurat") {
    stop("'srt' is not a Seurat object.")
  }
  library(scDblFinder)
  sce <- as.SingleCellExperiment(srt, assay = "RNA")
  sce <- scDblFinder(sce, dbr = db_rate, verbose = FALSE, ...)
  srt[["db.scDblFinder_score"]] <- sce[["scDblFinder.score"]]
  srt[["db.scDblFinder_class"]] <- sce[["scDblFinder.class"]]
  return(srt)
}

db_scds <- function(srt, db_rate = ncol(srt) / 1000 * 0.01, method = "hybrid", ...) {
  if (class(srt) != "Seurat") {
    stop("'srt' is not a Seurat object.")
  }
  library(scds)
  sce <- as.SingleCellExperiment(srt, assay = "RNA")
  sce <- cxds_bcds_hybrid(sce)
  srt[["db.cxds_score"]] <- sce[["cxds_score"]]
  srt[["db.bcds_score"]] <- sce[["bcds_score"]]
  srt[["db.hybrid_score"]] <- sce[["hybrid_score"]]
  ntop <- ceiling(db_rate * ncol(sce))
  db_out <- names(sort(srt[[paste0("db.", method, "_score"), drop = T]], decreasing = T)[1:ntop])
  srt[[paste0("db.scds_", method, "_class")]] <- "singlet"
  srt[[paste0("db.scds_", method, "_class")]][db_out, ] <- "doublet"
  return(srt)
}

db_Scrublet <- function(srt, db_rate = ncol(srt) / 1000 * 0.01, ...) {
  if (class(srt) != "Seurat") {
    stop("'srt' is not a Seurat object.")
  }
  library(reticulate)
  scr <- reticulate::import("scrublet")
  raw_counts <- r_to_py(t(as.matrix(GetAssayData(object = srt, assay = "RNA", slot = "counts"))))
  scrub <- scr$Scrublet(raw_counts, expected_doublet_rate = db_rate)
  res <- scrub$scrub_doublets()
  doublet_scores <- res[[1]]
  predicted_doublets <- res[[2]]

  srt[["db.Scrublet_score"]] <- doublet_scores
  srt[["db.Scrublet_class"]] <- sapply(predicted_doublets, function(i) {
    switch(as.character(i),
      "FALSE" = "singlet",
      "TRUE" = "doublet"
    )
  })
  return(srt)
}

db_DoubletDetection <- function(srt, db_rate = ncol(srt) / 1000 * 0.01, ...) {
  if (class(srt) != "Seurat") {
    stop("'srt' is not a Seurat object.")
  }
  library(reticulate)
  doubletdetection <- reticulate::import("doubletdetection")
  counts <- GetAssayData(object = srt, assay = "RNA", slot = "counts")
  clf <- doubletdetection$BoostClassifier(
    n_iters = as.integer(5),
    use_phenograph = FALSE,
    standard_scaling = TRUE
  )
  labels <- clf$fit(Matrix::t(counts))$predict()
  scores <- clf$doublet_score()

  srt[["db.DoubletDetection_score"]] <- scores
  srt[["db.DoubletDetection_class"]] <- sapply(labels, function(i) {
    switch(as.character(i),
      "0" = "singlet",
      "1" = "doublet"
    )
  })
  return(srt)
}

RunDoubletCalling <- function(srt, db_rate = ncol(srt) / 1000 * 0.01, db_method = "scDblFinder", ...) {
  if (class(srt) != "Seurat") {
    stop("'srt' is not a Seurat object.")
  }
  if (db_method %in% c("scDblFinder", "Scrublet", "DoubletDetection", "scds_cxds", "scds_bcds", "scds_hybrid")) {
    message("Run doublet-calling with ", db_method, " method...")
    methods <- unlist(strsplit(db_method, "_"))
    method1 <- methods[1]
    method2 <- methods[2]
    if (is.na(method2)) {
      args1 <- mget(names(formals()), sys.frame(sys.nframe()))
      args2 <- as.list(match.call())
    } else {
      args1 <- c(mget(names(formals()), sys.frame(sys.nframe())), method = method2)
      args2 <- c(as.list(match.call()), method = method2)
    }
    for (n in names(args2)) {
      args1[[n]] <- args2[[n]]
    }
    args1 <- args1[!names(args1) %in% c("db_method", "...")]
    tryCatch(expr = {
      srt <- base::do.call(
        what = paste0("db_", method1),
        args = args1
      )
    }, error = function(e) {
      message(e)
    })
    return(srt)
  } else {
    stop(paste(db_method, "is not a suppoted doublet-calling method!"))
  }
}

CellQC <- function(srt, db_method = "scDblFinder", UMI_threshold = 3000, gene_threshold = 1000,
                   mito_threshold = 20, mito_pattern = c("MT-", "Mt-", "mt-"),
                   ribo_threshold = 50, ribo_pattern = c("RP[SL]\\d+(\\w+)$", "Rp[sl]\\d+(\\w+)$", "rp[sl]\\d+(\\w+)$"),
                   species = NULL, species_gene_prefix = NULL,
                   species_percent = 5, species_nCount = 5000,
                   seed = 11) {
  set.seed(seed)

  if (length(species) != length(species_gene_prefix)) {
    stop("Parameter error! 'species_gene_prefix' must be the same length as 'species'.")
  }

  ntotal <- ncol(srt)

  db_out <- c()
  if (!is.null(db_method)) {
    for (dbm in db_method) {
      srt <- RunDoubletCalling(srt = srt, db_method = dbm)
      db_out <- unique(c(db_out, colnames(srt)[srt[[paste0("db.", dbm, "_class"), drop = T]] == "doublet"]))
    }
  }

  srt[["percent.top_20"]] <- pct_counts_in_top_20_features <- addPerCellQC(as.SingleCellExperiment(srt, assay = "RNA"), percent_top = c(20))$percent.top_20
  srt[["log10_nFeature_RNA"]] <- log10_nFeature_RNA <- log10(srt[["nFeature_RNA", drop = T]])
  srt[["log10_nCount_RNA"]] <- log10_nCount_RNA <- log10(srt[["nCount_RNA", drop = T]])
  mod <- loess(log10_nFeature_RNA ~ log10_nCount_RNA)
  pred <- predict(mod, newdata = data.frame(log10_nCount_RNA = log10_nCount_RNA))
  srt[["featcount_dist"]] <- featcount_dist <- log10_nFeature_RNA - pred

  filters <- c(
    "log10_nCount_RNA:higher:2.5",
    "log10_nCount_RNA:lower:5",
    "log10_nFeature_RNA:higher:2.5",
    "log10_nFeature_RNA:lower:5",
    "pct_counts_in_top_20_features:both:5",
    "featcount_dist:both:5"
  )
  qc_out <- lapply(strsplit(filters, ":"), function(f) {
    colnames(srt)[isOutlier(get(f[1]),
      log = FALSE,
      nmads = as.numeric(f[3]),
      type = f[2]
    )]
  })
  qc_out <- table(unlist(qc_out))
  qc_out <- names(qc_out)[qc_out >= 2]

  for (n in 1:length(species)) {
    if (n > 0) {
      sp <- species[n]
      prefix <- species_gene_prefix[n]
      sp_genes <- rownames(srt)[grep(pattern = paste0("^", prefix), x = rownames(srt))]
      srt[[paste0(c("nCount_RNA", sp), collapse = ".")]] <- colSums(srt@assays$RNA@counts[sp_genes, ])
      srt[[paste0(c("nFeature_RNA", sp), collapse = ".")]] <- colSums(srt@assays$RNA@counts[sp_genes, ] > 0)
      srt[[paste0(c("percent.mito", sp), collapse = ".")]] <- PercentageFeatureSet(object = srt, pattern = paste0("(", paste0("^", prefix, "-*", mito_pattern), ")", collapse = "|"))
      srt[[paste0(c("percent.ribo", sp), collapse = ".")]] <- PercentageFeatureSet(object = srt, pattern = paste0("(", paste0("^", prefix, "-*", ribo_pattern), ")", collapse = "|"))
      srt[[paste0(c("percent.genome", sp), collapse = ".")]] <- PercentageFeatureSet(object = srt, pattern = paste0("^", prefix))
    }
  }

  umi_out <- colnames(srt)[srt[[paste0(c("nCount_RNA", species[1]), collapse = "."), drop = TRUE]] <= UMI_threshold]
  gene_out <- colnames(srt)[srt[[paste0(c("nFeature_RNA", species[1]), collapse = "."), drop = TRUE]] <= gene_threshold]
  mito_out <- colnames(srt)[srt[[paste0(c("percent.mito", species[1]), collapse = "."), drop = TRUE]] >= mito_threshold]
  ribo_out <- colnames(srt)[srt[[paste0(c("percent.ribo", species[1]), collapse = "."), drop = TRUE]] >= ribo_threshold]

  species_out <- c()
  if (length(species) >= 2) {
    for (n in 2:length(species)) {
      sp <- species[n]
      species_out <- c(species_out, colnames(srt)[srt[[paste0(c("percent.genome", sp), collapse = "."), drop = T]] >= species_percent])
      species_out <- unique(species_out)
    }
  }

  total_out <- unique(c(db_out, qc_out, umi_out, gene_out, mito_out, ribo_out, species_out))
  srt[["CellFilterng"]] <- ifelse(colnames(srt) %in% total_out, "Discarded", "Kept")
  srt[["CellFilterng"]] <- factor(srt[["CellFilterng", drop = TRUE]], levels = c("Kept", "Discarded"))

  cat(">>>", "Total cells:", ntotal, "\n")
  cat(">>>", "Cells which are filtered out:", length(total_out), "\n")
  cat("...", length(db_out), "potential doublets", "\n")
  cat("...", length(qc_out), "unqualified cells", "\n")
  cat("...", length(umi_out), "low-UMI cells", "\n")
  cat("...", length(gene_out), "low-gene cells", "\n")
  cat("...", length(mito_out), "high-mito cells", "\n")
  cat("...", length(ribo_out), "high-ribo cells", "\n")
  cat("...", length(species_out), "species-contaminated cells", "\n")
  cat(">>>", "Remained cells after filtering:", ntotal - length(total_out), "\n")
  for (out in c("db_out", "qc_out", "umi_out", "gene_out", "mito_out", "ribo_out", "species_out")) {
    srt[[out]] <- ifelse(colnames(srt) %in% get(out), "Discarded", "Kept")
    srt[[out]] <- factor(srt[[out, drop = TRUE]], levels = c("Kept", "Discarded"))
  }
  return(srt)
}

GeneConvert <- function(geneID, geneID_from_IDtype, geneID_to_IDtype, species_from, species_to, Ensembl_version = "current_release") {
  library(biomaRt)
  library(dplyr)
  library(httr)
  httr::set_config(httr::config(ssl_verifypeer = FALSE))

  species_from_split <- unlist(strsplit(species_from, split = "_"))
  species_to_split <- unlist(strsplit(species_to, split = "_"))
  species_from_simp <- paste0(tolower(substring(species_from_split[1], 1, 1)), species_from_split[2])
  species_to_simp <- paste0(tolower(substring(species_to_split[1], 1, 1)), species_to_split[2])
  geneID_from_IDtype <- sapply(geneID_from_IDtype, tolower)
  geneID_to_IDtype <- sapply(geneID_to_IDtype, tolower)

  if (geneID_from_IDtype == "symbol") {
    geneID_from_IDtype <- c("ensembl_symbol", "entrez_symbol")
  }
  from_IDtype <- sapply(geneID_from_IDtype, function(x) {
    switch(x,
      "ensembl_symbol" = "external_gene_name",
      "ensembl_id" = "ensembl_gene_id",
      "entrez_symbol" = "entrezgene_accession",
      "entrez_id" = "entrezgene_id"
    )
  })
  to_IDtype <- switch(geneID_to_IDtype,
    "symbol" = "external_gene_name",
    "entrez_symbol" = "external_gene_name",
    "ensembl_symbol" = "external_gene_name",
    "ensembl_id" = "ensembl_gene_id",
    "entrez_id" = "entrezgene_id"
  )

  message("Connect to the Ensembl archives...")
  archives <- NULL
  ntry <- 0
  while (is.null(archives)) {
    ntry <- ntry + 1
    archives <- tryCatch(expr = {
      listEnsemblArchives()
    }, error = function(e) {
      message("Get errors when connecting with EnsemblArchives...\nRetrying...")
      Sys.sleep(1)
      return(NULL)
    })
    if (ntry > 5) {
      stop("Stop connecting...")
    }
  }

  if (Ensembl_version == "current_release") {
    url <- archives[which(archives$current_release == "*"), "url"]
  } else {
    url <- archives[which(archives$version == Ensembl_version), "url"]
  }

  message("Connect to the biomart...")
  mart <- NULL
  ntry <- 0
  while (is.null(mart)) {
    ntry <- ntry + 1
    mart <- tryCatch(
      expr = {
        useMart("ensembl")
      },
      error = function(e) {
        message("Get errors when connecting with ensembl mart...\nRetrying...")
        Sys.sleep(1)
        return(NULL)
      }
    )
    if (ntry > 5) {
      stop("Stop connecting...")
    }
  }
  Datasets <- listDatasets(mart)
  dataset <- paste0(species_from_simp, "_gene_ensembl")
  if (!dataset %in% Datasets$dataset) {
    warning(paste0("Can not find the dataset for the species: ", species_from, " (", dataset, ")"))
    return(list(geneID_res = NULL, Datasets = Datasets, Attributes = NULL))
  }

  message("Use the ", Ensembl_version, " version of biomart...")
  mart <- NULL
  ntry <- 0
  while (is.null(mart)) {
    ntry <- ntry + 1
    mart <- tryCatch(
      expr = {
        useMart(biomart = "ensembl", dataset = dataset, host = url)
      },
      error = function(e) {
        message("Get errors when connecting with ensembl mart...\nRetrying...")
        Sys.sleep(1)
        return(NULL)
      }
    )
    if (ntry > 5) {
      stop("Stop connecting...")
    }
  }

  if (species_from == species_to) {
    to_attr <- to_IDtype
  } else {
    to_IDtype <- switch(to_IDtype,
      "external_gene_name" = "associated_gene_name",
      "ensembl_gene_id" = "ensembl_gene"
    )
    to_attr <- paste(species_to_simp, to_IDtype, sep = "_homolog_")
  }

  Attributes <- listAttributes(mart)
  geneID_res_list <- list()
  total <- length(geneID)
  message("Match and convert the geneID...\n")
  if (to_attr %in% Attributes$name) {
    if (species_from != species_to & !all(from_IDtype %in% c("external_gene_name", "ensembl_gene_id"))) {
      for (from_attr in from_IDtype) {
        if (length(geneID) > 0) {
          geneID_res1 <- getBM(
            mart = mart,
            attributes = c(from_attr, "ensembl_gene_id"),
            filters = from_attr,
            values = list(geneID)
          )
          if (nrow(geneID_res1) == 0) {
            next
          }
          colnames(geneID_res1) <- c("from_geneID", "ensembl_gene_id")
          from_name <- geneID_from_IDtype[which(from_IDtype == from_attr)]
          geneID_res1[, "from_IDtype"] <- from_name

          geneID_res2 <- getBM(
            mart = mart,
            attributes = c("ensembl_gene_id", to_attr),
            filters = "ensembl_gene_id",
            values = list(geneID_res1[, "ensembl_gene_id"])
          )
          if (nrow(geneID_res2) == 0) {
            next
          }
          colnames(geneID_res2) <- c("ensembl_gene_id", "to_geneID")
          geneID_res2[, "to_IDtype"] <- geneID_to_IDtype
          geneID_res_list[[from_attr]] <- merge(x = geneID_res1, y = geneID_res2, by = "ensembl_gene_id")
          ismap <- geneID %in% geneID_res_list[[from_attr]][, "from_geneID"]
          message(paste(sum(ismap), "genes mapped with", from_name))
          geneID <- geneID[!ismap]
        }
      }
      message(
        paste0(
          paste0(rep("=", 30), collapse = ""), "\n",
          total - length(geneID), " genes mapped\n",
          length(geneID), " genes unmapped"
        ), "\n",
        paste0(rep("=", 30), collapse = ""), "\n"
      )
      geneID_res <- bind_rows(geneID_res_list)
    } else {
      for (from_attr in from_IDtype) {
        if (length(geneID) > 0) {
          geneID_res1 <- getBM(
            mart = mart,
            attributes = c("ensembl_gene_id", from_attr, to_attr),
            filters = from_attr,
            values = list(geneID)
          )
          if (nrow(geneID_res1) == 0) {
            next
          }
          colnames(geneID_res1) <- c("ensembl_gene_id", "from_geneID", "to_geneID")
          from_name <- geneID_from_IDtype[which(from_IDtype == from_attr)]
          geneID_res1[, "from_IDtype"] <- from_name
          geneID_res1[, "to_IDtype"] <- geneID_to_IDtype
          geneID_res1 <- geneID_res1[, c("ensembl_gene_id", "from_geneID", "from_IDtype", "to_geneID", "to_IDtype")]
          geneID_res_list[[from_attr]] <- geneID_res1
          ismap <- geneID %in% geneID_res_list[[from_attr]][, "from_geneID"]
          message(paste(sum(ismap), "genes mapped with", from_name))
          geneID <- geneID[!ismap]
        }
      }
      message(
        paste0(
          paste0(rep("=", 30), collapse = ""), "\n",
          total - length(geneID), " genes mapped\n",
          length(geneID), " genes unmapped"
        ), "\n",
        paste0(rep("=", 30), collapse = ""), "\n"
      )
      geneID_res <- bind_rows(geneID_res_list)
    }
  } else {
    warning(paste0("Can not find the attribute for the species: ", species_from, " (", to_attr, ")"))
    return(list(geneID_res = NULL, geneID_sim = NULL, Datasets = Datasets, Attributes = Attributes))
  }
  if (nrow(geneID_res) == 0) {
    warning(paste0("No gene mapped"))
    return(list(geneID_res = NULL, geneID_sim = NULL, Datasets = Datasets, Attributes = Attributes))
  }
  geneID_sim <- geneID_res %>%
    group_by(from_geneID) %>%
    mutate(from_geneID = unique(from_geneID), to_geneID = list(unique(to_geneID[!to_geneID %in% c("", NA)])))
  geneID_sim <- unique(as.data.frame(geneID_sim[, c("from_geneID", "to_geneID")]))
  geneID_sim <- geneID_sim[sapply(geneID_sim$to_geneID, length) > 0, ]
  rownames(geneID_sim) <- geneID_sim[, "from_geneID"]

  return(list(geneID_res = geneID_res, geneID_sim = geneID_sim, Datasets = Datasets, Attributes = Attributes))
}

CC_GenePrefetch <- function(species) {
  cc_S_genes <- Seurat::cc.genes.updated.2019$s.genes
  cc_G2M_genes <- Seurat::cc.genes.updated.2019$g2m.genes
  if (species != "Homo_sapiens") {
    res <- GeneConvert(
      geneID = c(cc_S_genes, cc_G2M_genes),
      geneID_from_IDtype = "symbol",
      geneID_to_IDtype = "ensembl_symbol",
      species_from = "Homo_sapiens",
      species_to = species
    )
    genes <- res[["geneID_sim"]]
    cc_S_genes <- unlist(genes[cc_S_genes[cc_S_genes %in% rownames(genes)], "to_geneID"])
    cc_G2M_genes <- unlist(genes[cc_G2M_genes[cc_G2M_genes %in% rownames(genes)], "to_geneID"])
  }
  return(list(
    cc_S_genes = cc_S_genes,
    cc_G2M_genes = cc_G2M_genes
  ))
}

Check_srtList <- function(srtList, batch = "orig.ident",
                          do_normalization = NULL, do_HVF_finding = TRUE,
                          normalization_method, HVF_source = "separate", nHVF = 3000, hvf = NULL,
                          vars_to_regress = NULL,
                          exogenous_genes = NULL, ...) {
  cat(paste0("[", Sys.time(), "]", " Checking srtList... ...\n"))
  library(Seurat)
  library(sctransform)
  library(glmGamPoi)

  if (class(srtList) != "list" | any(sapply(srtList, class) != "Seurat")) {
    stop("'srtList' is not a list of Seurat object.")
  }
  if (length(batch) != 1) {
    stop("'batch' must be a vector to specify the batch column in Seurat object!")
  }
  if (!all(sapply(srtList, function(x) {
    "orig.ident" %in% colnames(x@meta.data)
  }))) {
    stop(paste0("batch column('", batch, "') was not found in one or more object of the srtList!"))
  }
  if (!normalization_method %in% c("logCPM", "SCT")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT'")
  }
  if (!HVF_source %in% c("global", "separate")) {
    stop("'HVF_source' must be one of: 'global','separate'")
  }
  if (any(sapply(srtList, ncol) < 10)) {
    stop(paste0("Seurat objects in srtList contain less than 10 cells. srtList index: ", which(sapply(srtList, ncol) < 10)))
  }

  genelist <- lapply(srtList, function(x) {
    sort(rownames(GetAssayData(x, slot = "counts", assay = "RNA")))
  })
  if (length(unique(genelist)) != 1) {
    warning("'srtList' have different feature names! Will subset the common features for downstream analysis!")
    cf <- lapply(srtList, rownames) %>% Reduce(intersect, .)
    for (i in 1:length(srtList)) {
      srtList[[i]] <- subset(srtList[[i]], features = cf)
    }
  }

  celllist <- unlist(lapply(srtList, colnames))
  if (length(celllist) != length(unique(celllist))) {
    stop("'srtList' have duplicated cell names!")
  }

  for (i in 1:length(srtList)) {
    if (!is.factor(srtList[[i]][[batch, drop = T]])) {
      srtList[[i]][[batch, drop = T]] <- factor(srtList[[i]][[batch, drop = T]], levels = unique(srtList[[i]][[batch, drop = T]]))
    }
    if (!"RNA" %in% Seurat::Assays(srtList[[i]])) {
      stop(paste("srtList", i, "does not contain 'RNA' assay."))
    }
    DefaultAssay(srtList[[i]]) <- "RNA"
    if (isTRUE(do_normalization) | (is.null(do_normalization) & identical(
      x = GetAssayData(srtList[[i]], slot = "counts"),
      y = GetAssayData(srtList[[i]], slot = "data")
    ))) {
      cat("Perform NormalizeData(logCPM) on the data ", i, "/", length(srtList), " of the srtList...\n", sep = "")
      srtList[[i]] <- suppressWarnings(NormalizeData(object = srtList[[i]], normalization.method = "LogNormalize", verbose = FALSE))
    }
    if (is.null(hvf)) {
      if (isTRUE(do_HVF_finding) | (is.null(do_HVF_finding) & !"vst.variance.standardized" %in% colnames(srtList[[i]]@assays$RNA@meta.features))) {
        cat("Perform FindVariableFeatures on the data ", i, "/", length(srtList), " of the srtList...\n", sep = "")
        srtList[[i]] <- FindVariableFeatures(srtList[[i]], verbose = FALSE)
      }
      m <- GetAssayData(srtList[[i]], slot = "counts")
      gene_keep <- rownames(m)[Matrix::rowSums(m >= 1) >= 5]
      VariableFeatures(srtList[[i]]) <- srtList[[i]]@assays$RNA@meta.features %>%
        filter(vst.variance.standardized > 1 &
          (!rownames(.) %in% exogenous_genes) &
          rownames(.) %in% gene_keep) %>%
        dplyr::arrange(desc(vst.variance.standardized)) %>%
        rownames(.) %>%
        head(n = nHVF)
    }

    if (normalization_method %in% c("SCT")) {
      if (isTRUE(do_normalization) | isTRUE(do_HVF_finding) | !"SCT" %in% Seurat::Assays(srtList[[i]])) {
        cat("Perform SCTransform on the data", i, "of the srtList...\n")
        srtList[[i]] <- SCTransform(
          object = srtList[[i]],
          method = "glmGamPoi",
          variable.features.n = nHVF,
          vars.to.regress = vars_to_regress,
          return.only.var.genes = FALSE,
          min_cells = 5,
          assay = "RNA",
          verbose = FALSE
        )
      } else {
        DefaultAssay(srtList[[i]]) <- "SCT"
      }
      if (!"residual_variance" %in% colnames(srtList[[i]]@assays$SCT@meta.features)) {
        if (length(srtList[[i]]@assays$SCT@SCTModel.list) > 1) {
          index <- which(sapply(srtList[[i]]@assays$SCT@SCTModel.list, function(x) nrow(x@cell.attributes) == ncol(srtList[[i]])))
        } else {
          index <- 1
        }
        model <- srtList[[i]]@assays$SCT@SCTModel.list[[index]]
        feature.attr <- SCTResults(object = model, slot = "feature.attributes")
        nfeatures <- min(nHVF, nrow(x = feature.attr))
        top.features <- rownames(x = feature.attr)[order(feature.attr$residual_variance,
          decreasing = TRUE
        )[1:nHVF]]
        VariableFeatures(object = srtList[[i]]) <- top.features
        srtList[[i]]@assays$SCT@meta.features <- feature.attr
      }
      m <- GetAssayData(srtList[[i]], slot = "counts")
      gene_keep <- rownames(m)[Matrix::rowSums(m >= 1) >= 5]
      VariableFeatures(srtList[[i]]) <- srtList[[i]]@assays$SCT@meta.features %>%
        filter((!rownames(.) %in% exogenous_genes) &
          rownames(.) %in% gene_keep) %>%
        dplyr::arrange(desc(residual_variance)) %>%
        rownames(.) %>%
        head(n = nHVF)
    }
  }

  if (is.null(hvf)) {
    if (HVF_source == "global") {
      cat("Perform global HVF calculation on the merged datasets from the srtList...\n")
      gene_common <- lapply(srtList, function(x) {
        m <- GetAssayData(x, slot = "counts", assay = "RNA")
        gene_keep <- rownames(m)[Matrix::rowSums(m >= 1) >= 5 * length(srtList)]
        return(gene_keep)
      }) %>% Reduce(intersect, .)
      srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
      srtMerge <- FindVariableFeatures(srtMerge, assay = "RNA", verbose = FALSE)
      hvf <- srtMerge@assays$RNA@meta.features %>%
        filter(vst.variance.standardized > 1 &
          (!rownames(.) %in% exogenous_genes) &
          rownames(.) %in% gene_common) %>%
        dplyr::arrange(desc(vst.variance.standardized)) %>%
        rownames(.) %>%
        head(n = nHVF)
    }
    if (HVF_source == "separate") {
      cat("Perform separate HVF calculation with SelectIntegrationFeatures from the existed HVF in srtList...\n")
      hvf <- SelectIntegrationFeatures(object.list = srtList, nfeatures = nHVF, verbose = FALSE)
    }
  } else {
    cf <- lapply(srtList, function(x) {
      rownames(GetAssayData(x, slot = "counts"))
    }) %>% Reduce(intersect, .)
    hvf <- hvf[hvf %in% cf]
  }

  if (normalization_method %in% c("SCT")) {
    srtList <- PrepSCTIntegration(object.list = srtList, anchor.features = hvf, assay = "SCT", verbose = FALSE)
  }
  cat(paste0("[", Sys.time(), "]", " Finished checking.\n"))

  return(list(
    srtList = srtList,
    hvf = hvf
  ))
}

Check_srtMerge <- function(srtMerge, batch = "orig.ident", do_normalization = NULL, do_HVF_finding = TRUE, do_scaling = TRUE,
                           normalization_method, vars_to_regress = NULL,
                           HVF_source = "separate", nHVF = 3000, hvf = NULL,
                           exogenous_genes = NULL, ...) {
  if (class(srtMerge) != "Seurat") {
    stop("'srtMerge' is not a Seurat object.")
  }
  if (length(batch) != 1) {
    stop("'batch' must be a vector to specify the batch column in srtMerge object!")
  }
  if (!batch %in% colnames(srtMerge@meta.data)) {
    stop(paste0("No batch column('", batch, "') found in the srtMerge object!"))
  }
  if (!is.factor(srtMerge[[batch, drop = T]])) {
    srtMerge[[batch, drop = T]] <- factor(srtMerge[[batch, drop = T]], levels = unique(srtMerge[[batch, drop = T]]))
  }
  if (isTRUE(do_scaling) | (is.null(do_scaling) & nrow(GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")) != nrow(GetAssayData(srtMerge, slot = "data", assay = "RNA")))) {
    cat("Perform ScaleData on the data...\n")
    srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
  }
  scale_data <- GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")

  cat("Spliting srtMerge into srtList... ...\n")
  srtList <- SplitObject(object = srtMerge, split.by = batch)

  checked <- Check_srtList(
    srtList = srtList, batch = batch,
    do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
    normalization_method = normalization_method, vars_to_regress = vars_to_regress,
    HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
    exogenous_genes = exogenous_genes
  )
  srtList <- checked[["srtList"]]
  hvf <- checked[["hvf"]]
  srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
  srtMerge <- SetAssayData(srtMerge, slot = "scale.data", assay = "RNA", new.data = scale_data)
  VariableFeatures(srtMerge) <- hvf

  return(list(
    srtMerge = srtMerge,
    srtList = srtList,
    hvf = hvf
  ))
}

Check_srtIntegrated <- function(srtIntegrated, do_normalization = NULL, hvf, batch, vars_to_regress, ...) {
  raw_DefaultAssay <- DefaultAssay(object = srtIntegrated)

  DefaultAssay(object = srtIntegrated) <- "RNA"
  if (isTRUE(do_normalization) | (is.null(do_normalization) & identical(
    x = GetAssayData(srtIntegrated, slot = "counts"),
    y = GetAssayData(srtIntegrated, slot = "data")
  ))) {
    srtIntegrated <- NormalizeData(object = srtIntegrated)
  }
  if (length(VariableFeatures(srtIntegrated)) == 0) {
    hvf <- hvf[hvf %in% rownames(GetAssayData(srtIntegrated, slot = "counts"))]
    VariableFeatures(srtIntegrated) <- hvf
  }
  if (nrow(GetAssayData(srtIntegrated, slot = "scale.data")) != nrow(GetAssayData(srtIntegrated, slot = "data"))) {
    cat("Perform ScaleData on the data(Check_srtIntegrated)...\n")
    srtIntegrated <- Seurat::ScaleData(object = srtIntegrated, features = rownames(srtIntegrated), vars.to.regress = vars_to_regress)
  }

  DefaultAssay(object = srtIntegrated) <- raw_DefaultAssay
  srtIntegrated@project.name <- paste0(unique(srtIntegrated[[batch, drop = TRUE]]), collapse = ",")
  srtIntegrated[[batch]] <- factor(srtIntegrated[[batch, drop = TRUE]],
    levels = unique(srtIntegrated[[batch, drop = TRUE]])
  )
  return(srtIntegrated)
}

RenameFeatures <- function(srt, newnames = NULL) {
  for (i in Seurat::Assays(srt)) {
    assay <- GetAssay(srt, i)
    for (d in c("counts", "data", "scale.data", "meta.features")) {
      if (dim(slot(assay, d))[1] == length(newnames)) {
        rownames(slot(assay, d)) <- newnames
      } else {
        if (identical(dim(slot(assay, d)), as.integer(c(0, 0)))) {
          next
        } else {
          message(paste0("Slot ", d, " have a different number of features."))
        }
      }
    }
    srt@assays[[i]] <- assay
  }
  return(srt)
}

SrtReorder <- function(srt, features = NULL, reorder_by = NULL, slot = "data", assay = NULL) {
  if (is.null(features)) {
    features <- VariableFeatures(srt)
  }
  features <- intersect(x = features, y = rownames(x = srt))
  if (is.null(reorder_by)) {
    srt$ident <- Idents(srt)
  } else {
    srt$ident <- srt[[reorder_by, drop = TRUE]]
  }
  if (is.null(assay)) {
    assay <- DefaultAssay(srt)
  }
  data.avg <- AverageExpression(object = srt, features = features, slot = slot, assays = assay, group.by = "ident")[[1]]
  data.dist <- dist(x = t(x = data.avg[features, ]))
  hc <- hclust(d = data.dist)
  ident_new <- plyr::mapvalues(x = srt$ident, from = hc$labels[hc$order], to = 1:length(hc$labels))
  ident_new <- factor(ident_new, levels = 1:length(hc$labels))
  Idents(srt) <- srt$ident <- ident_new
  return(srt)
}

SrtAppend <- function(srt_raw, srt_append, append_pattern, show_warning = TRUE) {
  if (class(srt_raw) != "Seurat" | class(srt_append) != "Seurat") {
    stop("'srt_raw' or 'srt_append' is not a Seurat object.")
  }
  for (slot_nm in slotNames(srt_append)) {
    if (slot_nm %in% c("assays", "meta.data", "graphs", "neighbors", "reductions", "images", "misc", "tools")) {
      for (info in names(slot(srt_append, name = slot_nm))) {
        if (!info %in% names(slot(srt_raw, name = slot_nm)) | grepl(pattern = append_pattern, x = info)) {
          slot(srt_raw, name = slot_nm)[[info]] <- slot(srt_append, name = slot_nm)[[info]]
        }
      }
    } else {
      if (isTRUE(show_warning)) {
        warning(paste0(slot_nm, "is not appended."))
      }
    }
  }
  return(srt_raw)
}

RunDimReduction <- function(srt, prefix = NULL, features = NULL,
                            liner_reduction = NULL, liner_reduction_dims = NULL, liner_reduction_distance = "cosine",
                            nonliner_reduction = NULL, reduction_use = NULL, dims_use = NULL, nonliner_reduction_dims = NULL, nonliner_reduction_distance = "euclidean", verbose = TRUE) {
  library(SeuratWrappers)
  library(glmpca)
  library(intrinsicDimension)

  if (is.null(liner_reduction) + is.null(liner_reduction_dims) == 0) {
    if (!liner_reduction %in% c("pca", "ica", "nmf", "mds", "glmpca")) {
      stop("'liner_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
    }
    key_use <- switch(liner_reduction,
      "pca" = "PC_",
      "ica" = "IC_",
      "nmf" = "BE_",
      "mds" = "MDS_",
      "glmpca" = "GLMPC_"
    )
    if (liner_reduction == "pca") {
      srt <- RunPCA(
        object = srt, features = features, npcs = liner_reduction_dims,
        reduction.name = paste0(prefix, liner_reduction),
        reduction.key = paste0(prefix, key_use),
        verbose = verbose
      )
    }
    if (liner_reduction == "glmpca") {
      srt <- RunGLMPCA(
        object = srt, features = features, L = liner_reduction_dims,
        reduction.name = paste0(prefix, liner_reduction),
        reduction.key = paste0(prefix, key_use),
        verbose = verbose
      )
    }
    if (liner_reduction == "ica") {
      srt <- RunICA(
        object = srt, features = features, nics = liner_reduction_dims,
        reduction.name = paste0(prefix, liner_reduction),
        reduction.key = paste0(prefix, key_use),
        verbose = verbose
      )
    }
    if (liner_reduction == "nmf") {
      srt <- RunNMF(
        object = srt, features = features, nbes = liner_reduction_dims,
        reduction.name = paste0(prefix, liner_reduction),
        reduction.key = paste0(prefix, key_use),
        verbose = verbose
      )
    }
    if (liner_reduction == "mds") {
      srt <- RunMDS(
        object = srt, features = features, nmds = liner_reduction_dims,
        dist.method = liner_reduction_distance,
        reduction.name = paste0(prefix, liner_reduction),
        reduction.key = paste0(prefix, key_use),
        verbose = verbose
      )
    }
    if (liner_reduction %in% c("glmpca", "nmf")) {
      dims_estimate <- 1:liner_reduction_dims
    }
    dim_est <- maxLikGlobalDimEst(data = Embeddings(srt, reduction = paste0(prefix, liner_reduction)), k = 20, iterations = 100)[["dim.est"]]
    if (!is.na(dim_est)) {
      if (ceiling(dim_est) < 10) {
        dims_estimate <- 1:10
      } else {
        dims_estimate <- 1:ceiling(dim_est)
      }
    } else {
      dims_estimate <- 1:30
    }
    srt@reductions[[paste0(prefix, liner_reduction)]]@misc[["dims_estimate"]] <- dims_estimate
  } else if (is.null(nonliner_reduction) + is.null(nonliner_reduction_dims) + is.null(reduction_use) == 0) {
    if (!nonliner_reduction %in% c("umap", "tsne", "dm", "fdg", "isomap", "dbmap", "phate")) {
      stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
    }
    if (is.null(dims_use)) {
      dims_use <- 1:30
    }
    if (is.null(prefix)) {
      prefix <- reduction_use
    }
    if (nonliner_reduction == "umap") {
      srt <- RunUMAP(
        object = srt, reduction = reduction_use, features = features,
        dims = dims_use, n.components = nonliner_reduction_dims, umap.method = "uwot",
        reduction.name = paste0(prefix, "UMAP", nonliner_reduction_dims, "d"),
        reduction.key = paste0(prefix, "UMAP", nonliner_reduction_dims, "d_"),
        return.model = TRUE, verbose = verbose
      )
    }
    if (nonliner_reduction == "tsne") {
      srt <- RunTSNE(
        object = srt, reduction = reduction_use, features = features,
        dims = dims_use, dim.embed = nonliner_reduction_dims, tsne.method = "Rtsne",
        reduction.name = paste0(prefix, "TSNE", nonliner_reduction_dims, "d"),
        reduction.key = paste0(prefix, "TSNE", nonliner_reduction_dims, "d_"),
        num_threads = 0, verbose = verbose
      )
    }
    if (nonliner_reduction == "dm") {
      srt <- RunDM(
        object = srt, reduction = reduction_use, features = features,
        dims = dims_use, ndcs = nonliner_reduction_dims, dist.method = nonliner_reduction_distance,
        reduction.name = paste0(prefix, "DM", nonliner_reduction_dims, "d"),
        reduction.key = paste0(prefix, "DM", nonliner_reduction_dims, "d_"),
        verbose = verbose
      )
    }
    srt@reductions[[paste0(prefix, toupper(nonliner_reduction), nonliner_reduction_dims, "d")]]@misc[["dims_use"]] <- dims_use
    srt@reductions[[paste0(prefix, toupper(nonliner_reduction), nonliner_reduction_dims, "d")]]@misc[["reduction_use"]] <- reduction_use
  } else {
    stop("No reduction method provided.")
  }
  return(srt)
}

Uncorrected_integrate <- function(srtList = NULL, srtMerge = NULL, batch = "orig.ident", append = TRUE,
                                  do_normalization = NULL, normalization_method = "logCPM",
                                  do_HVF_finding = TRUE, HVF_source = "separate", nHVF = 3000, hvf = NULL,
                                  do_scaling = TRUE, vars_to_regress = NULL,
                                  liner_reduction = "pca", liner_reduction_dims = 50, liner_reduction_distance = "cosine", liner_reduction_dims_use = NULL,
                                  nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "euclidean",
                                  cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                                  exogenous_genes = NULL, seed = 11, ...) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (!liner_reduction %in% c("pca", "ica", "nmf", "mds", "glmpca")) {
    stop("'liner_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (length(liner_reduction) > 1) {
    warning("Only the first method in the 'liner_reduction' will be used.")
    liner_reduction <- liner_reduction[1]
  }
  if (!nonliner_reduction %in% c("umap", "tsne", "dm", "fdg", "isomap", "dbmap", "phate")) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm_index' must be one of 'louvain', 'slm', 'leiden'.")
  }

  set.seed(seed)
  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) & !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
    if (isTRUE(do_scaling) | (is.null(do_scaling) & nrow(GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")) != nrow(GetAssayData(srtMerge, slot = "data", assay = "RNA")))) {
      cat("Perform ScaleData on the srtMerge...\n")
      srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
    }
    scale_data <- GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
  }
  if (exists("scale_data")) {
    srtMerge <- SetAssayData(srtMerge, slot = "scale.data", assay = "RNA", new.data = scale_data)
  } else {
    cat("Perform ScaleData on the data...\n")
    srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
  }
  do_scaling <- NULL
  if (is.null(srtList) & !is.null(srtMerge)) {
    checked <- Check_srtMerge(
      srtMerge = srtMerge, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding, do_scaling = do_scaling,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
    do_scaling <- NULL
  }

  srtIntegrated <- Standard_SCP(
    srt = srtMerge, prefix = "Uncorrected",
    do_normalization = do_normalization, normalization_method = normalization_method,
    do_HVF_finding = do_HVF_finding, nHVF = nHVF, hvf = hvf,
    do_scaling = do_scaling, vars_to_regress = vars_to_regress,
    liner_reduction = liner_reduction, liner_reduction_dims = liner_reduction_dims, liner_reduction_distance = liner_reduction_distance, liner_reduction_dims_use = liner_reduction_dims_use,
    nonliner_reduction = nonliner_reduction, nonliner_reduction_dims = nonliner_reduction_dims, nonliner_reduction_distance = nonliner_reduction_distance,
    cluster_algorithm = cluster_algorithm, cluster_resolution = cluster_resolution, cluster_reorder = cluster_reorder,
    exogenous_genes = exogenous_genes, seed = seed
  )
  srtIntegrated[[paste0("Uncorrectedclusters")]] <- srtIntegrated[[paste0("Uncorrected", liner_reduction, "clusters")]]
  srtIntegrated[[paste0("Uncorrected", liner_reduction, "clusters")]] <- NULL
  for (nr in nonliner_reduction) {
    for (n in nonliner_reduction_dims) {
      reduc <- srtIntegrated@reductions[[paste0("Uncorrected", liner_reduction, toupper(nr), n, "d")]]
      srtIntegrated@reductions[[paste0("Uncorrected", liner_reduction, toupper(nr), n, "d")]] <- NULL
      srtIntegrated@reductions[[paste0("Uncorrected", toupper(nr), n, "d")]] <- reduc
    }
  }
  srtIntegrated <- Check_srtIntegrated(srtIntegrated, do_normalization = do_normalization, hvf = hvf, batch = batch, vars_to_regress = vars_to_regress)

  if (isTRUE(append) & !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, append_pattern = "Uncorrected", show_warning = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

Seurat_integrate <- function(srtList = NULL, srtMerge = NULL, batch = "orig.ident", append = TRUE,
                             do_normalization = NULL, normalization_method = "logCPM",
                             do_HVF_finding = TRUE, HVF_source = "separate", nHVF = 3000, hvf = NULL,
                             do_scaling = TRUE, vars_to_regress = NULL,
                             liner_reduction = "pca", liner_reduction_dims = 50, liner_reduction_dims_use = NULL, liner_reduction_distance = "cosine",
                             nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "euclidean",
                             cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                             exogenous_genes = NULL, seed = 11, ...) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (!liner_reduction %in% c("pca", "ica", "nmf", "mds", "glmpca")) {
    stop("'liner_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (length(liner_reduction) > 1) {
    warning("Only the first method in the 'liner_reduction' will be used.")
    liner_reduction <- liner_reduction[1]
  }
  if (!nonliner_reduction %in% c("umap", "tsne", "dm", "fdg", "isomap", "dbmap", "phate")) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm_index' must be one of 'louvain', 'slm', 'leiden'.")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) & !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
    if (isTRUE(do_scaling) | (is.null(do_scaling) & nrow(GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")) != nrow(GetAssayData(srtMerge, slot = "data", assay = "RNA")))) {
      cat("Perform ScaleData on the srtMerge...\n")
      srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
    }
    scale_data <- GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
  }
  if (is.null(srtList) & !is.null(srtMerge)) {
    srtList <- SplitObject(object = srtMerge, split.by = batch)
    checked <- Check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
  }

  srt_anchors <- FindIntegrationAnchors(
    object.list = srtList,
    normalization.method = switch(normalization_method,
      "logCPM" = "LogNormalize",
      "SCT" = "SCT"
    ),
    anchor.features = hvf,
    verbose = FALSE
  )
  srtIntegrated <- IntegrateData(
    anchorset = srt_anchors,
    new.assay.name = "Seurat",
    normalization.method = switch(normalization_method,
      "logCPM" = "LogNormalize",
      "SCT" = "SCT"
    ),
    features.to.integrate = c(
      hvf,
      Reduce(union, lapply(srtList, VariableFeatures)),
      Reduce(intersect, lapply(srtList, rownames))
    ), verbose = FALSE
  )
  DefaultAssay(srtIntegrated) <- "Seurat"

  if (exists("scale_data")) {
    srtIntegrated <- SetAssayData(srtIntegrated, slot = "scale.data", assay = "RNA", new.data = scale_data)
  }
  srtIntegrated <- Check_srtIntegrated(srtIntegrated, do_normalization = do_normalization, hvf = hvf, batch = batch, vars_to_regress = vars_to_regress)
  if (isTRUE(do_scaling) | (is.null(do_scaling) & nrow(GetAssayData(srtIntegrated, slot = "scale.data")) != nrow(GetAssayData(srtIntegrated, slot = "data")))) {
    cat("Perform ScaleData on the data(after integrated)...\n")
    srtIntegrated <- Seurat::ScaleData(object = srtIntegrated, features = hvf, vars.to.regress = vars_to_regress)
  }

  cat("Perform linear dimension reduction (", liner_reduction, ") on the data...\n", sep = "")
  srtIntegrated <- RunDimReduction(
    srt = srtIntegrated, prefix = "Seurat", features = hvf, liner_reduction_distance = liner_reduction_distance,
    liner_reduction = liner_reduction, liner_reduction_dims = liner_reduction_dims,
    verbose = FALSE
  )
  if (is.null(liner_reduction_dims_use)) {
    liner_reduction_dims_use <- srtIntegrated@reductions[[paste0("Seurat", liner_reduction)]]@misc[["dims_estimate"]]
  }

  cat("Perform FindClusters (", cluster_algorithm, ") on the data...\n", sep = "")
  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = paste0("Seurat", liner_reduction), dims = liner_reduction_dims_use, force.recalc = T, graph.name = paste0("Seurat", liner_reduction, "_", c("KNN", "SNN")), verbose = FALSE)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, graph.name = paste0("Seurat", liner_reduction, "_SNN"), verbose = FALSE)
  if (isTRUE(cluster_reorder)) {
    srtIntegrated <- SrtReorder(srtIntegrated, features = hvf, reorder_by = "seurat_clusters", slot = "data")
  }
  srtIntegrated[["seurat_clusters"]] <- NULL
  srtIntegrated[["Seuratclusters"]] <- Idents(srtIntegrated)

  for (nr in nonliner_reduction) {
    cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
    for (n in nonliner_reduction_dims) {
      srtIntegrated <- RunDimReduction(
        srt = srtIntegrated, prefix = "Seurat",
        reduction_use = paste0("Seurat", liner_reduction), dims_use = liner_reduction_dims_use,
        nonliner_reduction = nr, nonliner_reduction_dims = n,
        nonliner_reduction_distance = nonliner_reduction_distance,
        verbose = FALSE
      )
    }
  }

  DefaultAssay(srtIntegrated) <- "RNA"
  if (isTRUE(append) & !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, append_pattern = "Seurat", show_warning = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

fastMNN_integrate <- function(srtList = NULL, srtMerge = NULL, batch = "orig.ident", append = TRUE,
                              do_normalization = NULL, normalization_method = "logCPM",
                              do_HVF_finding = TRUE, HVF_source = "separate", nHVF = 3000, hvf = NULL,
                              do_scaling = TRUE, vars_to_regress = NULL,
                              fastMNN_dims_use = NULL,
                              nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "euclidean",
                              cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                              exogenous_genes = NULL, seed = 11, ...) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (!nonliner_reduction %in% c("umap", "tsne", "dm", "fdg", "isomap", "dbmap", "phate")) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm_index' must be one of 'louvain', 'slm', 'leiden'.")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  library(SeuratWrappers)
  library(BiocParallel)

  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) & !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
    if (isTRUE(do_scaling) | (is.null(do_scaling) & nrow(GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")) != nrow(GetAssayData(srtMerge, slot = "data", assay = "RNA")))) {
      cat("Perform ScaleData on the srtMerge...\n")
      srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
    }
    scale_data <- GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
  }
  if (is.null(srtList) & !is.null(srtMerge)) {
    srtList <- SplitObject(object = srtMerge, split.by = batch)
    checked <- Check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
  }

  srtIntegrated <- RunFastMNN(
    object.list = srtList,
    features = hvf,
    reduction.name = "fastMNN",
    reduction.key = "fastMNN_",
    assay = DefaultAssay(srtList[[1]]),
    assay.type = "logcounts",
    BPPARAM = MulticoreParam()
  )
  srtIntegrated@assays$mnn.reconstructed <- NULL

  if (exists("scale_data")) {
    srtIntegrated <- SetAssayData(srtIntegrated, slot = "scale.data", assay = "RNA", new.data = scale_data)
  }
  srtIntegrated <- Check_srtIntegrated(srtIntegrated, do_normalization = do_normalization, hvf = hvf, batch = batch, vars_to_regress = vars_to_regress)

  if (is.null(fastMNN_dims_use)) {
    dim_est <- maxLikGlobalDimEst(data = Embeddings(srtIntegrated, reduction = "fastMNN"), k = 20, iterations = 100)[["dim.est"]]
    if (!is.na(dim_est)) {
      if (ceiling(dim_est) < 10) {
        fastMNN_dims_use <- 1:10
      } else {
        fastMNN_dims_use <- 1:ceiling(dim_est)
      }
    } else {
      fastMNN_dims_use <- 1:30
    }
  }
  srtIntegrated@reductions[["fastMNN"]]@misc[["dims_estimate"]] <- fastMNN_dims_use

  cat("Perform FindClusters (", cluster_algorithm, ") on the data...\n", sep = "")
  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "fastMNN", dims = fastMNN_dims_use, force.recalc = T, graph.name = paste0("fastMNN", "_", c("KNN", "SNN")), verbose = FALSE)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, graph.name = paste0("fastMNN", "_SNN"), verbose = FALSE)
  if (isTRUE(cluster_reorder)) {
    srtIntegrated <- SrtReorder(srtIntegrated, features = hvf, reorder_by = "seurat_clusters", slot = "data")
  }
  srtIntegrated[["seurat_clusters"]] <- NULL
  srtIntegrated[["fastMNNclusters"]] <- Idents(srtIntegrated)

  for (nr in nonliner_reduction) {
    cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
    for (n in nonliner_reduction_dims) {
      srtIntegrated <- RunDimReduction(
        srt = srtIntegrated, prefix = "fastMNN",
        reduction_use = "fastMNN", dims_use = fastMNN_dims_use,
        nonliner_reduction = nr, nonliner_reduction_dims = n,
        nonliner_reduction_distance = nonliner_reduction_distance,
        verbose = FALSE
      )
    }
  }

  DefaultAssay(srtIntegrated) <- "RNA"
  if (isTRUE(append) & !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, append_pattern = "fastMNN", show_warning = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

Harmony_integrate <- function(srtList = NULL, srtMerge = NULL, batch = "orig.ident", append = TRUE,
                              do_normalization = NULL, normalization_method = "logCPM",
                              do_HVF_finding = TRUE, HVF_source = "separate", nHVF = 3000, hvf = NULL,
                              do_scaling = TRUE, vars_to_regress = NULL,
                              liner_reduction = "pca", liner_reduction_dims = 50, liner_reduction_dims_use = NULL, liner_reduction_distance = "cosine",
                              Harmony_dims_use = NULL,
                              nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "euclidean",
                              cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                              exogenous_genes = NULL, seed = 11, ...) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (!liner_reduction %in% c("pca", "ica", "nmf", "mds", "glmpca")) {
    stop("'liner_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (length(liner_reduction) > 1) {
    warning("Only the first method in the 'liner_reduction' will be used.")
    liner_reduction <- liner_reduction[1]
  }
  if (!nonliner_reduction %in% c("umap", "tsne", "dm", "fdg", "isomap", "dbmap", "phate")) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm_index' must be one of 'louvain', 'slm', 'leiden'.")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  library(SeuratWrappers)
  library(harmony)

  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) & !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
    if (isTRUE(do_scaling) | (is.null(do_scaling) & nrow(GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")) != nrow(GetAssayData(srtMerge, slot = "data", assay = "RNA")))) {
      cat("Perform ScaleData on the srtMerge...\n")
      srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
    }
    scale_data <- GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
  }
  if (exists("scale_data")) {
    srtMerge <- SetAssayData(srtMerge, slot = "scale.data", assay = "RNA", new.data = scale_data)
  } else {
    cat("Perform ScaleData on the data...\n")
    srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
  }
  do_scaling <- NULL
  if (is.null(srtList) & !is.null(srtMerge)) {
    checked <- Check_srtMerge(
      srtMerge = srtMerge, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding, do_scaling = do_scaling,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
    do_scaling <- NULL
  }

  cat("Perform linear dimension reduction (", liner_reduction, ") on the data...\n", sep = "")
  srtMerge <- RunDimReduction(
    srt = srtMerge, prefix = "Harmony", features = hvf, liner_reduction_distance = liner_reduction_distance,
    liner_reduction = liner_reduction, liner_reduction_dims = liner_reduction_dims,
    verbose = FALSE
  )
  if (is.null(liner_reduction_dims_use)) {
    liner_reduction_dims_use <- srtMerge@reductions[[paste0("Harmony", liner_reduction)]]@misc[["dims_estimate"]]
  }

  srtIntegrated <- RunHarmony(
    object = srtMerge,
    group.by.vars = batch,
    assay.use = DefaultAssay(srtMerge),
    reduction = paste0("Harmony", liner_reduction),
    dims.use = liner_reduction_dims_use,
    reduction.save = "Harmony",
    verbose = FALSE
  )
  srtMerge <- NULL

  srtIntegrated <- Check_srtIntegrated(srtIntegrated, do_normalization = do_normalization, hvf = hvf, batch = batch, vars_to_regress = vars_to_regress)

  if (is.null(Harmony_dims_use)) {
    dim_est <- maxLikGlobalDimEst(data = Embeddings(srtIntegrated, reduction = "Harmony"), k = 20, iterations = 100)[["dim.est"]]
    if (!is.na(dim_est)) {
      if (ceiling(dim_est) < 10) {
        Harmony_dims_use <- 1:10
      } else {
        Harmony_dims_use <- 1:ceiling(dim_est)
      }
    } else {
      Harmony_dims_use <- 1:30
    }
  }
  srtIntegrated@reductions[["Harmony"]]@misc[["dims_estimate"]] <- Harmony_dims_use

  cat("Perform FindClusters (", cluster_algorithm, ") on the data...\n", sep = "")
  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "Harmony", dims = Harmony_dims_use, force.recalc = T, graph.name = paste0("Harmony", "_", c("KNN", "SNN")), verbose = FALSE)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, graph.name = paste0("Harmony", "_SNN"), verbose = FALSE)
  if (isTRUE(cluster_reorder)) {
    srtIntegrated <- SrtReorder(srtIntegrated, features = hvf, reorder_by = "seurat_clusters", slot = "data")
  }
  srtIntegrated[["seurat_clusters"]] <- NULL
  srtIntegrated[["Harmonyclusters"]] <- Idents(srtIntegrated)

  for (nr in nonliner_reduction) {
    cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
    for (n in nonliner_reduction_dims) {
      srtIntegrated <- RunDimReduction(
        srt = srtIntegrated, prefix = "Harmony",
        reduction_use = "Harmony", dims_use = Harmony_dims_use,
        nonliner_reduction = nr, nonliner_reduction_dims = n,
        nonliner_reduction_distance = nonliner_reduction_distance,
        verbose = FALSE
      )
    }
  }

  DefaultAssay(srtIntegrated) <- "RNA"
  if (isTRUE(append) & !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, append_pattern = "Harmony", show_warning = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

Scanorama_integrate <- function(srtList = NULL, srtMerge = NULL, batch = "orig.ident", append = TRUE,
                                do_normalization = NULL, normalization_method = "logCPM",
                                do_HVF_finding = TRUE, HVF_source = "separate", nHVF = 3000, hvf = NULL,
                                do_scaling = TRUE, vars_to_regress = NULL,
                                liner_reduction = "pca", liner_reduction_dims = 50, liner_reduction_dims_use = NULL, liner_reduction_distance = "cosine",
                                nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "euclidean",
                                cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                                exogenous_genes = NULL, seed = 11, ...) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (!liner_reduction %in% c("pca", "ica", "nmf", "mds", "glmpca")) {
    stop("'liner_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (length(liner_reduction) > 1) {
    warning("Only the first method in the 'liner_reduction' will be used.")
    liner_reduction <- liner_reduction[1]
  }
  if (!nonliner_reduction %in% c("umap", "tsne", "dm", "fdg", "isomap", "dbmap", "phate")) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm_index' must be one of 'louvain', 'slm', 'leiden'.")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  library(reticulate)
  tryCatch(expr = {
    scanorama <- reticulate::import("scanorama")
  }, error = function(e) {
    stop("No python module named 'scanorama'")
  })

  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) & !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
    if (isTRUE(do_scaling) | (is.null(do_scaling) & nrow(GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")) != nrow(GetAssayData(srtMerge, slot = "data", assay = "RNA")))) {
      cat("Perform ScaleData on the srtMerge...\n")
      srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
    }
    scale_data <- GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
  }
  if (is.null(srtList) & !is.null(srtMerge)) {
    srtList <- SplitObject(object = srtMerge, split.by = batch)
    checked <- Check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
  }

  cat("Perform scanorama$correct... ...\n")
  assaylist <- list()
  genelist <- list()
  for (i in 1:length(srtList)) {
    assaylist[[i]] <- t(as.matrix(GetAssayData(object = srtList[[i]], slot = "data", assay = DefaultAssay(srtList[[1]]))))
    genelist[[i]] <- rownames(srtList[[i]])
  }
  integrated.corrected.data <- scanorama$correct(
    datasets_full = assaylist,
    genes_list = genelist,
    return_dimred = TRUE,
    return_dense = TRUE,
    union = FALSE
  )
  cor_value <- integrated.corrected.data[[2]] %>%
    plyr::rbind.fill.matrix() %>%
    t()
  rownames(cor_value) <- integrated.corrected.data[[3]]
  colnames(cor_value) <- unlist(sapply(assaylist, rownames))
  dim_reduction <- integrated.corrected.data[[1]] %>% plyr::rbind.fill.matrix()
  rownames(dim_reduction) <- unlist(sapply(assaylist, rownames))
  colnames(dim_reduction) <- paste0("Scanorama_", 1:100)
  stdevs <- apply(dim_reduction, MARGIN = 2, FUN = sd)

  srtIntegrated <- Reduce(function(x, y) merge(x, y), srtList)
  srtIntegrated[["Scanorama"]] <- CreateAssayObject(data = cor_value)
  srtIntegrated[["Scanorama_reduction"]] <- CreateDimReducObject(embeddings = dim_reduction, assay = "Scanorama", stdev = stdevs, key = "Scanorama_")
  DefaultAssay(srtIntegrated) <- "Scanorama"

  if (exists("scale_data")) {
    srtIntegrated <- SetAssayData(srtIntegrated, slot = "scale.data", assay = "RNA", new.data = scale_data)
  }
  srtIntegrated <- Check_srtIntegrated(srtIntegrated, do_normalization = do_normalization, hvf = hvf, batch = batch, vars_to_regress = vars_to_regress)

  cat("Perform ScaleData on the data(after integrated)...\n")
  srtIntegrated <- Seurat::ScaleData(srtIntegrated, features = hvf, vars.to.regress = vars_to_regress)

  cat("Perform linear dimension reduction (", liner_reduction, ") on the data...\n", sep = "")
  srtIntegrated <- RunDimReduction(
    srt = srtIntegrated, prefix = "Scanorama", features = hvf, liner_reduction_distance = liner_reduction_distance,
    liner_reduction = liner_reduction, liner_reduction_dims = liner_reduction_dims,
    verbose = FALSE
  )
  if (is.null(liner_reduction_dims_use)) {
    liner_reduction_dims_use <- srtIntegrated@reductions[[paste0("Scanorama", liner_reduction)]]@misc[["dims_estimate"]]
  }

  cat("Perform FindClusters (", cluster_algorithm, ") on the data...\n", sep = "")
  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = paste0("Scanorama", liner_reduction), dims = liner_reduction_dims_use, force.recalc = T, graph.name = paste0("Scanorama", liner_reduction, "_", c("KNN", "SNN")), verbose = FALSE)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, graph.name = paste0("Scanorama", liner_reduction, "_SNN"), verbose = FALSE)
  if (isTRUE(cluster_reorder)) {
    srtIntegrated <- SrtReorder(srtIntegrated, features = hvf, reorder_by = "seurat_clusters", slot = "data")
  }
  srtIntegrated[["seurat_clusters"]] <- NULL
  srtIntegrated[["Scanoramaclusters"]] <- Idents(srtIntegrated)

  for (nr in nonliner_reduction) {
    cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
    for (n in nonliner_reduction_dims) {
      srtIntegrated <- RunDimReduction(
        srt = srtIntegrated, prefix = "Scanorama",
        reduction_use = paste0("Scanorama", liner_reduction), dims_use = liner_reduction_dims_use,
        nonliner_reduction = nr, nonliner_reduction_dims = n,
        nonliner_reduction_distance = nonliner_reduction_distance,
        verbose = FALSE
      )
    }
  }

  DefaultAssay(srtIntegrated) <- "RNA"
  if (isTRUE(append) & !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, append_pattern = "Scanorama", show_warning = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

BBKNN_integrate <- function(srtList = NULL, srtMerge = NULL, batch = "orig.ident", append = TRUE,
                            do_normalization = NULL, normalization_method = "logCPM",
                            do_HVF_finding = TRUE, HVF_source = "separate", nHVF = 3000, hvf = NULL,
                            do_scaling = TRUE, vars_to_regress = NULL,
                            liner_reduction = "pca", liner_reduction_dims = 50, liner_reduction_dims_use = NULL, liner_reduction_distance = "cosine",
                            nonliner_reduction_dims = c(2, 3),
                            cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                            exogenous_genes = NULL, seed = 11, ...) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (!liner_reduction %in% c("pca", "ica", "nmf", "mds", "glmpca")) {
    stop("'liner_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (length(liner_reduction) > 1) {
    warning("Only the first method in the 'liner_reduction' will be used.")
    liner_reduction <- liner_reduction[1]
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm_index' must be one of 'louvain', 'slm', 'leiden'.")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  library(reticulate)
  bbknn <- reticulate::import("bbknn", convert = FALSE)

  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) & !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
    if (isTRUE(do_scaling) | (is.null(do_scaling) & nrow(GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")) != nrow(GetAssayData(srtMerge, slot = "data", assay = "RNA")))) {
      cat("Perform ScaleData on the srtMerge...\n")
      srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
    }
    scale_data <- GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
  }
  if (exists("scale_data")) {
    srtMerge <- SetAssayData(srtMerge, slot = "scale.data", assay = "RNA", new.data = scale_data)
  } else {
    cat("Perform ScaleData on the data...\n")
    srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
  }
  do_scaling <- NULL
  if (is.null(srtList) & !is.null(srtMerge)) {
    checked <- Check_srtMerge(
      srtMerge = srtMerge, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding, do_scaling = do_scaling,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
    do_scaling <- NULL
  }

  cat("Perform linear dimension reduction (", liner_reduction, ") on the data...\n", sep = "")
  srtMerge <- RunDimReduction(
    srt = srtMerge, prefix = "BBKNN", features = hvf, liner_reduction_distance = liner_reduction_distance,
    liner_reduction = liner_reduction, liner_reduction_dims = liner_reduction_dims,
    verbose = FALSE
  )
  if (is.null(liner_reduction_dims_use)) {
    liner_reduction_dims_use <- srtMerge@reductions[[paste0("BBKNN", liner_reduction)]]@misc[["dims_estimate"]]
  }

  pca <- reticulate::r_to_py(Embeddings(srtMerge, reduction = paste0("BBKNN", liner_reduction))[, liner_reduction_dims_use])
  bem <- bbknn$bbknn_matrix(pca, batch_list = srtMerge[[batch, drop = TRUE]])
  bem <- reticulate::py_to_r(bem)
  bbknn_graph <- as(as(bem[[2]], "CsparseMatrix"), "dgCMatrix")
  rownames(bbknn_graph) <- colnames(bbknn_graph) <- colnames(srtMerge)
  srtMerge@graphs$BBKNN <- as.Graph(bbknn_graph)
  srtIntegrated <- srtMerge
  srtMerge <- NULL
  srtIntegrated <- Check_srtIntegrated(srtIntegrated, do_normalization = do_normalization, hvf = hvf, batch = batch, vars_to_regress = vars_to_regress)

  cat("Perform FindClusters (", cluster_algorithm, ") on the data...\n", sep = "")
  srtIntegrated <- FindClusters(object = srtIntegrated, graph.name = "BBKNN", resolution = cluster_resolution, algorithm = cluster_algorithm_index, verbose = FALSE)
  if (isTRUE(cluster_reorder)) {
    srtIntegrated <- SrtReorder(srtIntegrated, features = hvf, reorder_by = "seurat_clusters", slot = "data")
  }
  srtIntegrated[["seurat_clusters"]] <- NULL
  srtIntegrated[["BBKNNclusters"]] <- Idents(srtIntegrated)

  cat("Perform nonlinear dimension reduction (umap) on the data...\n")
  for (n in nonliner_reduction_dims) {
    srtIntegrated <- RunUMAP(
      object = srtIntegrated, graph = "BBKNN",
      umap.method = "umap-learn", metric = "correlation",
      reduction.name = paste0("BBKNN", "UMAP", n, "d"),
      reduction.key = paste0("BBKNN", "UMAP", n, "d_"),
      n.components = as.integer(n), verbose = FALSE
    )
  }

  DefaultAssay(srtIntegrated) <- "RNA"
  if (isTRUE(append) & !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, append_pattern = "BBKNN", show_warning = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

CSS_integrate <- function(srtList = NULL, srtMerge = NULL, batch = "orig.ident", append = TRUE,
                          do_normalization = NULL, normalization_method = "logCPM",
                          do_HVF_finding = TRUE, HVF_source = "separate", nHVF = 3000, hvf = NULL,
                          do_scaling = TRUE, vars_to_regress = NULL,
                          liner_reduction = "pca", liner_reduction_dims = 50, liner_reduction_dims_use = NULL, liner_reduction_distance = "cosine",
                          nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "euclidean",
                          cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                          exogenous_genes = NULL, seed = 11, ...) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (!liner_reduction %in% c("pca", "ica", "nmf", "mds", "glmpca")) {
    stop("'liner_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (length(liner_reduction) > 1) {
    warning("Only the first method in the 'liner_reduction' will be used.")
    liner_reduction <- liner_reduction[1]
  }
  if (!nonliner_reduction %in% c("umap", "tsne", "dm", "fdg", "isomap", "dbmap", "phate")) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm_index' must be one of 'louvain', 'slm', 'leiden'.")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  if (!require(simspec)) {
    devtools::install_github("zh542370159/simspec@v1.0.3")
  }

  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) & !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
    if (isTRUE(do_scaling) | (is.null(do_scaling) & nrow(GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")) != nrow(GetAssayData(srtMerge, slot = "data", assay = "RNA")))) {
      cat("Perform ScaleData on the srtMerge...\n")
      srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
    }
    scale_data <- GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
  }
  if (exists("scale_data")) {
    srtMerge <- SetAssayData(srtMerge, slot = "scale.data", assay = "RNA", new.data = scale_data)
  } else {
    cat("Perform ScaleData on the data...\n")
    srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
  }
  do_scaling <- NULL
  if (is.null(srtList) & !is.null(srtMerge)) {
    checked <- Check_srtMerge(
      srtMerge = srtMerge, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding, do_scaling = do_scaling,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
    do_scaling <- NULL
  }

  cat("Perform linear dimension reduction (", liner_reduction, ") on the data...\n", sep = "")
  srtMerge <- RunDimReduction(
    srt = srtMerge, prefix = "CSS", features = hvf, liner_reduction_distance = liner_reduction_distance,
    liner_reduction = liner_reduction, liner_reduction_dims = liner_reduction_dims,
    verbose = FALSE
  )
  if (is.null(liner_reduction_dims_use)) {
    liner_reduction_dims_use <- srtMerge@reductions[[paste0("CSS", liner_reduction)]]@misc[["dims_estimate"]]
  }

  srtIntegrated <- cluster_sim_spectrum(
    object = srtMerge,
    use_dr = paste0("CSS", liner_reduction),
    dims_use = liner_reduction_dims_use,
    var_genes = hvf,
    label_tag = batch,
    reduction.name = "CSS",
    reduction.key = "CSS_",
    verbose = FALSE
  )
  srtMerge <- NULL

  srtIntegrated <- Check_srtIntegrated(srtIntegrated, do_normalization = do_normalization, hvf = hvf, batch = batch, vars_to_regress = vars_to_regress)
  CSS_dims_use <- 1:ncol(Embeddings(srtIntegrated, reduction = "CSS"))
  srtIntegrated@reductions[["CSS"]]@misc[["dims_estimate"]] <- CSS_dims_use

  cat("Perform FindClusters (", cluster_algorithm, ") on the data...\n", sep = "")
  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "CSS", dims = CSS_dims_use, force.recalc = T, graph.name = paste0("CSS", "_", c("KNN", "SNN")), verbose = FALSE)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, graph.name = paste0("CSS", "_SNN"), verbose = FALSE)
  if (isTRUE(cluster_reorder)) {
    srtIntegrated <- SrtReorder(srtIntegrated, features = hvf, reorder_by = "seurat_clusters", slot = "data")
  }
  srtIntegrated[["seurat_clusters"]] <- NULL
  srtIntegrated[["CSSclusters"]] <- Idents(srtIntegrated)

  for (nr in nonliner_reduction) {
    cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
    for (n in nonliner_reduction_dims) {
      srtIntegrated <- RunDimReduction(
        srt = srtIntegrated, prefix = "CSS",
        reduction_use = "CSS", dims_use = CSS_dims_use,
        nonliner_reduction = nr, nonliner_reduction_dims = n,
        nonliner_reduction_distance = nonliner_reduction_distance,
        verbose = FALSE
      )
    }
  }

  DefaultAssay(srtIntegrated) <- "RNA"
  if (isTRUE(append) & !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, append_pattern = "CSS", show_warning = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

LIGER_integrate <- function(srtList = NULL, srtMerge = NULL, batch = "orig.ident", append = TRUE,
                            do_normalization = NULL, normalization_method = "logCPM",
                            do_HVF_finding = TRUE, HVF_source = "separate", nHVF = 3000, hvf = NULL,
                            do_scaling = TRUE, vars_to_regress = NULL,
                            nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "euclidean",
                            cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                            exogenous_genes = NULL, seed = 11, ...) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (!nonliner_reduction %in% c("umap", "tsne", "dm", "fdg", "isomap", "dbmap", "phate")) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm_index' must be one of 'louvain', 'slm', 'leiden'.")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  library(SeuratWrappers)
  library(rliger)

  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) & !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
    if (isTRUE(do_scaling) | (is.null(do_scaling) & nrow(GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")) != nrow(GetAssayData(srtMerge, slot = "data", assay = "RNA")))) {
      cat("Perform ScaleData on the srtMerge...\n")
      srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
    }
    scale_data <- GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
  }
  if (exists("scale_data")) {
    srtMerge <- SetAssayData(srtMerge, slot = "scale.data", assay = "RNA", new.data = scale_data)
  } else {
    cat("Perform ScaleData on the data...\n")
    srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
  }
  do_scaling <- NULL
  if (is.null(srtList) & !is.null(srtMerge)) {
    checked <- Check_srtMerge(
      srtMerge = srtMerge, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding, do_scaling = do_scaling,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
    do_scaling <- NULL
  }

  cat("Perform ScaleData on the data(no center for LIGER)...\n")
  srtMerge_prep <- Seurat::ScaleData(object = srtMerge, features = hvf, split.by = batch, do.center = FALSE, vars.to.regress = vars_to_regress)
  srtMerge_prep <- RunOptimizeALS(srtMerge_prep,
    k = 20,
    lambda = 5,
    split.by = batch,
    verbose = FALSE
  )
  srtIntegrated <- RunQuantileNorm(srtMerge_prep,
    reduction.name = "LIGER",
    reduction.key = "LIGER_",
    split.by = batch,
    verbose = FALSE
  )
  srtMerge_prep <- NULL

  if (exists("scale_data")) {
    srtIntegrated <- SetAssayData(srtIntegrated, slot = "scale.data", assay = "RNA", new.data = scale_data)
  }
  srtIntegrated <- Check_srtIntegrated(srtIntegrated, do_normalization = do_normalization, hvf = hvf, batch = batch, vars_to_regress = vars_to_regress)
  LIGER_dims_use <- 1:ncol(Embeddings(srtIntegrated, reduction = "LIGER"))
  srtIntegrated@reductions[["LIGER"]]@misc[["dims_estimate"]] <- LIGER_dims_use

  cat("Perform FindClusters (", cluster_algorithm, ") on the data...\n", sep = "")
  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "LIGER", dims = LIGER_dims_use, force.recalc = T, graph.name = paste0("LIGER", "_", c("KNN", "SNN")), verbose = FALSE)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, graph.name = paste0("LIGER", "_SNN"), verbose = FALSE)
  if (isTRUE(cluster_reorder)) {
    srtIntegrated <- SrtReorder(srtIntegrated, features = hvf, reorder_by = "seurat_clusters", slot = "data")
  }
  srtIntegrated[["seurat_clusters"]] <- NULL
  srtIntegrated[["LIGERclusters"]] <- Idents(srtIntegrated)

  for (nr in nonliner_reduction) {
    cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
    for (n in nonliner_reduction_dims) {
      srtIntegrated <- RunDimReduction(
        srt = srtIntegrated, prefix = "LIGER",
        reduction_use = "LIGER", dims_use = LIGER_dims_use,
        nonliner_reduction = nr, nonliner_reduction_dims = n,
        nonliner_reduction_distance = nonliner_reduction_distance,
        verbose = FALSE
      )
    }
  }

  DefaultAssay(srtIntegrated) <- "RNA"
  if (isTRUE(append) & !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, append_pattern = "LIGER", show_warning = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

scMerge_integrate_deprecated <- function(srtList = NULL, srtMerge = NULL, batch = "orig.ident", append = TRUE,
                                         do_normalization = NULL, normalization_method = "logCPM",
                                         do_HVF_finding = TRUE, HVF_source = "separate", nHVF = 3000, hvf = NULL,
                                         do_scaling = TRUE, vars_to_regress = NULL,
                                         liner_reduction = "pca", liner_reduction_dims = 50, liner_reduction_dims_use = NULL, liner_reduction_distance = "cosine",
                                         nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "euclidean",
                                         cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                                         exogenous_genes = NULL, seed = 11, ...) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (!liner_reduction %in% c("pca", "ica", "nmf", "mds", "glmpca")) {
    stop("'liner_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!nonliner_reduction %in% c("umap", "tsne", "dm", "fdg", "isomap", "dbmap", "phate")) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm_index' must be one of 'louvain', 'slm', 'leiden'.")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  library(scMerge)

  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) & !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
    if (isTRUE(do_scaling) | (is.null(do_scaling) & nrow(GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")) != nrow(GetAssayData(srtMerge, slot = "data", assay = "RNA")))) {
      cat("Perform ScaleData on the srtMerge...\n")
      srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
    }
    scale_data <- GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
  }
  if (exists("scale_data")) {
    srtMerge <- SetAssayData(srtMerge, slot = "scale.data", assay = "RNA", new.data = scale_data)
  } else {
    cat("Perform ScaleData on the data...\n")
    srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
  }
  do_scaling <- NULL
  if (is.null(srtList) & !is.null(srtMerge)) {
    checked <- Check_srtMerge(
      srtMerge = srtMerge, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding, do_scaling = do_scaling,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
    do_scaling <- NULL
  }

  assay_use <- GetAssay(
    object = srtMerge,
    assay = switch(normalization_method,
      "logCPM" = "RNA",
      "SCT" = "SCT"
    )
  )
  sce <- SingleCellExperiment(
    assays = list(counts = assay_use@counts, logcounts = assay_use@data),
    colData = DataFrame(srtMerge@meta.data)
  )
  assay(sce, "counts") <- as(counts(sce), "dgeMatrix")
  assay(sce, "logcounts") <- as(logcounts(sce), "dgeMatrix")

  data("segList", package = "scMerge")
  scSEG <- segList$human$human_scSEG
  scSEG <- scSEG[scSEG %in% rownames(sce)]

  kmeansK <- sapply(srtList, function(x) {
    if (!"seurat_clusters" %in% colnames(x@meta.data)) {
      x <- x %>%
        RunPCA(npcs = maxPC, verbose = FALSE) %>%
        FindNeighbors() %>%
        FindClusters(resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
    }
    i <- nlevels(x[["seurat_clusters", drop = TRUE]])
  })

  sce_scMerge <- scMerge(
    sce_combine = sce,
    ctl = scSEG,
    kmeansK = kmeansK,
    batch_name = batch,
    assay_name = "scMerge",
    BSPARAM = IrlbaParam(),
    # BPPARAM = MulticoreParam(),
    plot_igraph = FALSE
  )
  assay(sce_scMerge, "counts") <- as(assay(sce_scMerge, "counts"), "dgCMatrix")
  assay(sce_scMerge, "scMerge") <- as(assay(sce_scMerge, "scMerge"), "dgCMatrix")
  srtIntegrated <- as.Seurat(
    x = sce_scMerge,
    counts = "counts",
    data = "scMerge"
  )
  srtIntegrated <- RenameAssays(srtIntegrated, originalexp = "scMerge")
  srtIntegrated@assays$RNA <- srtMerge@assays$RNA
  DefaultAssay(srtIntegrated) <- "scMerge"
  sce_scMerge <- sce <- NULL

  srtIntegrated <- Check_srtIntegrated(srtIntegrated, do_normalization = do_normalization, hvf = hvf, batch = batch, vars_to_regress = vars_to_regress)

  cat("Perform ScaleData on the data(after integrated)...\n")
  srtIntegrated <- Seurat::ScaleData(srtIntegrated, features = hvf, vars.to.regress = vars_to_regress)
  srtIntegrated <- RunPCA(
    object = srtIntegrated, npcs = maxPC, features = hvf,
    reduction.name = paste0("scMerge", liner_reduction),
    reduction.key = paste0("scMerge", liner_reduction),
    verbose = FALSE
  )
  if (is.null(dims)) {
    dims <- 1:ceiling(maxLikGlobalDimEst(data = Embeddings(srtIntegrated, reduction = paste0(reduction_prefix, "pca")), k = 20, iterations = 100)[["dim.est"]])
  }
  srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]] <- dims

  cat("Perform nonlinear dimension reduction on the data...\n")
  for (n in reduction_components) {
    if ("umap" %in% reduction) {
      srtIntegrated <- RunUMAP(
        object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"),
        reduction.name = paste0(reduction_prefix, "UMAP", n, "d"),
        reduction.key = paste0(reduction_prefix, "UMAP", n, "d_"),
        dims = dims, n.components = n, umap.method = "uwot",
        return.model = TRUE, verbose = FALSE
      )
    }
    if ("tsne" %in% reduction) {
      srtIntegrated <- RunTSNE(
        object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"),
        reduction.name = paste0(reduction_prefix, "TSNE", n, "d"),
        reduction.key = paste0(reduction_prefix, "TSNE", n, "d_"),
        dims = dims, dim.embed = n, tsne.method = "Rtsne",
        perplexity = max(ceiling(ncol(srtIntegrated) * 0.01), 30), max_iter = 2000,
        num_threads = 0, verbose = TRUE
      )
    }
  }

  cat("Perform FindClusters on the data...\n")
  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"), dims = dims, force.recalc = T, graph.name = paste0(reduction_prefix, "_", c("KNN", "SNN")), verbose = FALSE)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = resolution, algorithm = 1, graph.name = paste0(reduction_prefix, "_SNN"), verbose = FALSE)
  if (isTRUE(reorder)) {
    srtIntegrated <- SrtReorder(srtIntegrated, features = hvf, reorder_by = "seurat_clusters", slot = "data")
  }
  srtIntegrated[["seurat_clusters"]] <- NULL
  srtIntegrated[[paste0(reduction_prefix, "clusters")]] <- Idents(srtIntegrated)

  DefaultAssay(srtIntegrated) <- "RNA"
  if (isTRUE(append) & !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, append_pattern = "scMerge", show_warning = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

ZINBWaVE_integrate_deprecated <- function(srtList = NULL, srtMerge = NULL, batch = "orig.ident", append = TRUE,
                                          do_normalization = NULL, normalization_method = "logCPM",
                                          do_HVF_finding = TRUE, HVF_source = "separate", nHVF = 3000, hvf = NULL,
                                          do_scaling = TRUE, vars_to_regress = NULL,
                                          liner_reduction = "pca", liner_reduction_dims = 50, liner_reduction_dims_use = NULL, liner_reduction_distance = "cosine",
                                          nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "euclidean",
                                          cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                                          exogenous_genes = NULL, seed = 11, ...) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (!liner_reduction %in% c("pca", "ica", "nmf", "mds", "glmpca")) {
    stop("'liner_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!nonliner_reduction %in% c("umap", "tsne", "dm", "fdg", "isomap", "dbmap", "phate")) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm_index' must be one of 'louvain', 'slm', 'leiden'.")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  library(zinbwave)

  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) & !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
    if (isTRUE(do_scaling) | (is.null(do_scaling) & nrow(GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")) != nrow(GetAssayData(srtMerge, slot = "data", assay = "RNA")))) {
      cat("Perform ScaleData on the srtMerge...\n")
      srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
    }
    scale_data <- GetAssayData(srtMerge, slot = "scale.data", assay = "RNA")
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
  }
  if (exists("scale_data")) {
    srtMerge <- SetAssayData(srtMerge, slot = "scale.data", assay = "RNA", new.data = scale_data)
  } else {
    cat("Perform ScaleData on the data...\n")
    srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
  }
  do_scaling <- NULL
  if (is.null(srtList) & !is.null(srtMerge)) {
    checked <- Check_srtMerge(
      srtMerge = srtMerge, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding, do_scaling = do_scaling,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
    do_scaling <- NULL
  }

  sce <- as.SingleCellExperiment(srtMerge)
  assay(sce, "counts") <- as(counts(sce), "matrix")
  if (ncol(sce) < 10000) {
    sce_zinbwave <- zinbwave(
      Y = sce,
      K = 2,
      X = paste0("~", batch),
      which_assay = "counts",
      which_genes = hvf,
      epsilon = length(hvf),
      normalizedValues = TRUE,
      residuals = TRUE,
      BPPARAM = MulticoreParam()
    )
  } else {
    state <- 1
    while (state != 0) {
      tryCatch(expr = {
        sce_zinbwave <- zinbsurf(
          Y = sce,
          K = 2,
          X = paste0("~", batch),
          which_assay = "counts",
          which_genes = hvf,
          epsilon = length(hvf),
          prop_fit = 0.2,
          BPPARAM = MulticoreParam(workers = 8)
        )
        state <- 0
      }, error = function(error) {
        message(error)
        cat("Resampling the cells for zinbsurf ......\n")
        state <- state + 1
        if (state > 5) {
          stop("Resampling too many times. Stop the integration.\n")
        }
      })
    }
  }

  srtIntegrated <- as.Seurat(
    x = sce_zinbwave,
    counts = "counts",
    assay = "ZINBWaVE"
  )
  srtIntegrated@assays$RNA <- srtMerge@assays$RNA
  DefaultAssay(srtIntegrated) <- "ZINBWaVE"
  sce_zinbwave <- sce <- NULL

  srtIntegrated <- Check_srtIntegrated(srtIntegrated, do_normalization = do_normalization, hvf = hvf, batch = batch, vars_to_regress = vars_to_regress)

  cat("Perform NormalizeData on the data(after integrated)...\n")
  srtIntegrated <- NormalizeData(object = srtIntegrated, normalization.method = "LogNormalize")
  cat("Perform ScaleData on the data(after integrated)...\n")
  srtIntegrated <- Seurat::ScaleData(srtIntegrated, features = hvf, vars.to.regress = vars_to_regress)
  srtIntegrated <- RunPCA(
    object = srtIntegrated, npcs = maxPC, features = hvf,
    reduction.name = paste0(reduction_prefix, "pca"),
    reduction.key = paste0(reduction_prefix, "pca_"),
    verbose = FALSE
  )
  if (is.null(dims)) {
    dims <- 1:ceiling(maxLikGlobalDimEst(data = Embeddings(srtIntegrated, reduction = paste0(reduction_prefix, "pca")), k = 20, iterations = 100)[["dim.est"]])
  }
  srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]] <- dims

  cat("Perform nonlinear dimension reduction on the data...\n")
  for (n in reduction_components) {
    if ("umap" %in% reduction) {
      srtIntegrated <- RunUMAP(
        object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"),
        reduction.name = paste0(reduction_prefix, "UMAP", n, "d"),
        reduction.key = paste0(reduction_prefix, "UMAP", n, "d_"),
        dims = dims, n.components = n, umap.method = "uwot",
        return.model = TRUE, verbose = FALSE
      )
    }
    if ("tsne" %in% reduction) {
      srtIntegrated <- RunTSNE(
        object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"),
        reduction.name = paste0(reduction_prefix, "TSNE", n, "d"),
        reduction.key = paste0(reduction_prefix, "TSNE", n, "d_"),
        dims = dims, dim.embed = n, tsne.method = "Rtsne",
        perplexity = max(ceiling(ncol(srtIntegrated) * 0.01), 30), max_iter = 2000,
        num_threads = 0, verbose = TRUE
      )
    }
  }

  cat("Perform FindClusters on the data...\n")
  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"), dims = dims, force.recalc = T, graph.name = paste0(reduction_prefix, "_", c("KNN", "SNN")), verbose = FALSE)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = resolution, algorithm = 1, graph.name = paste0(reduction_prefix, "_SNN"), verbose = FALSE)
  if (isTRUE(reorder)) {
    srtIntegrated <- SrtReorder(srtIntegrated, features = hvf, reorder_by = "seurat_clusters", slot = "data")
  }
  srtIntegrated[["seurat_clusters"]] <- NULL
  srtIntegrated[[paste0(reduction_prefix, "clusters")]] <- Idents(srtIntegrated)

  DefaultAssay(srtIntegrated) <- "RNA"
  if (isTRUE(append) & !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, append_pattern = "ZINBWaVE", show_warning = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

Standard_SCP <- function(srt, prefix = "Standard",
                         do_normalization = NULL, normalization_method = "logCPM",
                         do_HVF_finding = TRUE, nHVF = 3000, hvf = NULL,
                         do_scaling = TRUE, vars_to_regress = NULL,
                         liner_reduction = "pca", liner_reduction_dims = 50, liner_reduction_distance = "cosine", liner_reduction_dims_use = NULL,
                         nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "euclidean",
                         cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                         exogenous_genes = NULL, seed = 11, ...) {
  if (class(srt) != "Seurat") {
    stop("'srt' is not a Seurat object.")
  }
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (!liner_reduction %in% c("pca", "ica", "nmf", "mds", "glmpca")) {
    stop("'liner_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!nonliner_reduction %in% c("umap", "tsne", "dm", "fdg", "isomap", "dbmap", "phate")) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm_index' must be one of 'louvain', 'slm', 'leiden'.")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  library(glmGamPoi)
  library(intrinsicDimension)
  library(dplyr)

  DefaultAssay(srt) <- "RNA"

  if (isTRUE(do_normalization) | (is.null(do_normalization) & identical(
    x = GetAssayData(srt, slot = "counts"),
    y = GetAssayData(srt, slot = "data")
  ))) {
    cat("Perform NormalizeData(logCPM) on the data...\n")
    srt <- NormalizeData(object = srt, normalization.method = "LogNormalize", verbose = FALSE)
  }
  if (is.null(hvf)) {
    if (isTRUE(do_HVF_finding) | (is.null(do_HVF_finding) & !"vst.variance.standardized" %in% colnames(srt@assays$RNA@meta.features))) {
      cat("Perform FindVariableFeatures on the data...\n")
      srt <- FindVariableFeatures(srt, verbose = FALSE)
    }
    m <- GetAssayData(srt, slot = "counts")
    gene_keep <- rownames(m)[Matrix::rowSums(m >= 1) >= 5]
    VariableFeatures(srt) <- hvf <- srt@assays$RNA@meta.features %>%
      filter(vst.variance.standardized > 1 &
        (!rownames(.) %in% exogenous_genes) &
        rownames(.) %in% gene_keep) %>%
      dplyr::arrange(desc(vst.variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
  } else {
    hvf <- hvf[hvf %in% rownames(GetAssayData(srt, slot = "counts"))]
    VariableFeatures(srt) <- hvf
  }
  if (isTRUE(do_scaling) | (is.null(do_scaling) & nrow(GetAssayData(srt, slot = "scale.data")) != nrow(GetAssayData(srt, slot = "data")))) {
    cat("Perform ScaleData on the data...\n")
    srt <- Seurat::ScaleData(object = srt, features = rownames(srt), vars.to.regress = vars_to_regress)
  }

  if (normalization_method %in% c("SCT")) {
    if (isTRUE(do_normalization) | isTRUE(do_HVF_finding) | (isTRUE(do_scaling)) | !"SCT" %in% Seurat::Assays(srt)) {
      cat("Perform SCTransform on the data...\n")
      srt <- SCTransform(
        object = srt,
        method = "glmGamPoi",
        variable.features.n = nHVF,
        vars.to.regress = vars_to_regress,
        return.only.var.genes = FALSE,
        min_cells = 5,
        assay = "RNA"
      )
    }
    DefaultAssay(srt) <- "SCT"
    if (is.null(hvf)) {
      if (!"residual_variance" %in% colnames(srt@assays$SCT@meta.features)) {
        feature.attr <- SCTResults(object = srt, slot = "feature.attributes")
        nfeatures <- min(nHVF, nrow(x = feature.attr))
        top.features <- rownames(x = feature.attr)[order(feature.attr$residual_variance,
          decreasing = TRUE
        )[1:nHVF]]
        VariableFeatures(object = srt) <- top.features
        srt@assays$SCT@meta.features <- feature.attr
      }
      m <- GetAssayData(srt, slot = "counts")
      gene_keep <- rownames(m)[Matrix::rowSums(m >= 1) >= 5]
      VariableFeatures(srt) <- hvf <- srt@assays$SCT@meta.features %>%
        filter((!rownames(.) %in% exogenous_genes) &
          rownames(.) %in% gene_keep) %>%
        dplyr::arrange(desc(residual_variance)) %>%
        rownames(.) %>%
        head(n = nHVF)
    } else {
      hvf <- hvf[hvf %in% rownames(GetAssayData(srt, slot = "counts"))]
      VariableFeatures(srt) <- hvf
    }
  }

  for (lr in liner_reduction) {
    cat("Perform linear dimension reduction (", lr, ") on the data...\n", sep = "")
    srt <- RunDimReduction(
      srt = srt, prefix = prefix, features = hvf, liner_reduction_distance = liner_reduction_distance,
      liner_reduction = lr, liner_reduction_dims = liner_reduction_dims,
      verbose = FALSE
    )
    if (is.null(liner_reduction_dims_use)) {
      liner_reduction_dims_use <- srt@reductions[[paste0(prefix, lr)]]@misc[["dims_estimate"]]
    }

    cat("Perform FindClusters (", cluster_algorithm, ") on the data...\n", sep = "")
    srt <- FindNeighbors(object = srt, reduction = paste0(prefix, lr), dims = liner_reduction_dims_use, force.recalc = T, graph.name = paste0(prefix, lr, "_", c("KNN", "SNN")), verbose = FALSE)
    srt <- FindClusters(object = srt, resolution = cluster_resolution, algorithm = cluster_algorithm_index, graph.name = paste0(prefix, lr, "_SNN"), verbose = FALSE)
    if (isTRUE(cluster_reorder)) {
      srt <- SrtReorder(srt, features = hvf, reorder_by = "seurat_clusters", slot = "data")
    }
    srt[["seurat_clusters"]] <- NULL
    srt[[paste0(prefix, lr, "clusters")]] <- Idents(srt)

    for (nr in nonliner_reduction) {
      cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
      for (n in nonliner_reduction_dims) {
        srt <- RunDimReduction(
          srt = srt, prefix = NULL,
          reduction_use = paste0(prefix, lr), dims_use = liner_reduction_dims_use,
          nonliner_reduction = nr, nonliner_reduction_dims = n,
          nonliner_reduction_distance = nonliner_reduction_distance,
          verbose = FALSE
        )
      }
    }
  }

  DefaultAssay(srt) <- "RNA"
  return(srt)
}

Integration_SCP <- function(srtList = NULL, srtMerge = NULL, append = TRUE, batch = "orig.ident",
                            integration_method = "Uncorrected",
                            do_normalization = NULL, normalization_method = "logCPM",
                            do_HVF_finding = TRUE, HVF_source = "separate", nHVF = 3000, hvf = NULL,
                            do_scaling = TRUE, vars_to_regress = NULL,
                            liner_reduction = "pca", liner_reduction_dims = 50, liner_reduction_distance = "cosine", liner_reduction_dims_use = NULL,
                            nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "euclidean",
                            cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                            exogenous_genes = NULL, seed = 11, ...) {
  if (is.null(srtList) & is.null(srtMerge)) {
    stop("Must be provided with one of the 'srtList' and 'srtMerge'")
  }
  if (length(integration_method) == 1 & integration_method %in% c("Uncorrected", "Seurat", "fastMNN", "Harmony", "Scanorama", "BBKNN", "CSS", "LIGER")) {
    args1 <- mget(names(formals()))
    args2 <- as.list(match.call())
    for (n in names(args2)) {
      args1[[n]] <- args2[[n]]
    }
    args1 <- args1[!names(args1) %in% c("integration_method", "...")]
    cat(paste0("[", Sys.time(), "] ", paste0("Start ", integration_method, "_integrate"), "...\n"))
    tryCatch(expr = {
      srtIntegrated <- base::do.call(
        what = paste0(integration_method, "_integrate"),
        args = args1
      )
    }, error = function(e) {
      message(e)
      message(cat("\n", paste0("[", Sys.time(), "]", " Stop the integration...\n")))
      stop()
    })
    cat(paste0("[", Sys.time(), "] ", paste0(integration_method, "_integrate done...\n")))
    if (exists("srtIntegrated")) {
      return(srtIntegrated)
    }
    if (!is.null(srtMerge)) {
      return(srtMerge)
    }
    return(NULL)
  } else {
    stop(paste(integration_method, "is not a suppoted integration method!"))
  }
}

RunDEtest <- function(srt, FindAllMarkers = TRUE, FindPairMarkers = FALSE,
                      group_by = NULL, cell_group1 = NULL, cell_group2 = NULL,
                      foldchange_threshold = 1.5, min_percent = 0.1,
                      BPPARAM = MulticoreParam(), force = FALSE, ...) {
  library(BiocParallel)
  library(dplyr)

  if (!is.null(cell_group1)) {
    if (!all(cell_group1) %in% colnames(srt)) {
      stop("cell_group1 has some cells not in the Seurat object.")
    }
    if (is.null(cell_group2)) {
      cell_group2 <- colnames(srt)[!colnames(srt) %in% cell_group1]
    }
    if (!all(cell_group2) %in% colnames(srt)) {
      stop("cell_group2 has some cells not in the Seurat object.")
    }
    cat("Perform FindMarkers(Wilcoxon) for custom cell groups...\n")
    markers <- FindMarkers(
      object = Seurat::Assays(srt, "RNA"), slot = "data",
      cells.1 = cell_group1,
      cells.2 = cell_group2,
      only.pos = T,
      logfc.threshold = log2(foldchange_threshold),
      min.pct = min_percent,
      test.use = "wilcox"
    )
    markers[, "gene"] <- rownames(markers)
    markers[, "group1"] <- "cell_group1"
    markers[, "group2"] <- "cell_group2"
    rownames(markers) <- NULL
    markers[, "group1"] <- factor(markers[, "group1"])
    markers[, "p_val_adj"] <- p.adjust(markers[, "p_val"], method = "BH")
    markers[, "DE_group_number"] <- as.integer(table(markers[["gene"]])[markers[, "gene"]])
    srt@misc[["DEtest_custom"]][["CellMarkers_Wilcoxon"]] <- markers
    srt@misc[["DEtest_custom"]][["cell_group1"]] <- cell_group1
    srt@misc[["DEtest_custom"]][["cell_group2"]] <- cell_group2
  }
  if (is.null(group_by)) {
    cell_group <- Idents(srt)
    group_by <- "active.ident"
  } else {
    cell_group <- srt[[group_by, drop = TRUE]]
  }
  if (!is.factor(cell_group)) {
    cell_group <- factor(cell_group, levels = unique(cell_group))
  }
  if (isTRUE(FindAllMarkers)) {
    cat("Perform FindAllMarkers(Wilcoxon)...\n")
    AllMarkers_Wilcoxon <- bplapply(levels(cell_group), FUN = function(group) {
      markers <- FindMarkers(
        object = Seurat::Assays(srt, "RNA"), slot = "data",
        cells.1 = names(cell_group)[cell_group == group],
        cells.2 = names(cell_group)[cell_group != group],
        only.pos = T,
        logfc.threshold = log2(foldchange_threshold),
        min.pct = min_percent,
        test.use = "wilcox"
      )
      if (nrow(markers) > 0) {
        markers[, "gene"] <- rownames(markers)
        markers[, "group1"] <- as.character(group)
        markers[, "group2"] <- "others"
        return(markers)
      } else {
        return(NULL)
      }
    }, BPPARAM = BPPARAM)
    AllMarkers_Wilcoxon <- dplyr::bind_rows(AllMarkers_Wilcoxon)
    rownames(AllMarkers_Wilcoxon) <- NULL
    AllMarkers_Wilcoxon[, "group1"] <- factor(AllMarkers_Wilcoxon[, "group1"], levels = levels(cell_group))
    AllMarkers_Wilcoxon[, "p_val_adj"] <- p.adjust(AllMarkers_Wilcoxon[, "p_val"], method = "BH")
    AllMarkers_Wilcoxon[, "DE_group_number"] <- as.integer(table(AllMarkers_Wilcoxon[["gene"]])[AllMarkers_Wilcoxon[, "gene"]])
    srt@misc[[paste0("DEtest_", group_by)]][["AllMarkers_Wilcoxon"]] <- AllMarkers_Wilcoxon
  }

  if (isTRUE(FindPairMarkers)) {
    if (nlevels(cell_group) > 30 & (!isTRUE(force))) {
      warning("Too many groups for FindPairMarkers function. If you want to force to run, set force=TRUE.")
      break
    }
    cat("Perform Wilcoxon FindPairMarkers...\n")
    pair <- expand.grid(x = levels(cell_group), y = levels(cell_group))
    pair <- pair[pair[, 1] != pair[, 2], ]
    PairMarkers_Wilcoxon <- bplapply(1:nrow(pair), function(i) {
      markers <- FindMarkers(
        object = Seurat::Assays(srt, "RNA"), slot = "data",
        cells.1 = names(cell_group)[cell_group == pair[i, 1]],
        cells.2 = names(cell_group)[cell_group == pair[i, 2]],
        only.pos = T,
        logfc.threshold = log2(foldchange_threshold),
        min.pct = min_percent,
        test.use = "wilcox"
      )
      if (nrow(markers) > 0) {
        markers[, "gene"] <- rownames(markers)
        markers[, "group1"] <- as.character(pair[i, 1])
        markers[, "group2"] <- as.character(pair[i, 2])
        return(markers)
      } else {
        return(NULL)
      }
    }, BPPARAM = BPPARAM)
    PairMarkers_Wilcoxon <- dplyr::bind_rows(PairMarkers_Wilcoxon)
    rownames(PairMarkers_Wilcoxon) <- NULL
    PairMarkers_Wilcoxon[, "group1"] <- factor(PairMarkers_Wilcoxon[, "group1"], levels = levels(cell_group))
    PairMarkers_Wilcoxon[, "p_val_adj"] <- p.adjust(PairMarkers_Wilcoxon[, "p_val"], method = "BH")
    PairMarkers_Wilcoxon[, "DE_group_number"] <- as.integer(table(PairMarkers_Wilcoxon[["gene"]])[PairMarkers_Wilcoxon[, "gene"]])
    PairMarkers_matrix <- as.data.frame.matrix(table(PairMarkers_Wilcoxon[, c("gene", "group1")]))
    PairMarkers_Wilcoxon[, "DE_group"] <- apply(PairMarkers_matrix, 1, function(x) {
      paste0(colnames(PairMarkers_matrix)[x > 0], collapse = ";")
    })[PairMarkers_Wilcoxon[, "gene"]]
    srt@misc[[paste0("DEtest_", group_by)]][["PairMarkers_Wilcoxon"]] <- PairMarkers_Wilcoxon
    srt@misc[[paste0("DEtest_", group_by)]][["PairMarkers_matrix"]] <- PairMarkers_matrix
  }
  return(srt)
}

srt_to_adata <- function(srt = NULL) {
  library(reticulate)
  library(Seurat)
  library(Matrix)
  if (!is.null(srt)) {
    sc <- import("scanpy", convert = FALSE)
    adata <- sc$AnnData(
      X = Matrix::t(GetAssayData(srt, assay = "RNA", slot = "counts")),
      obs = srt[[]],
      var = data.frame(features = rownames(srt))
    )

    layer_list <- list()
    for (assay in Seurat::Assays(srt)) {
      if (assay %in% c("spliced", "unspliced")) {
        layer_list[[assay]] <- t(GetAssayData(srt, assay = assay, slot = "counts"))
      }
    }
    adata$layers <- layer_list

    reduction_list <- list()
    for (reduction in Seurat::Reductions(srt)) {
      reduction_list[[paste0("X_", reduction)]] <- Seurat::Embeddings(srt, reduction = reduction)
    }
    adata$obsm <- reduction_list
    adata$var_names <- rownames(srt)
    return(adata)
  } else {
    stop("'srt' must be provided")
  }
}

kegg_get <- function(url) {
  temp <- tempfile()
  ntry <- 0
  status <- NULL
  while (is.null(status)) {
    ntry <- ntry + 1
    status <- tryCatch(expr = {
      download.file(url, destfile = temp, method = "wget", quiet = TRUE)
    }, error = function(e) {
      message("Get errors when connecting with KEGG...\nRetrying...")
      Sys.sleep(1)
      return(NULL)
    })
    if (ntry > 5) {
      stop("Stop connecting...")
    }
  }
  content <- readLines(temp) %>%
    strsplit(., "\t") %>%
    do.call("rbind", .)
  res <- data.frame(from = content[, 1], to = content[, 2])
  return(res)
}

RunEnrichment <- function(geneID, geneID_groups, IDtype = "symbol", species = "Homo_sapiens",
                          enrichment = "GO_BP", GO_simplify = TRUE, db_update = FALSE,
                          TERM2GENE = NULL, TERM2NAME = NULL,
                          show_plot = TRUE, topN = 6, pvalueCutoff = NULL, padjustCutoff = 0.05, palette = "Spectral",
                          BPPARAM = MulticoreParam()) {
  library(BiocParallel)
  library(stringr)
  library(taxize)
  library(clusterProfiler)
  library(GOSemSim)
  library(reactome.db)
  library(PFAM.db)
  library(GO.db)
  library(rWikiPathways)

  sp <- unlist(strsplit(species, split = "_"))
  org_sp <- paste0("org.", paste0(substring(sp, 1, 1), collapse = ""), ".eg.db")
  kegg_sp <- tolower(paste0(substring(sp[1], 1, 1), substring(sp[2], 1, 2), collapse = ""))
  complex_sp <- reactome_sp <- wiki_sp <- gsub(pattern = "_", replacement = " ", x = species)
  uid <- get_uid(complex_sp, messages = FALSE)
  common_name <- sci2comm(uid, db = "ncbi") %>% unlist()
  common_name <- paste(toupper(substr(common_name, 1, 1)), substr(common_name, 2, nchar(common_name)), sep = "")

  if (!require(org_sp, character.only = T)) {
    stop("No annoatation package found for the ", species)
  } else {
    library(org_sp, character.only = TRUE)
    orgdb <- get(org_sp)
    orgdbCHR <- get(paste0(gsub(pattern = ".db", "", org_sp), "CHR"))
  }

  if (is.factor(geneID_groups)) {
    geneID_groups <- factor(geneID_groups, levels = unique(geneID_groups))
  }
  names(geneID_groups) <- geneID

  ## GeneConvert
  if (IDtype != "entrez_id") {
    res <- GeneConvert(
      geneID = geneID,
      geneID_from_IDtype = IDtype,
      geneID_to_IDtype = "entrez_id",
      species_from = species,
      species_to = species
    )
    geneMap <- na.omit(res$geneID_res[, c("from_geneID", "to_geneID")])
    geneMap[, "geneID_groups"] <- geneID_groups[geneMap[, "from_geneID"]]
    geneID <- geneMap[, "to_geneID"]
    geneID_groups <- geneMap[, "geneID_groups"]
    names(geneID_groups) <- geneID
  }
  
  if (is.null(db_list)) {
    db_list <- list()
  }
  
  if (is.null(TERM2GENE)) {
    # TERM2GENE #####
    if (!enrichment%in%names(db_list)) {
      
    }
    ## GO ---------------------------------------------------------------------------
    if (enrichment %in% c("GO_BP", "GO_CC", "GO_MF")) {
      sub_enrichment <- unlist(strsplit(enrichment, split = "_"))[2]
      bg <- AnnotationDbi::select(orgdb, keys = keys(orgdb), columns = c("GOALL"))
      bg$EVIDENCEALL <- NULL
      bg <- unique(bg[!is.na(bg$GOALL), ])
      bg2 <- AnnotationDbi::select(GO.db, keys = keys(GO.db), columns = c("GOID", "TERM", "DEFINITION"))
      bg <- merge(x = bg, by.x = "GOALL", y = bg2, by.y = "GOID", all.x = T)
      TERM2GENE <- bg[which(bg$ONTOLOGYALL == sub_enrichment), c(1, 2)]
      TERM2NAME <- bg[which(bg$ONTOLOGYALL == sub_enrichment), c(1, 4)]
      semData <- godata(orgdb, ont = sub_enrichment)
    }

    ## KEGG ---------------------------------------------------------------------------
    if (enrichment %in% "KEGG") {
      kegg_db <- "pathway"
      kegg_pathwaygene_url <- paste0("http://rest.kegg.jp/link/", kegg_sp, "/", kegg_db, collapse = "")
      TERM2GENE <- kegg_get(kegg_pathwaygene_url)
      TERM2GENE[, 1] <- gsub(pattern = "[^:]+:", replacement = "", x = TERM2GENE[, 1])
      TERM2GENE[, 2] <- gsub(pattern = "[^:]+:", replacement = "", x = TERM2GENE[, 2])
      kegg_pathwayname_url <- paste0("http://rest.kegg.jp/list/", kegg_db, collapse = "")
      TERM2NAME <- kegg_get(kegg_pathwayname_url)
      TERM2NAME[, 1] <- gsub(pattern = "path:map", replacement = kegg_sp, TERM2NAME[, 1])
    }

    ## WikiPathwhay ---------------------------------------------------------------------------
    if (enrichment %in% "WikiPathwhay") {
      tempdir <- tempdir()
      gmt_files <- list.files(tempdir)[grep(".gmt", x = list.files(tempdir))]
      file.remove(gmt_files)
      downloadPathwayArchive(organism = wiki_sp, format = "gmt", date = "current", destpath = tempdir)
      wiki_gmt <- read.gmt(paste0(tempdir, "/", list.files(tempdir)[grep(".gmt", x = list.files(tempdir))]))
      wiki_gmt <- apply(wiki_gmt, 1, function(x) {
        wikiid <- str_split(x[["term"]], pattern = "%")[[1]][3]
        wikiterm <- str_split(x[["term"]], pattern = "%")[[1]][1]
        gmt <- x[["gene"]]
        data.frame(v0 = wikiid, v1 = gmt, v2 = wikiterm, stringsAsFactors = F)
      })
      bg <- bind_rows(wiki_gmt)
      TERM2GENE <- bg[, c(1, 2)]
      TERM2NAME <- bg[, c(1, 3)]
    }

    ## Reactome ---------------------------------------------------------------------------
    if (enrichment %in% "Reactome") {
      bg <- AnnotationDbi::select(reactome.db, keys = keys(reactome.db), columns = c("PATHID", "PATHNAME"))
      bg <- bg[str_detect(bg$PATHNAME, pattern = reactome_sp), ]
      bg <- na.omit(bg)
      bg$PATHNAME <- gsub(x = bg$PATHNAME, pattern = paste0("^", reactome_sp, ": "), replacement = "", perl = T)
      TERM2GENE <- bg[, c(2, 1)]
      TERM2NAME <- bg[, c(2, 3)]
    }

    ## Protein complex ---------------------------------------------------------------------------
    if (enrichment %in% c("Protein complex")) {
      temp <- tempfile()
      download.file("http://mips.helmholtz-muenchen.de/corum/download/coreComplexes.txt.zip", temp)
      df <- read.table(unz(temp, "coreComplexes.txt"), header = T, sep = "\t", stringsAsFactors = F, fill = T, quote = "")
      unlink(temp)
      df <- df[which(df$Organism == common_name), ]
      s <- strsplit(df$subunits.Entrez.IDs., split = ";")
      complex <- data.frame(V1 = rep(df$ComplexName, sapply(s, length)), V2 = unlist(s), V3 = rep(paste0("ComplexID:", df$ComplexID), sapply(s, length)))
      complex$V1 <- trimws(complex$V1) %>% gsub(pattern = "\\([^\\)]*\\)$", replacement = "", perl = FALSE)
      complex <- complex[!duplicated(complex), ]
      complex <- complex[which(complex$V2 != "None"), ]
      TERM2GENE <- complex[, c(3, 2)]
      TERM2NAME <- complex[, c(3, 1)]
    }

    ## PFAM ---------------------------------------------------------------------------
    if (enrichment %in% c("PFAM")) {
      bg <- AnnotationDbi::select(orgdb, keys = keys(orgdb), columns = "PFAM")
      bg <- unique(bg[!is.na(bg$PFAM), ])
      bg2 <- as.data.frame(PFAMDE2AC[mappedkeys(PFAMDE2AC)])
      rownames(bg2) <- bg2[["ac"]]
      bg[["PFAM_name"]] <- bg2[bg$PFAM, "de"]
      bg[is.na(bg[["PFAM_name"]]), "PFAM_name"] <- bg[is.na(bg[["PFAM_name"]]), "PFAM"]
      TERM2GENE <- complex[, c(2, 1)]
      TERM2NAME <- complex[, c(2, 3)]
    }

    ## Chromosome ---------------------------------------------------------------------------
    if (enrichment %in% c("Chromosome")) {
      chr <- as.data.frame(orgdbCHR[mappedkeys(orgdbCHR)])
      chr[, 2] <- paste0("chr", chr[, 2])
      TERM2GENE <- chr[, c(2, 1)]
      TERM2NAME <- chr[, c(2, 2)]
    }
  }

  res_list <- bplapply(unique(geneID_groups), FUN = function(group) {
    gene <- geneID[geneID_groups == group]
    enrich_res <- enricher(
      gene = gene,
      minGSSize = 10,
      maxGSSize = 500,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      TERM2GENE = TERM2GENE,
      TERM2NAME = TERM2NAME,
      universe = NULL
    )
    if (nrow(enrich_res@result > 0)) {
      if (isTRUE(GO_simplify) & enrichment %in% c("GO_BP", "GO_CC", "GO_MF")) {
        enrich_res@ontology <- sub_enrichment
        sim <- clusterProfiler::simplify(enrich_res,
          measure = "Rel",
          semData = semData
        )
        enrich_res@result <- sim@result
      }
      enrich_res <- setReadable(enrich_res, get(org_sp), keyType = "ENTREZID")
      enrich_res@result$Enrichment <- enrichment
      enrich_res@result$Groups <- group
      return(enrich_res@result)
    } else {
      return(NULL)
    }
  }, BPPARAM = BPPARAM)
  res <- bind_rows(res_list)
  p <- PlotEnrichment(
    x = res, topN = topN, pvalueCutoff = pvalueCutoff,
    padjustCutoff = padjustCutoff, palette = palette
  )
  if (isTRUE(show_plot)) {
    print(p)
  }
  return(list(enrichment = res, plot = p))
}

PlotEnrichment <- function(x, Enrichment = "Enrichment", Groups = "Groups", topN = 6,
                           pvalueCutoff = NULL, padjustCutoff = 0.05, palette = "Spectral") {
  if (is.null(pvalueCutoff) & is.null(padjustCutoff)) {
    stop("One of 'pvalueCutoff' or 'padjustCutoff' must be specified")
  }
  if (!is.factor(x[[Enrichment]])) {
    x[[Enrichment]] <- factor(x[[Enrichment]], levels = unique(x[[Enrichment]]))
  }
  df <- x %>%
    group_by(across(all_of(c(Enrichment, Groups)))) %>%
    filter(pvalue < pvalueCutoff || p.adjust < padjustCutoff) %>%
    top_n(-topN, wt = pvalue) %>%
    arrange(desc(pvalue)) %>%
    as.data.frame()
  y <- ifelse(is.null(padjustCutoff), "pvalue", "p.adjust")

  p <- ggplot(df, aes(
    x = .data[["Description"]], y = -log10(.data[[y]]),
    fill = .data[[Enrichment]],
    label = str_extract(GeneRatio, pattern = "\\d+(?=\\/)")
  )) +
    geom_bar(stat = "identity", color = "black") +
    geom_text(hjust = -0.3, size = 3.5) +
    scale_fill_manual(
      name = "Enrichment:",
      values = palette_zh(levels(df[[Enrichment]]), palette = palette),
      na.value = "grey90",
      guide = ifelse(length(unique(df[[Enrichment]])) > 1, "legend", "none")
    ) +
    scale_y_continuous(limits = c(0, 1.1 * max(-log10(df$p.adjust)))) +
    facet_wrap(formula(paste0(Enrichment, "~", Groups)), scales = "free", ncol = 2, dir = "v") +
    coord_flip() +
    theme_zh(aspect.ratio = 0.5)

  return(p)
}
