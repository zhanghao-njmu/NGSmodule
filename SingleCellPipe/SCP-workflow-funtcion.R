db_scDblFinder <- function(srt, db_rate = ncol(srt) / 1000 * 0.01, ...) {
  if (class(srt) != "Seurat") {
    stop("'srt' is not a Seurat object.",
      call. = FALSE
    )
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
    stop("'srt' is not a Seurat object.",
      call. = FALSE
    )
  }
  library(scds)
  sce <- as.SingleCellExperiment(srt, assay = "RNA")
  sce <- cxds_bcds_hybrid(sce)
  srt[["db.cxds_score"]] <- sce[["cxds_score"]]
  srt[["db.bcds_score"]] <- sce[["bcds_score"]]
  srt[["db.hybrid_score"]] <- sce[["hybrid_score"]]
  ntop <- ceiling(db_rate * ncol(sce))
  db_out <- names(sort(srt[[paste0("db.",method, "_score"), drop = T]], decreasing = T)[1:ntop])
  srt[[paste0("db.scds_", method, "_class")]] <- "singlet"
  srt[[paste0("db.scds_", method, "_class")]][db_out, ] <- "doublet"
  return(srt)
}

db_Scrublet <- function(srt, db_rate = ncol(srt) / 1000 * 0.01, ...) {
  if (class(srt) != "Seurat") {
    stop("'srt' is not a Seurat object.",
      call. = FALSE
    )
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
    stop("'srt' is not a Seurat object.",
      call. = FALSE
    )
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
    stop("'srt' is not a Seurat object.",
      call. = FALSE
    )
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
    stop(paste(db_method, "is not a suppoted doublet-calling method!"),
      call. = FALSE
    )
  }
}

Check_srtList <- function(srtList, do_normalization = NULL, normalization_method, vars_to_regress = NULL,
                          HVF_source = "separate", nHVF = 3000, hvf = NULL,
                          exogenous_genes = NULL, ...) {
  cat(paste0("[", Sys.time(), "]", " Checking srtList... ...\n"))
  library(Seurat)
  library(sctransform)
  library(glmGamPoi)

  if (class(srtList) != "list" | any(sapply(srtList, class) != "Seurat")) {
    stop("'srtList' is not a list of Seurat object.",
      call. = FALSE
    )
  }
  if (!normalization_method %in% c("logCPM", "SCT")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT'",
      call. = FALSE
    )
  }
  if (!HVF_source %in% c("global", "separate")) {
    stop("'HVF_source' must be one of: 'global','separate'",
      call. = FALSE
    )
  }

  genelist <- lapply(srtList, function(x) {
    sort(rownames(GetAssayData(x, slot = "counts", assay = "RNA")))
  })
  if (length(unique(genelist)) != 1) {
    warning("'srtList' have different feature names! Will subset the common features for downstream analysis!",
      call. = FALSE
    )
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
    if (!"RNA" %in% Seurat::Assays(srtList[[i]])) {
      stop(paste("srtList", i, "does not contain 'RNA' assay."),
        call. = FALSE
      )
    }
    DefaultAssay(srtList[[i]]) <- "RNA"
    if (isTRUE(do_normalization) | (is.null(do_normalization) & identical(
      x = GetAssayData(srtList[[i]], slot = "counts"),
      y = GetAssayData(srtList[[i]], slot = "data")
    ))) {
      cat("Perform NormalizeData(logCPM) on the data", i, "of the srtList...\n")
      srtList[[i]] <- NormalizeData(object = srtList[[i]], normalization.method = "LogNormalize", verbose = FALSE)
    }
    if (!"vst.variance.standardized" %in% colnames(srtList[[i]]@assays$RNA@meta.features)) {
      cat("Perform FindVariableFeatures on the data", i, "of the srtList...\n")
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
    DefaultAssay(srtList[[i]]) <- "RNA"

    if (normalization_method %in% c("SCT")) {
      if (!"SCT" %in% Seurat::Assays(srtList[[i]])) {
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

Check_srtMerge <- function(srtMerge, do_normalization = NULL, normalization_method, vars_to_regress = NULL, batch,
                           HVF_source = "separate", nHVF = 3000, hvf = NULL,
                           exogenous_genes = NULL, ...) {
  if (class(srtMerge) != "Seurat") {
    stop("'srtMerge' is not a Seurat object.",
      call. = FALSE
    )
  }
  if (length(batch) != 1) {
    stop("batch must be a vector to specify the batch column in srtMerge object!",
      call. = FALSE
    )
  }
  if (!batch %in% colnames(srtMerge@meta.data)) {
    stop(paste0("No batch column('", batch, "') found in the srtMerge object!"),
      call. = FALSE
    )
  }
  if (!is.factor(srtMerge[[batch, drop = T]])) {
    srtMerge[[batch, drop = T]] <- factor(srtMerge[[batch, drop = T]], levels = unique(srtMerge[[batch, drop = T]]))
  }
  cat("Spliting srtMerge into srtList... ...\n")
  srtList <- SplitObject(object = srtMerge, split.by = batch)

  checked <- Check_srtList(srtList,
    do_normalization = do_normalization,
    normalization_method = normalization_method, vars_to_regress = vars_to_regress,
    HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
    exogenous_genes = exogenous_genes
  )
  srtList <- checked[["srtList"]]
  hvf <- checked[["hvf"]]
  srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
  VariableFeatures(srtMerge) <- hvf

  return(list(
    srtMerge = srtMerge,
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
    srtIntegrated <- Seurat::ScaleData(object = srtIntegrated, features = rownames(srtIntegrated), vars.to.regress = vars_to_regress)
  }

  DefaultAssay(object = srtIntegrated) <- raw_DefaultAssay
  srtIntegrated@project.name <- paste0(unique(srtIntegrated[[batch, drop = TRUE]]), collapse = ",")
  srtIntegrated[[batch]] <- factor(srtIntegrated[[batch, drop = TRUE]],
    levels = unique(srtIntegrated[[batch, drop = TRUE]])
  )
  return(srtIntegrated)
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
  from_IDtype <- sapply(geneID_from_IDtype, function(x) {
    switch(x,
      "ensembl_symbol" = "external_gene_name",
      "ensembl_id" = "ensembl_gene_id",
      "entrez_symbol" = "entrezgene_accession",
      "entrez_id" = "entrezgene_id"
    )
  })
  to_IDtype <- switch(geneID_to_IDtype,
    "ensembl_symbol" = "external_gene_name",
    "ensembl_id" = "ensembl_gene_id"
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
  # web <- read_html(httr::RETRY("GET", "http://www.ensembl.org/info/website/archives/index.html?redirect=no", times = 1000, timeout(1000)))
  # urls <- web %>% html_nodes("ul") %>% html_nodes("strong") %>% html_nodes("a") %>% html_attr("href")
  # version <- web %>% html_nodes("ul") %>% html_nodes("strong") %>% html_nodes("a") %>% html_text(trim = TRUE) %>%
  #   gsub(pattern = "(Ensembl )|(:.*)",replacement = "",x = .,perl = T)
  # archives <- data.frame(version=version,url=urls,stringsAsFactors = F)

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
    sep <- ""
  } else {
    sep <- "_homolog_"
    to_IDtype <- switch(to_IDtype,
      "external_gene_name" = "associated_gene_name",
      "ensembl_gene_id" = "ensembl_gene"
    )
  }
  to_attr <- paste(species_to_simp, to_IDtype, sep = sep)

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
          colnames(geneID_res1) <- c("from_geneID", "ensembl_gene_id")
          from_name <- geneID_from_IDtype[which(from_IDtype == from_attr)]
          geneID_res1[, "from_IDtype"] <- from_name

          geneID_res2 <- getBM(
            mart = mart,
            attributes = c("ensembl_gene_id", to_attr),
            filters = "ensembl_gene_id",
            values = list(geneID_res1[, "ensembl_gene_id"])
          )
          colnames(geneID_res2) <- c("ensembl_gene_id", "to_geneID")
          geneID_res2[, "to_IDtype"] <- geneID_to_IDtype

          geneID_res_merge <- merge(x = geneID_res1, y = geneID_res2, by = "ensembl_gene_id")

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
        geneID_res1 <- getBM(
          mart = mart,
          attributes = c("ensembl_gene_id", from_attr, to_attr),
          filters = from_attr,
          values = list(geneID)
        )
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
  geneID_sim <- geneID_res %>%
    group_by(from_geneID) %>%
    mutate(from_geneID = unique(from_geneID), to_geneID = list(unique(to_geneID[to_geneID != ""])))
  geneID_sim <- unique(as.data.frame(geneID_sim[,c("from_geneID", "to_geneID")]))
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
      geneID_from_IDtype = c("ensembl_symbol", "entrez_symbol"),
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

Uncorrected_integrate <- function(srtList = NULL, srtMerge = NULL, append = FALSE,
                                  do_normalization = NULL, normalization_method = "logCPM", batch = "orig.ident",
                                  vars_to_regress = NULL,
                                  HVF_source = "separate", nHVF = 3000, hvf = NULL,
                                  maxPC = 100, dims = NULL, resolution = 0.8, reorder = TRUE,
                                  reduction = "umap", reduction_prefix = "Uncorrected", reduction_components = 2:3,
                                  exogenous_genes = NULL, seed = 11, ...) {
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
  if (isTRUE(append)) {
    if (is.null(srtMerge)) {
      stop("srtMerge must be provided when 'append = TRUE'.")
    } else {
      srtMerge_raw <- srtMerge
    }
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(srtList,
      do_normalization = do_normalization,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
  }
  if (is.null(srtList) & !is.null(srtMerge)) {
    checked <- Check_srtMerge(srtMerge,
      do_normalization = do_normalization,
      normalization_method = normalization_method, batch = batch, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = NULL
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
  }

  srtIntegrated <- Standard_SCP(
    srt = srtMerge, normalization_method = normalization_method, nHVF = nHVF, hvf = hvf, vars_to_regress = vars_to_regress,
    maxPC = maxPC, dims = dims, reduction_components = reduction_components, resolution = resolution, reduction = reduction, reduction_prefix = reduction_prefix, reorder = reorder,
    exogenous_genes = exogenous_genes, seed = seed
  )
  srtIntegrated@project.name <- paste0(unique(srtIntegrated[[batch, drop = TRUE]]), collapse = ",")

  if (isTRUE(append)) {
    DefaultAssay(srtMerge_raw) <- "RNA"
    srtMerge_raw@assays$RNA <- srtIntegrated@assays$RNA
    srtMerge_raw@graphs[paste0(reduction_prefix, "_", c("NN", "SNN"))] <- srtIntegrated@graphs[paste0(reduction_prefix, "_", c("NN", "SNN"))]
    srtMerge_raw[[paste0(reduction_prefix, "pca")]] <- srtIntegrated[[paste0(reduction_prefix, "pca")]]
    srtMerge_raw@misc[[paste0(reduction_prefix, "Dims")]] <- srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]]
    srtMerge_raw@misc[[paste0(reduction_prefix, "HVF")]] <- hvf
    srtMerge_raw[[paste0(reduction_prefix, "clusters")]] <- srtIntegrated[[paste0(reduction_prefix, "clusters")]]
    for (n in reduction_components) {
      for (i in reduction) {
        srtMerge_raw[[paste0(reduction_prefix, toupper(i), n, "d")]] <- srtIntegrated[[paste0(reduction_prefix, toupper(i), n, "d")]]
      }
    }
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

Seurat_integrate <- function(srtList = NULL, srtMerge = NULL, append = FALSE,
                             do_normalization = NULL, normalization_method = "logCPM", batch = "orig.ident",
                             vars_to_regress = NULL,
                             HVF_source = "separate", nHVF = 3000, hvf = NULL,
                             maxPC = 100, dims = NULL, resolution = 0.8, reorder = TRUE,
                             reduction = "umap", reduction_prefix = "Seurat", reduction_components = 2:3,
                             exogenous_genes = NULL, seed = 11, ...) {
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
  if (isTRUE(append)) {
    if (is.null(srtMerge)) {
      stop("srtMerge must be provided when 'append = TRUE'.")
    } else {
      srtMerge_raw <- srtMerge
    }
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(srtList,
      do_normalization = do_normalization,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
  }
  if (is.null(srtList) & !is.null(srtMerge)) {
    srtList <- SplitObject(object = srtMerge, split.by = batch)
    checked <- Check_srtList(srtList,
      do_normalization = do_normalization,
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
    anchor.features = hvf
  )
  srtIntegrated <- IntegrateData(
    anchorset = srt_anchors,
    normalization.method = switch(normalization_method,
      "logCPM" = "LogNormalize",
      "SCT" = "SCT"
    ),
    features.to.integrate = c(
      hvf,
      Reduce(union, lapply(srtList, VariableFeatures)),
      Reduce(intersect, c(
        lapply(srtList, rownames),
        list(c(cc_S_genes, cc_G2M_genes))
      ))
    )
  )
  srtIntegrated <- RenameAssays(object = srtIntegrated, integrated = "Seurat")
  DefaultAssay(srtIntegrated) <- "Seurat"

  srtIntegrated <- Check_srtIntegrated(srtIntegrated, do_normalization = do_normalization, hvf = hvf, batch = batch, vars_to_regress = vars_to_regress)

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

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"), dims = dims, force.recalc = T, graph.name = paste0(reduction_prefix, "_", c("NN", "SNN")), verbose = FALSE)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000, graph.name = paste0(reduction_prefix, "_SNN"), verbose = FALSE)
  if (isTRUE(reorder)) {
    srtIntegrated <- SrtReorder(srtIntegrated, features = hvf, reorder_by = "seurat_clusters", slot = "data")
  }
  srtIntegrated[["seurat_clusters"]] <- NULL
  srtIntegrated[[paste0(reduction_prefix, "clusters")]] <- Idents(srtIntegrated)

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

  DefaultAssay(srtIntegrated) <- "RNA"
  if (isTRUE(append)) {
    DefaultAssay(srtMerge_raw) <- "RNA"
    srtMerge_raw@assays$RNA <- srtIntegrated@assays$RNA
    srtMerge_raw@assays$Seurat <- srtIntegrated@assays$Seurat
    srtMerge_raw@graphs[paste0(reduction_prefix, "_", c("NN", "SNN"))] <- srtIntegrated@graphs[paste0(reduction_prefix, "_", c("NN", "SNN"))]
    srtMerge_raw[[paste0(reduction_prefix, "pca")]] <- srtIntegrated[[paste0(reduction_prefix, "pca")]]
    srtMerge_raw@misc[[paste0(reduction_prefix, "Dims")]] <- srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]]
    srtMerge_raw@misc[[paste0(reduction_prefix, "HVF")]] <- hvf
    srtMerge_raw[[paste0(reduction_prefix, "clusters")]] <- srtIntegrated[[paste0(reduction_prefix, "clusters")]]
    for (n in reduction_components) {
      for (i in reduction) {
        srtMerge_raw[[paste0(reduction_prefix, toupper(i), n, "d")]] <- srtIntegrated[[paste0(reduction_prefix, toupper(i), n, "d")]]
      }
    }
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

fastMNN_integrate <- function(srtList = NULL, srtMerge = NULL, append = FALSE,
                              do_normalization = NULL, normalization_method = "logCPM", batch = "orig.ident",
                              vars_to_regress = NULL,
                              HVF_source = "separate", nHVF = 3000, hvf = NULL,
                              maxPC = 100, dims = NULL, resolution = 0.8, reorder = TRUE,
                              reduction = "umap", reduction_prefix = "fastMNN", reduction_components = 2:3,
                              exogenous_genes = NULL, seed = 11, ...) {
  set.seed(seed)
  library(SeuratWrappers)

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
  if (isTRUE(append)) {
    if (is.null(srtMerge)) {
      stop("srtMerge must be provided when 'append = TRUE'.")
    } else {
      srtMerge_raw <- srtMerge
    }
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(srtList,
      do_normalization = do_normalization,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
  }
  if (is.null(srtList) & !is.null(srtMerge)) {
    srtList <- SplitObject(object = srtMerge, split.by = batch)
    checked <- Check_srtList(srtList,
      do_normalization = do_normalization,
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
    BPPARAM = MulticoreParam()
  )

  srtIntegrated <- Check_srtIntegrated(srtIntegrated, do_normalization = do_normalization, hvf = hvf, batch = batch, vars_to_regress = vars_to_regress)

  if (is.null(dims)) {
    dims <- 1:ceiling(maxLikGlobalDimEst(data = Embeddings(srtIntegrated, reduction = "fastMNN"), k = 20, iterations = 100)[["dim.est"]])
  }
  srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]] <- dims

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "fastMNN", dims = dims, force.recalc = T, graph.name = paste0(reduction_prefix, "_", c("NN", "SNN")), verbose = FALSE)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000, graph.name = paste0(reduction_prefix, "_SNN"), verbose = FALSE)
  if (isTRUE(reorder)) {
    srtIntegrated <- SrtReorder(srtIntegrated, features = hvf, reorder_by = "seurat_clusters", slot = "data")
  }
  srtIntegrated[["seurat_clusters"]] <- NULL
  srtIntegrated[[paste0(reduction_prefix, "clusters")]] <- Idents(srtIntegrated)

  cat("Perform nonlinear dimension reduction on the data...\n")
  for (n in reduction_components) {
    if ("umap" %in% reduction) {
      srtIntegrated <- RunUMAP(
        object = srtIntegrated, reduction = "fastMNN",
        reduction.name = paste0(reduction_prefix, "UMAP", n, "d"),
        reduction.key = paste0(reduction_prefix, "UMAP", n, "d_"),
        dims = dims, n.components = n, umap.method = "uwot",
        return.model = TRUE, verbose = FALSE
      )
    }
    if ("tsne" %in% reduction) {
      srtIntegrated <- RunTSNE(
        object = srtIntegrated, reduction = "fastMNN",
        reduction.name = paste0(reduction_prefix, "TSNE", n, "d"),
        reduction.key = paste0(reduction_prefix, "TSNE", n, "d_"),
        dims = dims, dim.embed = n, tsne.method = "Rtsne",
        perplexity = max(ceiling(ncol(srtIntegrated) * 0.01), 30), max_iter = 2000,
        num_threads = 0, verbose = TRUE
      )
    }
  }

  DefaultAssay(srtIntegrated) <- "RNA"
  if (isTRUE(append)) {
    DefaultAssay(srtMerge_raw) <- "RNA"
    srtMerge_raw@assays$RNA <- srtIntegrated@assays$RNA
    srtMerge_raw@reductions$fastMNN <- srtIntegrated@reductions$fastMNN
    srtMerge_raw@graphs[paste0(reduction_prefix, "_", c("NN", "SNN"))] <- srtIntegrated@graphs[paste0(reduction_prefix, "_", c("NN", "SNN"))]
    srtMerge_raw@misc[[paste0(reduction_prefix, "Dims")]] <- srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]]
    srtMerge_raw@misc[[paste0(reduction_prefix, "HVF")]] <- hvf
    srtMerge_raw[[paste0(reduction_prefix, "clusters")]] <- srtIntegrated[[paste0(reduction_prefix, "clusters")]]
    for (n in reduction_components) {
      for (i in reduction) {
        srtMerge_raw[[paste0(reduction_prefix, toupper(i), n, "d")]] <- srtIntegrated[[paste0(reduction_prefix, toupper(i), n, "d")]]
      }
    }
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

Harmony_integrate <- function(srtList = NULL, srtMerge = NULL, append = FALSE,
                              do_normalization = NULL, normalization_method = "logCPM", batch = "orig.ident",
                              vars_to_regress = NULL,
                              HVF_source = "separate", nHVF = 3000, hvf = NULL,
                              maxPC = 100, dims = NULL, resolution = 0.8, reorder = TRUE,
                              reduction = "umap", reduction_prefix = "Harmony", reduction_components = 2:3,
                              exogenous_genes = NULL, seed = 11, ...) {
  set.seed(seed)
  library(SeuratWrappers)

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
  if (isTRUE(append)) {
    if (is.null(srtMerge)) {
      stop("srtMerge must be provided when 'append = TRUE'.")
    } else {
      srtMerge_raw <- srtMerge
    }
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(srtList,
      do_normalization = do_normalization,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
  }
  if (is.null(srtList) & !is.null(srtMerge)) {
    checked <- Check_srtMerge(srtMerge,
      do_normalization = do_normalization,
      normalization_method = normalization_method, batch = batch, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = NULL
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
  }
  if (nrow(GetAssayData(srtMerge, slot = "scale.data")) != nrow(GetAssayData(srtMerge, slot = "data"))) {
    cat("Perform ScaleData on the data...\n")
    srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
  }
  srtMerge <- RunPCA(
    object = srtMerge, npcs = maxPC, features = hvf,
    reduction.name = paste0(reduction_prefix, "pca"),
    reduction.key = paste0(reduction_prefix, "pca_"),
    verbose = FALSE
  )

  srtIntegrated <- RunHarmony(
    object = srtMerge,
    group.by.vars = batch,
    reduction = paste0(reduction_prefix, "pca"),
    reduction.save = "Harmony",
    assay.use = DefaultAssay(srtMerge)
  )
  srtMerge <- NULL

  srtIntegrated <- Check_srtIntegrated(srtIntegrated, do_normalization = do_normalization, hvf = hvf, batch = batch, vars_to_regress = vars_to_regress)

  if (is.null(dims)) {
    dims <- 1:ceiling(maxLikGlobalDimEst(data = Embeddings(srtIntegrated, reduction = "Harmony"), k = 20, iterations = 100)[["dim.est"]])
  }
  srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]] <- dims

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "Harmony", dims = dims, force.recalc = T, graph.name = paste0(reduction_prefix, "_", c("NN", "SNN")), verbose = FALSE)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000, graph.name = paste0(reduction_prefix, "_SNN"), verbose = FALSE)
  if (isTRUE(reorder)) {
    srtIntegrated <- SrtReorder(srtIntegrated, features = hvf, reorder_by = "seurat_clusters", slot = "data")
  }
  srtIntegrated[["seurat_clusters"]] <- NULL
  srtIntegrated[[paste0(reduction_prefix, "clusters")]] <- Idents(srtIntegrated)

  cat("Perform nonlinear dimension reduction on the data...\n")
  for (n in reduction_components) {
    if ("umap" %in% reduction) {
      srtIntegrated <- RunUMAP(
        object = srtIntegrated, reduction = "Harmony",
        reduction.name = paste0(reduction_prefix, "UMAP", n, "d"),
        reduction.key = paste0(reduction_prefix, "UMAP", n, "d_"),
        dims = dims, n.components = n, umap.method = "uwot",
        return.model = TRUE, verbose = FALSE
      )
    }
    if ("tsne" %in% reduction) {
      srtIntegrated <- RunTSNE(
        object = srtIntegrated, reduction = "Harmony",
        reduction.name = paste0(reduction_prefix, "TSNE", n, "d"),
        reduction.key = paste0(reduction_prefix, "TSNE", n, "d_"),
        dims = dims, dim.embed = n, tsne.method = "Rtsne",
        perplexity = max(ceiling(ncol(srtIntegrated) * 0.01), 30), max_iter = 2000,
        num_threads = 0, verbose = TRUE
      )
    }
  }

  DefaultAssay(srtIntegrated) <- "RNA"
  if (isTRUE(append)) {
    DefaultAssay(srtMerge_raw) <- "RNA"
    srtMerge_raw@assays$RNA <- srtIntegrated@assays$RNA
    srtMerge_raw@reductions$Harmony <- srtIntegrated@reductions$Harmony
    srtMerge_raw@graphs[paste0(reduction_prefix, "_", c("NN", "SNN"))] <- srtIntegrated@graphs[paste0(reduction_prefix, "_", c("NN", "SNN"))]
    srtMerge_raw@misc[[paste0(reduction_prefix, "Dims")]] <- srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]]
    srtMerge_raw@misc[[paste0(reduction_prefix, "HVF")]] <- hvf
    srtMerge_raw[[paste0(reduction_prefix, "pca")]] <- srtIntegrated[[paste0(reduction_prefix, "pca")]]
    srtMerge_raw[[paste0(reduction_prefix, "clusters")]] <- srtIntegrated[[paste0(reduction_prefix, "clusters")]]
    for (n in reduction_components) {
      for (i in reduction) {
        srtMerge_raw[[paste0(reduction_prefix, toupper(i), n, "d")]] <- srtIntegrated[[paste0(reduction_prefix, toupper(i), n, "d")]]
      }
    }
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

Scanorama_integrate <- function(srtList = NULL, srtMerge = NULL, append = FALSE,
                                do_normalization = NULL, normalization_method = "logCPM", batch = "orig.ident",
                                vars_to_regress = NULL,
                                HVF_source = "separate", nHVF = 3000, hvf = NULL,
                                maxPC = 100, dims = NULL, resolution = 0.8, reorder = TRUE,
                                reduction = "umap", reduction_prefix = "Scanorama", reduction_components = 2:3,
                                exogenous_genes = NULL, seed = 11, ...) {
  set.seed(seed)
  library(reticulate)
  library(plyr)
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
  if (isTRUE(append)) {
    if (is.null(srtMerge)) {
      stop("srtMerge must be provided when 'append = TRUE'.")
    } else {
      srtMerge_raw <- srtMerge
    }
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(srtList,
      do_normalization = do_normalization,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
  }
  if (is.null(srtList) & !is.null(srtMerge)) {
    srtList <- SplitObject(object = srtMerge, split.by = batch)
    checked <- Check_srtList(srtList,
      do_normalization = do_normalization,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
  }

  assaylist <- list()
  genelist <- list()
  for (i in 1:length(srtList)) {
    assaylist[[i]] <- t(as.matrix(GetAssayData(object = srtList[[i]], slot = "data")))
    genelist[[i]] <- rownames(srtList[[i]])
  }

  cat("Perform scanorama$correct... ...\n")
  integrated.corrected.data <- scanorama$correct(
    datasets_full = assaylist,
    genes_list = genelist,
    return_dimred = TRUE,
    return_dense = TRUE,
    union = FALSE
  )

  cor_value <- integrated.corrected.data[[2]] %>%
    rbind.fill.matrix() %>%
    t()
  rownames(cor_value) <- integrated.corrected.data[[3]]
  colnames(cor_value) <- unlist(sapply(assaylist, rownames))

  dim_reduction <- integrated.corrected.data[[1]] %>% rbind.fill.matrix()
  rownames(dim_reduction) <- unlist(sapply(assaylist, rownames))
  colnames(dim_reduction) <- paste0("Scanorama_", 1:100)
  stdevs <- apply(dim_reduction, MARGIN = 2, FUN = sd)

  srtIntegrated <- Reduce(function(x, y) merge(x, y), srtList)
  srtIntegrated[["Scanorama"]] <- CreateAssayObject(data = cor_value)
  srtIntegrated[["Scanorama_reduction"]] <- CreateDimReducObject(embeddings = dim_reduction, assay = "Scanorama", stdev = stdevs, key = "Scanorama_")
  DefaultAssay(srtIntegrated) <- "Scanorama"

  srtIntegrated <- Check_srtIntegrated(srtIntegrated, do_normalization = do_normalization, hvf = hvf, batch = batch, vars_to_regress = vars_to_regress)

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

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"), dims = dims, force.recalc = T, graph.name = paste0(reduction_prefix, "_", c("NN", "SNN")), verbose = FALSE)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000, graph.name = paste0(reduction_prefix, "_SNN"), verbose = FALSE)
  if (isTRUE(reorder)) {
    srtIntegrated <- SrtReorder(srtIntegrated, features = hvf, reorder_by = "seurat_clusters", slot = "data")
  }
  srtIntegrated[["seurat_clusters"]] <- NULL
  srtIntegrated[[paste0(reduction_prefix, "clusters")]] <- Idents(srtIntegrated)

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

  DefaultAssay(srtIntegrated) <- "RNA"
  if (isTRUE(append)) {
    DefaultAssay(srtMerge_raw) <- "RNA"
    srtMerge_raw@assays$RNA <- srtIntegrated@assays$RNA
    srtMerge_raw@assays$Scanorama <- srtIntegrated@assays$Scanorama
    srtMerge_raw@graphs[paste0(reduction_prefix, "_", c("NN", "SNN"))] <- srtIntegrated@graphs[paste0(reduction_prefix, "_", c("NN", "SNN"))]
    srtMerge_raw@reductions$Scanorama_reduction <- srtIntegrated@reductions$Scanorama_reduction
    srtMerge_raw[[paste0(reduction_prefix, "pca")]] <- srtIntegrated[[paste0(reduction_prefix, "pca")]]
    srtMerge_raw@misc[[paste0(reduction_prefix, "Dims")]] <- srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]]
    srtMerge_raw@misc[[paste0(reduction_prefix, "HVF")]] <- hvf
    srtMerge_raw[[paste0(reduction_prefix, "clusters")]] <- srtIntegrated[[paste0(reduction_prefix, "clusters")]]
    for (n in reduction_components) {
      for (i in reduction) {
        srtMerge_raw[[paste0(reduction_prefix, toupper(i), n, "d")]] <- srtIntegrated[[paste0(reduction_prefix, toupper(i), n, "d")]]
      }
    }
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

BBKNN_integrate <- function(srtList = NULL, srtMerge = NULL, append = FALSE,
                            do_normalization = NULL, normalization_method = "logCPM", batch = "orig.ident",
                            vars_to_regress = NULL,
                            HVF_source = "separate", nHVF = 3000, hvf = NULL,
                            maxPC = 100, dims = NULL, resolution = 0.8, reorder = TRUE,
                            reduction_prefix = "BBKNN", reduction_components = 2:3,
                            exogenous_genes = NULL, seed = 11, ...) {
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
  if (isTRUE(append)) {
    if (is.null(srtMerge)) {
      stop("srtMerge must be provided when 'append = TRUE'.")
    } else {
      srtMerge_raw <- srtMerge
    }
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(srtList,
      do_normalization = do_normalization,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
  }
  if (is.null(srtList) & !is.null(srtMerge)) {
    checked <- Check_srtMerge(srtMerge,
      do_normalization = do_normalization,
      normalization_method = normalization_method, batch = batch, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = NULL
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
  }
  if (nrow(GetAssayData(srtMerge, slot = "scale.data")) != nrow(GetAssayData(srtMerge, slot = "data"))) {
    cat("Perform ScaleData on the data...\n")
    srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
  }
  srtMerge <- RunPCA(
    object = srtMerge, npcs = maxPC, features = hvf,
    reduction.name = paste0(reduction_prefix, "pca"),
    reduction.key = paste0(reduction_prefix, "pca_"),
    verbose = FALSE
  )
  if (is.null(dims)) {
    dims <- 1:ceiling(maxLikGlobalDimEst(data = Embeddings(srtMerge, reduction = paste0(reduction_prefix, "pca")), k = 20, iterations = 100)[["dim.est"]])
  }
  srtMerge@misc[[paste0(reduction_prefix, "Dims")]] <- dims

  pca <- reticulate::r_to_py(Embeddings(srtMerge, reduction = paste0(reduction_prefix, "pca"))[, dims])
  bem <- bbknn$bbknn_pca_matrix(pca, batch_list = srtMerge[[batch, drop = TRUE]])
  bem <- reticulate::py_to_r(bem)
  bbknn_graph <- as(as(bem[[2]], "CsparseMatrix"), "dgCMatrix")
  rownames(bbknn_graph) <- colnames(bbknn_graph) <- colnames(srtMerge)
  srtMerge@graphs$BBKNN <- as.Graph(bbknn_graph)
  srtIntegrated <- srtMerge
  srtMerge <- NULL

  srtIntegrated <- Check_srtIntegrated(srtIntegrated, do_normalization = do_normalization, hvf = hvf, batch = batch, vars_to_regress = vars_to_regress)

  srtIntegrated <- FindClusters(object = srtIntegrated, graph.name = "BBKNN", resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  if (isTRUE(reorder)) {
    srtIntegrated <- SrtReorder(srtIntegrated, features = hvf, reorder_by = "seurat_clusters", slot = "data")
  }
  srtIntegrated[["seurat_clusters"]] <- NULL
  srtIntegrated[[paste0(reduction_prefix, "clusters")]] <- Idents(srtIntegrated)

  cat("Perform nonlinear dimension reduction on the data...\n")
  reduction <- "umap"
  for (n in reduction_components) {
    srtIntegrated <- RunUMAP(
      object = srtIntegrated, graph = "BBKNN",
      umap.method = "umap-learn", metric = "correlation",
      reduction.name = paste0(reduction_prefix, "UMAP", n, "d"),
      reduction.key = paste0(reduction_prefix, "UMAP", n, "d_"),
      n.components = as.integer(n),
      return.model = TRUE, verbose = FALSE
    )
  }

  DefaultAssay(srtIntegrated) <- "RNA"
  if (isTRUE(append)) {
    DefaultAssay(srtMerge_raw) <- "RNA"
    srtMerge_raw@assays$RNA <- srtIntegrated@assays$RNA
    srtMerge_raw@graphs$BBKNN <- srtIntegrated@graphs$BBKNN
    srtMerge_raw@misc[[paste0(reduction_prefix, "Dims")]] <- srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]]
    srtMerge_raw@misc[[paste0(reduction_prefix, "HVF")]] <- hvf
    srtMerge_raw[[paste0(reduction_prefix, "pca")]] <- srtIntegrated[[paste0(reduction_prefix, "pca")]]
    srtMerge_raw[[paste0(reduction_prefix, "clusters")]] <- srtIntegrated[[paste0(reduction_prefix, "clusters")]]
    for (n in reduction_components) {
      for (i in reduction) {
        srtMerge_raw[[paste0(reduction_prefix, toupper(i), n, "d")]] <- srtIntegrated[[paste0(reduction_prefix, toupper(i), n, "d")]]
      }
    }
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

CSS_integrate <- function(srtList = NULL, srtMerge = NULL, append = FALSE,
                          do_normalization = NULL, normalization_method = "logCPM", batch = "orig.ident",
                          vars_to_regress = NULL,
                          HVF_source = "separate", nHVF = 3000, hvf = NULL,
                          maxPC = 100, dims = NULL, resolution = 0.8, reorder = TRUE,
                          reduction = "umap", reduction_prefix = "CSS", reduction_components = 2:3,
                          exogenous_genes = NULL, seed = 11, ...) {
  set.seed(seed)
  library(simspec) # devtools::install_github("zh542370159/simspec@v1.0.3")

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
  if (isTRUE(append)) {
    if (is.null(srtMerge)) {
      stop("srtMerge must be provided when 'append = TRUE'.")
    } else {
      srtMerge_raw <- srtMerge
    }
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(srtList,
      do_normalization = do_normalization,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
  }
  if (is.null(srtList) & !is.null(srtMerge)) {
    checked <- Check_srtMerge(srtMerge,
      do_normalization = do_normalization,
      normalization_method = normalization_method, batch = batch, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = NULL
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
  }
  if (nrow(GetAssayData(srtMerge, slot = "scale.data")) != nrow(GetAssayData(srtMerge, slot = "data"))) {
    cat("Perform ScaleData on the data...\n")
    srtMerge <- Seurat::ScaleData(object = srtMerge, features = rownames(srtMerge), vars.to.regress = vars_to_regress)
  }
  srtMerge <- RunPCA(
    object = srtMerge, npcs = maxPC, features = hvf,
    reduction.name = paste0(reduction_prefix, "pca"),
    reduction.key = paste0(reduction_prefix, "pca_"),
    verbose = FALSE
  )
  if (is.null(dims)) {
    dims <- 1:ceiling(maxLikGlobalDimEst(data = Embeddings(srtMerge, reduction = paste0(reduction_prefix, "pca")), k = 20, iterations = 100)[["dim.est"]])
  }

  srtIntegrated <- cluster_sim_spectrum(
    object = srtMerge,
    use_dr = paste0(reduction_prefix, "pca"),
    dims_use = dims,
    var_genes = hvf,
    label_tag = batch,
    reduction.name = "CSS",
    reduction.key = "CSS_"
  )
  srtMerge <- NULL

  srtIntegrated <- Check_srtIntegrated(srtIntegrated, do_normalization = do_normalization, hvf = hvf, batch = batch, vars_to_regress = vars_to_regress)

  CSSdims <- 1:ncol(Embeddings(srtIntegrated, reduction = "CSS"))
  srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]] <- CSSdims

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "CSS", dims = CSSdims, force.recalc = T, graph.name = paste0(reduction_prefix, "_", c("NN", "SNN")), verbose = FALSE)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000, graph.name = paste0(reduction_prefix, "_SNN"), verbose = FALSE)
  if (isTRUE(reorder)) {
    srtIntegrated <- SrtReorder(srtIntegrated, features = hvf, reorder_by = "seurat_clusters", slot = "data")
  }
  srtIntegrated[["seurat_clusters"]] <- NULL
  srtIntegrated[[paste0(reduction_prefix, "clusters")]] <- Idents(srtIntegrated)

  cat("Perform nonlinear dimension reduction on the data...\n")
  for (n in reduction_components) {
    if ("umap" %in% reduction) {
      srtIntegrated <- RunUMAP(
        object = srtIntegrated, reduction = "CSS",
        reduction.name = paste0(reduction_prefix, "UMAP", n, "d"),
        reduction.key = paste0(reduction_prefix, "UMAP", n, "d_"),
        dims = CSSdims, n.components = n, umap.method = "uwot",
        return.model = TRUE, verbose = FALSE
      )
    }
    if ("tsne" %in% reduction) {
      srtIntegrated <- RunTSNE(
        object = srtIntegrated, reduction = "CSS",
        reduction.name = paste0(reduction_prefix, "TSNE", n, "d"),
        reduction.key = paste0(reduction_prefix, "TSNE", n, "d_"),
        dims = dims, dim.embed = n, tsne.method = "Rtsne",
        perplexity = max(ceiling(ncol(srtIntegrated) * 0.01), 30), max_iter = 2000,
        num_threads = 0, verbose = TRUE
      )
    }
  }

  DefaultAssay(srtIntegrated) <- "RNA"
  if (isTRUE(append)) {
    DefaultAssay(srtMerge_raw) <- "RNA"
    srtMerge_raw@assays$RNA <- srtIntegrated@assays$RNA
    srtMerge_raw@reductions$CSS <- srtIntegrated@reductions$CSS
    srtMerge_raw@graphs[paste0(reduction_prefix, "_", c("NN", "SNN"))] <- srtIntegrated@graphs[paste0(reduction_prefix, "_", c("NN", "SNN"))]
    srtMerge_raw@misc[[paste0(reduction_prefix, "Dims")]] <- srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]]
    srtMerge_raw@misc[[paste0(reduction_prefix, "HVF")]] <- hvf
    srtMerge_raw[[paste0(reduction_prefix, "pca")]] <- srtIntegrated[[paste0(reduction_prefix, "pca")]]
    srtMerge_raw[[paste0(reduction_prefix, "clusters")]] <- srtIntegrated[[paste0(reduction_prefix, "clusters")]]
    for (n in reduction_components) {
      for (i in reduction) {
        srtMerge_raw[[paste0(reduction_prefix, toupper(i), n, "d")]] <- srtIntegrated[[paste0(reduction_prefix, toupper(i), n, "d")]]
      }
    }
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

LIGER_integrate <- function(srtList = NULL, srtMerge = NULL, append = FALSE,
                            do_normalization = NULL, normalization_method = "logCPM", batch = "orig.ident",
                            vars_to_regress = NULL,
                            HVF_source = "separate", nHVF = 3000, hvf = NULL,
                            maxPC = 100, dims = NULL, resolution = 0.8, reorder = TRUE,
                            reduction = "umap", reduction_prefix = "LIGER", reduction_components = 2:3,
                            exogenous_genes = NULL, seed = 11, ...) {
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
  if (isTRUE(append)) {
    if (is.null(srtMerge)) {
      stop("srtMerge must be provided when 'append = TRUE'.")
    } else {
      srtMerge_raw <- srtMerge
    }
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(srtList,
      do_normalization = do_normalization,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
  }
  if (is.null(srtList) & !is.null(srtMerge)) {
    checked <- Check_srtMerge(srtMerge,
      do_normalization = do_normalization,
      normalization_method = normalization_method, batch = batch, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = NULL
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
  }
  srtMerge <- Seurat::ScaleData(object = srtMerge, features = hvf, split.by = batch, do.center = FALSE, vars.to.regress = vars_to_regress)

  srtMerge <- RunOptimizeALS(srtMerge,
    k = 20,
    lambda = 5,
    split.by = batch
  )
  srtIntegrated <- RunQuantileNorm(srtMerge,
    reduction.name = "LIGER",
    reduction.key = "LIGER_",
    split.by = batch
  )
  srtMerge <- NULL

  srtIntegrated <- Check_srtIntegrated(srtIntegrated, do_normalization = do_normalization, hvf = hvf, batch = batch, vars_to_regress = vars_to_regress)

  dims <- 1:ncol(Embeddings(srtIntegrated, reduction = "LIGER"))
  srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]] <- dims

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "LIGER", dims = dims, force.recalc = T, graph.name = paste0(reduction_prefix, "_", c("NN", "SNN")), verbose = FALSE)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000, graph.name = paste0(reduction_prefix, "_SNN"), verbose = FALSE)
  if (isTRUE(reorder)) {
    srtIntegrated <- SrtReorder(srtIntegrated, features = hvf, reorder_by = "seurat_clusters", slot = "data")
  }
  srtIntegrated[["seurat_clusters"]] <- NULL
  srtIntegrated[[paste0(reduction_prefix, "clusters")]] <- Idents(srtIntegrated)

  cat("Perform nonlinear dimension reduction on the data...\n")
  for (n in reduction_components) {
    if ("umap" %in% reduction) {
      srtIntegrated <- RunUMAP(
        object = srtIntegrated, reduction = "LIGER",
        reduction.name = paste0(reduction_prefix, "UMAP", n, "d"),
        reduction.key = paste0(reduction_prefix, "UMAP", n, "d_"),
        dims = dims, n.components = n, umap.method = "uwot",
        return.model = TRUE, verbose = FALSE
      )
    }
    if ("tsne" %in% reduction) {
      srtIntegrated <- RunTSNE(
        object = srtIntegrated, reduction = "LIGER",
        reduction.name = paste0(reduction_prefix, "TSNE", n, "d"),
        reduction.key = paste0(reduction_prefix, "TSNE", n, "d_"),
        dims = dims, dim.embed = n, tsne.method = "Rtsne",
        perplexity = max(ceiling(ncol(srtIntegrated) * 0.01), 30), max_iter = 2000,
        num_threads = 0, verbose = TRUE
      )
    }
  }

  DefaultAssay(srtIntegrated) <- "RNA"
  if (isTRUE(append)) {
    DefaultAssay(srtMerge_raw) <- "RNA"
    srtMerge_raw@assays$RNA <- srtIntegrated@assays$RNA
    srtMerge_raw@reductions$LIGER <- srtIntegrated@reductions$LIGER
    srtMerge_raw@graphs[paste0(reduction_prefix, "_", c("NN", "SNN"))] <- srtIntegrated@graphs[paste0(reduction_prefix, "_", c("NN", "SNN"))]
    srtMerge_raw@misc[[paste0(reduction_prefix, "Dims")]] <- srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]]
    srtMerge_raw@misc[[paste0(reduction_prefix, "HVF")]] <- hvf
    srtMerge_raw[[paste0(reduction_prefix, "clusters")]] <- srtIntegrated[[paste0(reduction_prefix, "clusters")]]
    for (n in reduction_components) {
      for (i in reduction) {
        srtMerge_raw[[paste0(reduction_prefix, toupper(i), n, "d")]] <- srtIntegrated[[paste0(reduction_prefix, toupper(i), n, "d")]]
      }
    }
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

## slow
scMerge_integrate_deprecated <- function(srtList = NULL, srtMerge = NULL, append = FALSE,
                                         do_normalization = NULL, normalization_method = "logCPM", batch = "orig.ident",
                                         vars_to_regress = NULL,
                                         HVF_source = "separate", nHVF = 3000, hvf = NULL,
                                         maxPC = 100, dims = NULL, resolution = 0.8, reorder = TRUE,
                                         reduction = "umap", reduction_prefix = "scMerge", reduction_components = 2:3,
                                         exogenous_genes = NULL, seed = 11, ...) {
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
  if (isTRUE(append)) {
    if (is.null(srtMerge)) {
      stop("srtMerge must be provided when 'append = TRUE'.")
    } else {
      srtMerge_raw <- srtMerge
    }
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(srtList,
      do_normalization = do_normalization,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
  }
  if (is.null(srtList) & !is.null(srtMerge)) {
    checked <- Check_srtMerge(srtMerge,
      do_normalization = do_normalization,
      normalization_method = normalization_method, batch = batch, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = NULL
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
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

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"), dims = dims, force.recalc = T, graph.name = paste0(reduction_prefix, "_", c("NN", "SNN")), verbose = FALSE)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000, graph.name = paste0(reduction_prefix, "_SNN"), verbose = FALSE)
  if (isTRUE(reorder)) {
    srtIntegrated <- SrtReorder(srtIntegrated, features = hvf, reorder_by = "seurat_clusters", slot = "data")
  }
  srtIntegrated[["seurat_clusters"]] <- NULL
  srtIntegrated[[paste0(reduction_prefix, "clusters")]] <- Idents(srtIntegrated)

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

  DefaultAssay(srtIntegrated) <- "RNA"
  if (isTRUE(append)) {
    DefaultAssay(srtMerge_raw) <- "RNA"
    srtMerge_raw@assays$scMerge <- srtIntegrated@assays$scMerge
    srtMerge_raw@graphs[paste0(reduction_prefix, "_", c("NN", "SNN"))] <- srtIntegrated@graphs[paste0(reduction_prefix, "_", c("NN", "SNN"))]
    srtMerge_raw[[paste0(reduction_prefix, "pca")]] <- srtIntegrated[[paste0(reduction_prefix, "pca")]]
    srtMerge_raw@misc[[paste0(reduction_prefix, "Dims")]] <- srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]]
    srtMerge_raw@misc[[paste0(reduction_prefix, "HVF")]] <- hvf
    srtMerge_raw[[paste0(reduction_prefix, "clusters")]] <- srtIntegrated[[paste0(reduction_prefix, "clusters")]]
    for (n in reduction_components) {
      for (i in reduction) {
        srtMerge_raw[[paste0(reduction_prefix, toupper(i), n, "d")]] <- srtIntegrated[[paste0(reduction_prefix, toupper(i), n, "d")]]
      }
    }
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

ZINBWaVE_integrate_deprecated <- function(srtList = NULL, srtMerge = NULL, append = FALSE,
                                          do_normalization = NULL, normalization_method = "logCPM", batch = "orig.ident",
                                          vars_to_regress = NULL,
                                          HVF_source = "separate", nHVF = 3000, hvf = NULL,
                                          maxPC = 100, dims = NULL, resolution = 0.8, reorder = TRUE,
                                          reduction = "umap", reduction_prefix = "ZINBWaVE", reduction_components = 2:3,
                                          exogenous_genes = NULL, seed = 11, ...) {
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
  if (isTRUE(append)) {
    if (is.null(srtMerge)) {
      stop("srtMerge must be provided when 'append = TRUE'.")
    } else {
      srtMerge_raw <- srtMerge
    }
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(srtList,
      do_normalization = do_normalization,
      normalization_method = normalization_method, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
  }
  if (is.null(srtList) & !is.null(srtMerge)) {
    checked <- Check_srtMerge(srtMerge,
      do_normalization = do_normalization,
      normalization_method = normalization_method, batch = batch, vars_to_regress = vars_to_regress,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = NULL
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
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

  srtIntegrated <- NormalizeData(object = srtIntegrated, normalization.method = "LogNormalize")
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

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"), dims = dims, force.recalc = T, graph.name = paste0(reduction_prefix, "_", c("NN", "SNN")), verbose = FALSE)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000, graph.name = paste0(reduction_prefix, "_SNN"), verbose = FALSE)
  if (isTRUE(reorder)) {
    srtIntegrated <- SrtReorder(srtIntegrated, features = hvf, reorder_by = "seurat_clusters", slot = "data")
  }
  srtIntegrated[["seurat_clusters"]] <- NULL
  srtIntegrated[[paste0(reduction_prefix, "clusters")]] <- Idents(srtIntegrated)

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

  DefaultAssay(srtIntegrated) <- "RNA"
  if (isTRUE(append)) {
    DefaultAssay(srtMerge_raw) <- "RNA"
    srtMerge_raw@assays$ZINBWaVE <- srtIntegrated@assays$ZINBWaVE
    srtMerge_raw@graphs[paste0(reduction_prefix, "_", c("NN", "SNN"))] <- srtIntegrated@graphs[paste0(reduction_prefix, "_", c("NN", "SNN"))]
    srtMerge_raw[[paste0(reduction_prefix, "pca")]] <- srtIntegrated[[paste0(reduction_prefix, "pca")]]
    srtMerge_raw@misc[[paste0(reduction_prefix, "Dims")]] <- srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]]
    srtMerge_raw@misc[[paste0(reduction_prefix, "HVF")]] <- hvf
    srtMerge_raw[[paste0(reduction_prefix, "clusters")]] <- srtIntegrated[[paste0(reduction_prefix, "clusters")]]
    for (n in reduction_components) {
      for (i in reduction) {
        srtMerge_raw[[paste0(reduction_prefix, toupper(i), n, "d")]] <- srtIntegrated[[paste0(reduction_prefix, toupper(i), n, "d")]]
      }
    }
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}


Standard_SCP <- function(srt, do_normalization = NULL, normalization_method = "logCPM", nHVF = 3000, hvf = NULL,
                         vars_to_regress = NULL,
                         maxPC = 100, dims = NULL, resolution = 0.8, reorder = TRUE,
                         reduction = "umap", reduction_prefix = "Standard", reduction_components = 2:3,
                         exogenous_genes = NULL, seed = 11) {
  set.seed(seed)
  library(glmGamPoi)
  library(intrinsicDimension)
  library(dplyr)

  if (class(srt) != "Seurat") {
    stop("'srt' is not a Seurat object.",
      call. = FALSE
    )
  }
  if (!normalization_method %in% c("logCPM", "SCT")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT'",
      call. = FALSE
    )
  }

  DefaultAssay(srt) <- "RNA"

  if (isTRUE(do_normalization) | (is.null(do_normalization) & identical(
    x = GetAssayData(srt, slot = "counts"),
    y = GetAssayData(srt, slot = "data")
  ))) {
    cat("Perform NormalizeData(logCPM) on the data...\n")
    srt <- NormalizeData(object = srt, normalization.method = "LogNormalize", verbose = FALSE)
  }
  if (is.null(hvf)) {
    if (!"vst.variance.standardized" %in% colnames(srt@assays$RNA@meta.features)) {
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
  srt@misc[[paste0(reduction_prefix, "HVF")]] <- hvf
  if (nrow(GetAssayData(srt, slot = "scale.data")) != nrow(GetAssayData(srt, slot = "data"))) {
    cat("Perform ScaleData on the data...\n")
    srt <- Seurat::ScaleData(object = srt, features = rownames(srt), vars.to.regress = vars_to_regress)
  }

  if (normalization_method %in% c("SCT")) {
    if (!"SCT" %in% Seurat::Assays(srt)) {
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

  cat("Perform PCA on the data...\n")
  srt <- RunPCA(
    object = srt, npcs = maxPC, features = hvf,
    reduction.name = paste0(reduction_prefix, "pca"),
    reduction.key = paste0(reduction_prefix, "pca_"),
    verbose = FALSE
  )

  if (is.null(dims)) {
    dim_est <- maxLikGlobalDimEst(data = Embeddings(srt, reduction = paste0(reduction_prefix, "pca")), k = 20, iterations = 100)[["dim.est"]]
    if (!is.na(dim_est)) {
      dims <- 1:ceiling(dim_est)
    } else {
      dims <- 1:20
    }
  }
  srt@misc[[paste0(reduction_prefix, "Dims")]] <- dims

  srt <- FindNeighbors(object = srt, reduction = paste0(reduction_prefix, "pca"), dims = dims, force.recalc = T, graph.name = paste0(reduction_prefix, "_", c("NN", "SNN")), verbose = FALSE)
  srt <- FindClusters(object = srt, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000, graph.name = paste0(reduction_prefix, "_SNN"), verbose = FALSE)
  if (isTRUE(reorder)) {
    srt <- SrtReorder(srt, features = hvf, reorder_by = "seurat_clusters", slot = "data")
  }
  srt[["seurat_clusters"]] <- NULL
  srt[[paste0(reduction_prefix, "clusters")]] <- Idents(srt)

  cat("Perform nonlinear dimension reduction on the data...\n")
  for (n in reduction_components) {
    if ("umap" %in% reduction) {
      srt <- RunUMAP(
        object = srt, reduction = paste0(reduction_prefix, "pca"),
        reduction.name = paste0(reduction_prefix, "UMAP", n, "d"),
        reduction.key = paste0(reduction_prefix, "UMAP", n, "d_"),
        dims = dims, n.components = n, umap.method = "uwot",
        return.model = TRUE, verbose = FALSE
      )
    }
    if ("tsne" %in% reduction) {
      srt <- RunTSNE(
        object = srt, reduction = paste0(reduction_prefix, "pca"),
        reduction.name = paste0(reduction_prefix, "TSNE", n, "d"),
        reduction.key = paste0(reduction_prefix, "TSNE", n, "d_"),
        dims = dims, dim.embed = n, tsne.method = "Rtsne",
        perplexity = max(ceiling(ncol(srt) * 0.01), 30), max_iter = 2000,
        num_threads = 0, verbose = TRUE
      )
    }
  }

  DefaultAssay(srt) <- "RNA"
  return(srt)
}

Integration_SCP <- function(srtList = NULL, srtMerge = NULL, append = FALSE,
                            do_normalization = NULL, integration_method = "Uncorrected", batch = "orig.ident",
                            normalization_method = "logCPM", vars_to_regress = NULL,
                            HVF_source = "separate", nHVF = 3000, hvf = NULL,
                            maxPC = 100, dims = NULL, resolution = 0.8, reorder = TRUE,
                            reduction = "umap", reduction_components = 2:3,
                            exogenous_genes = NULL, seed = 11, ...) {
  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (length(integration_method) == 1 & integration_method %in% c("Uncorrected", "Seurat", "fastMNN", "Harmony", "Scanorama", "BBKNN", "CSS", "LIGER")) {
    args1 <- c(mget(names(formals())), reduction_prefix = integration_method)
    args2 <- c(as.list(match.call()), reduction_prefix = integration_method)
    for (n in names(args2)) {
      args1[[n]] <- args2[[n]]
    }
    args1 <- args1[!names(args1) %in% c("integration_method", "...")]
    tryCatch(expr = {
      srtIntegrated <- base::do.call(
        what = paste0(integration_method, "_integrate"),
        args = args1
      )
    }, error = function(e) {
      message(e)
      message(cat("\n", paste0("[", Sys.time(), "]", " Stop the integration...\n")))
      stop(call. = FALSE)
    })
    if (exists("srtIntegrated")) {
      return(srtIntegrated)
    }
    if (!is.null(srtMerge)) {
      return(srtMerge)
    }
    return(NULL)
  } else {
    stop(paste(integration_method, "is not a suppoted integration method!"),
      call. = FALSE
    )
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
    srt@tools[["DEtest_custom"]][["CellMarkers_Wilcoxon"]] <- markers
    srt@tools[["DEtest_custom"]][["cell_group1"]] <- cell_group1
    srt@tools[["DEtest_custom"]][["cell_group2"]] <- cell_group2
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
      } else {
        return(NULL)
      }
      return(markers)
    }, BPPARAM = BPPARAM)
    AllMarkers_Wilcoxon <- dplyr::bind_rows(AllMarkers_Wilcoxon)
    rownames(AllMarkers_Wilcoxon) <- NULL
    AllMarkers_Wilcoxon[, "group1"] <- factor(AllMarkers_Wilcoxon[, "group1"], levels = levels(cell_group))
    AllMarkers_Wilcoxon[, "p_val_adj"] <- p.adjust(AllMarkers_Wilcoxon[, "p_val"], method = "BH")
    # AllMarkers_Wilcoxon <- AllMarkers_Wilcoxon %>%
    #   dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
    #   dplyr::group_by(gene) %>%
    #   dplyr::mutate(groupMaxExp = dplyr::first(group1)) %>%
    #   dplyr::ungroup() %>%
    #   dplyr::arrange(groupMaxExp) %>%
    #   as.data.frame()
    AllMarkers_Wilcoxon[, "DE_group_number"] <- as.integer(table(AllMarkers_Wilcoxon[["gene"]])[AllMarkers_Wilcoxon[, "gene"]])
    srt@tools[[paste0("DEtest_", group_by)]][["AllMarkers_Wilcoxon"]] <- AllMarkers_Wilcoxon
  }

  if (isTRUE(FindPairMarkers)) {
    if (nlevels(cell_group) > 30 & (!isTRUE(force))) {
      warning("Too many groups for FindPairMarkers function. If you want to force to run, set force=TRUE.")
      break
    }
    cat("Perform Wilcoxon FindPairMarkers...\n")
    # pair <- t(combn(nlevels(x = cell_group),2))
    # pair <- apply(pair, c(1,2), function(i){levels(cell_group)[i]})
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
      markers[, "gene"] <- rownames(markers)
      markers[, "group1"] <- as.character(pair[i, 1])
      markers[, "group2"] <- as.character(pair[i, 2])
      return(markers)
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
    srt@tools[[paste0("DEtest_", group_by)]][["PairMarkers_Wilcoxon"]] <- PairMarkers_Wilcoxon
    srt@tools[[paste0("DEtest_", group_by)]][["PairMarkers_matrix"]] <- PairMarkers_matrix
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
      if (assay %in% c("spliced", "unspliced", "ambiguous")) {
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


RunEnrichment <- function(gene, group, enrichment = "GO") {

}
