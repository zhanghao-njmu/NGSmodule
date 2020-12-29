Check_srtList <- function(srtList, normalization_method,
                          HVF_source = "separate", nHVF = 3000, hvf = NULL,
                          exogenous_genes = NULL, ...) {
  require(Seurat)
  require(sctransform)
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
    rownames(GetAssayData(x, slot = "counts", assay = "RNA"))
  })
  if (length(unique(genelist)) != 1) {
    stop("'srtList' must have identical feature names!",
      call. = FALSE
    )
  }

  for (i in 1:length(srtList)) {
    DefaultAssay(srtList[[i]]) <- "RNA"
    if (identical(
      x = GetAssayData(srtList[[i]], slot = "counts"),
      y = GetAssayData(srtList[[i]], slot = "data")
    )) {
      srtList[[i]] <- NormalizeData(object = srtList[[i]], normalization.method = "LogNormalize", verbose = FALSE)
    }
    if (!"variance.standardized" %in% colnames(srtList[[i]]@assays$RNA@meta.features)) {
      srtList[[i]] <- FindVariableFeatures(srtList[[i]], verbose = FALSE)
    }
    m <- GetAssayData(srtList[[i]], slot = "counts")
    gene_keep <- rownames(m)[Matrix::rowSums(m > 5) > 5]
    VariableFeatures(srtList[[i]]) <- HVFInfo(srtList[[i]]) %>%
      filter(variance.standardized > 1 &
        (!rownames(.) %in% exogenous_genes) &
        rownames(.) %in% gene_keep) %>%
      dplyr::arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
    # if (nrow(GetAssayData(srtList[[i]], slot = "scale.data")) != nrow(GetAssayData(srtList[[i]], slot = "data"))) {
    #   srtList[[i]] <- ScaleData(object = srtList[[i]], features = rownames(srtList[[i]]), verbose = FALSE)
    # }
    DefaultAssay(srtList[[i]]) <- "RNA"

    if (normalization_method %in% c("SCT")) {
      if (!"SCT" %in% Seurat::Assays(srtList[[i]])) {
        srtList[[i]] <- SCTransform(
          object = srtList[[i]],
          variable.features.n = nHVF,
          return.only.var.genes = FALSE,
          min_cells = 5,
          assay = "RNA",
          verbose = FALSE
        )
      } else {
        DefaultAssay(srtList[[i]]) <- "SCT"
      }
      if (!"variance.standardized" %in% colnames(srtList[[i]]@assays$SCT@meta.features)) {
        srtList[[i]] <- FindVariableFeatures(srtList[[i]], verbose = FALSE)
      }
      m <- GetAssayData(srtList[[i]], slot = "counts")
      gene_keep <- rownames(m)[Matrix::rowSums(m > 5) > 5]
      VariableFeatures(srtList[[i]]) <- HVFInfo(srtList[[i]]) %>%
        filter(variance.standardized > 1 &
          (!rownames(.) %in% exogenous_genes) &
          rownames(.) %in% gene_keep) %>%
        dplyr::arrange(desc(variance.standardized)) %>%
        rownames(.) %>%
        head(n = nHVF)
    }
  }

  if (is.null(hvf)) {
    if (HVF_source == "global") {
      gene_common <- lapply(srtList, function(x) {
        m <- GetAssayData(x, slot = "counts")
        gene_keep <- rownames(m)[Matrix::rowSums(m > 5) > 5 * length(srtList)]
        return(gene_keep)
      }) %>% Reduce(intersect, .)
      srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
      hvf <- FindVariableFeatures(srtMerge, verbose = FALSE) %>%
        HVFInfo(.) %>%
        filter(variance.standardized > 1 &
          (!rownames(.) %in% exogenous_genes) &
          rownames(.) %in% gene_common) %>%
        dplyr::arrange(desc(variance.standardized)) %>%
        rownames(.) %>%
        head(n = nHVF)
    }
    if (HVF_source == "separate") {
      hvf <- SelectIntegrationFeatures(object.list = srtList, nfeatures = nHVF, verbose = FALSE)
    }
  } else {
    hvf <- hvf[hvf %in% rownames(GetAssayData(srtList[[1]], slot = "counts"))]
  }

  if (normalization_method %in% c("SCT")) {
    srtList <- PrepSCTIntegration(object.list = srtList, anchor.features = hvf, verbose = FALSE)
  }

  return(list(
    srtList = srtList,
    hvf = hvf
  ))
}

Check_srtMerge <- function(srtMerge, normalization_method, batch,
                           HVF_source = "separate", nHVF = 3000, hvf = NULL,
                           exogenous_genes = NULL, ...) {
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

  srtList <- SplitObject(object = srtMerge, split.by = batch)

  checked <- Check_srtList(srtList,
    normalization_method = normalization_method,
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

Check_srtIntegrated <- function(srtIntegrated, hvf, batch, ...) {
  raw_DefaultAssay <- DefaultAssay(object = srtIntegrated)

  DefaultAssay(object = srtIntegrated) <- "RNA"
  if (identical(
    x = GetAssayData(srtIntegrated, slot = "counts"),
    y = GetAssayData(srtIntegrated, slot = "data")
  )) {
    srtIntegrated <- NormalizeData(object = srtIntegrated)
  }
  if (length(VariableFeatures(srtIntegrated)) == 0) {
    hvf <- hvf[hvf %in% rownames(GetAssayData(srtIntegrated, slot = "counts"))]
    VariableFeatures(srtIntegrated) <- hvf
  }
  if (nrow(GetAssayData(srtIntegrated, slot = "scale.data")) != nrow(GetAssayData(srtIntegrated, slot = "data"))) {
    srtIntegrated <- ScaleData(object = srtIntegrated, features = rownames(srtIntegrated))
  }

  DefaultAssay(object = srtIntegrated) <- raw_DefaultAssay
  srtIntegrated@project.name <- paste0(unique(srtIntegrated[[batch, drop = TRUE]]), collapse = ",")
  srtIntegrated[[batch]] <- factor(srtIntegrated[[batch, drop = TRUE]],
    levels = unique(srtIntegrated[[batch, drop = TRUE]])
  )
  return(srtIntegrated)
}

CC_GenePrefetch <- function(species) {
  if (species == "Homo_sapiens") {
    cc_S_genes <- Seurat::cc.genes.updated.2019$s.genes
    cc_G2M_genes <- Seurat::cc.genes.updated.2019$g2m.genes
  } else {
    species_split <- unlist(strsplit(species, split = "_"))
    species_homolog <- paste0(tolower(substring(species_split[1], 1, 1)), species_split[2], "_homolog_associated_gene_name")

    archives <- listEnsemblArchives()
    # web <- read_html(httr::RETRY("GET", "http://www.ensembl.org/info/website/archives/index.html?redirect=no", times = 1000, timeout(1000)))
    # urls <- web %>% html_nodes("ul") %>% html_nodes("strong") %>% html_nodes("a") %>% html_attr("href")
    # version <- web %>% html_nodes("ul") %>% html_nodes("strong") %>% html_nodes("a") %>% html_text(trim = TRUE) %>%
    #   gsub(pattern = "(Ensembl )|(:.*)",replacement = "",x = .,perl = T)
    # archives <- data.frame(version=version,url=urls,stringsAsFactors = F)
    url <- archives[which(archives$version == Ensembl_version), "url"]

    mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = url)
    homolog <- listAttributes(mart)$name

    if (species_homolog %in% homolog) {
      cc_S_genes <- getBM(
        mart = mart,
        attributes = c(species_homolog),
        filters = c("external_gene_name"),
        values = list(Seurat::cc.genes.updated.2019$s.genes)
      )[[1]]
      cc_G2M_genes <- getBM(
        mart = mart,
        attributes = c(species_homolog),
        filters = c("external_gene_name"),
        values = list(Seurat::cc.genes.updated.2019$g2m.genes)
      )[[1]]

      if (length(cc_S_genes) < 3 | length(cc_G2M_genes) < 3) {
        warning(paste0("number of cell-cycle homolog genes is too small. CellCycleScoring will not performed."))
      }
    } else {
      warning(paste0("Can not find the homolog attributes for the species: ", species, " (", species_homolog, ")"))
      cc_S_genes <- cc_S_genes <- NA
    }
  }
  return(list(
    cc_S_genes = cc_S_genes,
    cc_G2M_genes = cc_G2M_genes
  ))
}

CC_module <- function(srt, cc_S_genes, cc_G2M_genes, ...) {
  if (length(cc_S_genes) >= 3 & length(cc_G2M_genes) >= 3) {
    srt <- CellCycleScoring(
      object = srt,
      s.features = cc_S_genes,
      g2m.features = cc_G2M_genes,
      set.ident = FALSE
    )
    srt[["CC.Difference"]] <- srt[["S.Score"]] - srt[["G2M.Score"]]
    srt[["Phase"]] <- factor(srt[["Phase", drop = TRUE]], levels = c("G1", "S", "G2M"))
  }
  return(srt)
}

Uncorrected_integrate <- function(srtList = NULL, srtMerge = NULL, append = TRUE,
                                  normalization_method = "logCPM", batch = "orig.ident",
                                  HVF_source = "separate", nHVF = 3000, hvf = NULL,
                                  maxPC = 100, resolution = 0.8,
                                  reduction = c("tsne", "umap"), reduction_prefix = "Uncorrected_",
                                  cc_S_genes = NULL, cc_G2M_genes = NULL,
                                  exogenous_genes = NULL, ...) {
  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
    checked <- Check_srtMerge(srtMerge,
      normalization_method = normalization_method, batch = batch,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = NULL
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
  } else {
    checked <- Check_srtList(srtList,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
  }

  srtIntegrated <- Standard_SCP(
    srt = srtMerge, normalization_method = normalization_method, nHVF = nHVF, hvf = hvf,
    maxPC = maxPC, resolution = resolution, reduction = reduction, reduction_prefix = reduction_prefix,
    cc_S_genes = cc_S_genes, cc_G2M_genes = cc_G2M_genes,
    exogenous_genes = exogenous_genes
  )
  srtIntegrated@project.name <- paste0(unique(srtIntegrated[[batch, drop = TRUE]]), collapse = ",")

  if (!is.null(srtMerge) & isTRUE(append)) {
    srtMerge_raw[[paste0(reduction_prefix, "pca")]] <- srtIntegrated[[paste0(reduction_prefix, "pca")]]
    srtMerge_raw@misc[[paste0(reduction_prefix, "Dims")]] <- srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]]
    srtMerge_raw[[paste0(reduction_prefix, "clusters")]] <- srtIntegrated[[paste0(reduction_prefix, "clusters")]]
    for (i in reduction) {
      srtMerge_raw[[paste0(reduction_prefix, i)]] <- srtIntegrated[[paste0(reduction_prefix, i)]]
    }
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

Seurat_integrate <- function(srtList = NULL, srtMerge = NULL, append = TRUE,
                             normalization_method = "logCPM", batch = "orig.ident",
                             HVF_source = "separate", nHVF = 3000, hvf = NULL,
                             maxPC = 100, resolution = 0.8,
                             reduction = c("tsne", "umap"), reduction_prefix = "Seurat_",
                             cc_S_genes = NULL, cc_G2M_genes = NULL,
                             exogenous_genes = NULL, ...) {
  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
    srtList <- SplitObject(object = srtMerge, split.by = batch)
    checked <- Check_srtList(srtList,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(srtList,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
  }

  srt_anchors <- FindIntegrationAnchors(
    object.list = srtList,
    normalization.method = switch(normalization_method,
      "logCPM" = "LogNormalize", "SCT" = "SCT"
    ),
    anchor.features = hvf,
    dims = 1:30
  )
  srtIntegrated <- IntegrateData(
    anchorset = srt_anchors,
    normalization.method = switch(normalization_method,
      "logCPM" = "LogNormalize", "SCT" = "SCT"
    ),
    dims = 1:30,
    features.to.integrate = c(
      hvf,
      Reduce(union, lapply(srtList, VariableFeatures)),
      Reduce(intersect, c(
        lapply(srtList, rownames),
        list(c(cc_S_genes, cc_G2M_genes))
      ))
    )
  )
  RenameAssays(object = srtIntegrated, integrated = "Seurat")
  DefaultAssay(srtIntegrated) <- "Seurat"

  srtIntegrated <- Check_srtIntegrated(srtIntegrated, hvf = hvf, batch = batch)

  srtIntegrated <- ScaleData(srtIntegrated, features = hvf)
  srtIntegrated <- RunPCA(object = srtIntegrated, npcs = maxPC, features = hvf, reduction.name = paste0(reduction_prefix, "pca"))
  dims <- 1:ceiling(maxLikGlobalDimEst(data = Embeddings(srtIntegrated, reduction = paste0(reduction_prefix, "pca")), k = 20, iterations = 100)[["dim.est"]])
  srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]] <- dims

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"), dims = dims, force.recalc = T)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srtIntegrated <- BuildClusterTree(srtIntegrated, features = hvf, slot = "data", reorder = T, reorder.numeric = T)
  srtIntegrated[[paste0(reduction_prefix, "clusters")]] <- Idents(srtIntegrated)

  if ("umap" %in% reduction) {
    srtIntegrated <- RunUMAP(object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"), dims = dims, n.components = 2, umap.method = "uwot-learn", reduction.name = paste0(reduction_prefix, "umap"))
  }
  if ("tsne" %in% reduction) {
    srtIntegrated <- RunTSNE(
      object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"), dims = dims, dim.embed = 2, tsne.method = "Rtsne", reduction.name = paste0(reduction_prefix, "tsne"),
      perplexity = max(ceiling(ncol(srtIntegrated) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  srtIntegrated <- CC_module(srtIntegrated, cc_S_genes, cc_G2M_genes)

  DefaultAssay(srtIntegrated) <- "RNA"
  if (!is.null(srtMerge) & isTRUE(append)) {
    srtMerge_raw@assays$Seurat <- srtIntegrated@assays$Seurat
    srtMerge_raw[[paste0(reduction_prefix, "pca")]] <- srtIntegrated[[paste0(reduction_prefix, "pca")]]
    srtMerge_raw@misc[[paste0(reduction_prefix, "Dims")]] <- srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]]
    srtMerge_raw[[paste0(reduction_prefix, "clusters")]] <- srtIntegrated[[paste0(reduction_prefix, "clusters")]]
    for (i in reduction) {
      srtMerge_raw[[paste0(reduction_prefix, i)]] <- srtIntegrated[[paste0(reduction_prefix, i)]]
    }
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

fastMNN_integrate <- function(srtList = NULL, srtMerge = NULL, append = TRUE,
                              normalization_method = "logCPM", batch = "orig.ident",
                              HVF_source = "separate", nHVF = 3000, hvf = NULL,
                              maxPC = 100, resolution = 0.8,
                              reduction = c("tsne", "umap"), reduction_prefix = "fastMNN_",
                              cc_S_genes = NULL, cc_G2M_genes = NULL,
                              exogenous_genes = NULL, ...) {
  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
    srtList <- SplitObject(object = srtMerge, split.by = batch)
    checked <- Check_srtList(srtList,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(srtList,
      normalization_method = normalization_method,
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

  srtIntegrated <- Check_srtIntegrated(srtIntegrated, hvf = hvf, batch = batch)

  dims <- 1:ceiling(maxLikGlobalDimEst(data = Embeddings(srtIntegrated, reduction = "fastMNN"), k = 20, iterations = 100)[["dim.est"]])
  srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]] <- dims

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "fastMNN", dims = dims, force.recalc = T)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srtIntegrated <- BuildClusterTree(srtIntegrated, features = hvf, slot = "data", reorder = T, reorder.numeric = T)
  srtIntegrated[[paste0(reduction_prefix, "clusters")]] <- Idents(srtIntegrated)

  if ("umap" %in% reduction) {
    srtIntegrated <- RunUMAP(object = srtIntegrated, reduction = "fastMNN", dims = dims, n.components = 2, umap.method = "uwot-learn", reduction.name = paste0(reduction_prefix, "umap"))
  }
  if ("tsne" %in% reduction) {
    srtIntegrated <- RunTSNE(
      object = srtIntegrated, reduction = "fastMNN", dims = dims, dim.embed = 2, tsne.method = "Rtsne", reduction.name = paste0(reduction_prefix, "tsne"),
      perplexity = max(ceiling(ncol(srtIntegrated) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  srtIntegrated <- CC_module(srtIntegrated, cc_S_genes, cc_G2M_genes)

  DefaultAssay(srtIntegrated) <- "RNA"
  if (!is.null(srtMerge) & isTRUE(append)) {
    srtMerge_raw@reductions$fastMNN <- srtIntegrated@reductions$fastMNN
    srtMerge_raw@misc[[paste0(reduction_prefix, "Dims")]] <- srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]]
    srtMerge_raw[[paste0(reduction_prefix, "clusters")]] <- srtIntegrated[[paste0(reduction_prefix, "clusters")]]
    for (i in reduction) {
      srtMerge_raw[[paste0(reduction_prefix, i)]] <- srtIntegrated[[paste0(reduction_prefix, i)]]
    }
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

Harmony_integrate <- function(srtList = NULL, srtMerge = NULL, append = TRUE,
                              normalization_method = "logCPM", batch = "orig.ident",
                              HVF_source = "separate", nHVF = 3000, hvf = NULL,
                              maxPC = 100, resolution = 0.8,
                              reduction = c("tsne", "umap"), reduction_prefix = "Harmony_",
                              cc_S_genes = NULL, cc_G2M_genes = NULL,
                              exogenous_genes = NULL, ...) {
  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
    checked <- Check_srtMerge(srtMerge,
      normalization_method = normalization_method, batch = batch,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = NULL
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(srtList,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
  }

  srtMerge <- ScaleData(object = srtMerge, features = rownames(srtMerge))
  srtMerge <- RunPCA(object = srtMerge, npcs = maxPC, features = hvf)
  dims <- 1:ceiling(maxLikGlobalDimEst(data = Embeddings(srtMerge, reduction = "pca"), k = 20, iterations = 100)[["dim.est"]])

  srtIntegrated <- RunHarmony(
    object = srtMerge,
    group.by.vars = batch,
    reduction = "pca",
    dims.use = dims,
    reduction.save = "Harmony",
    assay.use = DefaultAssay(srtMerge)
  )
  srtMerge <- NULL

  srtIntegrated <- Check_srtIntegrated(srtIntegrated, hvf = hvf, batch = batch)

  dims <- 1:ceiling(maxLikGlobalDimEst(data = Embeddings(srtIntegrated, reduction = "Harmony"), k = 20, iterations = 100)[["dim.est"]])
  srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]] <- dims

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "Harmony", dims = dims, force.recalc = T)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srtIntegrated <- BuildClusterTree(srtIntegrated, features = hvf, slot = "data", reorder = T, reorder.numeric = T)
  srtIntegrated[[paste0(reduction_prefix, "clusters")]] <- Idents(srtIntegrated)

  if ("umap" %in% reduction) {
    srtIntegrated <- RunUMAP(object = srtIntegrated, reduction = "Harmony", dims = dims, n.components = 2, umap.method = "uwot-learn", reduction.name = paste0(reduction_prefix, "umap"))
  }
  if ("tsne" %in% reduction) {
    srtIntegrated <- RunTSNE(
      object = srtIntegrated, reduction = "Harmony", dims = dims, dim.embed = 2, tsne.method = "Rtsne", reduction.name = paste0(reduction_prefix, "tsne"),
      perplexity = max(ceiling(ncol(srtIntegrated) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  srtIntegrated <- CC_module(srtIntegrated, cc_S_genes, cc_G2M_genes)

  DefaultAssay(srtIntegrated) <- "RNA"
  if (!is.null(srtMerge) & isTRUE(append)) {
    srtMerge_raw@reductions$Harmony <- srtIntegrated@reductions$Harmony
    srtMerge_raw@misc[[paste0(reduction_prefix, "Dims")]] <- srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]]
    srtMerge_raw[[paste0(reduction_prefix, "clusters")]] <- srtIntegrated[[paste0(reduction_prefix, "clusters")]]
    for (i in reduction) {
      srtMerge_raw[[paste0(reduction_prefix, i)]] <- srtIntegrated[[paste0(reduction_prefix, i)]]
    }
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

Scanorama_integrate <- function(srtList = NULL, srtMerge = NULL, append = TRUE,
                                normalization_method = "logCPM", batch = "orig.ident",
                                HVF_source = "separate", nHVF = 3000, hvf = NULL,
                                maxPC = 100, resolution = 0.8,
                                reduction = c("tsne", "umap"), reduction_prefix = "Scanorama_",
                                cc_S_genes = NULL, cc_G2M_genes = NULL,
                                exogenous_genes = NULL, ...) {
  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
    srtList <- SplitObject(object = srtMerge, split.by = batch)
    checked <- Check_srtList(srtList,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(srtList,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
  }

  require(reticulate)
  scanorama <- reticulate::import("scanorama")

  assaylist <- list()
  genelist <- list()
  for (i in 1:length(srtList)) {
    assaylist[[i]] <- t(as.matrix(GetAssayData(object = srtList[[i]], slot = "data")))
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
    rbind.fill.matrix() %>%
    t()
  rownames(cor_value) <- integrated.corrected.data[[3]]
  colnames(cor_value) <- unlist(sapply(assaylist, rownames))

  dim_reduction <- integrated.corrected.data[[1]] %>% rbind.fill.matrix()
  rownames(dim_reduction) <- unlist(sapply(assaylist, rownames))
  colnames(dim_reduction) <- paste0("Scanorama_", 1:100)
  stdevs <- apply(dim_reduction, MARGIN = 2, FUN = sd)

  srtIntegrated <- Reduce(function(x, y) merge(x, y), srtList)
  srtIntegrated@assay$Scanorama <- CreateAssayObject(data = cor_value)
  srtIntegrated@reductions$Scanorama <- CreateDimReducObject(embeddings = dim_reduction, assay = "Scanorama", stdev = stdevs, key = "Scanorama_")
  DefaultAssay(srtIntegrated) <- "Scanorama"

  srtIntegrated <- Check_srtIntegrated(srtIntegrated, hvf = hvf, batch = batch)

  srtIntegrated <- ScaleData(srtIntegrated, features = hvf)
  srtIntegrated <- RunPCA(object = srtIntegrated, npcs = maxPC, features = hvf, reduction.name = paste0(reduction_prefix, "pca"))
  dims <- 1:ceiling(maxLikGlobalDimEst(data = Embeddings(srtIntegrated, reduction = paste0(reduction_prefix, "pca")), k = 20, iterations = 100)[["dim.est"]])
  srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]] <- dims

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"), dims = dims, force.recalc = T)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srtIntegrated <- BuildClusterTree(srtIntegrated, features = hvf, slot = "data", reorder = T, reorder.numeric = T)
  srtIntegrated[[paste0(reduction_prefix, "clusters")]] <- Idents(srtIntegrated)

  if ("umap" %in% reduction) {
    srtIntegrated <- RunUMAP(object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"), dims = dims, n.components = 2, umap.method = "uwot-learn", reduction.name = paste0(reduction_prefix, "umap"))
  }
  if ("tsne" %in% reduction) {
    srtIntegrated <- RunTSNE(
      object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"), dims = dims, dim.embed = 2, tsne.method = "Rtsne", reduction.name = paste0(reduction_prefix, "tsne"),
      perplexity = max(ceiling(ncol(srtIntegrated) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  srtIntegrated <- CC_module(srtIntegrated, cc_S_genes, cc_G2M_genes)

  DefaultAssay(srtIntegrated) <- "RNA"
  if (!is.null(srtMerge) & isTRUE(append)) {
    srtMerge_raw@assays$Scanorama <- srtIntegrated@assays$Scanorama
    srtMerge_raw@reductions$Scanorama <- srtIntegrated@reductions$Scanorama
    srtMerge_raw[[paste0(reduction_prefix, "pca")]] <- srtIntegrated[[paste0(reduction_prefix, "pca")]]
    srtMerge_raw@misc[[paste0(reduction_prefix, "Dims")]] <- srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]]
    srtMerge_raw[[paste0(reduction_prefix, "clusters")]] <- srtIntegrated[[paste0(reduction_prefix, "clusters")]]
    for (i in reduction) {
      srtMerge_raw[[paste0(reduction_prefix, i)]] <- srtIntegrated[[paste0(reduction_prefix, i)]]
    }
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

BBKNN_integrate <- function(srtList = NULL, srtMerge = NULL, append = TRUE,
                            normalization_method = "logCPM", batch = "orig.ident",
                            HVF_source = "separate", nHVF = 3000, hvf = NULL,
                            maxPC = 100, resolution = 0.8,
                            reduction_prefix = "BBKNN_",
                            cc_S_genes = NULL, cc_G2M_genes = NULL,
                            exogenous_genes = NULL, ...) {
  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
    checked <- Check_srtMerge(srtMerge,
      normalization_method = normalization_method, batch = batch,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = NULL
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(srtList,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
  }

  srtMerge <- ScaleData(object = srtMerge, features = rownames(srtMerge))
  srtMerge <- RunPCA(object = srtMerge, npcs = maxPC, features = hvf)
  dims <- 1:ceiling(maxLikGlobalDimEst(data = Embeddings(srtMerge, reduction = "pca"), k = 20, iterations = 100)[["dim.est"]])
  srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]] <- dims

  bbknn <- reticulate::import("bbknn", convert = FALSE)
  pca <- reticulate::r_to_py(Embeddings(srtMerge, reduction = "pca")[, dims])

  bem <- bbknn$bbknn_pca_matrix(pca, batch_list = srtMerge[[batch, drop = TRUE]])
  bem <- reticulate::py_to_r(bem)
  bbknn_graph <- as.matrix(bem[[2]])
  rownames(bbknn_graph) <- colnames(bbknn_graph) <- colnames(srtMerge)
  srtMerge@graphs$BBKNN <- as.Graph(bbknn_graph)
  srtIntegrated <- srtMerge
  srtMerge <- NULL

  srtIntegrated <- Check_srtIntegrated(srtIntegrated, hvf = hvf, batch = batch)

  srtMerge <- FindClusters(object = srtMerge, graph.name = "BBKNN", resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srtIntegrated <- BuildClusterTree(srtIntegrated, features = hvf, slot = "data", reorder = T, reorder.numeric = T)
  srtIntegrated[[paste0(reduction_prefix, "clusters")]] <- Idents(srtIntegrated)

  srtIntegrated <- RunUMAP(object = srtIntegrated, graph = "BBKNN", umap.method = "umap-learn", reduction.name = paste0(reduction_prefix, "umap"))

  srtIntegrated <- CC_module(srtIntegrated, cc_S_genes, cc_G2M_genes)

  DefaultAssay(srtIntegrated) <- "RNA"
  if (!is.null(srtMerge) & isTRUE(append)) {
    srtMerge_raw@graphs$BBKNN <- srtIntegrated@graphs$BBKNN
    srtMerge_raw@misc[[paste0(reduction_prefix, "Dims")]] <- srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]]
    srtMerge_raw[[paste0(reduction_prefix, "clusters")]] <- srtIntegrated[[paste0(reduction_prefix, "clusters")]]
    srtMerge_raw[[paste0(reduction_prefix, "umap")]] <- srtIntegrated[[paste0(reduction_prefix, "umap")]]
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

CSS_integrate <- function(srtList = NULL, srtMerge = NULL, append = TRUE,
                          normalization_method = "logCPM", batch = "orig.ident",
                          HVF_source = "separate", nHVF = 3000, hvf = NULL,
                          maxPC = 100, resolution = 0.8,
                          reduction = c("tsne", "umap"), reduction_prefix = "CSS_",
                          cc_S_genes = NULL, cc_G2M_genes = NULL,
                          exogenous_genes = NULL, ...) {
  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
    checked <- Check_srtMerge(srtMerge,
      normalization_method = normalization_method, batch = batch,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = NULL
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(srtList,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
  }

  srtMerge <- ScaleData(object = srtMerge, features = rownames(srtMerge))
  srtMerge <- RunPCA(object = srtMerge, npcs = maxPC, features = hvf)
  dims <- 1:ceiling(maxLikGlobalDimEst(data = Embeddings(srtMerge, reduction = "pca"), k = 20, iterations = 100)[["dim.est"]])

  srtIntegrated <- cluster_sim_spectrum(
    object = srtMerge,
    use_dr = "pca",
    dims_use = dims,
    var_genes = hvf,
    label_tag = batch,
    cluster_resolution = 0.4,
    # corr_method = "pearson",
    # spectrum_type = "corr_kernel"
    reduction.name = "CSS",
    reduction.key = "CSS_"
  )
  srtMerge <- NULL

  srtIntegrated <- Check_srtIntegrated(srtIntegrated, hvf = hvf, batch = batch)

  dims <- 1:ncol(Embeddings(srtIntegrated, reduction = "CSS"))
  srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]] <- dims

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "CSS", dims = dims, force.recalc = T)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srtIntegrated <- BuildClusterTree(srtIntegrated, features = hvf, slot = "data", reorder = T, reorder.numeric = T)
  srtIntegrated[[paste0(reduction_prefix, "clusters")]] <- Idents(srtIntegrated)

  if ("umap" %in% reduction) {
    srtIntegrated <- RunUMAP(object = srtIntegrated, reduction = "CSS", dims = dims, n.components = 2, umap.method = "uwot-learn", reduction.name = paste0(reduction_prefix, "umap"))
  }
  if ("tsne" %in% reduction) {
    srtIntegrated <- RunTSNE(
      object = srtIntegrated, reduction = "CSS", dims = dims, dim.embed = 2, tsne.method = "Rtsne", reduction.name = paste0(reduction_prefix, "tsne"),
      perplexity = max(ceiling(ncol(srtIntegrated) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  srtIntegrated <- CC_module(srtIntegrated, cc_S_genes, cc_G2M_genes)

  DefaultAssay(srtIntegrated) <- "RNA"
  if (!is.null(srtMerge) & isTRUE(append)) {
    srtMerge_raw@reductions$CSS <- srtIntegrated@reductions$CSS
    srtMerge_raw@misc[[paste0(reduction_prefix, "Dims")]] <- srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]]
    srtMerge_raw[[paste0(reduction_prefix, "clusters")]] <- srtIntegrated[[paste0(reduction_prefix, "clusters")]]
    for (i in reduction) {
      srtMerge_raw[[paste0(reduction_prefix, i)]] <- srtIntegrated[[paste0(reduction_prefix, i)]]
    }
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

LIGER_integrate <- function(srtList = NULL, srtMerge = NULL, append = TRUE,
                            normalization_method = "logCPM", batch = "orig.ident",
                            HVF_source = "separate", nHVF = 3000, hvf = NULL,
                            maxPC = 100, resolution = 0.8,
                            reduction = c("tsne", "umap"), reduction_prefix = "LIGER_",
                            cc_S_genes = NULL, cc_G2M_genes = NULL,
                            exogenous_genes = NULL, ...) {
  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
    checked <- Check_srtMerge(srtMerge,
      normalization_method = normalization_method, batch = batch,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = NULL
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(srtList,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
  }

  srtMerge <- ScaleData(object = srtMerge, features = hvf, split.by = batch, do.center = FALSE)

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

  srtIntegrated <- Check_srtIntegrated(srtIntegrated, hvf = hvf, batch = batch)

  dims <- 1:ncol(Embeddings(srtIntegrated, reduction = "LIGER"))
  srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]] <- dims

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "LIGER", dims = dims, force.recalc = T)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srtIntegrated <- BuildClusterTree(srtIntegrated, features = hvf, slot = "data", reorder = T, reorder.numeric = T)
  srtIntegrated[[paste0(reduction_prefix, "clusters")]] <- Idents(srtIntegrated)

  if ("umap" %in% reduction) {
    srtIntegrated <- RunUMAP(object = srtIntegrated, reduction = "LIGER", dims = dims, n.components = 2, umap.method = "uwot-learn", reduction.name = paste0(reduction_prefix, "umap"))
  }
  if ("tsne" %in% reduction) {
    srtIntegrated <- RunTSNE(
      object = srtIntegrated, reduction = "LIGER", dims = dims, dim.embed = 2, tsne.method = "Rtsne", reduction.name = paste0(reduction_prefix, "tsne"),
      perplexity = max(ceiling(ncol(srtIntegrated) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  srtIntegrated <- CC_module(srtIntegrated, cc_S_genes, cc_G2M_genes)

  DefaultAssay(srtIntegrated) <- "RNA"
  if (!is.null(srtMerge) & isTRUE(append)) {
    srtMerge_raw@reductions$LIGER <- srtIntegrated@reductions$LIGER
    srtMerge_raw@misc[[paste0(reduction_prefix, "Dims")]] <- srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]]
    srtMerge_raw[[paste0(reduction_prefix, "clusters")]] <- srtIntegrated[[paste0(reduction_prefix, "clusters")]]
    for (i in reduction) {
      srtMerge_raw[[paste0(reduction_prefix, i)]] <- srtIntegrated[[paste0(reduction_prefix, i)]]
    }
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

scMerge_integrate <- function(srtList = NULL, srtMerge = NULL, append = TRUE,
                              normalization_method = "logCPM", batch = "orig.ident",
                              HVF_source = "separate", nHVF = 3000, hvf = NULL,
                              maxPC = 100, resolution = 0.8,
                              reduction = c("tsne", "umap"), reduction_prefix = "scMerge_",
                              cc_S_genes = NULL, cc_G2M_genes = NULL,
                              exogenous_genes = NULL, ...) {
  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
    checked <- Check_srtMerge(srtMerge,
      normalization_method = normalization_method, batch = batch,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = NULL
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(srtList,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
  }

  sce <- as.SingleCellExperiment(srtMerge)
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
  srtIntegrated <- as.Seurat(sce_scMerge,
    counts = "counts",
    data = "scMerge",
    assay = "scMerge"
  )
  srtIntegrated@assays$RNA <- srtMerge@assays$RNA
  DefaultAssay(srtIntegrated) <- "scMerge"
  sce_scMerge <- sce <- NULL

  srtIntegrated <- Check_srtIntegrated(srtIntegrated, hvf = hvf, batch = batch)

  srtIntegrated <- ScaleData(srtIntegrated, features = hvf)
  srtIntegrated <- RunPCA(object = srtIntegrated, npcs = maxPC, features = hvf, reduction.name = paste0(reduction_prefix, "pca"))
  dims <- 1:ceiling(maxLikGlobalDimEst(data = Embeddings(srtIntegrated, reduction = paste0(reduction_prefix, "pca")), k = 20, iterations = 100)[["dim.est"]])
  srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]] <- dims

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"), dims = dims, force.recalc = T)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srtIntegrated <- BuildClusterTree(srtIntegrated, features = hvf, slot = "data", reorder = T, reorder.numeric = T)
  srtIntegrated[[paste0(reduction_prefix, "clusters")]] <- Idents(srtIntegrated)

  if ("umap" %in% reduction) {
    srtIntegrated <- RunUMAP(object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"), dims = dims, n.components = 2, umap.method = "uwot-learn", reduction.name = paste0(reduction_prefix, "umap"))
  }
  if ("tsne" %in% reduction) {
    srtIntegrated <- RunTSNE(
      object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"), dims = dims, dim.embed = 2, tsne.method = "Rtsne", reduction.name = paste0(reduction_prefix, "tsne"),
      perplexity = max(ceiling(ncol(srtIntegrated) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  srtIntegrated <- CC_module(srtIntegrated, cc_S_genes, cc_G2M_genes)

  DefaultAssay(srtIntegrated) <- "RNA"
  if (!is.null(srtMerge) & isTRUE(append)) {
    srtMerge_raw@assays$scMerge <- srtIntegrated@assays$scMerge
    srtMerge_raw[[paste0(reduction_prefix, "pca")]] <- srtIntegrated[[paste0(reduction_prefix, "pca")]]
    srtMerge_raw@misc[[paste0(reduction_prefix, "Dims")]] <- srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]]
    srtMerge_raw[[paste0(reduction_prefix, "clusters")]] <- srtIntegrated[[paste0(reduction_prefix, "clusters")]]
    for (i in reduction) {
      srtMerge_raw[[paste0(reduction_prefix, i)]] <- srtIntegrated[[paste0(reduction_prefix, i)]]
    }
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

ZINBWaVE_integrate <- function(srtList = NULL, srtMerge = NULL, append = TRUE,
                               normalization_method = "logCPM", batch = "orig.ident",
                               HVF_source = "separate", nHVF = 3000, hvf = NULL,
                               maxPC = 100, resolution = 0.8,
                               reduction = c("tsne", "umap"), reduction_prefix = "ZINBWaVE_",
                               cc_S_genes = NULL, cc_G2M_genes = NULL,
                               exogenous_genes = NULL, ...) {
  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
    checked <- Check_srtMerge(srtMerge,
      normalization_method = normalization_method, batch = batch,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = NULL
    )
    srtMerge <- checked[["srtMerge"]]
    hvf <- checked[["hvf"]]
  }
  if (!is.null(srtList)) {
    checked <- Check_srtList(srtList,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
      exogenous_genes = exogenous_genes
    )
    srtList <- checked[["srtList"]]
    hvf <- checked[["hvf"]]
    srtMerge <- Reduce(function(x, y) merge(x, y), srtList)
    VariableFeatures(srtMerge) <- hvf
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
        state <- 1
      }, error = function(error) {
        message(error)
        cat("Resampling the cells for zinbsurf ......\n")
        state <- 0
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

  srtIntegrated <- Check_srtIntegrated(srtIntegrated, hvf = hvf, batch = batch)

  srtIntegrated <- NormalizeData(object = srtIntegrated, normalization.method = "LogNormalize")
  srtIntegrated <- ScaleData(srtIntegrated, features = hvf)
  srtIntegrated <- RunPCA(object = srtIntegrated, npcs = maxPC, features = hvf, reduction.name = paste0(reduction_prefix, "pca"))
  dims <- 1:ceiling(maxLikGlobalDimEst(data = Embeddings(srtIntegrated, reduction = paste0(reduction_prefix, "pca")), k = 20, iterations = 100)[["dim.est"]])
  srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]] <- dims

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"), dims = dims, force.recalc = T)
  srtIntegrated <- FindClusters(object = srtIntegrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srtIntegrated <- BuildClusterTree(srtIntegrated, features = hvf, slot = "data", reorder = T, reorder.numeric = T)
  srtIntegrated[[paste0(reduction_prefix, "clusters")]] <- Idents(srtIntegrated)

  if ("umap" %in% reduction) {
    srtIntegrated <- RunUMAP(object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"), dims = dims, n.components = 2, umap.method = "uwot-learn", reduction.name = paste0(reduction_prefix, "umap"))
  }
  if ("tsne" %in% reduction) {
    srtIntegrated <- RunTSNE(
      object = srtIntegrated, reduction = paste0(reduction_prefix, "pca"), dims = dims, dim.embed = 2, tsne.method = "Rtsne", reduction.name = paste0(reduction_prefix, "tsne"),
      perplexity = max(ceiling(ncol(srtIntegrated) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  srtIntegrated <- CC_module(srtIntegrated, cc_S_genes, cc_G2M_genes)

  DefaultAssay(srtIntegrated) <- "RNA"
  if (!is.null(srtMerge) & isTRUE(append)) {
    srtMerge_raw@assays$ZINBWaVE <- srtIntegrated@assays$ZINBWaVE
    srtMerge_raw[[paste0(reduction_prefix, "pca")]] <- srtIntegrated[[paste0(reduction_prefix, "pca")]]
    srtMerge_raw@misc[[paste0(reduction_prefix, "Dims")]] <- srtIntegrated@misc[[paste0(reduction_prefix, "Dims")]]
    srtMerge_raw[[paste0(reduction_prefix, "clusters")]] <- srtIntegrated[[paste0(reduction_prefix, "clusters")]]
    for (i in reduction) {
      srtMerge_raw[[paste0(reduction_prefix, i)]] <- srtIntegrated[[paste0(reduction_prefix, i)]]
    }
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}


DEtest <- function(srt, FindAllMarkers = TRUE, FindPairMarkers = TRUE,
                   foldchange_threshold = 1.5, pvalue_threshold = 0.05, roc_threshold = 0.4,
                   BPPARAM = MulticoreParam(), ...) {
  if (isTRUE(FindAllMarkers)) {
    DefaultAssay(srt) <- "RNA"
    AllMarkers_Wilcoxon <- FindAllMarkers(
      object = srt, only.pos = T, logfc.threshold = log2(foldchange),
      test.use = "wilcox", return.thresh = pvalue_threshold # , latent.vars = "orig.ident"
    )
    AllMarkers_ROC <- FindAllMarkers(
      object = srt, only.pos = T, logfc.threshold = log2(foldchange),
      test.use = "roc", return.thresh = roc_threshold
    )
    srt@tools$FindAllMarkers <- setNames(
      object = list(AllMarkers_Wilcoxon, AllMarkers_ROC),
      nm = c("AllMarkers_Wilcoxon", "AllMarkers_ROC")
    )
  }

  if (isTRUE(FindPairMarkers)) {
    DefaultAssay(srt) <- "RNA"
    pair <- expand.grid(x = levels(Idents(srt)), y = levels(Idents(srt)))
    pair <- pair[pair[, 1] != pair[, 2], ]
    PairMarkers_Wilcoxon <- bplapply(1:nrow(pair), function(i) {
      res <- FindMarkers(
        ident.1 = as.character(pair[i, 1]), ident.2 = as.character(pair[i, 2]),
        object = srt, only.pos = T, logfc.threshold = log2(foldchange),
        test.use = "wilcox" # , latent.vars = "orig.ident"
      )
      res[, "ident.1"] <- as.character(pair[i, 1])
      res[, "ident.2"] <- as.character(pair[i, 2])
      res[, "gene"] <- rownames(res)
      res <- res[res[["p_val"]] < pvalue_threshold, ]

      return(res)
    }, BPPARAM = BPPARAM)
    PairMarkers_Wilcoxon <- bind_rows(PairMarkers_Wilcoxon)
    PairMarkers_Wilcoxon[, "DEnumber"] <- table(PairMarkers_Wilcoxon[["gene"]])

    srt@tools$FindPairMarkers <- setNames(
      object = PairMarkers_Wilcoxon,
      nm = "PairMarkers_Wilcoxon"
    )

    # exp1<- GetAssayData(object = srt,assay = "RNA",slot = "data")["NEUROD1",WhichCells(srt,idents = as.character(pair[1,1]))]
    # exp2<- GetAssayData(object = srt,assay = "RNA",slot = "data")["NEUROD1",WhichCells(srt,idents = as.character(pair[2,1]))]
    # sum(exp1!=0)/length(exp1)
    # sum(exp2!=0)/length(exp2)
    # log2(mean(exp(exp1))/mean(exp(exp2)))
  }

  return(srt)
}

Standard_SCP <- function(srt, normalization_method = "logCPM", nHVF = 3000, hvf = NULL,
                         maxPC = 100, resolution = 0.8,
                         reduction = c("tsne", "umap"), reduction_prefix = "",
                         cc_S_genes = NULL, cc_G2M_genes = NULL,
                         exogenous_genes = NULL, ...) {
  DefaultAssay(srt) <- "RNA"

  if (!normalization_method %in% c("logCPM", "SCT")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT'",
      call. = FALSE
    )
  }

  if (identical(
    x = GetAssayData(srt, slot = "counts"),
    y = GetAssayData(srt, slot = "data")
  )) {
    srt <- NormalizeData(object = srt, normalization.method = "LogNormalize")
  }
  if (is.null(hvf)) {
    if (!"variance.standardized" %in% colnames(srt@assays$RNA@meta.features)) {
      srt <- FindVariableFeatures(srt)
    }
    m <- GetAssayData(srt, slot = "counts")
    gene_keep <- rownames(m)[Matrix::rowSums(m > 5) > 5]
    VariableFeatures(srt) <- hvf <- HVFInfo(srt) %>%
      filter(variance.standardized > 1 &
        (!rownames(.) %in% exogenous_genes) &
        rownames(.) %in% gene_keep) %>%
      dplyr::arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
  } else {
    hvf <- hvf[hvf %in% rownames(GetAssayData(srt, slot = "counts"))]
    VariableFeatures(srt) <- hvf
  }
  if (nrow(GetAssayData(srt, slot = "scale.data")) != nrow(GetAssayData(srt, slot = "data"))) {
    srt <- ScaleData(object = srt, features = rownames(srt))
  }
  DefaultAssay(srt) <- "RNA"

  if (normalization_method %in% c("SCT")) {
    if (!"SCT" %in% Seurat::Assays(srt)) {
      srt <- SCTransform(
        object = srt,
        variable.features.n = nHVF,
        return.only.var.genes = FALSE,
        min_cells = 5,
        assay = "RNA"
      )
    }
    DefaultAssay(srt) <- "SCT"
    if (is.null(hvf)) {
      if (!"variance.standardized" %in% colnames(srt@assays$SCT@meta.features)) {
        srt <- FindVariableFeatures(srt)
      }
      m <- GetAssayData(srt, slot = "counts")
      gene_keep <- rownames(m)[Matrix::rowSums(m > 5) > 5]
      VariableFeatures(srt) <- hvf <- HVFInfo(srt) %>%
        filter(variance.standardized > 1 &
          (!rownames(.) %in% exogenous_genes) &
          rownames(.) %in% gene_keep) %>%
        dplyr::arrange(desc(variance.standardized)) %>%
        rownames(.) %>%
        head(n = nHVF)
    } else {
      hvf <- hvf[hvf %in% rownames(GetAssayData(srt, slot = "counts"))]
      VariableFeatures(srt) <- hvf
    }
  }

  srt <- RunPCA(object = srt, npcs = maxPC, features = hvf, reduction.name = paste0(reduction_prefix, "pca"))
  dims <- 1:ceiling(maxLikGlobalDimEst(data = Embeddings(srt, reduction = paste0(reduction_prefix, "pca")), k = 20, iterations = 100)[["dim.est"]])
  srt@misc[[paste0(reduction_prefix, "Dims")]] <- dims

  srt <- FindNeighbors(object = srt, reduction = paste0(reduction_prefix, "pca"), dims = dims, force.recalc = T)
  srt <- FindClusters(object = srt, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srt <- BuildClusterTree(srt, features = hvf, slot = "data", reorder = T, reorder.numeric = T)
  srt[[paste0(reduction_prefix, "clusters")]] <- Idents(srt)

  if ("umap" %in% reduction) {
    srt <- RunUMAP(object = srt, reduction = paste0(reduction_prefix, "pca"), dims = dims, n.components = 2, umap.method = "uwot-learn", reduction.name = paste0(reduction_prefix, "umap"))
  }
  if ("tsne" %in% reduction) {
    srt <- RunTSNE(
      object = srt, reduction = paste0(reduction_prefix, "pca"), dims = dims, dim.embed = 2, tsne.method = "Rtsne", reduction.name = paste0(reduction_prefix, "tsne"),
      perplexity = max(ceiling(ncol(srt) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  srt <- CC_module(srt, cc_S_genes, cc_G2M_genes)

  DefaultAssay(srt) <- "RNA"
  return(srt)
}

Integration_SCP <- function(srtList = NULL, srtMerge = NULL, append = TRUE,
                            integration_method = "Seurat", batch = "orig.ident",
                            normalization_method = "logCPM",
                            HVF_source = "separate", nHVF = 3000, hvf = NULL,
                            maxPC = 100, resolution = 0.8, reduction = c("tsne", "umap"),
                            cc_S_genes = NULL, cc_G2M_genes = NULL,
                            exogenous_genes = NULL, ...) {
  if (is.null(srtList) & is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (length(integration_method) == 1 & integration_method %in% c("Uncorrected", "Seurat", "fastMNN", "Harmony", "Scanorama", "BBKNN", "CSS", "LIGER", "scMerge", "ZINBWaVE")) {
    srtIntegrated <- base::do.call(
      what = paste0(integration_method, "_integrate"),
      args = as.list(match.call())
    )
    return(srtIntegrated)
  } else {
    stop(paste("Error!", integration_method, "is not a suppoted integration method!"),
      call. = FALSE
    )
  }
}

