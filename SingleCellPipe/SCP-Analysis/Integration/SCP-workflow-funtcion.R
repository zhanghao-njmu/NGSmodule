Check_scList <- function(sc_list, normalization_method = "logCPM",
                         HVF_source = "separate", nHVF = 3000, hvf = NULL,
                         exogenous_genes = NULL) {
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

  for (i in 1:length(sc_list)) {
    DefaultAssay(sc_list[[i]]) <- "RNA"
    if (identical(
      x = GetAssayData(sc_list[[i]], slot = "counts"),
      y = GetAssayData(sc_list[[i]], slot = "data")
    )) {
      sc_list[[i]] <- NormalizeData(object = sc_list[[i]], normalization.method = "LogNormalize", verbose = FALSE)
    }
    if (length(VariableFeatures(sc_list[[i]])) == 0) {
      sc_list[[i]] <- FindVariableFeatures(sc_list[[i]], verbose = FALSE)
    }
    VariableFeatures(sc_list[[i]]) <- HVFInfo(sc_list[[i]]) %>%
      filter(variance.standardized > 1 &
        (!rownames(.) %in% exogenous_genes)) %>%
      dplyr::arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
    if (nrow(GetAssayData(sc_list[[i]], slot = "scale.data")) == 0) {
      sc_list[[i]] <- ScaleData(object = sc_list[[i]], features = rownames(sc_list[[i]]), verbose = FALSE)
    }
    DefaultAssay(sc_list[[i]]) <- "RNA"

    if (normalization_method %in% c("SCT")) {
      if (!"SCT" %in% Seurat::Assays(sc_list[[i]])) {
        sc_list[[i]] <- SCTransform(
          object = sc_list[[i]],
          variable.features.n = nHVF,
          return.only.var.genes = FALSE,
          assay = "RNA",
          verbose = FALSE
        )
      } else {
        DefaultAssay(sc_list[[i]]) <- "SCT"
      }
      if (length(VariableFeatures(sc_list[[i]])) == 0) {
        sc_list[[i]] <- FindVariableFeatures(sc_list[[i]], verbose = FALSE)
        VariableFeatures(sc_list[[i]]) <- HVFInfo(sc_list[[i]]) %>%
          filter(variance.standardized > 1 &
            (!rownames(.) %in% exogenous_genes)) %>%
          dplyr::arrange(desc(variance.standardized)) %>%
          rownames(.) %>%
          head(n = nHVF)
      } else {
        VariableFeatures(sc_list[[i]]) <- HVFInfo(sc_list[[i]], selection.method = "sctransform") %>%
          filter((!rownames(.) %in% exogenous_genes)) %>%
          dplyr::arrange(desc(residual_variance)) %>%
          rownames(.) %>%
          head(n = nHVF)
      }
    }
  }

  if (is.null(hvf)) {
    if (HVF_source == "global") {
      gene_common <- lapply(sc_list, function(x) {
        m <- GetAssayData(x, slot = "counts")
        gene_keep <- rownames(m)[Matrix::rowSums(m > 0) >= 5]
        return(gene_keep)
      }) %>% Reduce(intersect, .)
      sc_merge <- Reduce(function(x, y) merge(x, y), sc_list)
      hvf <- FindVariableFeatures(sc_merge, verbose = FALSE) %>%
        HVFInfo(.) %>%
        filter(variance.standardized > 1 &
          (!rownames(.) %in% exogenous_genes) &
          rownames(.) %in% gene_common) %>%
        dplyr::arrange(desc(variance.standardized)) %>%
        rownames(.) %>%
        head(n = nHVF)
    }
    if (HVF_source == "separate") {
      hvf <- SelectIntegrationFeatures(object.list = sc_list, nfeatures = nHVF, verbose = FALSE)
    }
  }

  if (normalization_method %in% c("SCT")) {
    sc_list <- PrepSCTIntegration(object.list = sc_list, anchor.features = hvf, verbose = FALSE)
  }

  return(list(sc_list = sc_list, hvf = hvf))
}

Check_srtIntegrated <- function(srt_integrated, hvf) {
  raw_DefaultAssay <- DefaultAssay(object = srt_integrated)

  DefaultAssay(object = srt_integrated) <- "RNA"
  if (identical(
    x = GetAssayData(srt_integrated, slot = "counts"),
    y = GetAssayData(srt_integrated, slot = "data")
  )) {
    srt_integrated <- NormalizeData(object = srt_integrated)
  }
  if (length(VariableFeatures(srt_integrated)) == 0) {
    VariableFeatures(srt_integrated) <- hvf
  }
  if (nrow(GetAssayData(srt_integrated, slot = "scale.data")) != nrow(GetAssayData(srt_integrated, slot = "data"))) {
    srt_integrated <- ScaleData(object = srt_integrated, features = rownames(srt_integrated))
  }

  DefaultAssay(object = srt_integrated) <- raw_DefaultAssay
  srt_integrated@project.name <- paste0(unique(srt_integrated[["orig.ident", drop = TRUE]]), collapse = ",")
  srt_integrated[["orig.ident"]] <- factor(srt_integrated[["orig.ident", drop = TRUE]],
    levels = unique(srt_integrated[["orig.ident", drop = TRUE]])
  )
  return(srt_integrated)
}

CC_module <- function(sc, cc_S_genes, cc_G2M_genes) {
  if (length(cc_S_genes) >= 3 & length(cc_G2M_genes) >= 3) {
    sc <- CellCycleScoring(
      object = sc,
      s.features = cc_S_genes,
      g2m.features = cc_G2M_genes,
      set.ident = FALSE
    )
    sc[["CC.Difference"]] <- sc[["S.Score"]] - sc[["G2M.Score"]]
    sc[["Phase"]] <- factor(sc[["Phase", drop = TRUE]], levels = c("G1", "S", "G2M"))
  }
  return(sc)
}

Standard_SCP <- function(sc, normalization_method = "logCPM", nHVF = 3000, hvf = NULL,
                         maxPC = 100, resolution = 0.8, reduction = c("tsne", "umap"),
                         cc_S_genes = Seurat::cc.genes.updated.2019$s.genes, cc_G2M_genes = Seurat::cc.genes.updated.2019$g2m.genes,
                         exogenous_genes = NULL) {
  DefaultAssay(sc) <- "RNA"
  sc@project.name <- paste0(unique(sc[["orig.ident", drop = TRUE]]), collapse = ",")
  if (!normalization_method %in% c("logCPM", "SCT")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT'",
      call. = FALSE
    )
  }

  if (identical(
    x = GetAssayData(sc, slot = "counts"),
    y = GetAssayData(sc, slot = "data")
  )) {
    sc <- NormalizeData(object = sc, normalization.method = "LogNormalize")
  }
  if (is.null(hvf)) {
    if (length(VariableFeatures(sc)) == 0) {
      sc <- FindVariableFeatures(sc)
    }
    VariableFeatures(sc) <- hvf <- HVFInfo(sc) %>%
      filter(variance.standardized > 1 &
        (!rownames(.) %in% exogenous_genes)) %>%
      dplyr::arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
  } else {
    VariableFeatures(sc) <- hvf
  }
  if (nrow(GetAssayData(sc, slot = "scale.data")) == 0) {
    sc <- ScaleData(object = sc, features = rownames(sc))
  }
  DefaultAssay(sc) <- "RNA"

  if (normalization_method %in% c("SCT")) {
    if (!"SCT" %in% Seurat::Assays(sc)) {
      sc <- SCTransform(
        object = sc,
        variable.features.n = nHVF,
        return.only.var.genes = FALSE,
        assay = "RNA"
      )
    }
    if (is.null(hvf)) {
      if (length(VariableFeatures(sc)) == 0) {
        sc <- FindVariableFeatures(sc)
        VariableFeatures(sc) <- hvf <- HVFInfo(sc) %>%
          filter(variance.standardized > 1 &
            (!rownames(.) %in% exogenous_genes)) %>%
          dplyr::arrange(desc(variance.standardized)) %>%
          rownames(.) %>%
          head(n = nHVF)
      } else {
        VariableFeatures(sc) <- hvf <- HVFInfo(sc, selection.method = "sctransform") %>%
          filter((!rownames(.) %in% exogenous_genes)) %>%
          dplyr::arrange(desc(residual_variance)) %>%
          rownames(.) %>%
          head(n = nHVF)
      }
    } else {
      VariableFeatures(sc) <- hvf
    }
    DefaultAssay(sc) <- "SCT"
  }

  sc <- RunPCA(object = sc, npcs = maxPC, features = hvf)
  PC_use <- ceiling(maxLikGlobalDimEst(data = Embeddings(sc, reduction = "pca"), k = 20, iterations = 100)[["dim.est"]])
  sc@misc$PC_use <- PC_use

  sc <- FindNeighbors(object = sc, reduction = "pca", dims = 1:PC_use, force.recalc = T)
  sc <- FindClusters(object = sc, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  sc <- BuildClusterTree(sc, features = hvf, slot = "data", reorder = T, reorder.numeric = T)
  sc$seurat_clusters <- Idents(sc)

  if ("umap" %in% reduction) {
    sc <- RunUMAP(object = sc, reduction = "pca", dims = 1:PC_use, n.components = 2, umap.method = "uwot-learn")
  }
  if ("tsne" %in% reduction) {
    sc <- RunTSNE(
      object = sc, reduction = "pca", dims = 1:PC_use, dim.embed = 2, tsne.method = "Rtsne",
      perplexity = max(ceiling(ncol(sc) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  sc <- CC_module(sc, cc_S_genes, cc_G2M_genes)

  DefaultAssay(sc) <- "RNA"
  return(sc)
}

Seurat_integrate <- function(sc_list, normalization_method = "logCPM",
                             HVF_source = "separate", nHVF = 3000, hvf = NULL,
                             maxPC = 100, resolution = 0.8, reduction = c("tsne", "umap"),
                             cc_S_genes = Seurat::cc.genes.updated.2019$s.genes, cc_G2M_genes = Seurat::cc.genes.updated.2019$g2m.genes,
                             exogenous_genes = NULL) {
  checked <- Check_scList(sc_list,
    normalization_method = normalization_method,
    HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
    exogenous_genes = exogenous_genes
  )
  sc_list <- checked[["sc_list"]]
  hvf <- checked[["hvf"]]

  srt_anchors <- FindIntegrationAnchors(
    object.list = sc_list,
    normalization.method = switch(normalization_method,
      "logCPM" = "LogNormalize", "SCT" = "SCT"
    ),
    anchor.features = hvf,
    dims = 1:30
  )
  srt_integrated <- IntegrateData(
    anchorset = srt_anchors,
    normalization.method = switch(normalization_method,
      "logCPM" = "LogNormalize", "SCT" = "SCT"
    ),
    dims = 1:30,
    features.to.integrate = c(
      hvf,
      Reduce(union, lapply(sc_list, VariableFeatures)),
      Reduce(intersect, c(
        lapply(sc_list, rownames),
        list(c(cc_S_genes, cc_G2M_genes))
      ))
    )
  )

  srt_integrated <- Check_srtIntegrated(srt_integrated, hvf)

  srt_integrated <- ScaleData(srt_integrated, features = hvf)
  srt_integrated <- RunPCA(object = srt_integrated, npcs = maxPC, features = hvf)
  PC_use <- ceiling(maxLikGlobalDimEst(data = Embeddings(srt_integrated, reduction = "pca"), k = 20, iterations = 100)[["dim.est"]])
  srt_integrated@misc$PC_use <- PC_use

  srt_integrated <- FindNeighbors(object = srt_integrated, reduction = "pca", dims = 1:PC_use, force.recalc = T)
  srt_integrated <- FindClusters(object = srt_integrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srt_integrated <- BuildClusterTree(srt_integrated, features = hvf, slot = "data", reorder = T, reorder.numeric = T)
  srt_integrated$seurat_clusters <- Idents(srt_integrated)

  if ("umap" %in% reduction) {
    srt_integrated <- RunUMAP(object = srt_integrated, reduction = "pca", dims = 1:PC_use, n.components = 2, umap.method = "uwot-learn")
  }
  if ("tsne" %in% reduction) {
    srt_integrated <- RunTSNE(
      object = srt_integrated, reduction = "pca", dims = 1:PC_use, dim.embed = 2, tsne.method = "Rtsne",
      perplexity = max(ceiling(ncol(srt_integrated) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  srt_integrated <- CC_module(srt_integrated, cc_S_genes, cc_G2M_genes)

  DefaultAssay(srt_integrated) <- "RNA"
  return(srt_integrated)
}

fastMNN_integrate <- function(sc_list, normalization_method = "logCPM",
                              HVF_source = "separate", nHVF = 3000, hvf = NULL,
                              maxPC = 100, resolution = 0.8, reduction = c("tsne", "umap"),
                              cc_S_genes = Seurat::cc.genes.updated.2019$s.genes, cc_G2M_genes = Seurat::cc.genes.updated.2019$g2m.genes,
                              exogenous_genes = NULL) {
  checked <- Check_scList(sc_list,
    normalization_method = normalization_method,
    HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
    exogenous_genes = exogenous_genes
  )
  sc_list <- checked[["sc_list"]]
  hvf <- checked[["hvf"]]

  srt_integrated <- RunFastMNN(
    object.list = sc_list,
    features = hvf,
    d = maxPC,
    BPPARAM = MulticoreParam()
  )
  raw_DefaultAssay <- DefaultAssay(object = srt_integrated)

  srt_integrated <- Check_srtIntegrated(srt_integrated, hvf)

  PC_use <- ceiling(maxLikGlobalDimEst(data = Embeddings(srt_integrated, reduction = "mnn"), k = 20, iterations = 100)[["dim.est"]])
  srt_integrated@misc$PC_use <- PC_use

  srt_integrated <- FindNeighbors(object = srt_integrated, reduction = "mnn", dims = 1:PC_use, force.recalc = T)
  srt_integrated <- FindClusters(object = srt_integrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srt_integrated <- BuildClusterTree(srt_integrated, features = hvf, slot = "data", reorder = T, reorder.numeric = T)
  srt_integrated$seurat_clusters <- Idents(srt_integrated)

  if ("umap" %in% reduction) {
    srt_integrated <- RunUMAP(object = srt_integrated, reduction = "mnn", dims = 1:PC_use, n.components = 2, umap.method = "uwot-learn")
  }
  if ("tsne" %in% reduction) {
    srt_integrated <- RunTSNE(
      object = srt_integrated, reduction = "mnn", dims = 1:PC_use, dim.embed = 2, tsne.method = "Rtsne",
      perplexity = max(ceiling(ncol(srt_integrated) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  srt_integrated <- CC_module(srt_integrated, cc_S_genes, cc_G2M_genes)

  DefaultAssay(srt_integrated) <- "RNA"
  return(srt_integrated)
}

Harmony_integrate <- function(sc_list, normalization_method = "logCPM",
                              HVF_source = "separate", nHVF = 3000, hvf = NULL,
                              maxPC = 100, resolution = 0.8, reduction = c("tsne", "umap"),
                              cc_S_genes = Seurat::cc.genes.updated.2019$s.genes, cc_G2M_genes = Seurat::cc.genes.updated.2019$g2m.genes,
                              exogenous_genes = NULL) {
  checked <- Check_scList(sc_list,
    normalization_method = normalization_method,
    HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
    exogenous_genes = exogenous_genes
  )
  sc_list <- checked[["sc_list"]]
  hvf <- checked[["hvf"]]

  sc_merge <- Reduce(function(x, y) merge(x, y), sc_list)
  VariableFeatures(sc_merge) <- hvf
  sc_merge <- ScaleData(object = sc_merge, features = rownames(sc_merge))
  sc_merge <- RunPCA(object = sc_merge, npcs = maxPC, features = hvf)
  PC_use <- ceiling(maxLikGlobalDimEst(data = Embeddings(sc_merge, reduction = "pca"), k = 20, iterations = 100)[["dim.est"]])

  srt_integrated <- RunHarmony(
    object = sc_merge,
    group.by.vars = "orig.ident",
    reduction = "pca",
    dims.use = 1:PC_use,
    block.size = 0.01,
    max.iter.harmony = 100,
    max.iter.cluster = 200,
    epsilon.cluster = 1e-16,
    epsilon.harmony = 1e-16,
    assay.use = DefaultAssay(sc_merge)
  )
  raw_DefaultAssay <- DefaultAssay(object = srt_integrated)
  sc_merge <- NULL

  srt_integrated <- Check_srtIntegrated(srt_integrated, hvf)

  PC_use <- ceiling(maxLikGlobalDimEst(data = Embeddings(srt_integrated, reduction = "harmony"), k = 20, iterations = 100)[["dim.est"]])
  srt_integrated@misc$PC_use <- PC_use

  srt_integrated <- FindNeighbors(object = srt_integrated, reduction = "harmony", dims = 1:PC_use, force.recalc = T)
  srt_integrated <- FindClusters(object = srt_integrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srt_integrated <- BuildClusterTree(srt_integrated, features = hvf, slot = "data", reorder = T, reorder.numeric = T)
  srt_integrated$seurat_clusters <- Idents(srt_integrated)

  if ("umap" %in% reduction) {
    srt_integrated <- RunUMAP(object = srt_integrated, reduction = "harmony", dims = 1:PC_use, n.components = 2, umap.method = "uwot-learn")
  }
  if ("tsne" %in% reduction) {
    srt_integrated <- RunTSNE(
      object = srt_integrated, reduction = "harmony", dims = 1:PC_use, dim.embed = 2, tsne.method = "Rtsne",
      perplexity = max(ceiling(ncol(srt_integrated) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  srt_integrated <- CC_module(srt_integrated, cc_S_genes, cc_G2M_genes)

  DefaultAssay(srt_integrated) <- "RNA"
  return(srt_integrated)
}

Scanorama_integrate <- function(sc_list, normalization_method = "logCPM",
                                HVF_source = "separate", nHVF = 3000, hvf = NULL,
                                maxPC = 100, resolution = 0.8, reduction = c("tsne", "umap"),
                                cc_S_genes = Seurat::cc.genes.updated.2019$s.genes, cc_G2M_genes = Seurat::cc.genes.updated.2019$g2m.genes,
                                exogenous_genes = NULL) {
  checked <- Check_scList(sc_list,
    normalization_method = normalization_method,
    HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
    exogenous_genes = exogenous_genes
  )
  sc_list <- checked[["sc_list"]]
  hvf <- checked[["hvf"]]

  require(reticulate)
  scanorama <- reticulate::import("scanorama")

  assaylist <- list()
  genelist <- list()
  for (i in 1:length(sc_list)) {
    assaylist[[i]] <- t(as.matrix(GetAssayData(object = sc_list[[i]], slot = "data")))
    genelist[[i]] <- rownames(sc_list[[i]])
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
  colnames(dim_reduction) <- paste0("PC_", 1:100)
  stdevs <- apply(dim_reduction, MARGIN = 2, FUN = sd)

  srt_integrated <- Reduce(function(x, y) merge(x, y), sc_list)
  srt_integrated[["integrated"]] <- CreateAssayObject(data = cor_value)
  srt_integrated[["scanorama"]] <- CreateDimReducObject(embeddings = dim_reduction, assay = "integrated", stdev = stdevs, key = "scanorama_")

  srt_integrated <- Check_srtIntegrated(srt_integrated, hvf)

  srt_integrated <- ScaleData(srt_integrated, features = hvf)
  srt_integrated <- RunPCA(object = srt_integrated, npcs = maxPC, features = hvf)
  PC_use <- ceiling(maxLikGlobalDimEst(data = Embeddings(srt_integrated, reduction = "pca"), k = 20, iterations = 100)[["dim.est"]])
  srt_integrated@misc$PC_use <- PC_use

  srt_integrated <- FindNeighbors(object = srt_integrated, reduction = "pca", dims = 1:PC_use, force.recalc = T)
  srt_integrated <- FindClusters(object = srt_integrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srt_integrated <- BuildClusterTree(srt_integrated, features = hvf, slot = "data", reorder = T, reorder.numeric = T)
  srt_integrated$seurat_clusters <- Idents(srt_integrated)

  if ("umap" %in% reduction) {
    srt_integrated <- RunUMAP(object = srt_integrated, reduction = "pca", dims = 1:PC_use, n.components = 2, umap.method = "uwot-learn")
  }
  if ("tsne" %in% reduction) {
    srt_integrated <- RunTSNE(
      object = srt_integrated, reduction = "pca", dims = 1:PC_use, dim.embed = 2, tsne.method = "Rtsne",
      perplexity = max(ceiling(ncol(srt_integrated) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  srt_integrated <- CC_module(srt_integrated, cc_S_genes, cc_G2M_genes)

  DefaultAssay(srt_integrated) <- "RNA"
  return(srt_integrated)
}

BBKNN_integrate <- function(sc_list, normalization_method = "logCPM",
                            HVF_source = "separate", nHVF = 3000, hvf = NULL,
                            maxPC = 100, resolution = 0.8,
                            cc_S_genes = Seurat::cc.genes.updated.2019$s.genes, cc_G2M_genes = Seurat::cc.genes.updated.2019$g2m.genes,
                            exogenous_genes = NULL) {
  checked <- Check_scList(sc_list,
    normalization_method = normalization_method,
    HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
    exogenous_genes = exogenous_genes
  )
  sc_list <- checked[["sc_list"]]
  hvf <- checked[["hvf"]]

  sc_merge <- Reduce(function(x, y) merge(x, y), sc_list)
  VariableFeatures(sc_merge) <- hvf
  sc_merge <- ScaleData(object = sc_merge, features = rownames(sc_merge))
  sc_merge <- RunPCA(object = sc_merge, npcs = maxPC, features = hvf)

  bbknn <- reticulate::import("bbknn", convert = FALSE)
  pca <- reticulate::r_to_py(Embeddings(sc_merge, reduction = "pca"))

  bem <- bbknn$bbknn_pca_matrix(pca, batch_list = sc_merge[["orig.ident", drop = TRUE]])
  bem <- reticulate::py_to_r(bem)
  bbknn_graph <- as.matrix(bem[[2]])
  rownames(bbknn_graph) <- colnames(bbknn_graph) <- colnames(sc_merge)
  sc_merge@graphs$bbknn <- as.Graph(bbknn_graph)
  sc_merge <- FindClusters(object = sc_merge, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000, graph.name = "bbknn")
  srt_integrated <- sc_merge

  raw_DefaultAssay <- DefaultAssay(object = srt_integrated)
  sc_merge <- NULL

  srt_integrated <- Check_srtIntegrated(srt_integrated, hvf)

  srt_integrated <- BuildClusterTree(srt_integrated, features = hvf, slot = "data", reorder = T, reorder.numeric = T)
  srt_integrated$seurat_clusters <- Idents(srt_integrated)

  srt_integrated <- RunUMAP(object = srt_integrated, graph = "bbknn", umap.method = "umap-learn")

  srt_integrated <- CC_module(srt_integrated, cc_S_genes, cc_G2M_genes)

  DefaultAssay(srt_integrated) <- "RNA"
  return(srt_integrated)
}

CSS_integrate <- function(sc_list, normalization_method = "logCPM",
                          HVF_source = "separate", nHVF = 3000, hvf = NULL,
                          maxPC = 100, resolution = 0.8, reduction = c("tsne", "umap"),
                          cc_S_genes = Seurat::cc.genes.updated.2019$s.genes, cc_G2M_genes = Seurat::cc.genes.updated.2019$g2m.genes,
                          exogenous_genes = NULL) {
  checked <- Check_scList(sc_list,
    normalization_method = normalization_method,
    HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
    exogenous_genes = exogenous_genes
  )
  sc_list <- checked[["sc_list"]]
  hvf <- checked[["hvf"]]

  sc_merge <- Reduce(function(x, y) merge(x, y), sc_list)
  VariableFeatures(sc_merge) <- hvf
  sc_merge <- ScaleData(object = sc_merge, features = rownames(sc_merge))
  sc_merge <- RunPCA(object = sc_merge, npcs = maxPC, features = hvf)

  srt_integrated <- cluster_sim_spectrum(
    object = sc_merge,
    dims_use = 1:20,
    var_genes = hvf,
    label_tag = "orig.ident",
    cluster_resolution = 0.4,
    corr_method = "pearson",
    spectrum_type = "corr_kernel"
  )
  raw_DefaultAssay <- DefaultAssay(object = srt_integrated)
  sc_merge <- NULL

  srt_integrated <- Check_srtIntegrated(srt_integrated, hvf)

  PC_use <- ncol(Embeddings(srt_integrated, reduction = "css"))
  srt_integrated@misc$PC_use <- PC_use

  srt_integrated <- FindNeighbors(object = srt_integrated, reduction = "css", dims = 1:PC_use, force.recalc = T)
  srt_integrated <- FindClusters(object = srt_integrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srt_integrated <- BuildClusterTree(srt_integrated, features = hvf, slot = "data", reorder = T, reorder.numeric = T)
  srt_integrated$seurat_clusters <- Idents(srt_integrated)

  if ("umap" %in% reduction) {
    srt_integrated <- RunUMAP(object = srt_integrated, reduction = "css", dims = 1:PC_use, n.components = 2, umap.method = "uwot-learn")
  }
  if ("tsne" %in% reduction) {
    srt_integrated <- RunTSNE(
      object = srt_integrated, reduction = "css", dims = 1:PC_use, dim.embed = 2, tsne.method = "Rtsne",
      perplexity = max(ceiling(ncol(srt_integrated) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  srt_integrated <- CC_module(srt_integrated, cc_S_genes, cc_G2M_genes)

  DefaultAssay(srt_integrated) <- "RNA"
  return(srt_integrated)
}

LIGER_integrate <- function(sc_list, normalization_method = "logCPM",
                            HVF_source = "separate", nHVF = 3000, hvf = NULL,
                            maxPC = 100, resolution = 0.8, reduction = c("tsne", "umap"),
                            cc_S_genes = Seurat::cc.genes.updated.2019$s.genes, cc_G2M_genes = Seurat::cc.genes.updated.2019$g2m.genes,
                            exogenous_genes = NULL) {
  checked <- Check_scList(sc_list,
    normalization_method = normalization_method,
    HVF_source = HVF_source, nHVF = nHVF, hvf = hvf,
    exogenous_genes = exogenous_genes
  )
  sc_list <- checked[["sc_list"]]
  hvf <- checked[["hvf"]]

  sc_merge <- Reduce(function(x, y) merge(x, y), sc_list)
  VariableFeatures(sc_merge) <- hvf
  sc_merge <- ScaleData(object = sc_merge, features = hvf, split.by = "orig.ident", do.center = FALSE)

  sc_merge <- RunOptimizeALS(sc_merge, k = 20, lambda = 5, split.by = "orig.ident")
  srt_integrated <- RunQuantileNorm(sc_merge, split.by = "orig.ident")
  raw_DefaultAssay <- DefaultAssay(object = srt_integrated)
  sc_merge <- NULL

  srt_integrated <- Check_srtIntegrated(srt_integrated, hvf)

  PC_use <- srt_integrated@misc$PC_use <- ncol(srt_integrated[["iNMF"]])

  srt_integrated <- FindNeighbors(object = srt_integrated, reduction = "iNMF", dims = 1:PC_use, force.recalc = T)
  srt_integrated <- FindClusters(object = srt_integrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srt_integrated <- BuildClusterTree(srt_integrated, features = hvf, slot = "data", reorder = T, reorder.numeric = T)
  srt_integrated$seurat_clusters <- Idents(srt_integrated)

  if ("umap" %in% reduction) {
    srt_integrated <- RunUMAP(object = srt_integrated, reduction = "iNMF", dims = 1:PC_use, n.components = 2, umap.method = "uwot-learn")
  }
  if ("tsne" %in% reduction) {
    srt_integrated <- RunTSNE(
      object = srt_integrated, reduction = "iNMF", dims = 1:PC_use, dim.embed = 2, tsne.method = "Rtsne",
      perplexity = max(ceiling(ncol(srt_integrated) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  srt_integrated <- CC_module(srt_integrated, cc_S_genes, cc_G2M_genes)

  DefaultAssay(srt_integrated) <- "RNA"
  return(srt_integrated)
}




DEtest <- function(srt, FindAllMarkers = TRUE, FindPairMarkers = TRUE,
                   foldchange_threshold = 1.5, pvalue_threshold = 0.05, roc_threshold = 0.4,
                   BPPARAM = MulticoreParam()) {
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
