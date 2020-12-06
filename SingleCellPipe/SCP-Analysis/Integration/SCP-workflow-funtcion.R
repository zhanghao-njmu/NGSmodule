Standard_SCP <- function(sc, nHVF = 3000,
                         maxPC = 100, resolution = 0.8, reduction = c("tsne", "umap"),
                         cc_S_genes = Seurat::cc.genes.updated.2019$s.genes, cc_G2M_genes = Seurat::cc.genes.updated.2019$g2m.genes,
                         exogenous_genes = NULL, assay = "RNA") {
  DefaultAssay(sc) <- assay
  sc@project.name <- paste0(unique(sc[["orig.ident", drop = TRUE]]), collapse = ",")
  if (identical(
    x = GetAssayData(sc, slot = "counts"),
    y = GetAssayData(sc, slot = "data")
  )) {
    sc <- NormalizeData(object = sc)
  }
  if (length(VariableFeatures(sc)) == 0) {
    sc <- FindVariableFeatures(sc)
  }
  VariableFeatures(sc) <- hvf <- HVFInfo(sc) %>%
    filter(variance.standardized > 1 &
      (!rownames(.) %in% exogenous_genes)) %>%
    dplyr::arrange(desc(variance.standardized)) %>%
    rownames(.) %>%
    head(n = nHVF)
  if (nrow(GetAssayData(sc, slot = "scale.data")) == 0) {
    sc <- ScaleData(object = sc, features = rownames(sc))
  }

  sc <- RunPCA(object = sc, npcs = maxPC, features = VariableFeatures(sc))
  PC_use <- ceiling(maxLikGlobalDimEst(data = Embeddings(sc, reduction = "pca"), k = 20, iterations = 100)[["dim.est"]])
  sc@misc$PC_use <- PC_use

  sc <- FindNeighbors(object = sc, reduction = "pca", dims = 1:PC_use, force.recalc = T)
  sc <- FindClusters(object = sc, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  sc <- BuildClusterTree(sc, features = hvf, slot = "scale.data", reorder = T, reorder.numeric = T)
  sc$seurat_clusters <- Idents(sc)

  if ("umap" %in% reduction) {
    sc <- RunUMAP(object = sc, reduction = "pca", dims = 1:PC_use, n.components = 2, umap.method = "uwot")
  }
  if ("tsne" %in% reduction) {
    sc <- RunTSNE(
      object = sc, reduction = "pca", dims = 1:PC_use, dim.embed = 2, tsne.method = "Rtsne",
      perplexity = max(ceiling(ncol(sc) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

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

  DefaultAssay(sc) <- "RNA"
  return(sc)
}

SCTransform_SCP <- function(sc, nHVF = 3000,
                            maxPC = 100, resolution = 0.8, reduction = c("tsne", "umap"),
                            cc_S_genes = Seurat::cc.genes.updated.2019$s.genes, cc_G2M_genes = Seurat::cc.genes.updated.2019$g2m.genes,
                            exogenous_genes = NULL, assay = "RNA") {
  DefaultAssay(sc) <- assay
  sc@project.name <- paste0(unique(sc[["orig.ident", drop = TRUE]]), collapse = ",")
  if (!"SCT" %in% Seurat::Assays(sc)) {
    sc <- SCTransform(
      object = sc,
      variable.features.n = nHVF,
      return.only.var.genes = FALSE,
      assay = "RNA"
    )
  }
  DefaultAssay(sc) <- "SCT"

  VariableFeatures(sc) <- hvf <- HVFInfo(sc) %>%
    filter((!rownames(.) %in% exogenous_genes)) %>%
    dplyr::arrange(desc(residual_variance)) %>%
    rownames(.) %>%
    head(n = nHVF)
  sc <- RunPCA(object = sc, npcs = maxPC, features = VariableFeatures(sc))
  PC_use <- ceiling(maxLikGlobalDimEst(data = Embeddings(sc, reduction = "pca"), k = 20, iterations = 100)[["dim.est"]])
  sc@misc$PC_use <- PC_use

  sc <- FindNeighbors(object = sc, reduction = "pca", dims = 1:PC_use, force.recalc = T)
  sc <- FindClusters(object = sc, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  sc <- BuildClusterTree(sc, slot = "scale.data", features = hvf, reorder = T, reorder.numeric = T)
  sc$seurat_clusters <- Idents(sc)

  if ("umap" %in% reduction) {
    sc <- RunUMAP(object = sc, reduction = "pca", dims = 1:PC_use, n.components = 2, umap.method = "uwot")
  }
  if ("tsne" %in% reduction) {
    sc <- RunTSNE(
      object = sc, reduction = "pca", dims = 1:PC_use, dim.embed = 2, tsne.method = "Rtsne",
      perplexity = max(ceiling(ncol(sc) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

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

  DefaultAssay(sc) <- "RNA"
  return(sc)
}

Standard_integrate <- function(sc_list, HVF_source = "separate", nHVF = 3000,
                               anchor_dims = 1:30, integrate_dims = 1:30,
                               maxPC = 100, resolution = 0.8, reduction = c("tsne", "umap"),
                               cc_S_genes = Seurat::cc.genes.updated.2019$s.genes, cc_G2M_genes = Seurat::cc.genes.updated.2019$g2m.genes,
                               exogenous_genes = NULL) {
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
      sc_list[[i]] <- NormalizeData(object = sc_list[[i]])
    }
    if (length(VariableFeatures(sc_list[[i]])) == 0) {
      sc_list[[i]] <- FindVariableFeatures(sc_list[[i]])
    }
    VariableFeatures(sc_list[[i]]) <- HVFInfo(sc_list[[i]]) %>%
      filter(variance.standardized > 1 &
        (!rownames(.) %in% exogenous_genes)) %>%
      dplyr::arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
  }
  if (HVF_source == "global") {
    # gene_common <- lapply(sc_list, rownames) %>% Reduce(intersect, .)
    gene_common <- lapply(sc_list, function(x) {
      m <- GetAssayData(x, assay = "RNA", slot = "counts")
      gene_keep <- rownames(m)[Matrix::rowSums(m > 0) >= 5]
      return(gene_keep)
    }) %>% Reduce(intersect, .)
    sc_merge <- Reduce(function(x, y) merge(x, y), sc_list)
    DefaultAssay(sc_merge) <- "RNA"
    hvf <- NormalizeData(object = sc_merge) %>%
      FindVariableFeatures(.) %>%
      HVFInfo(.) %>%
      filter(variance.standardized > 1 &
        (!rownames(.) %in% exogenous_genes) &
        rownames(.) %in% gene_common) %>%
      dplyr::arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
    sc_merge <- NULL
  }
  if (HVF_source == "separate") {
    hvf <- SelectIntegrationFeatures(object.list = sc_list, nfeatures = nHVF)
  }

  srt_anchors <- FindIntegrationAnchors(
    object.list = sc_list,
    normalization.method = "LogNormalize",
    anchor.features = hvf,
    dims = anchor_dims
  )
  srt_integrated <- IntegrateData(
    anchorset = srt_anchors,
    normalization.method = "LogNormalize",
    dims = integrate_dims,
    features.to.integrate = c(
      hvf,
      Reduce(union, lapply(sc_list, VariableFeatures)),
      Reduce(intersect, c(
        lapply(sc_list, rownames),
        list(c(cc_S_genes, cc_G2M_genes))
      ))
    )
  )

  DefaultAssay(object = srt_integrated) <- "RNA"
  if (identical(
    x = GetAssayData(srt_integrated, slot = "counts"),
    y = GetAssayData(srt_integrated, slot = "data")
  )) {
    srt_integrated <- NormalizeData(object = srt_integrated)
  }
  if (nrow(GetAssayData(srt_integrated, slot = "scale.data")) == 0) {
    srt_integrated <- ScaleData(object = srt_integrated, features = rownames(srt_integrated))
  }

  DefaultAssay(object = srt_integrated) <- "integrated"
  srt_integrated@project.name <- paste0(unique(srt_integrated[["orig.ident", drop = TRUE]]), collapse = ",")
  srt_integrated[["orig.ident"]] <- factor(srt_integrated[["orig.ident", drop = TRUE]],
    levels = unique(srt_integrated[["orig.ident", drop = TRUE]])
  )

  srt_integrated <- ScaleData(srt_integrated, features = hvf)
  srt_integrated <- RunPCA(object = srt_integrated, npcs = maxPC, features = hvf)
  PC_use <- ceiling(maxLikGlobalDimEst(data = Embeddings(srt_integrated, reduction = "pca"), k = 20, iterations = 100)[["dim.est"]])
  srt_integrated@misc$PC_use <- PC_use

  srt_integrated <- FindNeighbors(object = srt_integrated, reduction = "pca", dims = 1:PC_use, force.recalc = T)
  srt_integrated <- FindClusters(object = srt_integrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srt_integrated <- BuildClusterTree(srt_integrated, features = hvf, slot = "scale.data", reorder = T, reorder.numeric = T)
  srt_integrated$seurat_clusters <- Idents(srt_integrated)

  if ("umap" %in% reduction) {
    srt_integrated <- RunUMAP(object = srt_integrated, reduction = "pca", dims = 1:PC_use, n.components = 2, umap.method = "uwot")
  }
  if ("tsne" %in% reduction) {
    srt_integrated <- RunTSNE(
      object = srt_integrated, reduction = "pca", dims = 1:PC_use, dim.embed = 2, tsne.method = "Rtsne",
      perplexity = max(ceiling(ncol(srt_integrated) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  if (length(cc_S_genes) >= 3 & length(cc_G2M_genes) >= 3) {
    srt_integrated <- CellCycleScoring(
      object = srt_integrated,
      s.features = cc_S_genes,
      g2m.features = cc_G2M_genes,
      set.ident = FALSE
    )
    srt_integrated[["CC.Difference"]] <- srt_integrated[["S.Score"]] - srt_integrated[["G2M.Score"]]
    srt_integrated[["Phase"]] <- factor(srt_integrated[["Phase", drop = TRUE]], levels = c("G1", "S", "G2M"))
  }

  DefaultAssay(srt_integrated) <- "RNA"
  return(srt_integrated)
}

SCTransform_integrate <- function(sc_list, HVF_source = "separate", nHVF = 3000,
                                  anchor_dims = 1:30, integrate_dims = 1:30,
                                  maxPC = 100, resolution = 0.8, reduction = c("tsne", "umap"),
                                  cc_S_genes = Seurat::cc.genes.updated.2019$s.genes, cc_G2M_genes = Seurat::cc.genes.updated.2019$g2m.genes,
                                  exogenous_genes = NULL) {
  if (!HVF_source %in% c("global", "separate")) {
    stop("'HVF_source' must be one of: 'global','separate'",
      call. = FALSE
    )
  }

  for (i in 1:length(sc_list)) {
    if (!"SCT" %in% Seurat::Assays(sc_list[[i]])) {
      sc_list[[i]] <- SCTransform(
        object = sc_list[[i]],
        variable.features.n = nHVF,
        return.only.var.genes = TRUE,
        assay = "RNA"
      )
    }
    DefaultAssay(sc_list[[i]]) <- "SCT"
    VariableFeatures(sc_list[[i]]) <- HVFInfo(sc_list[[i]]) %>%
      filter((!rownames(.) %in% exogenous_genes)) %>%
      dplyr::arrange(desc(residual_variance)) %>%
      rownames(.) %>%
      head(n = nHVF)
  }
  if (HVF_source == "global") {
    gene_common <- lapply(sc_list, rownames) %>% Reduce(intersect, .)
    sc_merge <- Reduce(function(x, y) merge(x, y), sc_list)
    DefaultAssay(sc_merge) <- "RNA"
    hvf <- NormalizeData(object = sc_merge) %>%
      FindVariableFeatures(.) %>%
      HVFInfo(.) %>%
      filter(variance.standardized > 1 &
        (!rownames(.) %in% exogenous_genes) &
        rownames(.) %in% gene_common) %>%
      dplyr::arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
    sc_merge <- NULL
  }
  if (HVF_source == "separate") {
    hvf <- SelectIntegrationFeatures(object.list = sc_list, nfeatures = nHVF)
  }

  sc_list <- PrepSCTIntegration(object.list = sc_list, anchor.features = hvf)
  srt_anchors <- FindIntegrationAnchors(
    object.list = sc_list,
    normalization.method = "SCT",
    anchor.features = hvf,
    dims = anchor_dims
  )
  srt_integrated <- IntegrateData(
    anchorset = srt_anchors,
    normalization.method = "SCT",
    dims = integrate_dims,
    features.to.integrate = c(
      hvf,
      Reduce(union, lapply(sc_list, VariableFeatures)),
      Reduce(intersect, c(
        lapply(sc_list, rownames),
        list(c(cc_S_genes, cc_G2M_genes))
      ))
    )
  )

  DefaultAssay(object = srt_integrated) <- "RNA"
  if (identical(
    x = GetAssayData(srt_integrated, slot = "counts"),
    y = GetAssayData(srt_integrated, slot = "data")
  )) {
    srt_integrated <- NormalizeData(object = srt_integrated)
  }
  if (nrow(GetAssayData(srt_integrated, slot = "scale.data")) == 0) {
    srt_integrated <- ScaleData(object = srt_integrated, features = rownames(srt_integrated))
  }

  DefaultAssay(object = srt_integrated) <- "integrated"
  srt_integrated@project.name <- paste0(unique(srt_integrated[["orig.ident", drop = TRUE]]), collapse = ",")
  srt_integrated[["orig.ident"]] <- factor(srt_integrated[["orig.ident", drop = TRUE]],
    levels = unique(srt_integrated[["orig.ident", drop = TRUE]])
  )

  srt_integrated <- ScaleData(srt_integrated, features = hvf)
  srt_integrated <- RunPCA(object = srt_integrated, npcs = maxPC, features = hvf)
  PC_use <- ceiling(maxLikGlobalDimEst(data = Embeddings(srt_integrated, reduction = "pca"), k = 20, iterations = 100)[["dim.est"]])
  srt_integrated@misc$PC_use <- PC_use

  srt_integrated <- FindNeighbors(object = srt_integrated, reduction = "pca", dims = 1:PC_use, force.recalc = T)
  srt_integrated <- FindClusters(object = srt_integrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srt_integrated <- BuildClusterTree(srt_integrated, features = hvf, slot = "scale.data", reorder = T, reorder.numeric = T)
  srt_integrated$seurat_clusters <- Idents(srt_integrated)

  if ("umap" %in% reduction) {
    srt_integrated <- RunUMAP(object = srt_integrated, reduction = "pca", dims = 1:PC_use, n.components = 2, umap.method = "uwot")
  }
  if ("tsne" %in% reduction) {
    srt_integrated <- RunTSNE(
      object = srt_integrated, reduction = "pca", dims = 1:PC_use, dim.embed = 2, tsne.method = "Rtsne",
      perplexity = max(ceiling(ncol(srt_integrated) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  if (length(cc_S_genes) >= 3 & length(cc_G2M_genes) >= 3) {
    srt_integrated <- CellCycleScoring(
      object = srt_integrated,
      s.features = cc_S_genes,
      g2m.features = cc_G2M_genes,
      set.ident = FALSE
    )
    srt_integrated[["CC.Difference"]] <- srt_integrated[["S.Score"]] - srt_integrated[["G2M.Score"]]
    srt_integrated[["Phase"]] <- factor(srt_integrated[["Phase", drop = TRUE]], levels = c("G1", "S", "G2M"))
  }

  DefaultAssay(srt_integrated) <- "RNA"
  return(srt_integrated)
}

fastMNN_integrate <- function(sc_list, HVF_source = "separate", nHVF = 3000,
                              maxPC = 100, resolution = 0.8, reduction = c("tsne", "umap"),
                              cc_S_genes = Seurat::cc.genes.updated.2019$s.genes, cc_G2M_genes = Seurat::cc.genes.updated.2019$g2m.genes,
                              exogenous_genes = NULL) {
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
      sc_list[[i]] <- NormalizeData(object = sc_list[[i]])
    }
    if (length(VariableFeatures(sc_list[[i]])) == 0) {
      sc_list[[i]] <- FindVariableFeatures(sc_list[[i]])
    }
    VariableFeatures(sc_list[[i]]) <- HVFInfo(sc_list[[i]]) %>%
      filter(variance.standardized > 1 &
        (!rownames(.) %in% exogenous_genes)) %>%
      dplyr::arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
  }
  if (HVF_source == "global") {
    # gene_common <- lapply(sc_list, rownames) %>% Reduce(intersect, .)
    gene_common <- lapply(sc_list, function(x) {
      m <- GetAssayData(x, assay = "RNA", slot = "counts")
      gene_keep <- rownames(m)[Matrix::rowSums(m > 0) >= 5]
      return(gene_keep)
    }) %>% Reduce(intersect, .)
    sc_merge <- Reduce(function(x, y) merge(x, y), sc_list)
    DefaultAssay(sc_merge) <- "RNA"
    hvf <- NormalizeData(object = sc_merge) %>%
      FindVariableFeatures(.) %>%
      HVFInfo(.) %>%
      filter(variance.standardized > 1 &
        (!rownames(.) %in% exogenous_genes) &
        rownames(.) %in% gene_common) %>%
      dplyr::arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
    sc_merge <- NULL
  }
  if (HVF_source == "separate") {
    hvf <- SelectIntegrationFeatures(object.list = sc_list, nfeatures = nHVF)
  }

  srt_integrated <- RunFastMNN(
    object.list = sc_list,
    features = hvf,
    d = maxPC,
    BPPARAM = MulticoreParam()
  )

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
  if (nrow(GetAssayData(srt_integrated, slot = "scale.data")) == 0) {
    srt_integrated <- ScaleData(object = srt_integrated, features = rownames(srt_integrated))
  }

  srt_integrated@project.name <- paste0(unique(srt_integrated[["orig.ident", drop = TRUE]]), collapse = ",")
  srt_integrated[["orig.ident"]] <- factor(srt_integrated[["orig.ident", drop = TRUE]],
    levels = unique(srt_integrated[["orig.ident", drop = TRUE]])
  )

  PC_use <- ceiling(maxLikGlobalDimEst(data = Embeddings(srt_integrated, reduction = "mnn"), k = 20, iterations = 100)[["dim.est"]])
  srt_integrated@misc$PC_use <- PC_use

  srt_integrated <- FindNeighbors(object = srt_integrated, reduction = "mnn", dims = 1:PC_use, force.recalc = T)
  srt_integrated <- FindClusters(object = srt_integrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srt_integrated <- BuildClusterTree(srt_integrated, features = hvf, slot = "scale.data", reorder = T, reorder.numeric = T)
  srt_integrated$seurat_clusters <- Idents(srt_integrated)

  if ("umap" %in% reduction) {
    srt_integrated <- RunUMAP(object = srt_integrated, reduction = "mnn", dims = 1:PC_use, n.components = 2, umap.method = "uwot")
  }
  if ("tsne" %in% reduction) {
    srt_integrated <- RunTSNE(
      object = srt_integrated, reduction = "mnn", dims = 1:PC_use, dim.embed = 2, tsne.method = "Rtsne",
      perplexity = max(ceiling(ncol(srt_integrated) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  if (length(cc_S_genes) >= 3 & length(cc_G2M_genes) >= 3) {
    srt_integrated <- CellCycleScoring(
      object = srt_integrated,
      s.features = cc_S_genes,
      g2m.features = cc_G2M_genes,
      set.ident = FALSE
    )
    srt_integrated[["CC.Difference"]] <- srt_integrated[["S.Score"]] - srt_integrated[["G2M.Score"]]
    srt_integrated[["Phase"]] <- factor(srt_integrated[["Phase", drop = TRUE]], levels = c("G1", "S", "G2M"))
  }

  DefaultAssay(srt_integrated) <- "RNA"
  return(srt_integrated)
}

Harmony_integrate <- function(sc_list, HVF_source = "separate", nHVF = 3000,
                              maxPC = 100, resolution = 0.8, reduction = c("tsne", "umap"),
                              cc_S_genes = Seurat::cc.genes.updated.2019$s.genes, cc_G2M_genes = Seurat::cc.genes.updated.2019$g2m.genes,
                              exogenous_genes = NULL) {
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
      sc_list[[i]] <- NormalizeData(object = sc_list[[i]])
    }
    if (length(VariableFeatures(sc_list[[i]])) == 0) {
      sc_list[[i]] <- FindVariableFeatures(sc_list[[i]])
    }
    VariableFeatures(sc_list[[i]]) <- HVFInfo(sc_list[[i]]) %>%
      filter(variance.standardized > 1 &
        (!rownames(.) %in% exogenous_genes)) %>%
      dplyr::arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
  }

  sc_merge <- Reduce(function(x, y) merge(x, y), sc_list)
  DefaultAssay(sc_merge) <- "RNA"
  sc_merge <- NormalizeData(object = sc_merge) %>%
    FindVariableFeatures(.)

  if (HVF_source == "global") {
    # gene_common <- lapply(sc_list, rownames) %>% Reduce(intersect, .)
    gene_common <- lapply(sc_list, function(x) {
      m <- GetAssayData(x, assay = "RNA", slot = "counts")
      gene_keep <- rownames(m)[Matrix::rowSums(m > 0) >= 5]
      return(gene_keep)
    }) %>% Reduce(intersect, .)
    hvf <- HVFInfo(sc_merge) %>%
      filter(variance.standardized > 1 &
        (!rownames(.) %in% exogenous_genes) &
        rownames(.) %in% gene_common) %>%
      dplyr::arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
  }
  if (HVF_source == "separate") {
    hvf <- SelectIntegrationFeatures(object.list = sc_list, nfeatures = nHVF)
  }

  sc_merge <- ScaleData(object = sc_merge, features = rownames(sc_merge))
  sc_merge <- RunPCA(object = sc_merge, npcs = maxPC, features = hvf)

  srt_integrated <- RunHarmony(
    object = sc_merge,
    group.by.vars = "orig.ident",
    reduction = "pca",
    max.iter.harmony = 100,
    max.iter.cluster = 200,
    epsilon.cluster = 1e-16,
    epsilon.harmony = 1e-16
  )
  sc_merge <- NULL

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
  if (nrow(GetAssayData(srt_integrated, slot = "scale.data")) == 0) {
    srt_integrated <- ScaleData(object = srt_integrated, features = rownames(srt_integrated))
  }

  srt_integrated@project.name <- paste0(unique(srt_integrated[["orig.ident", drop = TRUE]]), collapse = ",")
  srt_integrated[["orig.ident"]] <- factor(srt_integrated[["orig.ident", drop = TRUE]],
    levels = unique(srt_integrated[["orig.ident", drop = TRUE]])
  )

  PC_use <- ceiling(maxLikGlobalDimEst(data = Embeddings(srt_integrated, reduction = "harmony"), k = 20, iterations = 100)[["dim.est"]])
  srt_integrated@misc$PC_use <- PC_use

  srt_integrated <- FindNeighbors(object = srt_integrated, reduction = "harmony", dims = 1:PC_use, force.recalc = T)
  srt_integrated <- FindClusters(object = srt_integrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srt_integrated <- BuildClusterTree(srt_integrated, features = hvf, slot = "scale.data", reorder = T, reorder.numeric = T)
  srt_integrated$seurat_clusters <- Idents(srt_integrated)

  if ("umap" %in% reduction) {
    srt_integrated <- RunUMAP(object = srt_integrated, reduction = "harmony", dims = 1:PC_use, n.components = 2, umap.method = "uwot")
  }
  if ("tsne" %in% reduction) {
    srt_integrated <- RunTSNE(
      object = srt_integrated, reduction = "harmony", dims = 1:PC_use, dim.embed = 2, tsne.method = "Rtsne",
      perplexity = max(ceiling(ncol(srt_integrated) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  if (length(cc_S_genes) >= 3 & length(cc_G2M_genes) >= 3) {
    srt_integrated <- CellCycleScoring(
      object = srt_integrated,
      s.features = cc_S_genes,
      g2m.features = cc_G2M_genes,
      set.ident = FALSE
    )
    srt_integrated[["CC.Difference"]] <- srt_integrated[["S.Score"]] - srt_integrated[["G2M.Score"]]
    srt_integrated[["Phase"]] <- factor(srt_integrated[["Phase", drop = TRUE]], levels = c("G1", "S", "G2M"))
  }

  DefaultAssay(srt_integrated) <- "RNA"
  return(srt_integrated)
}

Scanorama_integrate <- function(sc_list, HVF_source = "separate", nHVF = 3000,
                                maxPC = 100, resolution = 0.8, reduction = c("tsne", "umap"),
                                cc_S_genes = Seurat::cc.genes.updated.2019$s.genes, cc_G2M_genes = Seurat::cc.genes.updated.2019$g2m.genes,
                                exogenous_genes = NULL) {
  if (!HVF_source %in% c("global", "separate")) {
    stop("'HVF_source' must be one of: 'global','separate'",
      call. = FALSE
    )
  }

  require(reticulate)
  scanorama <- import("scanorama")

  for (i in 1:length(sc_list)) {
    DefaultAssay(sc_list[[i]]) <- "RNA"
    if (identical(
      x = GetAssayData(sc_list[[i]], slot = "counts"),
      y = GetAssayData(sc_list[[i]], slot = "data")
    )) {
      sc_list[[i]] <- NormalizeData(object = sc_list[[i]])
    }
    if (length(VariableFeatures(sc_list[[i]])) == 0) {
      sc_list[[i]] <- FindVariableFeatures(sc_list[[i]])
    }
    VariableFeatures(sc_list[[i]]) <- HVFInfo(sc_list[[i]]) %>%
      filter(variance.standardized > 1 &
        (!rownames(.) %in% exogenous_genes)) %>%
      dplyr::arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
  }
  if (HVF_source == "global") {
    # gene_common <- lapply(sc_list, rownames) %>% Reduce(intersect, .)
    gene_common <- lapply(sc_list, function(x) {
      m <- GetAssayData(x, assay = "RNA", slot = "counts")
      gene_keep <- rownames(m)[Matrix::rowSums(m > 0) >= 5]
      return(gene_keep)
    }) %>% Reduce(intersect, .)
    sc_merge <- Reduce(function(x, y) merge(x, y), sc_list)
    DefaultAssay(sc_merge) <- "RNA"
    hvf <- NormalizeData(object = sc_merge) %>%
      FindVariableFeatures(.) %>%
      HVFInfo(.) %>%
      filter(variance.standardized > 1 &
        (!rownames(.) %in% exogenous_genes) &
        rownames(.) %in% gene_common) %>%
      dplyr::arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
    sc_merge <- NULL
  }
  if (HVF_source == "separate") {
    hvf <- SelectIntegrationFeatures(object.list = sc_list, nfeatures = nHVF)
  }
  assaylist <- list()
  genelist <- list()
  for (i in 1:length(sc_list)) {
    assaylist[[i]] <- t(as.matrix(GetAssayData(sc_list[[i]], "data")))
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
  if (nrow(GetAssayData(srt_integrated, slot = "scale.data")) == 0) {
    srt_integrated <- ScaleData(object = srt_integrated, features = rownames(srt_integrated))
  }

  DefaultAssay(object = srt_integrated) <- "integrated"
  srt_integrated@project.name <- paste0(unique(srt_integrated[["orig.ident", drop = TRUE]]), collapse = ",")
  srt_integrated[["orig.ident"]] <- factor(srt_integrated[["orig.ident", drop = TRUE]],
    levels = unique(srt_integrated[["orig.ident", drop = TRUE]])
  )

  srt_integrated <- ScaleData(srt_integrated, features = hvf)
  srt_integrated <- RunPCA(object = srt_integrated, npcs = maxPC, features = hvf)
  PC_use <- ceiling(maxLikGlobalDimEst(data = Embeddings(srt_integrated, reduction = "pca"), k = 20, iterations = 100)[["dim.est"]])
  srt_integrated@misc$PC_use <- PC_use

  srt_integrated <- FindNeighbors(object = srt_integrated, reduction = "pca", dims = 1:PC_use, force.recalc = T)
  srt_integrated <- FindClusters(object = srt_integrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srt_integrated <- BuildClusterTree(srt_integrated, features = hvf, slot = "scale.data", reorder = T, reorder.numeric = T)
  srt_integrated$seurat_clusters <- Idents(srt_integrated)

  if ("umap" %in% reduction) {
    srt_integrated <- RunUMAP(object = srt_integrated, reduction = "pca", dims = 1:PC_use, n.components = 2, umap.method = "uwot")
  }
  if ("tsne" %in% reduction) {
    srt_integrated <- RunTSNE(
      object = srt_integrated, reduction = "pca", dims = 1:PC_use, dim.embed = 2, tsne.method = "Rtsne",
      perplexity = max(ceiling(ncol(srt_integrated) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  if (length(cc_S_genes) >= 3 & length(cc_G2M_genes) >= 3) {
    srt_integrated <- CellCycleScoring(
      object = srt_integrated,
      s.features = cc_S_genes,
      g2m.features = cc_G2M_genes,
      set.ident = FALSE
    )
    srt_integrated[["CC.Difference"]] <- srt_integrated[["S.Score"]] - srt_integrated[["G2M.Score"]]
    srt_integrated[["Phase"]] <- factor(srt_integrated[["Phase", drop = TRUE]], levels = c("G1", "S", "G2M"))
  }

  DefaultAssay(srt_integrated) <- "RNA"
  return(srt_integrated)
}

CSS_integrate <- function(sc_list, HVF_source = "separate", nHVF = 3000,
                          maxPC = 100, resolution = 0.8, reduction = c("tsne", "umap"),
                          cc_S_genes = Seurat::cc.genes.updated.2019$s.genes, cc_G2M_genes = Seurat::cc.genes.updated.2019$g2m.genes,
                          exogenous_genes = NULL) {
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
      sc_list[[i]] <- NormalizeData(object = sc_list[[i]])
    }
    if (length(VariableFeatures(sc_list[[i]])) == 0) {
      sc_list[[i]] <- FindVariableFeatures(sc_list[[i]])
    }
    VariableFeatures(sc_list[[i]]) <- HVFInfo(sc_list[[i]]) %>%
      filter(variance.standardized > 1 &
        (!rownames(.) %in% exogenous_genes)) %>%
      dplyr::arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
  }

  sc_merge <- Reduce(function(x, y) merge(x, y), sc_list)
  DefaultAssay(sc_merge) <- "RNA"
  sc_merge <- NormalizeData(object = sc_merge) %>%
    FindVariableFeatures(.)

  if (HVF_source == "global") {
    # gene_common <- lapply(sc_list, rownames) %>% Reduce(intersect, .)
    gene_common <- lapply(sc_list, function(x) {
      m <- GetAssayData(x, assay = "RNA", slot = "counts")
      gene_keep <- rownames(m)[Matrix::rowSums(m > 0) >= 5]
      return(gene_keep)
    }) %>% Reduce(intersect, .)
    hvf <- HVFInfo(sc_merge) %>%
      filter(variance.standardized > 1 &
        (!rownames(.) %in% exogenous_genes) &
        rownames(.) %in% gene_common) %>%
      dplyr::arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
  }
  if (HVF_source == "separate") {
    hvf <- SelectIntegrationFeatures(object.list = sc_list, nfeatures = nHVF)
  }

  sc_merge <- ScaleData(object = sc_merge, features = rownames(sc_merge))
  sc_merge <- RunPCA(object = sc_merge, npcs = maxPC, features = hvf)

  srt_integrated <- cluster_sim_spectrum(
    object = sc_merge,
    label_tag = "orig.ident",
    cluster_resolution = resolution,
    corr_method = "pearson",
    spectrum_type = "corr_kernel"
  )
  sc_merge <- NULL

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
  if (nrow(GetAssayData(srt_integrated, slot = "scale.data")) == 0) {
    srt_integrated <- ScaleData(object = srt_integrated, features = rownames(srt_integrated))
  }

  srt_integrated@project.name <- paste0(unique(srt_integrated[["orig.ident", drop = TRUE]]), collapse = ",")
  srt_integrated[["orig.ident"]] <- factor(srt_integrated[["orig.ident", drop = TRUE]],
    levels = unique(srt_integrated[["orig.ident", drop = TRUE]])
  )

  PC_use <- ncol(Embeddings(srt_integrated, reduction = "css"))
  srt_integrated@misc$PC_use <- PC_use

  srt_integrated <- FindNeighbors(object = srt_integrated, reduction = "css", dims = 1:PC_use, force.recalc = T)
  srt_integrated <- FindClusters(object = srt_integrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srt_integrated <- BuildClusterTree(srt_integrated, features = hvf, slot = "scale.data", reorder = T, reorder.numeric = T)
  srt_integrated$seurat_clusters <- Idents(srt_integrated)

  if ("umap" %in% reduction) {
    srt_integrated <- RunUMAP(object = srt_integrated, reduction = "css", dims = 1:PC_use, n.components = 2, umap.method = "uwot")
  }
  if ("tsne" %in% reduction) {
    srt_integrated <- RunTSNE(
      object = srt_integrated, reduction = "css", dims = 1:PC_use, dim.embed = 2, tsne.method = "Rtsne",
      perplexity = max(ceiling(ncol(srt_integrated) * 0.01), 30), max_iter = 2000, num_threads = 0, verbose = T
    )
  }

  if (length(cc_S_genes) >= 3 & length(cc_G2M_genes) >= 3) {
    srt_integrated <- CellCycleScoring(
      object = srt_integrated,
      s.features = cc_S_genes,
      g2m.features = cc_G2M_genes,
      set.ident = FALSE
    )
    srt_integrated[["CC.Difference"]] <- srt_integrated[["S.Score"]] - srt_integrated[["G2M.Score"]]
    srt_integrated[["Phase"]] <- factor(srt_integrated[["Phase", drop = TRUE]], levels = c("G1", "S", "G2M"))
  }

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
