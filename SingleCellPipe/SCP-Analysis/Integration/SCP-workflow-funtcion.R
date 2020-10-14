Standard_SCP <- function(sc, nHVF = 3000, maxPC = 100, resolution = 0.8,
                         cc_S_genes = Seurat::cc.genes.updated.2019$s.genes, cc_G2M_genes = Seurat::cc.genes.updated.2019$g2m.genes,
                         exogenous_genes = NULL, assay = "RNA") {
  DefaultAssay(sc) <- assay
  if (identical(
    x = GetAssayData(sc, slot = "counts"),
    y = GetAssayData(sc, slot = "data")
  )) {
    sc <- NormalizeData(object = sc)
  }
  if (length(VariableFeatures(sc)) == 0) {
    sc <- FindVariableFeatures(sc)
  }
  VariableFeatures(sc) <-
    HVFInfo(sc) %>%
    filter(variance.standardized > 1 &
      (!rownames(.) %in% exogenous_genes)) %>%
    arrange(desc(variance.standardized)) %>%
    rownames(.) %>%
    head(n = nHVF)
  if (nrow(GetAssayData(sc, slot = "scale.data")) == 0) {
    sc <- ScaleData(object = sc, features = rownames(sc))
  }

  sc <- RunPCA(object = sc, npcs = maxPC, features = VariableFeatures(sc))
  PC_use <- ceiling(maxLikGlobalDimEst(data = Embeddings(sc, reduction = "pca"), k = 20, iterations = 100)[["dim.est"]])
  sc <- RunUMAP(object = sc, reduction = "pca", dims = 1:PC_use, n.components = 3, umap.method = "uwot")
  sc <- FindNeighbors(object = sc, reduction = "pca", dims = 1:PC_use, force.recalc = T)
  sc <- FindClusters(object = sc, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  sc <- BuildClusterTree(sc, reorder = T)
  Idents(sc) <- sc[["seurat_clusters"]] <-
    mapvalues(Idents(sc),
      from = levels(Idents(sc)),
      to = 1:length(levels(Idents(sc)))
    )
  sc@tools$BuildClusterTree$tip.label <-
    mapvalues(sc@tools$BuildClusterTree$tip.label,
      from = sc@tools$BuildClusterTree$tip.label,
      to = 1:length(levels(Idents(sc)))
    )

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
  Markers_MAST <- FindAllMarkers(
    object = sc, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25,
    test.use = "MAST", latent.vars = "orig.ident"
  )
  Markers_ROC <- FindAllMarkers(
    object = sc, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25,
    test.use = "roc"
  )
  sc@tools$FindAllMarkers <- setNames(
    object = list(Markers_MAST, Markers_ROC),
    nm = c("Markers_MAST", "Markers_ROC")
  )

  return(sc)
}

SCTransform_SCP <- function(sc, nHVF = 3000, maxPC = 100, resolution = 0.8,
                            cc_S_genes = Seurat::cc.genes.updated.2019$s.genes, cc_G2M_genes = Seurat::cc.genes.updated.2019$g2m.genes,
                            exogenous_genes = NULL, assay = "RNA") {
  DefaultAssay(sc) <- assay
  if (!"SCT" %in% Seurat::Assays(sc)) {
    sc <- SCTransform(
      object = sc,
      variable.features.n = nHVF,
      return.only.var.genes = FALSE,
    )
  }
  DefaultAssay(sc) <- "SCT"

  VariableFeatures(sc) <- HVFInfo(sc) %>%
    filter((!rownames(.) %in% exogenous_genes)) %>%
    arrange(desc(residual_variance)) %>%
    rownames(.) %>%
    head(n = nHVF)
  sc <- RunPCA(object = sc, npcs = maxPC, features = VariableFeatures(sc))
  PC_use <- ceiling(maxLikGlobalDimEst(data = Embeddings(sc, reduction = "pca"), k = 20, iterations = 100)[["dim.est"]])
  sc <- RunUMAP(object = sc, reduction = "pca", dims = 1:PC_use, n.components = 3, umap.method = "uwot")
  sc <- FindNeighbors(object = sc, reduction = "pca", dims = 1:PC_use, force.recalc = T)
  sc <- FindClusters(object = sc, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  sc <- BuildClusterTree(sc, reorder = T)
  Idents(sc) <- sc[["seurat_clusters"]] <-
    mapvalues(Idents(sc),
      from = levels(Idents(sc)),
      to = 1:length(levels(Idents(sc)))
    )
  sc@tools$BuildClusterTree$tip.label <-
    mapvalues(sc@tools$BuildClusterTree$tip.label,
      from = sc@tools$BuildClusterTree$tip.label,
      to = 1:length(levels(Idents(sc)))
    )

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
  Markers_MAST <- FindAllMarkers(
    object = sc, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25,
    test.use = "MAST", latent.vars = "orig.ident"
  )
  Markers_ROC <- FindAllMarkers(
    object = sc, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25,
    test.use = "roc"
  )
  sc@tools$FindAllMarkers <- setNames(
    object = list(Markers_MAST, Markers_ROC),
    nm = c("Markers_MAST", "Markers_ROC")
  )

  return(sc)
}

Standard_integrate <- function(sc_list, nHVF = 3000, anchor_dims = 1:30, integrate_dims = 1:30, maxPC = 100, resolution = 0.8,
                               HVF_source = "separate",
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
      arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
  }
  if (HVF_source == "global") {
    gene_common <- lapply(sc_list, rownames) %>% Reduce(intersect, .)
    sc_merge <- Reduce(function(x, y) merge(x, y), sc_list)
    hvf <- NormalizeData(object = sc_merge) %>%
      FindVariableFeatures(.) %>%
      HVFInfo(.) %>%
      filter(variance.standardized > 1 &
        (!rownames(.) %in% exogenous_genes) &
        rownames(.) %in% gene_common) %>%
      arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
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
  srt_integrated <- RunUMAP(object = srt_integrated, reduction = "pca", dims = 1:PC_use, n.components = 3, umap.method = "uwot")
  srt_integrated <- FindNeighbors(object = srt_integrated, reduction = "pca", dims = 1:PC_use, force.recalc = T)
  srt_integrated <- FindClusters(object = srt_integrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srt_integrated <- BuildClusterTree(srt_integrated, reorder = T)
  Idents(srt_integrated) <- srt_integrated[["seurat_clusters"]] <-
    mapvalues(Idents(srt_integrated),
      from = levels(Idents(srt_integrated)),
      to = 1:length(levels(Idents(srt_integrated)))
    )
  srt_integrated@tools$BuildClusterTree$tip.label <-
    mapvalues(srt_integrated@tools$BuildClusterTree$tip.label,
      from = srt_integrated@tools$BuildClusterTree$tip.label,
      to = 1:length(levels(Idents(srt_integrated)))
    )

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
  Markers_MAST <- FindAllMarkers(
    object = srt_integrated, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25,
    test.use = "MAST", latent.vars = "orig.ident"
  )
  Markers_ROC <- FindAllMarkers(
    object = srt_integrated, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25,
    test.use = "roc"
  )
  srt_integrated@tools$FindAllMarkers <- setNames(
    object = list(Markers_MAST, Markers_ROC),
    nm = c("Markers_MAST", "Markers_ROC")
  )

  return(srt_integrated)
}

SCTransform_integrate <- function(sc_list, nHVF = 3000, anchor_dims = 1:30, integrate_dims = 1:30, maxPC = 100, resolution = 0.8,
                                  HVF_source = "separate",
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
        return.only.var.genes = FALSE,
        assay = "RNA"
      )
    }
    DefaultAssay(sc_list[[i]]) <- "SCT"
    VariableFeatures(sc_list[[i]]) <- HVFInfo(sc_list[[i]]) %>%
      filter((!rownames(.) %in% exogenous_genes)) %>%
      arrange(desc(residual_variance)) %>%
      rownames(.) %>%
      head(n = nHVF)
  }
  if (HVF_source == "global") {
    gene_common <- lapply(sc_list, rownames) %>% Reduce(intersect, .)
    sc_merge <- Reduce(function(x, y) merge(x, y), sc_list)
    hvf <- NormalizeData(object = sc_merge) %>%
      FindVariableFeatures(.) %>%
      HVFInfo(.) %>%
      filter((!rownames(.) %in% exogenous_genes) &
        rownames(.) %in% gene_common) %>%
      arrange(desc(residual_variance)) %>%
      rownames(.) %>%
      head(n = nHVF)
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
  srt_integrated <- RunUMAP(object = srt_integrated, reduction = "pca", dims = 1:PC_use, n.components = 3, umap.method = "uwot")
  srt_integrated <- FindNeighbors(object = srt_integrated, reduction = "pca", dims = 1:PC_use, force.recalc = T)
  srt_integrated <- FindClusters(object = srt_integrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srt_integrated <- BuildClusterTree(srt_integrated, reorder = T)
  Idents(srt_integrated) <- srt_integrated[["seurat_clusters"]] <-
    mapvalues(Idents(srt_integrated),
      from = levels(Idents(srt_integrated)),
      to = 1:length(levels(Idents(srt_integrated)))
    )
  srt_integrated@tools$BuildClusterTree$tip.label <-
    mapvalues(srt_integrated@tools$BuildClusterTree$tip.label,
      from = srt_integrated@tools$BuildClusterTree$tip.label,
      to = 1:length(levels(Idents(srt_integrated)))
    )

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
  Markers_MAST <- FindAllMarkers(
    object = srt_integrated, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25,
    test.use = "MAST", latent.vars = "orig.ident"
  )
  Markers_ROC <- FindAllMarkers(
    object = srt_integrated, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25,
    test.use = "roc"
  )
  srt_integrated@tools$FindAllMarkers <- setNames(
    object = list(Markers_MAST, Markers_ROC),
    nm = c("Markers_MAST", "Markers_ROC")
  )

  return(srt_integrated)
}

fastMNN_integrate <- function(sc_list, nHVF = 3000, maxPC = 100, resolution = 0.8,
                              HVF_source = "separate",
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
      arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
  }
  if (HVF_source == "global") {
    gene_common <- lapply(sc_list, rownames) %>% Reduce(intersect, .)
    sc_merge <- Reduce(function(x, y) merge(x, y), sc_list)
    hvf <- NormalizeData(object = sc_merge) %>%
      FindVariableFeatures(.) %>%
      HVFInfo(.) %>%
      filter(variance.standardized > 1 &
        (!rownames(.) %in% exogenous_genes) &
        rownames(.) %in% gene_common) %>%
      arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
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
  srt_integrated <- RunUMAP(object = srt_integrated, reduction = "mnn", dims = 1:PC_use, n.components = 3, umap.method = "uwot")
  srt_integrated <- FindNeighbors(object = srt_integrated, reduction = "mnn", dims = 1:PC_use, force.recalc = T)
  srt_integrated <- FindClusters(object = srt_integrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srt_integrated <- BuildClusterTree(srt_integrated, reorder = T)
  Idents(srt_integrated) <- srt_integrated[["seurat_clusters"]] <-
    mapvalues(Idents(srt_integrated),
      from = levels(Idents(srt_integrated)),
      to = 1:length(levels(Idents(srt_integrated)))
    )
  srt_integrated@tools$BuildClusterTree$tip.label <-
    mapvalues(srt_integrated@tools$BuildClusterTree$tip.label,
      from = srt_integrated@tools$BuildClusterTree$tip.label,
      to = 1:length(levels(Idents(srt_integrated)))
    )

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
  Markers_MAST <- FindAllMarkers(
    object = srt_integrated, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25,
    test.use = "MAST", latent.vars = "orig.ident"
  )
  Markers_ROC <- FindAllMarkers(
    object = srt_integrated, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25,
    test.use = "roc"
  )
  srt_integrated@tools$FindAllMarkers <- setNames(
    object = list(Markers_MAST, Markers_ROC),
    nm = c("Markers_MAST", "Markers_ROC")
  )

  return(srt_integrated)
}

Harmony_integrate <- function(sc_list, nHVF = 3000, maxPC = 100, resolution = 0.8,
                              HVF_source = "separate",
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
      arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
  }

  sc_merge <- Reduce(function(x, y) merge(x, y), sc_list)
  sc_merge <- NormalizeData(object = sc_merge)
  sc_merge <- FindVariableFeatures(object = sc_merge)
  sc_merge <- ScaleData(object = sc_merge, features = rownames(sc_merge))

  if (HVF_source == "global") {
    gene_common <- lapply(sc_list, rownames) %>% Reduce(intersect, .)
    hvf <- HVFInfo(sc_merge) %>%
      filter(variance.standardized > 1 &
        (!rownames(.) %in% exogenous_genes) &
        rownames(.) %in% gene_common) %>%
      arrange(desc(variance.standardized)) %>%
      rownames(.) %>%
      head(n = nHVF)
  }
  if (HVF_source == "separate") {
    hvf <- SelectIntegrationFeatures(object.list = sc_list, nfeatures = nHVF)
  }

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
  srt_integrated <- RunUMAP(object = srt_integrated, reduction = "harmony", dims = 1:PC_use, n.components = 3, umap.method = "uwot")
  srt_integrated <- FindNeighbors(object = srt_integrated, reduction = "harmony", dims = 1:PC_use, force.recalc = T)
  srt_integrated <- FindClusters(object = srt_integrated, resolution = resolution, algorithm = 1, n.start = 100, n.iter = 10000)
  srt_integrated <- BuildClusterTree(srt_integrated, reorder = T)
  Idents(srt_integrated) <- srt_integrated[["seurat_clusters"]] <-
    mapvalues(Idents(srt_integrated),
      from = levels(Idents(srt_integrated)),
      to = 1:length(levels(Idents(srt_integrated)))
    )
  srt_integrated@tools$BuildClusterTree$tip.label <-
    mapvalues(srt_integrated@tools$BuildClusterTree$tip.label,
      from = srt_integrated@tools$BuildClusterTree$tip.label,
      to = 1:length(levels(Idents(srt_integrated)))
    )
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
  Markers_MAST <- FindAllMarkers(
    object = srt_integrated, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25,
    test.use = "MAST", latent.vars = "orig.ident"
  )
  Markers_ROC <- FindAllMarkers(
    object = srt_integrated, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25,
    test.use = "roc"
  )
  srt_integrated@tools$FindAllMarkers <- setNames(
    object = list(Markers_MAST, Markers_ROC),
    nm = c("Markers_MAST", "Markers_ROC")
  )

  return(srt_integrated)
}
