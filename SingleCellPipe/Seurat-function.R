#' Run Non-negative matrix factorization(NMF)
#'
#' Run a NMF dimensionality reduction.
#'
#' @param object An object
#' @param ... Arguments passed to \code{\link[RcppML]{nmf}}.
#'
#' @return Returns Seurat object with the NMF calculation stored in the reductions slot.
#'
#' @export
#'
#' @rdname RunNMF
#' @export RunNMF
#'
RunNMF <- function(object, ...) {
  UseMethod(generic = "RunNMF", object = object)
}

#' @param assay Name of Assay nmf is being run on.
#' @param slot Name of slot nmf is being run on.
#' @param nbes Total Number of BEs("basis experiment", aka "metagene") to compute and store (50 by default).
#' @param rev.nmf By default computes the nmf on the cell x gene matrix. Setting
#' to true will compute it on gene x cell matrix.
#' @param verbose Print the top genes associated with high/low loadings for
#' the BEs.
#' @param ndims.print BEs to print genes for.
#' @param nfeatures.print Number of genes to print for each BE.
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. "BE_" by default.
#' @param seed.use Set a random seed. By default, sets the seed to 42. Setting
#' NULL will not set a seed.
#'
#' @importFrom RcppML nmf mse
#' @importFrom NMF nmf .basis .coef
#' @importFrom stats prcomp
#' @importFrom utils capture.output
#'
#' @rdname RunNMF
#' @concept dimensional_reduction
#' @export
#'
RunNMF.default <- function(object,
                           assay = NULL,
                           slot = "data",
                           nbes = 50,
                           nmf.method = "RcppML",
                           tol = 1e-5,
                           maxit = 100,
                           rev.nmf = FALSE,
                           verbose = TRUE,
                           ndims.print = 1:5,
                           nfeatures.print = 30,
                           reduction.key = "BE_",
                           seed.use = 42,
                           ...) {
  require("RcppML")
  require("NMF")
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (rev.nmf) {
    object <- Matrix::t(object)
  }
  nbes <- min(nbes, nrow(x = object) - 1)
  if (nmf.method == "RcppML") {
    nmf.results <- RcppML::nmf(
      A = t(object), k = nbes, tol = tol, maxit = maxit,
      seed = seed.use, verbose = verbose, ...
    )
    cell.embeddings <- nmf.results$w
    feature.loadings <- t(nmf.results$h)
    d <- nmf.results$d
    iter <- nmf.results$iter
    mse <- RcppML::mse(t(object), cell.embeddings, nmf.results$d, t(feature.loadings))
  }
  if (nmf.method == "NMF") {
    nmf.results <- NMF::nmf(x = as.matrix(t(object)), rank = nbes, seed = seed.use)
    cell.embeddings <- nmf.results@fit@W
    feature.loadings <- t(nmf.results@fit@H)
    d <- iter <- tol <- NULL
    mse <- RcppML::mse(t(object), cell.embeddings, NULL, t(feature.loadings))
  }

  rownames(x = feature.loadings) <- rownames(x = object)
  colnames(x = feature.loadings) <- paste0(reduction.key, 1:nbes)
  rownames(x = cell.embeddings) <- colnames(x = object)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = assay,
    key = reduction.key,
    misc = list(slot = slot, mse = mse, d = d, tol = tol, iter = iter)
  )
  if (verbose) {
    msg <- capture.output(print(
      x = reduction.data,
      dims = ndims.print,
      nfeatures = nfeatures.print
    ))
    message(paste(msg, collapse = "\n"))
  }
  return(reduction.data)
}

#' @param features Features to compute nmf on. If features=NULL, nmf will be run
#' using the variable features for the Assay.
#'
#' @rdname RunNMF
#' @concept dimensional_reduction
#' @export
#' @method RunNMF Assay
#'
RunNMF.Assay <- function(object,
                         assay = NULL,
                         slot = "data",
                         features = NULL,
                         nbes = 50,
                         nmf.method = "RcppML",
                         tol = 1e-5,
                         maxit = 100,
                         rev.nmf = FALSE,
                         verbose = TRUE,
                         ndims.print = 1:5,
                         nfeatures.print = 30,
                         reduction.key = "BE_",
                         seed.use = 42,
                         ...) {
  features <- features %||% VariableFeatures(object = object)
  data.use <- GetAssayData(object = object, slot = slot)
  features.var <- apply(
    X = data.use[features, ], MARGIN = 1,
    FUN = var
  )
  features.keep <- features[features.var > 0]
  data.use <- data.use[features.keep, ]
  reduction.data <- RunNMF(
    object = data.use,
    assay = assay,
    slot = slot,
    nbes = nbes,
    nmf.method = nmf.method,
    tol = tol,
    maxit = maxit,
    rev.nmf = rev.nmf,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  return(reduction.data)
}

#' @param reduction.name dimensional reduction name,  nmf by default
#'
#' @rdname RunNMF
#' @concept dimensional_reduction
#' @export
#' @method RunNMF Seurat
#'
RunNMF.Seurat <- function(object,
                          assay = NULL,
                          slot = "data",
                          features = NULL,
                          nbes = 50,
                          nmf.method = "RcppML",
                          tol = 1e-5,
                          maxit = 100,
                          rev.nmf = FALSE,
                          verbose = TRUE,
                          ndims.print = 1:5,
                          nfeatures.print = 30,
                          reduction.name = "nmf",
                          reduction.key = "BE_",
                          seed.use = 42,
                          ...) {
  features <- features %||% VariableFeatures(object = object)
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  reduction.data <- RunNMF(
    object = assay.data,
    assay = assay,
    slot = slot,
    features = features,
    nbes = nbes,
    nmf.method = nmf.method,
    tol = tol,
    maxit = maxit,
    rev.nmf = rev.nmf,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' Run multi-dimensional scaling(MDS)
#'
#' Run a MDS dimensionality reduction.
#'
#' @param object An object
#'
#' @return Returns Seurat object with the MDS calculation stored in the reductions slot.
#'
#' @export
#'
#' @rdname RunMDS
#' @export RunMDS
#'
RunMDS <- function(object, ...) {
  UseMethod(generic = "RunMDS", object = object)
}

#' @param assay Name of Assay mds is being run on.
#' @param slot Name of slot mds is being run on.
#' @param nmds Total Number of MDS to compute and store (50 by default).
#' @param rev.mds By default computes the mds on the cell x gene matrix. Setting
#' to true will compute it on gene x cell matrix.
#' @param verbose Print the top genes associated with high/low loadings for
#' the MDS.
#' @param ndims.print MDS to print genes for.
#' @param nfeatures.print Number of genes to print for each MDS.
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. "MDS_" by default.
#' @param seed.use Set a random seed. By default, sets the seed to 42. Setting
#' NULL will not set a seed.
#'
#' @importFrom parallelDist parDist
#' @importFrom stats prcomp
#' @importFrom utils capture.output
#'
#' @rdname RunMDS
#' @concept dimensional_reduction
#' @export
#'
RunMDS.default <- function(object,
                           assay = NULL,
                           slot = "scale.data",
                           nmds = 50,
                           dist.method = "cosine",
                           mds.method = "cmdscale",
                           rev.mds = FALSE,
                           ndims.print = 1:5,
                           nfeatures.print = 30,
                           reduction.key = "MDS_",
                           seed.use = 42,
                           verbose = TRUE,
                           ...) {
  require("parallelDist")
  require("MASS")
  message("Use ", assay, "/", slot, " to measure '", dist.method, "' distance.")
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (rev.mds) {
    object <- Matrix::t(object)
  }
  nmds <- min(nmds, nrow(x = object) - 1)
  x <- t(as.matrix(object))
  dist.method.keep <- dist.method
  if (dist.method %in% c("pearson", "spearman")) {
    if (dist.method == "spearman") {
      x <- t(apply(x, 1, rank))
    }
    x <- t(apply(x, 1, function(x) x - mean(x)))
    dist.method <- "cosine"
  }
  cell.dist <- parDist(x = x, method = dist.method)
  if (mds.method == "cmdscale") {
    mds.results <- cmdscale(cell.dist, k = nmds, eig = TRUE)
  }
  if (mds.method == "isoMDS") {
    mds.results2 <- isoMDS(cell.dist, k = nmds)
  }
  if (mds.method == "sammon") {
    mds.results3 <- sammon(cell.dist, k = nmds)
  }
  cell.embeddings <- mds.results$points

  rownames(x = cell.embeddings) <- colnames(x = object)
  colnames(x = cell.embeddings) <- paste0(reduction.key, 1:nmds)
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    assay = assay,
    key = reduction.key,
    misc = list(slot = slot, dist.method = dist.method.keep, mds.method = mds.method)
  )
  if (verbose) {
    msg <- capture.output(print(
      x = reduction.data,
      dims = ndims.print,
      nfeatures = nfeatures.print
    ))
    message(paste(msg, collapse = "\n"))
  }
  return(reduction.data)
}

#' @param features Features to compute mds on. If features=NULL, mds will be run
#' using the variable features for the Assay.
#'
#' @rdname RunMDS
#' @concept dimensional_reduction
#' @export
#' @method RunMDS Assay
#'
RunMDS.Assay <- function(object,
                         assay = NULL,
                         slot = "scale.data",
                         features = NULL,
                         nmds = 50,
                         dist.method = "cosine",
                         mds.method = "cmdscale",
                         rev.mds = FALSE,
                         verbose = TRUE,
                         ndims.print = 1:5,
                         nfeatures.print = 30,
                         reduction.key = "MDS_",
                         seed.use = 42,
                         ...) {
  features <- features %||% VariableFeatures(object = object)
  data.use <- GetAssayData(object = object, slot = slot)
  features.var <- apply(
    X = data.use[features, ], MARGIN = 1,
    FUN = var
  )
  features.keep <- features[features.var > 0]
  data.use <- data.use[features.keep, ]
  reduction.data <- RunMDS(
    object = data.use,
    assay = assay,
    slot = slot,
    nmds = nmds,
    dist.method = dist.method,
    mds.method = mds.method,
    rev.mds = rev.mds,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  return(reduction.data)
}

#' @param reduction.name dimensional reduction name, mds by default
#'
#' @rdname RunMDS
#' @concept dimensional_reduction
#' @export
#' @method RunMDS Seurat
#'
RunMDS.Seurat <- function(object,
                          assay = NULL,
                          slot = "scale.data",
                          features = NULL,
                          nmds = 50,
                          dist.method = "cosine",
                          mds.method = "cmdscale",
                          rev.mds = FALSE,
                          verbose = TRUE,
                          ndims.print = 1:5,
                          nfeatures.print = 30,
                          reduction.name = "mds",
                          reduction.key = "MDS_",
                          seed.use = 42,
                          ...) {
  features <- features %||% VariableFeatures(object = object)
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  reduction.data <- RunMDS(
    object = assay.data,
    assay = assay,
    slot = slot,
    features = features,
    nmds = nmds,
    dist.method = dist.method,
    mds.method = mds.method,
    rev.mds = rev.mds,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' Run DiffusionMap(DM)
#'
#' Run a DM dimensionality reduction.
#'
#' @param object An object
#' @param ... Arguments passed to \code{\link[diffusionMap]{diffuse}}.
#'
#' @return Returns Seurat object with the DM calculation stored in the reductions slot.
#'
#' @export
#'
#' @rdname RunDM
#' @export RunDM
#'
RunDM <- function(object, ...) {
  UseMethod(generic = "RunDM", object = object)
}

#' @param ndcs Total Number of DM to compute and store (2 by default).
#' @param assay Name of Assay dm is being run on.
#' @param slot Name of slot dm is being run on.
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. "DM_" by default.
#' @param seed.use Set a random seed. By default, sets the seed to 42. Setting
#' NULL will not set a seed.
#' @param verbose Show a progressbar and other progress information
#'
#' @importFrom parallelDist parDist
#' @importFrom stats prcomp
#' @importFrom utils capture.output
#'
#' @rdname RunDM
#' @concept dimensional_reduction
#' @export
#'
RunDM.dist <- function(object,
                       exp,
                       ndcs = 2,
                       assay = NULL,
                       slot = "scale.data",
                       reduction.key = "DM_",
                       seed.use = 42,
                       verbose = TRUE,
                       ...) {
  require("destiny")
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  dm.results <- DiffusionMap(distance = object, n_eigs = ndcs, verbose = verbose)
  # gene.relevance <- gene_relevance(coords=dm.results@eigenvectors, exprs=exp,verbose=verbose)
  cell.embeddings <- dm.results@eigenvectors
  rownames(x = cell.embeddings) <- attr(object, "Labels")
  colnames(x = cell.embeddings) <- paste0(reduction.key, 1:ndcs)
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    assay = assay,
    key = reduction.key,
    misc = list(slot = slot, dm.results = dm.results)
  )
  return(reduction.data)
}


#' @rdname RunDM
#' @concept dimensional_reduction
#' @export
#' @method RunDM Assay
#'
RunDM.matrix <- function(object,
                         ndcs = 2,
                         dist.method = "euclidean",
                         assay = NULL,
                         slot = "scale.data",
                         reduction.key = "DM_",
                         seed.use = 42,
                         verbose = TRUE,
                         ...) {
  require("parallelDist")
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  ndcs <- min(ndcs, nrow(x = object) - 1)
  x <- as.matrix(object)
  if (dist.method %in% c("pearson", "spearman")) {
    if (dist.method == "spearman") {
      x <- t(apply(x, 1, rank))
    }
    x <- t(apply(x, 1, function(x) x - mean(x)))
    dist.method <- "cosine"
  }
  cell.dist <- parDist(x = x, method = dist.method)
  reduction.data <- RunDM(
    object = cell.dist,
    exp = x,
    ndcs = ndcs,
    assay = assay,
    slot = slot,
    reduction.key = reduction.key,
    seed.use = seed.use,
    verbose = verbose,
    ...
  )
  return(reduction.data)
}

#' @param reduction.name dimensional reduction name, dm by default
#'
#' @rdname RunDM
#' @concept dimensional_reduction
#' @export
#' @method RunDM Seurat
#'
RunDM.Seurat <- function(object,
                         reduction = "pca",
                         dims = 1:30,
                         assay = NULL,
                         slot = "scale.data",
                         features = NULL,
                         ndcs = 2,
                         dist.method = "euclidean",
                         reduction.name = "dm",
                         reduction.key = "DM_",
                         seed.use = 42,
                         verbose = TRUE,
                         ...) {
  if (!is.null(x = features)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- as.matrix(x = t(x = GetAssayData(object = object,slot = slot, assay = assay)[features, , drop = FALSE]))
    if (ncol(x = data.use) < ndcs) {
      stop("Please provide as many or more features than ndcs: ",
        length(x = features), " features provided, ",
        ndcs, " Diffusion components requested",
        call. = FALSE
      )
    }
  } else if (!is.null(x = dims)) {
    data.use <- Embeddings(object[[reduction]])[, dims]
    assay <- DefaultAssay(object = object[[reduction]])
    if (length(x = dims) < ndcs) {
      stop("Please provide as many or more dims than ndcs: ",
        length(x = dims), " dims provided, ", ndcs,
        " DiffusionMap components requested",
        call. = FALSE
      )
    }
  } else {
    stop("Please specify one of dims, features")
  }
  reduction.data <- RunDM(
    object = data.use,
    reduction = reduction,
    dims = dims,
    assay = assay,
    slot = slot,
    features = features,
    ndcs = ndcs,
    dist.method = dist.method,
    reduction.key = reduction.key,
    seed.use = seed.use,
    verbose = verbose,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object <- LogSeuratCommand(object = object)
  return(object)
}
