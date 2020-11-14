#' Read bismark coverage file as a methylKit object
#'
#' Bismark aligner can output methylation information per base in
#' multiple different formats. This function reads coverage files,
#' which have chr,start,end, number of cytosines (methylated bases)
#' and number of thymines (unmethylated bases).
#'
#' @param location a list or vector of file paths to coverage files
#'
#' @param sample.id a list or vector of sample ids
#' @param assembly a string for genome assembly. Any string would work.
#' @param treatment if there are multiple files to be read a treatment
#'                  vector should be supplied.
#' @param context a string for context of methylation such as: "CpG" or "CHG"
#' @param min.cov a numeric value for minimum coverage. Bases that have coverage
#' below this value will be removed.
#'
#' @return methylRaw or methylRawList objects
readBismarkCoverage <- function(location, sample.id, assembly = "unknown", treatment,
                                context = "CpG", min.cov = 10) {
  if (length(location) > 1) {
    stopifnot(
      length(location) == length(sample.id),
      length(location) == length(treatment)
    )
  }

  result <- list()
  for (i in 1:length(location)) {
    cat(paste0("[",i,"/",length(location),"]"," ++++++"), sample.id[i], "++++++", "\n")
    df <- fread(location[[i]], data.table = FALSE)

    # remove low coverage stuff
    df <- df[(df[, 5] + df[, 6]) >= min.cov, ]

    # make the object (arrange columns of df), put it in a list
    result[[i]] <- new("methylRaw", data.frame(
      chr = df[, 1], start = df[, 2], end = df[, 3],
      strand = "*", coverage = (df[, 5] + df[, 6]),
      numCs = df[, 5], numTs = df[, 6]
    ),
    sample.id = sample.id[[i]],
    assembly = assembly, context = context, resolution = "base"
    )
  }

  if (length(result) == 1) {
    return(result[[1]])
  } else {
    new("methylRawList", result, treatment = treatment)
  }
}

#' Read bismark cytosine report file as a methylKit object
#'
#' Bismark aligner can output methylation information per base in
#' multiple different formats. This function reads cytosine report files,
#' which have chr,start, strand, number of cytosines (methylated bases)
#' and number of thymines (unmethylated bases),context, trinucletide context.
#'
#' @param location a list or vector of file paths to coverage files
#'
#' @param sample.id a list or vector of sample ids
#' @param assembly a string for genome assembly. Any string would work.
#' @param treatment if there are multiple files to be read a treatment
#'                  vector should be supplied.
#' @param context a string for context of methylation such as: "CpG" or "CHG"
#' @param min.cov a numeric value for minimum coverage. Bases that have coverage
#' below this value will be removed.
#'
#' @return methylRaw or methylRawList objects

readBismarkCytosineReport <- function(location, sample.id, assembly = "unknown", treatment,
                                      context = "CpG", min.cov = 10) {
  if (length(location) > 1) {
    stopifnot(
      length(location) == length(sample.id),
      length(location) == length(treatment)
    )
  }

  result <- list()
  for (i in 1:length(location)) {
    cat(paste0("[",i,"/",length(location),"]"," ++++++"), sample.id[i], "++++++", "\n")
    df <- fread(location[[i]], data.table = FALSE)

    # remove low coverage stuff
    df <- df[(df[, 4] + df[, 5]) >= min.cov, ]

    # make the object (arrange columns of df), put it in a list
    result[[i]] <- new("methylRaw",
      data.frame(
        chr = df[, 1], start = df[, 2], end = df[, 2],
        strand = df[, 3], coverage = (df[, 4] + df[, 5]),
        numCs = df[, 4], numTs = df[, 5]
      ),
      sample.id = sample.id[[i]],
      assembly = assembly, context = context, resolution = "base"
    )
  }

  if (length(result) == 1) {
    return(result[[1]])
  } else {
    new("methylRawList", result, treatment = treatment)
  }
}


# ---------------------------------------------------------------------------- #
#' Run limma for differential methylation calling
#'
#' @param sim.methylBase a methylBase object from the methylKit library
#' @param transform a boolean indicating whether percent of methylation values
#'                  will be converted using a logistic transformation (Default: TRUE)
#' @return returns a methylDiff object from the methylKit library
limma.meth <- function(methylBase.obj, transform = TRUE) {
  require(limma)
  require(qvalue)

  group <- methylBase.obj@treatment
  design <- model.matrix(~group)

  # do the test in limma
  p.meth <- percMethylation(methylBase.obj) # get percent methylation values
  if (transform) {
    p.meth2 <- log2((p.meth + 1) / (100 - p.meth + 1))
    fit <- lmFit(p.meth2, design = design)
    fit2 <- eBayes(fit)
  } else {
    fit <- lmFit(p.meth, design = design)
    fit2 <- eBayes(fit)
  }
  # make the data for methylKit object
  df <- cbind(getData(methylBase.obj)[, 1:4],
    pvalue = fit2$p.value[, 2],
    qvalue = qvalue(fit2$p.value[, 2])$qvalues,
    meth.diff = rowMeans(p.meth[, methylBase.obj@treatment == 1]) - rowMeans(p.meth[, methylBase.obj@treatment == 0])
  )

  # create a new methylDiff object
  obj <- new("methylDiff", df,
    sample.ids = methylBase.obj@sample.ids,
    assembly = methylBase.obj@assembly,
    context = methylBase.obj@context, treatment = methylBase.obj@treatment,
    destranded = methylBase.obj@destranded, resolution = methylBase.obj@resolution
  )
}



# ---------------------------------------------------------------------------- #
# DSS functions

run.DSS <- function(methylBase.obj, difference = 5, cores = 1) {
  print(cores)

  dss.qvalue <- calculateDiffMethDSS(methylBase.obj,
    adjust = "qvalue",
    mc.cores = cores
  )
  dss.fdr <- calculateDiffMethDSS(methylBase.obj,
    adjust = "fdr",
    mc.cores = cores
  )
  dss.slim <- calculateDiffMethDSS(methylBase.obj,
    adjust = "SLIM",
    mc.cores = cores
  )
  mylist.dss <- list(dss.qvalue, dss.fdr, dss.slim)
  names(mylist.dss) <- c("qvalue", "fdr", "slim")

  mylist.dss.diff.gr <- mclapply(mylist.dss,
    function(x) {
      as(getMethylDiff(x, difference = difference), "GRanges")
    },
    mc.cores = 3
  )
  names(mylist.dss.diff.gr) <- c("qvalue", "fdr", "slim")

  return(list(mylist.dss, mylist.dss.diff.gr))
}

# ---------------------------------------------------------------------------- #
# methylKit functions

run.methylkit <- function(sim.methylBase, cores = 1,
                          difference = 5) {
  require("methylKit")
  adjusts <- c("SLIM", "fdr", "qvalue")
  overdispersion <- c(
    "none",
    "MN",
    "shrinkMN"
  )
  tests <- c("F", "Chisq")

  combined <- expand.grid(
    test = tests,
    adjust = adjusts,
    overd = overdispersion,
    stringsAsFactors = FALSE
  )
  combined <- cbind(combined,
    name = with(combined, paste("methylKit",
      test,
      adjust,
      overd,
      sep = "."
    ))
  )
  combined$name <- as.character(combined$name)

  methylKit.list <- list()
  for (i in 1:nrow(combined)) {
    co <- combined[i, ]
    methylkit.obj <- calculateDiffMeth(sim.methylBase,
      overdispersion = co$overd,
      adjust = co$adjust,
      test = co$test,
      mc.cores = cores
    )
    methylkit.obj.diff <- getMethylDiff(methylkit.obj, difference = difference)
    methylKit.list[[i]] <- as(methylkit.obj.diff, "GRanges")
  }
  names(methylKit.list) <- combined$name

  return(methylKit.list)
}

# ---------------------------------------------------------------------------- #
# BSmooth functions

smoothed.dmrFinder <- function(smoothed) {
  x <- smoothed

  if (class(x@stats[, "tstat.corrected"]) == "character") {
    x@stats <- apply(x@stats, 2, as.numeric)
  }
  if (any(is.na(x@stats))) {
    # b = apply(x@stats, 1, function(y) any(is.na(y)))
    b <- is.na(x@stats[, "tstat.corrected"])
    x@stats <- x@stats[-which(b), ]
  }

  bsmooth.cutoff4.5 <- function.error(try(my_dmrFinder(x,
    cutoff = c(-4.5, 4.5),
    qcutoff = NULL
  )))
  bsmooth.cutoff4.5 <- subset(
    bsmooth.cutoff4.5,
    n >= 3 & abs(meanDiff) >= 0.1
  )

  bsmooth.qcutoff1 <- function.error(try(my_dmrFinder(x,
    cutoff = NULL,
    qcutoff = c(0.01, 0.99)
  )))
  bsmooth.qcutoff1 <- subset(
    bsmooth.qcutoff1,
    n >= 3 & abs(meanDiff) >= 0.1
  )

  bsmooth.qcutoff10 <- function.error(try(my_dmrFinder(x,
    cutoff = NULL,
    qcutoff = c(0.1, 0.9)
  )))
  bsmooth.qcutoff10 <- subset(
    bsmooth.qcutoff10,
    n >= 3 & abs(meanDiff) >= 0.1
  )

  mylist <- list(bsmooth.cutoff4.5, bsmooth.qcutoff1, bsmooth.qcutoff10)
  names(mylist) <- c("bsmooth.cutoff4.5", "bsmooth.qcutoff1", "bsmooth.qcutoff10")
  return(mylist)
}

function.error <- function(dmrs0) {
  if (class(dmrs0) == "try-error") {
    print(dmrs0)
    return(NA)
  }
  if (is.null(dmrs0)) {
    return(GRanges())
  }
  return(dmrs0)
}

smoothed.dmrFinder.paper <- function(smoothed) {
  x <- smoothed

  if (class(x@stats[, "tstat.corrected"]) == "character") {
    x@stats <- apply(x@stats, 2, as.numeric)
  }
  if (any(is.na(x@stats))) {
    # b = apply(x@stats, 1, function(y) any(is.na(y)))
    b <- is.na(x@stats[, "tstat.corrected"])
    x@stats <- x@stats[-which(b), ]
  }

  bsmooth.qcutoff <- function.error(try(my_dmrFinder(x,
    cutoff = c(0.025, 0.975),
    qcutoff = NULL
  )))
  bsmooth.qcutoff.default_sub <- subset(
    bsmooth.qcutoff,
    n >= 3 & abs(meanDiff) >= 0.1
  )

  list(bsmooth.same.default = bsmooth.qcutoff.default_sub)
}


run.BSmooth <- function(methylBase.obj.list, cores = 1,
                        difference = 5) {
  ctcf.bsmooth.smooth.same <- lapply(methylBase.obj.list, function(x) {
    methylBase.obj <- rm.small.chr.methylBase(x)
    runBSmooth.smooth(methylBase.obj,
      estimate.var = "same",
      cores = 60
    )
  })
  ctcf.bsmooth.smooth.group2 <- lapply(methylBase.obj.list, function(x) {
    methylBase.obj <- rm.small.chr.methylBase(x)
    runBSmooth.smooth(methylBase.obj,
      estimate.var = "group2",
      cores = 60
    )
  })
}



require("methylKit")
require("bsseq") #BSmooth


convert.methylBase2BSseq.obj = function(methylBase.obj){
  require(methylKit)
  require(bsseq)
  
  # assumtion: there is the same number of treated and non treated samples
  nbr.samples = length(methylBase.obj@treatment)
  n.t = nbr.samples/2
  
  methylBase.d = getData(methylBase.obj)
  Covs <- as.matrix(methylBase.d[, methylBase.obj@coverage.index ])
  Meths <- as.matrix(methylBase.d[, methylBase.obj@numCs.index ])
  colData = data.frame(Type=methylBase.obj@treatment,
                       Pair=c(paste0("sample", 1:n.t),paste0("sample", 1:n.t)))
  rownames(colData) <- c(paste0(rep("C",n.t), 1:n.t), paste0(rep("T", n.t), 1:n.t))
  BS.obj <- BSseq(chr = methylBase.d$chr, pos = methylBase.d$start,
                  M = Meths, Cov = Covs, sampleNames = rownames(colData),
                  pData=colData)
  return(BS.obj)
}

runBSmooth.smooth = function(methylBase.obj, 
                             estimate.var="group2",
                             cores=10){
  
  nmbr.treat = sum(methylBase.obj@treatment==1)
  nmbr.notreat = sum(methylBase.obj@treatment==0)
  BSseq.obj = convert.methylBase2BSseq.obj(methylBase.obj)
  
  BS1.fit <- BSmooth(BSseq.obj, 
                     verbose = TRUE)
  BS1.tstat <- try(BSmooth.tstat(BS1.fit,
                                 group1 = paste0(rep("T", nmbr.treat), 1:nmbr.treat),
                                 group2 = paste0(rep("C",nmbr.notreat), 1:nmbr.notreat),
                                 estimate.var = estimate.var,
                                 local.correct = TRUE,
                                 verbose = TRUE,
                                 mc.cores=cores))
  if(class( BS1.tstat)=="try-error"){
    print(BS1.tstat)
    return(NA)
  }
  return(BS1.tstat)
}


runBSmooth = function(methylBase.obj, 
                      estimate.var="group2",
                      cutoff=c(-4.5, 4.5),
                      qcutoff=NULL,
                      cores=10){
  
  nmbr.treat = sum(methylBase.obj@treatment==1)
  nmbr.notreat = sum(methylBase.obj@treatment==0)
  BSseq.obj = convert.methylBase2BSseq.obj(methylBase.obj)
  
  BS1.fit <- BSmooth(BSseq.obj, 
                     verbose = TRUE)
  
  BS1.tstat <- try(BSmooth.tstat(BS1.fit,
                                 group1 = paste0(rep("T", nmbr.treat), 1:nmbr.treat),
                                 group2 = paste0(rep("C",nmbr.notreat), 1:nmbr.notreat),
                                 estimate.var = estimate.var,
                                 local.correct = TRUE,
                                 verbose = TRUE))
  if(class( BS1.tstat)=="try-error"){
    print(BS1.tstat)
    return(NA)
  }
  
  dmrs0 <- try(dmrFinder(BS1.tstat, 
                         cutoff = cutoff,
                         qcutoff = qcutoff)) #cutoff using these quantiles of the t-statistic.
  # Error in quantile.default(dmrStat, qcutoff) : 
  #   missing values and NaN's not allowed if 'na.rm' is FALSE
  if(class( dmrs0)=="try-error"){
    print(dmrs0)
    return(NA)
  }
  if(is.null(dmrs0)){
    return(GRanges())
  }
  
  dmrs.gr = makeGRangesFromDataFrame(dmrs0,
                                     keep.extra.columns=TRUE)
  return(dmrs.gr)
}

