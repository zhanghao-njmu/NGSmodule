findmarkergenes <- function(object, species = NULL, cluster = "All", ref = NULL, match_ref = FALSE,
                            cancer = NULL, tissue = NULL, cell_min_pct = 0.25, logfc = 0.25,
                            pvalue = 0.05, BPPARAM = SerialParam()) {
  if (is.null(ref)) {
    stop("ref must be specified.")
  }
  ndata <- object[["RNA"]]@data
  clu_info <- Seurat::Idents(object = object)
  clu_info <- data.frame(
    cell = names(clu_info), cluster = as.character(clu_info),
    stringsAsFactors = F
  )
  clu_info[, 1] <- as.character(clu_info[, 1])
  clu_info[, 2] <- as.character(clu_info[, 2])
  clu_num <- unique(clu_info[, 2])
  clu_num1 <- clu_num[nchar(clu_num) == 1]
  clu_num2 <- clu_num[nchar(clu_num) == 2]
  clu_num1 <- clu_num1[order(clu_num1)]
  clu_num2 <- clu_num2[order(clu_num2)]
  clu_num <- c(clu_num1, clu_num2)
  rm(object)
  cat(
    "Note: the raw data matrix includes", ncol(ndata), "cells and",
    nrow(ndata), "genes.", "\n"
  )
  if (length(clu_num) == 1) {
    stop("There is only one cluster, please do clustering analysis or select more clusters for Seurat object!")
  }
  Sys.sleep(2)
  clu_num1 <- clu_num
  if (is.null(species)) {
    cat("\n")
    stop("Please define the species! 'Human' or 'Mouse'.")
  }
  geneinfo <- geneinfo[geneinfo$species == species, ]
  cat("\n")
  cat(
    "---Revising gene symbols according to NCBI Gene symbols (updated in June 19, 2020, https://www.ncbi.nlm.nih.gov/gene) and no matched genes and duplicated genes will be removed.",
    "\n"
  )
  Sys.sleep(2)
  genename <- rownames(ndata)
  genename1 <- genename[genename %in% geneinfo$Symbol]
  genename2 <- genename[!genename %in% geneinfo$Symbol]
  genename3 <- genename2[genename2 %in% geneinfo$Synonyms]
  genename4 <- rep("NA", length(genename3))
  for (i in 1:length(genename3)) {
    d1 <- geneinfo[geneinfo$Synonyms == genename3[i], ]$Symbol
    if (length(d1) == 1) {
      genename4[i] <- d1
    }
  }
  genename3 <- c(genename1, genename3)
  genename4 <- c(genename1, genename4)
  genedata <- data.frame(
    raw_name = genename3, new_name = genename4,
    stringsAsFactors = F
  )
  genedata <- genedata[!genedata$new_name == "NA", ]
  genedata1 <- as.data.frame(table(genedata$new_name), stringsAsFactors = F)
  genedata1 <- genedata1[genedata1$Freq == 1, ]
  genedata <- genedata[genedata$new_name %in% genedata1$Var1, ]
  ndata <- ndata[genedata$raw_name, ]
  rownames(ndata) <- genedata$new_name
  cat("\n")
  cat(
    "Note: the new data matrix includes", ncol(ndata), "cells and",
    nrow(ndata), "genes.", "\n"
  )
  Sys.sleep(2)
  if (match_ref) {
    ndata3 <- ndata
    if (is.null(cancer)) {
      ref <- ref[ref$speciesType ==
        species & ref$cancerType == "Normal", ]
      if (is.null(tissue)) {
        cat("\n")
        stop("Please define the origin tissue of cells! Select one or more related tissue types.")
      }
      if (tissue[1] == "All") {
        tissue <- unique(ref$tissueType)
      }
      if (!all(tissue %in% ref$tissueType)) {
        cat("\n")
        stop(paste(paste0(tissue[!tissue %in% ref$tissueType], collapse = ","), ", not matched with the tissue types in ref database! Please select one or more related tissue types.",
          sep = ""
        ))
      }

      ref <- ref[ref$tissueType %in%
        tissue, ]
      cellmarker_num <- ref$geneSymbol
      tissue1 <- tissue[1]
      if (length(tissue) > 1) {
        for (i in 2:length(tissue)) {
          tissue1 <- paste(tissue1, tissue[i], sep = ", ")
        }
      }
      if (length(cellmarker_num) < 100) {
        cat("\n")
        cat(
          paste("Warning: there are only ", length(cellmarker_num),
            " potential marker genes in ref database for ",
            species, " on ", tissue1, "!",
            sep = ""
          ),
          "\n"
        )
      }
      if (length(cellmarker_num) >= 100) {
        cat("\n")
        cat(
          paste("Note: there are ", length(cellmarker_num),
            " potential marker genes in ref database for ",
            species, " on ", tissue1, ".",
            sep = ""
          ),
          "\n"
        )
      }
      ndata <- ndata[rownames(ndata) %in% cellmarker_num, ]
    }
    if (!is.null(cancer)) {
      ref <- ref[ref$speciesType ==
        species, ]
      if (cancer[1] == "All") {
        cancer <- unique(ref$cancerType)
      }
      if (!all(cancer %in% ref$cancerType)) {
        cat("\n")
        stop(paste(paste0(cancer[!cancer %in% ref$cancerType], collapse = ","), ", not matched with the cancer types in ref database! Please select one or more related cancer types.",
          sep = ""
        ))
      }
      ref <- ref[ref$cancerType %in%
        cancer, ]

      if (is.null(tissue)) {
        cat("\n")
        stop("Please define the origin tissue of cells! Select one or more related tissue types.")
      }
      if (tissue[1] == "All") {
        tissue <- unique(ref$tissueType)
      }
      if (!all(tissue %in% ref$tissueType)) {
        cat("\n")
        stop(paste(paste0(tissue[!tissue %in% ref$tissueType], collapse = ","), ", not matched with the tissue types in ref database! Please select one or more related tissue types.",
          sep = ""
        ))
      }
      ref <- ref[ref$tissueType %in%
        tissue, ]

      cellmarker_num <- ref$geneSymbol
      cancer1 <- cancer[1]
      if (length(cancer) > 1) {
        for (i in 2:length(cancer)) {
          cancer1 <- paste(cancer1, cancer[i], sep = ", ")
        }
      }
      tissue1 <- tissue[1]
      if (length(tissue) > 1) {
        for (i in 2:length(tissue)) {
          tissue1 <- paste(tissue1, tissue[i], sep = ", ")
        }
      }
      if (length(cellmarker_num) < 100) {
        cat("\n")
        cat(paste("Warning: there are only ", length(cellmarker_num),
          " potential marker genes in ref database for ",
          species, " ", cancer1, " on ", tissue1, "!",
          sep = ""
        ), "\n")
      }
      if (length(cellmarker_num) >= 100) {
        cat("\n")
        cat(paste("Note: there are ", length(cellmarker_num),
          " potential marker genes in ref database for ",
          species, " ", cancer1, " on ", tissue1, ".",
          sep = ""
        ), "\n")
      }
      ndata <- ndata[rownames(ndata) %in% cellmarker_num, ]
    }
  }
  if (nrow(ndata) == 0) {
    cat("\n")
    stop("There is no matched potential marker genes in the matrix! Please try again to find differential expressed genes by setting match_ref as FALSE!")
  }
  clu_pair <- NULL
  for (i in 1:length(clu_num)) {
    d1 <- data.frame(cluster1 = rep(clu_num[i], (length(clu_num) -
      1)), cluster2 = clu_num[-i], stringsAsFactors = F)
    clu_pair <- rbind(clu_pair, d1)
  }
  if (cluster[1] != "All") {
    cluster_match <- cluster %in% clu_num
    cluster_match <- which(cluster_match == "FALSE")
    if (length(cluster_match) > 0) {
      cat("\n")
      stop(paste(cluster[cluster_match], ", not matched with the cell clusters! Please select one or more related clusters.",
        sep = ""
      ))
    }
    clu_pair <- clu_pair[clu_pair$cluster1 %in% cluster, ]
    clu_num1 <- clu_num[clu_num %in% cluster]
  }
  cat(
    "Note: There are", length(clu_num), "clusters in Seurat object.",
    "\n"
  )
  Sys.sleep(2)
  cat(
    "Preparing to find potential marker genes for", cluster,
    "cluster(s).", "\n"
  )
  clu_marker <- NULL
  Sys.sleep(2)
  # for (i in 1:length(clu_num1)) {
  clu_marker_list <- bplapply(1:length(clu_num1), function(i) {
    cat(paste("Finding potential marker genes for cluster",
      clu_num1[i],
      sep = " "
    ), "\n")
    clu_pair1 <- clu_pair[clu_pair$cluster1 == clu_num1[i], ]
    clu_marker1 <- NULL
    # pb <- progress_bar$new(
    #   format = "[:bar] Finished::percent Remaining::eta",
    #   total = nrow(clu_pair1) * nrow(ndata), clear = FALSE,
    #   width = 60, complete = "+", incomplete = "-"
    # )
    for (j in 1:nrow(clu_pair1)) {
      clu_marker2 <- as.data.frame(matrix(
        data = 0, nrow = nrow(ndata),
        ncol = 6
      ))
      colnames(clu_marker2) <- c(
        "cluster", "gene", "pct",
        "comp_cluster", "avg_logfc", "pvalue"
      )
      clu_marker2[, 2] <- rownames(ndata)
      ndata1 <- ndata[, clu_info[clu_info$cluster == clu_pair1$cluster1[j], ]$cell]
      ndata2 <- ndata[, clu_info[clu_info$cluster == clu_pair1$cluster2[j], ]$cell]
      clu_marker_pct <- rep(0, nrow(clu_marker2))
      clu_marker_pvalue <- rep(1, nrow(clu_marker2))
      clu_marker2$cluster <- clu_pair1$cluster1[j]
      clu_marker2$comp_cluster <- clu_pair1$cluster2[j]
      clu_marker2$avg_logfc <- Matrix::rowMeans(ndata1) -
        Matrix::rowMeans(ndata2)
      for (k in 1:nrow(clu_marker2)) {
        genedata1 <- as.numeric(ndata1[k, ])
        if (length(genedata1[genedata1 > 0]) >= (length(genedata1) *
          cell_min_pct)) {
          clu_marker_pct[k] <- (length(genedata1[genedata1 >
            0])) / (length(genedata1))
          genedata2 <- as.numeric(ndata2[k, ])
          clu_marker_pvalue[k] <- stats::wilcox.test(
            genedata1,
            genedata2,
          )$p.value
        }
        # pb$tick()
      }
      clu_marker2$pct <- clu_marker_pct
      clu_marker2$pvalue <- clu_marker_pvalue
      clu_marker2 <- clu_marker2[clu_marker2$pct >= cell_min_pct &
        clu_marker2$avg_logfc >= logfc & clu_marker2$pvalue <
        pvalue, ]
      clu_marker1 <- rbind(clu_marker1, clu_marker2)
    }
    if (nrow(clu_marker1) > 0) {
      d1 <- as.data.frame(table(clu_marker1$gene), stringsAsFactors = F)
      d1 <- d1[d1$Freq == (length(clu_num) - 1), ]
      if (nrow(d1) > 1) {
        clu_marker1 <- clu_marker1[clu_marker1$gene %in%
          d1$Var1, ]
        clu_marker_gene <- unique(clu_marker1$gene)
        clu_marker_gene <- clu_marker_gene[order(clu_marker_gene)]
        avg_logfc <- NULL
        for (j in 1:length(clu_marker_gene)) {
          clu_marker_avg_logfc <- clu_marker1[clu_marker1$gene ==
            clu_marker_gene[j], ]
          avg_logfc[j] <- mean(clu_marker_avg_logfc$avg_logfc)
        }
        clu_marker1 <- unique(clu_marker1[, c(
          "cluster",
          "gene", "pct"
        )])
        clu_marker1 <- clu_marker1[order(clu_marker1$gene), ]
        clu_marker1$avg_logfc <- avg_logfc
        return(clu_marker1)
        # clu_marker <- rbind(clu_marker, clu_marker1)
      }
    }
    cat("---Done---", "\n")
    # Sys.sleep(2)
    # }
  }, BPPARAM = BPPARAM)
  clu_marker <- bind_rows(clu_marker_list)

  # rm(ndata1, ndata2)
  res <- list()
  if (match_ref) {
    res[[1]] <- ndata3
    rm(ndata3)
  }
  if (!match_ref) {
    res[[1]] <- ndata
    rm(ndata)
  }
  res[[2]] <- "NA"
  names(res) <- c("new_data_matrix", "clu_markers")
  if (is.null(clu_marker)) {
    cat("\n")
    cat("Warning: there is no potential marker genes found in the matrix! You may adjust the cell clusters by applying other clustering algorithm and try again.")
    return(res)
    stop()
  }
  if (nrow(clu_marker) == 0) {
    cat("\n")
    cat("Warning: there is no potential marker genes found in the matrix! You may adjust the cell clusters by applying other clustering algorithm and try again.")
    return(res)
    stop()
  }
  if (nrow(clu_marker) > 0) {
    res[[2]] <- clu_marker
    return(res)
  }
}


scCATCH <- function(object, species = NULL, ref = NULL, cancer = NULL, tissue = NULL) {
  if (is.null(species)) {
    stop("Please define the species! 'Human' or 'Mouse'.")
  }
  if (is.null(ref)) {
    stop("ref must be specified.")
  }
  ref <- ref[ref$speciesType == species, ]
  if (!is.data.frame(object)) {
    if (names(object)[2] != "clu_markers") {
      stop("Please select the clu_clusters as the input!")
    }
    object <- object$clu_markers
  }
  object$genesymbol <- object$gene
  if (is.null(cancer)) {
    ref <- ref[ref$cancerType == "Normal", ]
    if (is.null(tissue)) {
      stop("Please define the origin tissue of cells! Select one or more related tissue types.")
    }
    if (tissue[1] == "All") {
      tissue <- unique(ref$tissueType)
    }
    if (!all(tissue %in% ref$tissueType)) {
      cat("\n")
      stop(paste(paste0(tissue[!tissue %in% ref$tissueType], collapse = ","), ", not matched with the tissue types in ref database! Please select one or more related tissue types.",
        sep = ""
      ))
    }
    ref <- ref[ref$tissueType %in% tissue, ]
    cellmarker_num <- ref$geneSymbol
    tissue1 <- tissue[1]
    if (length(tissue) > 1) {
      for (i in 2:length(tissue)) {
        tissue1 <- paste(tissue1, tissue[i], sep = ", ")
      }
    }
    if (length(cellmarker_num) < 100) {
      cat(paste("Warning: there are only ", length(cellmarker_num),
        " potential marker genes in ref database for ",
        species, " on ", tissue1, "!",
        sep = ""
      ), "\n")
    }
    if (length(cellmarker_num) >= 100) {
      cat(paste("Note: there are ", length(cellmarker_num),
        " potential marker genes in ref database for ",
        species, " on ", tissue1, "!",
        sep = ""
      ), "\n")
    }
    Sys.sleep(2)
    cellmarker_num <- unique(ref$cellName)
    if (length(cellmarker_num) < 10) {
      cat(paste("Warning: there are only ", length(cellmarker_num),
        " cell types in ref database for ", species,
        " on ", tissue1, "!",
        sep = ""
      ), "\n")
    }
    if (length(cellmarker_num) >= 10) {
      cat(paste("Note: there are ", length(cellmarker_num),
        " cell types in ref database for ", species,
        " on ", tissue1, "!",
        sep = ""
      ), "\n")
    }
  }
  if (!is.null(cancer)) {
    ref <- ref[ref$speciesType ==
      species, ]
    if (cancer[1] == "All") {
      cancer <- unique(ref$cancerType)
    }
    if (!all(cancer %in% ref$cancerType)) {
      cat("\n")
      stop(paste(paste0(cancer[!cancer %in% ref$cancerType], collapse = ","), ", not matched with the cancer types in ref database! Please select one or more related cancer types.",
        sep = ""
      ))
    }
    ref <- ref[ref$cancerType %in% cancer, ]
    if (is.null(tissue)) {
      stop("Please define the origin tissue of cells! Select one or more related tissue types.")
    }
    tissue_match <- NULL
    tissue_match <- tissue %in% ref$tissueType
    tissue_match <- which(tissue_match == "FALSE")
    if (length(tissue_match) > 0) {
      stop(paste(tissue[tissue_match], ", not matched with the tissue types in ref database! Please select one or more related tissue types.",
        sep = ""
      ))
    }
    ref <- ref[ref$tissueType %in% tissue, ]
    cellmarker_num <- ref$geneSymbol
    cancer1 <- cancer[1]
    if (length(cancer) > 1) {
      for (i in 2:length(cancer)) {
        cancer1 <- paste(cancer1, cancer[i], sep = ", ")
      }
    }
    tissue1 <- tissue[1]
    if (length(tissue) > 1) {
      for (i in 2:length(tissue)) {
        tissue1 <- paste(tissue1, tissue[i], sep = ", ")
      }
    }
    if (length(cellmarker_num) < 100) {
      cat(paste("Warning: there are only ", length(cellmarker_num),
        " potential marker genes in ref database for ",
        species, " ", cancer1, " on ", tissue1, "!",
        sep = ""
      ), "\n")
    }
    if (length(cellmarker_num) >= 100) {
      cat(paste("Note: there are ", length(cellmarker_num),
        " potential marker genes in ref database for ",
        species, " ", cancer1, " on ", tissue1, "!",
        sep = ""
      ), "\n")
    }
    Sys.sleep(2)
    cellmarker_num <- unique(ref$cellName)
    if (length(cellmarker_num) < 10) {
      cat(paste("Warning: there are only ", length(cellmarker_num),
        " cell types in ref database for ", species,
        " on ", tissue1, "!",
        sep = ""
      ), "\n")
    }
    if (length(cellmarker_num) >= 10) {
      cat(paste("Note: there are ", length(cellmarker_num),
        " cell types in ref database for ", species,
        " on ", tissue1, "!",
        sep = ""
      ), "\n")
    }
  }
  clu.num <- as.character(unique(object$cluster))
  clu.num1 <- clu.num[nchar(clu.num) == 1]
  clu.num2 <- clu.num[nchar(clu.num) == 2]
  clu.num1 <- clu.num1[order(clu.num1)]
  clu.num2 <- clu.num2[order(clu.num2)]
  clu.num <- c(clu.num1, clu.num2)
  Sys.sleep(2)
  cat("Beginning evidence-based scoring and annotation", "\n")
  Sys.sleep(2)
  # pb <- progress_bar$new(
  #   format = "[:bar] Finished::percent Remaining::eta",
  #   total = length(clu.num), clear = FALSE, width = 60,
  #   complete = "+", incomplete = "-"
  # )
  clu_ann_res <- NULL
  for (i in 1:length(clu.num)) {
    Sys.sleep(1)
    rescluster_marker1 <- object[object$cluster == clu.num[i], ]$genesymbol
    cellsubtype1 <- "NA"
    cellsubtype2 <- "NA"
    cellsubtype3 <- "NA"
    celltype <- "NA"
    celltype_score <- "NA"
    clu_marker <- "NA"
    PMID <- "NA"
    res <- "I"
    if (length(rescluster_marker1) > 0) {
      rescluster_marker2 <- rescluster_marker1[rescluster_marker1 %in%
        ref$geneSymbol]
      res <- "II"
      if (length(rescluster_marker2) > 0) {
        clu_ann <- ref[ref$geneSymbol %in%
          rescluster_marker2, ]
        clu_ann_cellname <- as.data.frame(table(clu_ann$shortname),
          stringsAsFactors = F
        )
        m <- clu_ann_cellname$Freq + 1
        clu_ann_cellname$Freq <- clu_ann_cellname$Freq / m
        clu_ann_cellname <- clu_ann_cellname[order(clu_ann_cellname$Var1), ]
        clu_ann_article <- unique(clu_ann[, c(
          "PMID",
          "shortname"
        )])
        clu_ann_article <- as.data.frame(table(clu_ann_article$shortname),
          stringsAsFactors = F
        )
        m <- clu_ann_article$Freq + 1
        clu_ann_article$Freq <- clu_ann_article$Freq / m
        clu_ann_article <- clu_ann_article[order(clu_ann_article$Var1), ]
        clu_ann_cellname$Freq <- sqrt(clu_ann_article$Freq *
          clu_ann_cellname$Freq)
        clu_ann_cellname <- clu_ann_cellname[clu_ann_cellname$Freq ==
          max(clu_ann_cellname$Freq), ]
        if (nrow(clu_ann_cellname) > 1) {
          celltype <- clu_ann_cellname$Var1
          celltype_score <- round(clu_ann_cellname$Freq,
            digits = 2
          )
          clu_marker <- unique(clu_ann[clu_ann$shortname %in%
            clu_ann_cellname$Var1, ]$geneSymbol)
          PMID <- unique(clu_ann[clu_ann$shortname %in%
            clu_ann_cellname$Var1, ]$PMID)
          res <- "III"
        }
        if (nrow(clu_ann_cellname) == 1) {
          celltype <- clu_ann_cellname$Var1
          celltype_score <- round(clu_ann_cellname$Freq,
            digits = 2
          )
          clu_marker <- unique(clu_ann[clu_ann$shortname %in%
            clu_ann_cellname$Var1, ]$geneSymbol)
          PMID <- unique(clu_ann[clu_ann$shortname %in%
            clu_ann_cellname$Var1, ]$PMID)
          res <- "IV"
          clu_ann1 <- clu_ann[clu_ann$shortname == clu_ann_cellname$Var1, ]
          clu_ann1 <- clu_ann1[, c(
            "geneSymbol", "PMID",
            "third", "second", "first"
          )]
          clu_ann1$row <- 1:nrow(clu_ann1)
          clu_ann1$row <- as.character(clu_ann1$row)
          clu_ann1 <- clu_ann1[, c(6, 1:5)]
          clu_ann1 <- data.table::as.data.table(clu_ann1)
          clu_ann1 <- reshape2::melt(
            data = clu_ann1, id.vars = 1:3,
            measure.vars = c("third", "second", "first")
          )
          clu_ann1 <- as.data.frame(clu_ann1)
          clu_ann1$variable <- as.character(clu_ann1$variable)
          clu_ann1 <- clu_ann1[!is.na(clu_ann1$value), ]
          if (nrow(clu_ann1) > 0) {
            clu_ann_subcellname <- as.data.frame(table(clu_ann1$value),
              stringsAsFactors = F
            )
            m <- clu_ann_subcellname$Freq + 1
            clu_ann_subcellname$Freq <- clu_ann_subcellname$Freq / m
            clu_ann_subcellname <- clu_ann_subcellname[order(clu_ann_subcellname$Var1), ]
            clu_ann_subarticle <- unique(clu_ann1[
              ,
              c("PMID", "value")
            ])
            clu_ann_subarticle <- as.data.frame(table(clu_ann_subarticle$value),
              stringsAsFactors = F
            )
            m <- clu_ann_subarticle$Freq + 1
            clu_ann_subarticle$Freq <- clu_ann_subarticle$Freq / m
            clu_ann_subarticle <- clu_ann_subarticle[order(clu_ann_subarticle$Var1), ]
            clu_ann_subcellname$Freq <- sqrt(clu_ann_subcellname$Freq *
              clu_ann_subarticle$Freq)
            clu_ann_subcellname <- clu_ann_subcellname[clu_ann_subcellname$Freq ==
              max(clu_ann_subcellname$Freq), ]
            clu_ann_subcellname <- clu_ann_subcellname[clu_ann_subcellname$Freq >
              0.5, ]
            if (nrow(clu_ann_subcellname) == 1) {
              clu_ann_for_det <- clu_ann1[clu_ann1$value ==
                clu_ann_subcellname$Var1, ]
              if (unique(clu_ann_for_det$variable ==
                "first")) {
                cellsubtype1 <- unique(clu_ann_for_det$value)
                clu_marker <- unique(clu_ann_for_det[clu_ann_for_det$value %in%
                  cellsubtype1, ]$geneSymbol)
                PMID <- unique(clu_ann_for_det[clu_ann_for_det$value %in%
                  cellsubtype1, ]$PMID)
                clu_ann_for_second <- clu_ann1[clu_ann1$variable ==
                  "second", ]
                if (nrow(clu_ann_for_second) > 0) {
                  clu_second <- clu_ann1[clu_ann1$variable ==
                    "second", ]
                  clu_second_name <- as.data.frame(table(clu_second$value),
                    stringsAsFactors = F
                  )
                  m <- clu_second_name$Freq + 1
                  clu_second_name$Freq <- clu_second_name$Freq / m
                  clu_second_name <- clu_second_name[order(clu_second_name$Var1), ]
                  clu_second_article <- unique(clu_second[
                    ,
                    c("PMID", "value")
                  ])
                  clu_second_article <- as.data.frame(table(clu_second_article$value),
                    stringsAsFactors = F
                  )
                  m <- clu_second_article$Freq + 1
                  clu_second_article$Freq <- clu_second_article$Freq / m
                  clu_second_article <- clu_second_article[order(clu_second_article$Var1), ]
                  clu_second_name$Freq <- sqrt(clu_second_name$Freq *
                    clu_second_article$Freq)
                  clu_ann_second <- clu_ann[clu_ann$first ==
                    cellsubtype1, ]
                  clu_ann_second <- clu_ann_second[!is.na(clu_ann_second$first), ]
                  clu_ann_secondname <- as.data.frame(table(clu_ann_second$second),
                    stringsAsFactors = F
                  )
                  if (nrow(clu_ann_secondname) > 0) {
                    m <- clu_ann_secondname$Freq + 1
                    clu_ann_secondname$Freq <- clu_ann_secondname$Freq / m
                    clu_ann_secondname <- clu_ann_secondname[order(clu_ann_secondname$Var1), ]
                    clu_ann_secondarticle <- unique(clu_ann_second[
                      ,
                      c("PMID", "second")
                    ])
                    clu_ann_secondarticle <- as.data.frame(table(clu_ann_secondarticle$second),
                      stringsAsFactors = F
                    )
                    m <- clu_ann_secondarticle$Freq +
                      1
                    clu_ann_secondarticle$Freq <- clu_ann_secondarticle$Freq / m
                    clu_ann_secondarticle <- clu_ann_secondarticle[order(clu_ann_secondarticle$Var1), ]
                    clu_ann_secondname$Freq <- sqrt(clu_ann_secondname$Freq *
                      clu_ann_secondarticle$Freq)
                    clu_ann_secondname <- clu_ann_secondname[clu_ann_secondname$Freq ==
                      max(clu_ann_secondname$Freq), ]
                    clu_ann_secondname <- clu_ann_secondname[(clu_ann_secondname$Freq >
                      0.5) & (clu_ann_secondname$Freq >=
                      max(clu_second_name$Freq)), ]
                    if (nrow(clu_ann_secondname) ==
                      1) {
                      cellsubtype2 <- clu_ann_secondname$Var1
                      clu_marker <- unique(clu_ann1[clu_ann1$value %in%
                        cellsubtype2, ]$geneSymbol)
                      PMID <- unique(clu_ann1[clu_ann1$value %in%
                        cellsubtype2, ]$PMID)
                      clu_ann_for_third <- clu_ann1[clu_ann1$variable ==
                        "third", ]
                      if (nrow(clu_ann_for_third) >
                        0) {
                        clu_third <- clu_ann1[clu_ann1$variable ==
                          "third", ]
                        clu_third_name <- as.data.frame(table(clu_third$value),
                          stringsAsFactors = F
                        )
                        m <- clu_third_name$Freq + 1
                        clu_third_name$Freq <- clu_third_name$Freq / m
                        clu_third_name <- clu_third_name[order(clu_third_name$Var1), ]
                        clu_third_article <- unique(clu_third[
                          ,
                          c("PMID", "value")
                        ])
                        clu_third_article <- as.data.frame(table(clu_third_article$value),
                          stringsAsFactors = F
                        )
                        m <- clu_third_article$Freq +
                          1
                        clu_third_article$Freq <- clu_third_article$Freq / m
                        clu_third_article <- clu_third_article[order(clu_third_article$Var1), ]
                        clu_third_name$Freq <- sqrt(clu_third_name$Freq *
                          clu_third_article$Freq)
                        clu_ann_third <- clu_ann_second[clu_ann_second$second ==
                          cellsubtype2, ]
                        clu_ann_third <- clu_ann_third[!is.na(clu_ann_third$second), ]
                        clu_ann_thirdname <- as.data.frame(table(clu_ann_third$third),
                          stringsAsFactors = F
                        )
                        if (nrow(clu_ann_thirdname) >
                          0) {
                          m <- clu_ann_thirdname$Freq +
                            1
                          clu_ann_thirdname$Freq <- clu_ann_thirdname$Freq / m
                          clu_ann_thirdname <- clu_ann_thirdname[order(clu_ann_thirdname$Var1), ]
                          clu_ann_thirdarticle <- unique(clu_ann_third[
                            ,
                            c("PMID", "third")
                          ])
                          clu_ann_thirdarticle <- as.data.frame(table(clu_ann_thirdarticle$third),
                            stringsAsFactors = F
                          )
                          m <- clu_ann_thirdarticle$Freq +
                            1
                          clu_ann_thirdarticle$Freq <- clu_ann_thirdarticle$Freq / m
                          clu_ann_thirdarticle <- clu_ann_thirdarticle[order(clu_ann_thirdarticle$Var1), ]
                          clu_ann_thirdname$Freq <- sqrt(clu_ann_thirdname$Freq *
                            clu_ann_thirdarticle$Freq)
                          clu_ann_thirdname <- clu_ann_thirdname[clu_ann_thirdname$Freq ==
                            max(clu_ann_thirdname$Freq), ]
                          clu_ann_thirdname <- clu_ann_thirdname[(clu_ann_thirdname$Freq >
                            0.5) & (clu_ann_thirdname$Freq >=
                            max(clu_third_name$Freq)), ]
                          if (nrow(clu_ann_thirdname) >=
                            1) {
                            cellsubtype3 <- clu_ann_thirdname$Var1
                            clu_marker <- unique(clu_ann1[clu_ann1$value %in%
                              cellsubtype3, ]$geneSymbol)
                            PMID <- unique(clu_ann1[clu_ann1$value %in%
                              cellsubtype3, ]$PMID)
                          }
                        }
                      }
                    }
                  }
                }
              }
              if (unique(clu_ann_for_det$variable ==
                "second")) {
                cellsubtype2 <- unique(clu_ann_for_det$value)
                clu_marker <- unique(clu_ann_for_det[clu_ann_for_det$value %in%
                  cellsubtype2, ]$geneSymbol)
                PMID <- unique(clu_ann_for_det[clu_ann_for_det$value %in%
                  cellsubtype2, ]$PMID)
                clu_ann_for_first <- clu_ann1[clu_ann1$variable ==
                  "first", ]
                if (nrow(clu_ann_for_first) > 0) {
                  clu_first <- clu_ann1[clu_ann1$variable ==
                    "first", ]
                  clu_first_name <- as.data.frame(table(clu_first$value),
                    stringsAsFactors = F
                  )
                  m <- clu_first_name$Freq + 1
                  clu_first_name$Freq <- clu_first_name$Freq / m
                  clu_first_name <- clu_first_name[order(clu_first_name$Var1), ]
                  clu_first_article <- unique(clu_first[
                    ,
                    c("PMID", "value")
                  ])
                  clu_first_article <- as.data.frame(table(clu_first_article$value),
                    stringsAsFactors = F
                  )
                  m <- clu_first_article$Freq + 1
                  clu_first_article$Freq <- clu_first_article$Freq / m
                  clu_first_article <- clu_first_article[order(clu_first_article$Var1), ]
                  clu_first_name$Freq <- sqrt(clu_first_name$Freq *
                    clu_first_article$Freq)
                  clu_ann_first <- clu_ann[clu_ann$second ==
                    cellsubtype2, ]
                  clu_ann_first <- clu_ann_first[!is.na(clu_ann_first$second), ]
                  clu_ann_firstname <- as.data.frame(table(clu_ann_first$first),
                    stringsAsFactors = F
                  )
                  if (nrow(clu_ann_firstname) > 0) {
                    m <- clu_ann_firstname$Freq + 1
                    clu_ann_firstname$Freq <- clu_ann_firstname$Freq / m
                    clu_ann_firstname <- clu_ann_firstname[order(clu_ann_firstname$Var1), ]
                    clu_ann_firstarticle <- unique(clu_ann_first[
                      ,
                      c("PMID", "first")
                    ])
                    clu_ann_firstarticle <- as.data.frame(table(clu_ann_firstarticle$first),
                      stringsAsFactors = F
                    )
                    m <- clu_ann_firstarticle$Freq +
                      1
                    clu_ann_firstarticle$Freq <- clu_ann_firstarticle$Freq / m
                    clu_ann_firstarticle <- clu_ann_firstarticle[order(clu_ann_firstarticle$Var1), ]
                    clu_ann_firstname$Freq <- sqrt(clu_ann_firstname$Freq *
                      clu_ann_firstarticle$Freq)
                    clu_ann_firstname <- clu_ann_firstname[clu_ann_firstname$Freq ==
                      max(clu_ann_firstname$Freq), ]
                    clu_ann_firstname <- clu_ann_firstname[(clu_ann_firstname$Freq >
                      0.5) & (clu_ann_firstname$Freq >=
                      max(clu_ann_firstname$Freq)), ]
                    if (nrow(clu_ann_firstname) == 1) {
                      cellsubtype1 <- clu_ann_firstname$Var1
                      clu_marker <- unique(clu_ann1[clu_ann1$value %in%
                        cellsubtype1, ]$geneSymbol)
                      PMID <- unique(clu_ann1[clu_ann1$value %in%
                        cellsubtype1, ]$PMID)
                      clu_ann_for_third <- clu_ann1[clu_ann1$variable ==
                        "third", ]
                      if (nrow(clu_ann_for_third) >
                        0) {
                        clu_third <- clu_ann1[clu_ann1$variable ==
                          "third", ]
                        clu_third_name <- as.data.frame(table(clu_third$value),
                          stringsAsFactors = F
                        )
                        m <- clu_third_name$Freq + 1
                        clu_third_name$Freq <- clu_third_name$Freq / m
                        clu_third_name <- clu_third_name[order(clu_third_name$Var1), ]
                        clu_third_article <- unique(clu_third[
                          ,
                          c("PMID", "value")
                        ])
                        clu_third_article <- as.data.frame(table(clu_third_article$value),
                          stringsAsFactors = F
                        )
                        m <- clu_third_article$Freq +
                          1
                        clu_third_article$Freq <- clu_third_article$Freq / m
                        clu_third_article <- clu_third_article[order(clu_third_article$Var1), ]
                        clu_third_name$Freq <- sqrt(clu_third_name$Freq *
                          clu_third_article$Freq)
                        clu_ann_third <- clu_ann_first[clu_ann_first$first ==
                          cellsubtype1, ]
                        clu_ann_third <- clu_ann_third[!is.na(clu_ann_third$first), ]
                        clu_ann_thirdname <- as.data.frame(table(clu_ann_third$third),
                          stringsAsFactors = F
                        )
                        if (nrow(clu_ann_thirdname) >
                          0) {
                          m <- clu_ann_thirdname$Freq +
                            1
                          clu_ann_thirdname$Freq <- clu_ann_thirdname$Freq / m
                          clu_ann_thirdname <- clu_ann_thirdname[order(clu_ann_thirdname$Var1), ]
                          clu_ann_thirdarticle <- unique(clu_ann_third[
                            ,
                            c("PMID", "third")
                          ])
                          clu_ann_thirdarticle <- as.data.frame(table(clu_ann_thirdarticle$third),
                            stringsAsFactors = F
                          )
                          m <- clu_ann_thirdarticle$Freq +
                            1
                          clu_ann_thirdarticle$Freq <- clu_ann_thirdarticle$Freq / m
                          clu_ann_thirdarticle <- clu_ann_thirdarticle[order(clu_ann_thirdarticle$Var1), ]
                          clu_ann_thirdname$Freq <- sqrt(clu_ann_thirdname$Freq *
                            clu_ann_thirdarticle$Freq)
                          clu_ann_thirdname <- clu_ann_thirdname[clu_ann_thirdname$Freq ==
                            max(clu_ann_thirdname$Freq), ]
                          clu_ann_thirdname <- clu_ann_thirdname[(clu_ann_thirdname$Freq >
                            0.5) & (clu_ann_thirdname$Freq >=
                            max(clu_third_name$Freq)), ]
                          if (nrow(clu_ann_thirdname) >=
                            1) {
                            cellsubtype3 <- clu_ann_thirdname$Var1
                            clu_marker <- unique(clu_ann1[clu_ann1$value %in%
                              cellsubtype3, ]$geneSymbol)
                            PMID <- unique(clu_ann1[clu_ann1$value %in%
                              cellsubtype3, ]$PMID)
                          }
                        }
                      }
                    }
                  }
                }
              }
              if (unique(clu_ann_for_det$variable ==
                "third")) {
                cellsubtype3 <- unique(clu_ann_for_det$value)
                clu_marker <- unique(clu_ann_for_det[clu_ann_for_det$value %in%
                  cellsubtype3, ]$geneSymbol)
                PMID <- unique(clu_ann_for_det[clu_ann_for_det$value %in%
                  cellsubtype3, ]$PMID)
                clu_ann_for_second <- clu_ann1[clu_ann1$variable ==
                  "second", ]
                if (nrow(clu_ann_for_second) > 0) {
                  clu_second <- clu_ann1[clu_ann1$variable ==
                    "second", ]
                  clu_second_name <- as.data.frame(table(clu_second$value),
                    stringsAsFactors = F
                  )
                  m <- clu_second_name$Freq + 1
                  clu_second_name$Freq <- clu_second_name$Freq / m
                  clu_second_name <- clu_second_name[order(clu_second_name$Var1), ]
                  clu_second_article <- unique(clu_second[
                    ,
                    c("PMID", "value")
                  ])
                  clu_second_article <- as.data.frame(table(clu_second_article$value),
                    stringsAsFactors = F
                  )
                  m <- clu_second_article$Freq + 1
                  clu_second_article$Freq <- clu_second_article$Freq / m
                  clu_second_article <- clu_second_article[order(clu_second_article$Var1), ]
                  clu_second_name$Freq <- sqrt(clu_second_name$Freq *
                    clu_second_article$Freq)
                  clu_ann_second <- clu_ann[clu_ann$third ==
                    cellsubtype3, ]
                  clu_ann_second <- clu_ann_second[!is.na(clu_ann_second$third), ]
                  clu_ann_secondname <- as.data.frame(table(clu_ann_second$second),
                    stringsAsFactors = F
                  )
                  if (nrow(clu_ann_secondname) > 0) {
                    m <- clu_ann_secondname$Freq + 1
                    clu_ann_secondname$Freq <- clu_ann_secondname$Freq / m
                    clu_ann_secondname <- clu_ann_secondname[order(clu_ann_secondname$Var1), ]
                    clu_ann_secondarticle <- unique(clu_ann_second[
                      ,
                      c("PMID", "second")
                    ])
                    clu_ann_secondarticle <- as.data.frame(table(clu_ann_secondarticle$second),
                      stringsAsFactors = F
                    )
                    m <- clu_ann_secondarticle$Freq +
                      1
                    clu_ann_secondarticle$Freq <- clu_ann_secondarticle$Freq / m
                    clu_ann_secondarticle <- clu_ann_secondarticle[order(clu_ann_secondarticle$Var1), ]
                    clu_ann_secondname$Freq <- sqrt(clu_ann_secondname$Freq *
                      clu_ann_secondarticle$Freq)
                    clu_ann_secondname <- clu_ann_secondname[clu_ann_secondname$Freq ==
                      max(clu_ann_secondname$Freq), ]
                    clu_ann_secondname <- clu_ann_secondname[(clu_ann_secondname$Freq >
                      0.5) & (clu_ann_secondname$Freq >=
                      max(clu_second_name$Freq)), ]
                    if (nrow(clu_ann_secondname) ==
                      1) {
                      cellsubtype2 <- clu_ann_secondname$Var1
                      clu_marker <- unique(clu_ann1[clu_ann1$value %in%
                        cellsubtype2, ]$geneSymbol)
                      PMID <- unique(clu_ann1[clu_ann1$value %in%
                        cellsubtype2, ]$PMID)
                      clu_ann_for_first <- clu_ann1[clu_ann1$variable ==
                        "first", ]
                      if (nrow(clu_ann_for_first) >
                        0) {
                        clu_first <- clu_ann1[clu_ann1$variable ==
                          "first", ]
                        clu_first_name <- as.data.frame(table(clu_first$value),
                          stringsAsFactors = F
                        )
                        m <- clu_first_name$Freq + 1
                        clu_first_name$Freq <- clu_first_name$Freq / m
                        clu_first_name <- clu_first_name[order(clu_first_name$Var1), ]
                        clu_first_article <- unique(clu_first[
                          ,
                          c("PMID", "value")
                        ])
                        clu_first_article <- as.data.frame(table(clu_first_article$value),
                          stringsAsFactors = F
                        )
                        m <- clu_first_article$Freq +
                          1
                        clu_first_article$Freq <- clu_first_article$Freq / m
                        clu_first_article <- clu_first_article[order(clu_first_article$Var1), ]
                        clu_first_name$Freq <- sqrt(clu_first_name$Freq *
                          clu_first_article$Freq)
                        clu_ann_first <- clu_ann_second[clu_ann_second$second ==
                          cellsubtype2, ]
                        clu_ann_first <- clu_ann_first[!is.na(clu_ann_first$second), ]
                        clu_ann_firstname <- as.data.frame(table(clu_ann_first$first),
                          stringsAsFactors = F
                        )
                        if (nrow(clu_ann_firstname) >
                          0) {
                          m <- clu_ann_firstname$Freq +
                            1
                          clu_ann_firstname$Freq <- clu_ann_firstname$Freq / m
                          clu_ann_firstname <- clu_ann_firstname[order(clu_ann_firstname$Var1), ]
                          clu_ann_firstarticle <- unique(clu_ann_first[
                            ,
                            c("PMID", "first")
                          ])
                          clu_ann_firstarticle <- as.data.frame(table(clu_ann_firstarticle$first),
                            stringsAsFactors = F
                          )
                          m <- clu_ann_firstarticle$Freq +
                            1
                          clu_ann_firstarticle$Freq <- clu_ann_firstarticle$Freq / m
                          clu_ann_firstarticle <- clu_ann_firstarticle[order(clu_ann_firstarticle$Var1), ]
                          clu_ann_firstname$Freq <- sqrt(clu_ann_firstname$Freq *
                            clu_ann_firstarticle$Freq)
                          clu_ann_firstname <- clu_ann_firstname[clu_ann_firstname$Freq ==
                            max(clu_ann_firstname$Freq), ]
                          clu_ann_firstname <- clu_ann_firstname[(clu_ann_firstname$Freq >
                            0.5) & (clu_ann_firstname$Freq >=
                            max(clu_first_name$Freq)), ]
                          if (nrow(clu_ann_firstname) >=
                            1) {
                            cellsubtype1 <- clu_ann_firstname$Var1
                            clu_marker <- unique(clu_ann1[clu_ann1$value %in%
                              cellsubtype1, ]$geneSymbol)
                            PMID <- unique(clu_ann1[clu_ann1$value %in%
                              cellsubtype1, ]$PMID)
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    if (length(rescluster_marker1) == 0) {
      rescluster_marker1 <- "NA"
    }
    d1 <- rescluster_marker1[1]
    if (length(rescluster_marker1) > 1) {
      for (j in 2:length(rescluster_marker1)) {
        d1 <- paste(d1, rescluster_marker1[j], sep = ", ")
      }
    }
    rescluster_marker1 <- d1
    d1 <- celltype[1]
    if (length(celltype) > 1) {
      for (j in 2:length(celltype)) {
        d1 <- paste(d1, celltype[j], sep = ", ")
      }
    }
    celltype <- d1
    d1 <- celltype_score[1]
    if (length(celltype_score) > 1) {
      for (j in 2:length(celltype_score)) {
        d1 <- paste(d1, celltype_score[j], sep = ", ")
      }
    }
    celltype_score <- d1
    d1 <- cellsubtype1[1]
    if (length(cellsubtype1) > 1) {
      for (j in 2:length(cellsubtype1)) {
        d1 <- paste(d1, cellsubtype1[j], sep = "; ")
      }
    }
    cellsubtype1 <- d1
    d1 <- cellsubtype2[1]
    if (length(cellsubtype2) > 1) {
      for (j in 2:length(cellsubtype2)) {
        d1 <- paste(d1, cellsubtype2[j], sep = "; ")
      }
    }
    cellsubtype2 <- d1
    d1 <- cellsubtype3[1]
    if (length(cellsubtype3) > 1) {
      for (j in 2:length(cellsubtype3)) {
        d1 <- paste(d1, cellsubtype3[j], sep = "; ")
      }
    }
    cellsubtype3 <- d1
    d1 <- clu_marker[1]
    if (length(clu_marker) > 1) {
      for (j in 2:length(clu_marker)) {
        d1 <- paste(d1, clu_marker[j], sep = ", ")
      }
    }
    clu_marker <- d1
    d1 <- PMID[1]
    if (length(PMID) > 1) {
      for (j in 2:length(PMID)) {
        d1 <- paste(d1, PMID[j], sep = ", ")
      }
    }
    PMID <- d1
    clu_ann <- data.frame(
      cluster = clu.num[i], cluster_marker = rescluster_marker1,
      cellsubtype3 = cellsubtype3, cellsubtype2 = cellsubtype2,
      cellsubtype1 = cellsubtype1, cell_type = celltype,
      celltype_score = celltype_score, celltype_related_marker = clu_marker,
      PMID = PMID, stringsAsFactors = F
    )
    clu_ann_res <- rbind(clu_ann_res, clu_ann)
    # pb$tick()
  }
  for (i in 1:nrow(clu_ann_res)) {
    d1 <- as.character(clu_ann_res[i, c(
      "cellsubtype3",
      "cellsubtype2", "cellsubtype1", "cell_type"
    )])
    d1 <- d1[!is.na(d1)]
    d1 <- d1[!d1 == "NA"]
    d2 <- d1[1]
    if (length(d1) > 1) {
      for (j in 2:length(d1)) {
        d2 <- paste(d2, d1[j], sep = " ")
      }
    }
    clu_ann_res[i, "cell_type"] <- d2
  }
  clu_ann_res <- clu_ann_res[, -c(3, 4, 5)]
  cat("---Done---", "\n")
  Sys.sleep(2)
  return(clu_ann_res)
}
