theme_zh <- function(aspect.ratio = 1, ...) {
  library(ggplot2)
  args1 <- list(
    aspect.ratio = aspect.ratio,
    text = element_text(size = 14, color = "black"),
    title = element_text(size = 14, colour = "black"),
    axis.line = element_blank(),
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 14, colour = "black"),
    strip.text = element_text(size = 12, colour = "black", face = "italic", hjust = 0.5),
    strip.background = element_rect(fill = "white", linetype = 0),
    legend.text = element_text(size = 12, colour = "black"),
    legend.title = element_text(size = 14, colour = "black", hjust = 0),
    legend.key = element_rect(fill = "transparent", color = "transparent"),
    legend.background = element_blank(),
    plot.subtitle = element_text(size = 12, hjust = 0),
    panel.background = element_blank(),
    panel.border = element_rect(fill = "transparent", colour = "black", size = 1),
    complete = TRUE
  )
  args2 <- c(as.list(match.call()))
  call.envir <- parent.frame(1)
  args2 <- lapply(args2, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call.envir)
    } else {
      arg
    }
  })
  for (n in names(args2)) {
    args1[[n]] <- args2[[n]]
  }
  out <- base::do.call(
    what = theme,
    args = args1
  )
  return(ggplot2::theme_classic() + out)
}

palette_zh <- function(x, palette = "Paired", type = "auto", matched = FALSE, reverse = FALSE, NA_color = "grey90") {
  library(RColorBrewer)
  library(ggsci)
  brewer.pal.info <- RColorBrewer::brewer.pal.info
  ggsci_db <- ggsci:::ggsci_db
  if (!palette %in% c(rownames(brewer.pal.info), names(ggsci_db))) {
    stop(paste0("Invalid palette Must be one of ", paste0(c(rownames(brewer.pal.info), names(ggsci_db)), collapse = ",")))
  }
  if (palette %in% rownames(brewer.pal.info)) {
    pal_n <- brewer.pal.info[palette, "maxcolors"]
    pal_category <- brewer.pal.info[palette, "category"]
    if (pal_category == "div") {
      pal_color <- rev(brewer.pal(name = palette, n = pal_n))
    } else {
      if (palette == "Paired") {
        pal_color <- brewer.pal(12, "Paired")[c(1:4, 7, 8, 5, 6, 9, 10, 11, 12)]
      } else {
        pal_color <- brewer.pal(name = palette, n = pal_n)
      }
    }
  } else {
    pal_color <- ggsci_db[[palette]][[1]]
    pal_n <- length(pal_color)
  }
  if (isTRUE(reverse)) {
    pal_color <- rev(pal_color)
  }

  if (!type %in% c("auto", "discrete", "continuous")) {
    stop("'type' must be one of 'auto','discrete' and 'continuous'.")
  }
  if (type == "auto") {
    if (is.numeric(x)) {
      type <- "continuous"
    } else {
      type <- "discrete"
    }
  }

  if (type == "discrete") {
    if (!is.factor(x)) {
      x <- factor(x, levels = unique(as.character(x)))
    }
    n <- nlevels(x)
    color <- ifelse(rep(n, n) <= pal_n,
      pal_color[1:n],
      colorRampPalette(pal_color)(n)
    )
    names(color) <- levels(x)
    if (isTRUE(matched)) {
      color <- color[x]
    }
  } else {
    if (!is.numeric(x)) {
      stop("'x' must be type of numeric when use continuous color palettes.")
    }
    if (all(is.na(x))) {
      values <- rep(0, 100)
    } else {
      values <- cut(x, breaks = seq(min(x, na.rm = T), max(x, na.rm = T) + 1e-10, length.out = 100), include.lowest = T)
    }
    n <- nlevels(values)
    color <- ifelse(rep(n, n) <= pal_n,
      pal_color[1:n],
      colorRampPalette(pal_color)(n)
    )
    if (isTRUE(matched)) {
      color <- color[values]
    }
  }
  color[is.na(color)] <- NA_color
  return(color)
}

ClassDimPlot <- function(srt, group.by = "orig.ident", split.by = NULL, reduction = NULL,
                         label = FALSE, cells.highlight = NULL, cols.highlight = "black", sizes.highlight = 1,
                         palette = "Paired", legend.position = "right", combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE) {
  library(shadowtext)
  library(dplyr)
  library(ggplot2)
  if (is.null(split.by)) {
    split.by <- "All_cells"
    srt[[split.by]] <- factor("All_cells")
  }
  for (i in c(group.by, split.by)) {
    if (!i %in% colnames(srt@meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (!is.factor(srt[[i, drop = TRUE]])) {
      srt[[i, drop = TRUE]] <- factor(srt[[i, drop = TRUE]], levels = unique(as.character(srt[[i, drop = TRUE]])))
    }
  }
  if (is.null(reduction)) {
    reduction <- Seurat:::DefaultDimReduc(srt)
  }
  if (!reduction %in% Reductions(srt)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }

  dat_meta <- FetchData(object = srt, vars = c(group.by, split.by))
  dat_dim <- Embeddings(srt, reduction = reduction)
  dat_use <- cbind(dat_dim, dat_meta)
  ptsize <- min(1583 / nrow(dat_use), 1)
  plist <- list()
  for (g in group.by) {
    pal_color <- palette_zh(levels(dat_use[[g]]), palette = palette)
    for (s in levels(dat_use[[split.by]])) {
      dat <- dat_use
      cells_mask <- dat[[split.by]] != s
      dat[[g]][cells_mask] <- NA
      labels_tb <- table(dat[[g]])
      labels_tb <- labels_tb[labels_tb != 0]
      dat[, "x"] <- dat[, paste0(reduction, "_1")]
      dat[, "y"] <- dat[, paste0(reduction, "_2")]
      dat[, "g"] <- dat[, g]
      dat[, "Split"] <- s
      dat <- dat[order(dat[, "g"], decreasing = F, na.last = FALSE), ]
      naindex <- which(is.na(dat[, "g"]))
      naindex <- ifelse(length(naindex) > 0, max(naindex), 1)
      dat <- dat[c(1:naindex, sample((naindex + 1):nrow(dat))), ]
      subtitle <- paste0(s, " (nCell:", sum(!is.na(dat[[g]])), ")")
      p <- ggplot(dat, aes(x = x, y = y)) +
        geom_point(aes(color = g), size = ptsize) +
        labs(subtitle = subtitle, x = paste0(reduction, "_1"), y = paste0(reduction, "_2")) +
        scale_x_continuous(limits = c(min(dat_dim[, paste0(reduction, "_1")]), max(dat_dim[, paste0(reduction, "_1")]))) +
        scale_y_continuous(limits = c(min(dat_dim[, paste0(reduction, "_2")]), max(dat_dim[, paste0(reduction, "_2")]))) +
        facet_grid(Split ~ .) +
        theme_zh(aspect.ratio = 1, legend.position = legend.position)
      if (!is.null(cells.highlight)) {
        p$data[, "cells.highlight"] <- (!is.na(p$data[, g])) & rownames(p$data) %in% cells.highlight
        cell_df <- subset(p$data, cells.highlight == TRUE)
        cell_df[, "x"] <- cell_df[, paste0(reduction, "_1")]
        cell_df[, "y"] <- cell_df[, paste0(reduction, "_2")]
        cell_df[, "g"] <- cell_df[, g]
        if (nrow(cell_df) > 0) {
          point_size <- p$layers[[1]]$aes_params$size
          p <- p + geom_point(
            data = cell_df,
            aes(x = x, y = y, fill = g),
            shape = 21, color = cols.highlight, size = sizes.highlight, show.legend = FALSE
          )
        }
      }
      if (isTRUE(label)) {
        p$data[, "g"] <- p$data[, g]
        p$data[, "x"] <- p$data[, paste0(reduction, "_1")]
        p$data[, "y"] <- p$data[, paste0(reduction, "_2")]
        label_df <- p$data %>%
          dplyr::group_by(g) %>%
          dplyr::summarize(x = median(x), y = median(y)) %>%
          as.data.frame()
        colnames(label_df)[1] <- "label"
        label_df <- label_df[!is.na(label_df[, "label"]), ]
        label_df[, "rank"] <- 1:nrow(label_df)
        p <- p + geom_shadowtext(
          data = label_df, aes(x = x, y = y, label = rank),
          color = "white", bg.colour = "black", size = 4.5, inherit.aes = F
        ) + scale_color_manual(
          name = "Groups:",
          values = pal_color[names(labels_tb)],
          labels = paste0(1:length(labels_tb), ": ", names(labels_tb), "(", labels_tb, ")"),
          na.value = "grey90"
        ) + scale_fill_manual(
          name = "Groups:",
          values = pal_color[names(labels_tb)],
          labels = paste0(1:length(labels_tb), ": ", names(labels_tb), "(", labels_tb, ")"),
          na.value = "grey90"
        )
      } else {
        p <- p + scale_color_manual(
          name = "Groups:",
          values = pal_color[names(labels_tb)],
          labels = paste0(names(labels_tb), "(", labels_tb, ")"),
          na.value = "grey90"
        ) + scale_fill_manual(
          name = "Groups:",
          values = pal_color[names(labels_tb)],
          labels = paste0(names(labels_tb), "(", labels_tb, ")"),
          na.value = "grey90"
        )
      }
      p <- p + guides(color = guide_legend(
        title.hjust = 0,
        keywidth = 0,
        keyheight = 0,
        default.unit = "inch",
        override.aes = list(size = 4.5)
      ))
      plist[[paste0(s, ":", g)]] <- p
    }
  }
  if (isTRUE(combine)) {
    return(cowplot::plot_grid(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow, align = "hv", axis = "tblr"))
  } else {
    return(plist)
  }
}

ExpDimPlot <- function(srt, features = NULL, split.by = NULL, reduction = NULL, keep.scale = NULL,
                       cells.highlight = NULL, palette = "Spectral", nrow = NULL, ncol = NULL, byrow = TRUE, combine = TRUE) {
  library(shadowtext)
  if (is.null(split.by)) {
    split.by <- "All_cells"
    srt[[split.by]] <- factor("All_cells")
  }
  for (i in c(split.by)) {
    if (!i %in% colnames(srt@meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (!is.factor(srt[[i, drop = TRUE]])) {
      srt[[i, drop = TRUE]] <- factor(srt[[i, drop = TRUE]], levels = unique(as.character(srt[[i, drop = TRUE]])))
    }
  }
  if (is.null(reduction)) {
    reduction <- Seurat:::DefaultDimReduc(srt)
  }
  if (!reduction %in% Reductions(srt)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }

  nafeatures <- features[!features %in% c(rownames(srt), colnames(srt@meta.data))]
  if (length(nafeatures) > 0) {
    warning(paste0(nafeatures, collapse = ","), " are not in the features of srt.")
    features <- features[features != nafeatures]
  }
  dat_exp <- FetchData(object = srt, vars = features)
  if (!all(sapply(dat_exp, class) == "numeric")) {
    stop("'features' must be type of numeric variable.")
  }
  dat_sp <- FetchData(object = srt, vars = split.by)
  dat_dim <- Embeddings(srt, reduction = reduction)
  dat_use <- cbind(dat_dim, dat_sp)
  ptsize <- min(1583 / nrow(dat_use), 1)
  plist <- list()
  for (f in features) {
    for (s in levels(dat_sp[[split.by]])) {
      dat <- cbind(dat_use, dat_exp[, f, drop = F])
      dat[[f]][dat[, f] == 0] <- NA
      dat[, "x"] <- dat[, paste0(reduction, "_1")]
      dat[, "y"] <- dat[, paste0(reduction, "_2")]
      dat[, "f"] <- dat[, f]
      dat[, "Feature"] <- f
      dat[, "Split"] <- s
      cells_keep <- dat[[split.by]] == s
      dat <- dat[cells_keep, ]
      dat <- dat[order(dat[, "f"], decreasing = F, na.last = FALSE), ]
      colors <- palette_zh(dat[, f], type = "continuous", palette = palette, matched = FALSE)
      if (all(is.na(dat[, f]))) {
        colors_value <- rep(0, 100)
      } else {
        if (is.null(keep.scale)) {
          colors_value <- seq(min(dat[, f], na.rm = T), quantile(dat[, f], 0.99, na.rm = T) + 0.001, length.out = 100)
        } else {
          if (keep.scale == "feature") {
            colors_value <- seq(min(dat_exp[, f], na.rm = T), quantile(dat_exp[, f], 0.99, na.rm = T) + 0.001, length.out = 100)
          }
          if (keep.scale == "all") {
            colors_value <- seq(min(dat_exp[, features], na.rm = T), quantile(dat_exp[, features], 0.99, na.rm = T) + 0.001, length.out = 100)
          }
        }
      }
      dat[which(dat[, "f"] > max(colors_value)), "f"] <- max(colors_value)
      subtitle <- paste0(s, " (nPos:", sum(dat$f > 0, na.rm = T), ", ", round(sum(dat$f > 0, na.rm = T) / nrow(dat) * 100, 2), "%)")
      p <- ggplot(dat, aes(x = x, y = y)) +
        geom_point(aes(color = f), size = ptsize) +
        labs(subtitle = subtitle, x = paste0(reduction, "_1"), y = paste0(reduction, "_2")) +
        scale_x_continuous(limits = c(min(dat_use[, paste0(reduction, "_1")]), max(dat_use[, paste0(reduction, "_1")]))) +
        scale_y_continuous(limits = c(min(dat_use[, paste0(reduction, "_2")]), max(dat_use[, paste0(reduction, "_2")]))) +
        facet_grid(Split ~ Feature) +
        theme_zh(aspect.ratio = 1)

      if (!is.null(cells.highlight)) {
        p$data[, "cells.highlight"] <- rownames(p$data) %in% cells.highlight
        cell_df <- subset(p$data, cells.highlight == TRUE)
        if (nrow(cell_df) > 0) {
          point_size <- p$layers[[1]]$aes_params$size
          p <- p + geom_point(
            data = cell_df,
            aes(x = x, y = y, fill = f),
            shape = 21, color = "red", size = point_size + 1, show.legend = FALSE
          )
        }
      }
      if (all(is.na(dat[, f]))) {
        p <- p + scale_colour_gradient(
          name = "Value", na.value = "grey90",
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
        ) + scale_fill_gradient(
          name = "Value", na.value = "grey90",
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
        )
      } else {
        p <- p + scale_colour_gradientn(
          name = "Value", colours = colors, values = scales::rescale(colors_value), na.value = "grey90",
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
        ) + scale_fill_gradientn(
          name = "Value", colours = colors, values = scales::rescale(colors_value), na.value = "grey90",
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
        )
      }
      plist[[paste0(s, ":", f)]] <- p
    }
  }
  if (isTRUE(combine)) {
    return(cowplot::plot_grid(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow, align = "hv", axis = "tblr"))
  } else {
    return(plist)
  }
}

ExpDotPlot <- function(srt, genes = NULL, gene_groups = NULL, columns = c("orig.ident", "Uncorrectedclusters"), assay = "RNA",
                       cluster_rows = TRUE, cluster_columns = FALSE, grid_size = 0.5) {
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(dplyr)

  if (is.null(genes)) {
    genes <- VariableFeatures(srt)
  }
  genes <- intersect(x = genes, y = rownames(x = srt))
  if (is.null(columns)) {
    columns <- "orig.ident"
  }
  columns <- columns[columns %in% colnames(srt@meta.data)]
  if (length(columns) == 0) {
    stop("Stop plot! 'columns' is invalid!")
  }
  if (is.null(assay)) {
    assay <- DefaultAssay(srt)
  }
  if (is.null(slot)) {
    slot <- "data"
  }

  color_palette <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
  colors <- colorRamp2(seq(-2, 2, length = 100), color_palette)
  legend_point <- Legend(
    labels = paste0(seq(20, 100, length.out = 5), "%"), labels_gp = gpar(fontsize = 12),
    title = "Percent", title_position = "topcenter",
    title_gp = gpar(fontsize = 12, fontfamily = "sans"),
    type = "points", pch = 21,
    size = unit(grid_size * seq(0.2, 1, length.out = 5) + 0.2, "cm"),
    grid_height = unit(grid_size, "cm"),
    grid_width = unit(grid_size, "cm"),
    legend_gp = gpar(fill = "grey90"), border = FALSE,
    background = "transparent"
  )
  legend_color <- Legend(
    col_fun = colors,
    labels_gp = gpar(fontsize = 12),
    title = "Z-score", title_position = "topcenter",
    title_gp = gpar(fontsize = 12, fontfamily = "sans"),
    type = "grid",
    size = unit(grid_size * seq(0.2, 1, length.out = 5) + 0.2, "cm"),
    grid_height = unit(grid_size, "cm"),
    grid_width = unit(grid_size, "cm"),
    border = TRUE,
    background = "transparent"
  )

  dat <- FetchData(object = srt, vars = c(columns, genes), slot = "data")
  dotHT_list <- list()
  for (i in columns) {
    dat_exp <- aggregate(formula(paste0(".~", i)), dat[, c(i, genes)], FUN = function(x) {
      mean(expm1(x))
    }) %>% t()
    colnames(dat_exp) <- dat_exp[i, ]
    dat_exp <- dat_exp[-which(rownames(dat_exp) == i), ]
    dat_exp <- apply(dat_exp, c(1, 2), as.numeric)
    dat_exp <- t(scale(t(dat_exp)))
    dat_exp[is.na(dat_exp)] <- 0

    dat_pec <- aggregate(formula(paste0(".~", i)), dat, FUN = function(x) {
      sum(x > 0) / length(x)
    }) %>% t()
    colnames(dat_pec) <- dat_pec[i, ]
    dat_pec <- dat_pec[-which(rownames(dat_pec) == i), ]
    dat_pec <- apply(dat_pec, c(1, 2), as.numeric)
    dat_pec <- dat_pec[rownames(dat_exp), ]

    if (!is.null(gene_groups)) {
      row_color <- colorRampPalette(brewer.pal(n = 9, name = "Set1"))(length(unique(gene_groups)))
      names(row_color) <- unique(gene_groups)
      row_split <- factor(gene_groups, levels = unique(gene_groups))
      names(row_split) <- row_split
      df_left_annotation <- HeatmapAnnotation(
        foo1 = anno_block(
          gp = gpar(fill = row_color)
        ),
        width = unit(0.5, "cm"), which = "row"
      )
    } else {
      df_left_annotation <- row_split <- NULL
    }
    dotHT_list[[i]] <- Heatmap(dat_exp[genes, ],
      col = colors,
      show_row_names = TRUE,
      row_names_side = "left",
      cluster_columns = cluster_columns,
      cluster_rows = cluster_rows,
      row_split = row_split,
      column_split = NULL,
      cluster_row_slices = FALSE,
      cluster_column_slices = FALSE,
      row_title_rot = 0,
      cell_fun = function(j, i, x, y, w, h, fill) {
        p <- dat_pec[genes, ][i, j]
        grid.rect(x, y,
          width = w, height = h,
          gp = gpar(col = "white", lwd = 1, fill = "white")
        )
        grid.rect(x, y,
          width = w, height = h,
          gp = gpar(col = fill, lwd = 1, fill = alpha(fill, 0.5))
        )
        grid.circle(x, y,
          r = 1 * h * p / 2,
          gp = gpar(col = "black", lwd = 1, fill = fill)
        )
      },
      left_annotation = df_left_annotation,
      width = unit(ncol(dat_exp) * grid_size, "cm"),
      height = unit(length(genes) * grid_size, "cm"),
      show_heatmap_legend = FALSE
    )
  }
  ht_list <- NULL
  for (ht in dotHT_list) {
    ht_list <- ht_list + ht
  }
  grob <- grid.grabExpr({
    ComplexHeatmap::draw(ht_list, heatmap_legend_list = list(legend_color, legend_point))
  })
  p <- as_ggplot(grob)

  return(p)
}

ExpHeatmapPlot <- function(srt, genes = NULL, gene_groups = NULL, samplesize = 100,
                           columns = c("orig.ident", "Uncorrectedclusters"), assay = "RNA") {
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(dplyr)

  if (!is.null(gene_groups)) {
    if (length(gene_groups) != length(genes)) {
      stop("length(gene_groups)!=length(genes)")
    }
  }
  color_palette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(100)
  colors <- colorRamp2(seq(-2, 2, length = 100), color_palette)
  legend_color <- Legend(
    col_fun = colors,
    labels_gp = gpar(fontsize = 12),
    title = "Z-score", title_position = "topcenter",
    title_gp = gpar(fontsize = 12, fontfamily = "sans"),
    type = "grid",
    size = unit(0.5 * seq(0.2, 1, length.out = 5) + 0.2, "cm"),
    grid_height = unit(0.5, "cm"),
    grid_width = unit(0.5, "cm"),
    border = TRUE,
    background = "transparent"
  )
  dotHT_list <- list()
  for (i in columns) {
    if (!is.factor(srt[[i, drop = TRUE]])) {
      srt[[i, drop = TRUE]] <- factor(srt[[i, drop = TRUE]], levels = unique(as.character(srt[[i, drop = TRUE]])))
    }
    cell <- lapply(levels(srt[[i, drop = TRUE]]), function(x) {
      cells <- colnames(srt)[srt[[i, drop = TRUE]] == x]
      size <- ifelse(length(cells) > samplesize, samplesize, length(cells))
      cells_sample <- sample(cells, size)
      out <- setNames(rep(x, size), cells_sample)
      return(out)
    }) %>% unlist(use.names = T)
    cell <- factor(cell, levels = levels(srt[[i, drop = TRUE]]))

    if (nrow(GetAssayData(srt, slot = "scale.data", assay = assay)) != nrow(GetAssayData(srt, slot = "data", assay = assay))) {
      cat("Perform ScaleData on the data...\n")
      srt <- Seurat::ScaleData(object = srt, features = rownames(srt))
    }
    mat <- GetAssayData(srt, slot = "scale.data", assay = assay)[genes, names(cell)]
    df_left_annotation <- HeatmapAnnotation(
      Cell = anno_block(
        gp = gpar(fill = palette_zh(levels(cell))),
        show_name = FALSE
      ),
      width = unit(0.5, "cm"), which = "row", border = TRUE
    )
    df_top_annotation <- HeatmapAnnotation(
      Cell = anno_block(
        gp = gpar(fill = palette_zh(levels(cell), palette = "jco")),
        show_name = FALSE
      ),
      height = unit(0.5, "cm"), which = "column", border = TRUE
    )
    dotHT_list[[i]] <- Heatmap(
      matrix = mat,
      col = colors,
      row_split = gene_groups,
      row_title = table(gene_groups),
      row_title_rot = 0,
      column_split = cell,
      column_title_rot = 90,
      cluster_row_slices = F,
      cluster_column_slices = F,
      show_row_names = F,
      show_column_names = F,
      left_annotation = df_left_annotation,
      top_annotation = df_top_annotation,
      show_heatmap_legend = FALSE
    )
  }
  ht_list <- NULL
  for (ht in dotHT_list) {
    ht_list <- ht_list + ht
  }
  grob <- grid.grabExpr({
    ComplexHeatmap::draw(ht_list, heatmap_legend_list = list(legend_color))
  })
  p <- as_ggplot(grob)
  return(p)
}

SankeyPlot <- function(source, target) {
  library(networkD3)
  library(tidyverse)
  library(plotly)
  data <- cbind(source = source, target = target) %>%
    as.data.frame() %>%
    table() %>%
    as.data.frame.matrix()
  data_long <- data %>%
    rownames_to_column() %>%
    gather(key = "key", value = "value", -rowname) %>%
    filter(value > 0)
  colnames(data_long) <- c("source", "target", "value")

  nodes <- data.frame(name = c(as.character(data_long$source), as.character(data_long$target)) %>% unique())
  data_long$IDsource <- match(data_long$source, nodes$name) - 1
  data_long$IDtarget <- match(data_long$target, nodes$name) - 1

  p <- sankeyNetwork(
    Links = data_long, Nodes = nodes,
    Source = "IDsource", Target = "IDtarget",
    Value = "value", NodeID = "name",
    sinksRight = FALSE, nodeWidth = 40, fontSize = 15, nodePadding = 20
  )
  return(p)
}

SummaryPlot <- function(srt, groups = "orig.ident", clusters = "Standardclusters",
                        features = c("POU5F1", "percent.mito"),
                        reduction = "StandardUMAP2d",
                        Class_palette = "Paired", Exp_palette = "Spectral") {
  for (i in c(clusters, groups)) {
    if (!i %in% colnames(srt@meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (!is.factor(srt[[i, drop = TRUE]])) {
      srt[[i, drop = TRUE]] <- factor(srt[[i, drop = TRUE]], levels = unique(as.character(srt[[i, drop = TRUE]])))
    }
  }
  nafeatures <- features[!features %in% c(rownames(srt), colnames(srt@meta.data))]
  if (length(nafeatures) > 0) {
    warning(paste0(nafeatures, collapse = ","), " are not in the features of srt.")
    features <- features[features != nafeatures]
  }
  class_features <- c()
  exp_features <- c()
  for (i in features) {
    if (i %in% colnames(srt@meta.data)) {
      if (!is.numeric(srt@meta.data[[i]])) {
        class_features <- c(class_features, i)
      } else {
        exp_features <- c(exp_features, i)
      }
    } else {
      exp_features <- c(exp_features, i)
    }
  }
  if (length(features) == 0) {
    stop("No valid features to inspect.")
  }
  if (is.null(reduction)) {
    reduction <- Seurat:::DefaultDimReduc(srt)
  }
  if (!reduction %in% Reductions(srt)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }

  p1 <- ClassDimPlot(srt,
    reduction = reduction,
    group.by = groups,
    palette = Class_palette
  )

  p2 <- ClassDimPlot(srt,
    reduction = reduction,
    group.by = clusters,
    palette = Class_palette,
    label = T
  )
  p_a <- cowplot::plot_grid(plotlist = list(p1, p2), align = "hv", axis = "tblr", nrow = 2, byrow = T)

  p3_list1 <- ClassDimPlot(
    srt = srt,
    reduction = reduction,
    group.by = class_features,
    palette = Class_palette,
    combine = FALSE,
  )
  p3_list2 <- ExpDimPlot(
    srt = srt,
    reduction = reduction,
    features = exp_features,
    palette = Exp_palette,
    combine = FALSE,
  )
  p3_list <- c(p3_list1, p3_list2)

  p_b <- cowplot::plot_grid(plotlist = p3_list, align = "hv", axis = "tblr", nrow = 2, byrow = T)

  ncols <- ceiling(length(features) / 2)
  p_ab <- cowplot::plot_grid(p_a, p_b, align = "hv", axis = "tblr", nrow = 1, byrow = T, rel_widths = c(2, ncols))

  p4_list <- ClassDimPlot(srt,
    reduction = reduction,
    group.by = groups,
    split.by = groups,
    palette = Class_palette,
    combine = F
  )
  p_c <- cowplot::plot_grid(plotlist = p4_list, align = "hv", axis = "tblr", nrow = 1, byrow = T)

  p <- cowplot::plot_grid(p_ab, p_c, align = "hv", axis = "tblr", nrow = 2, byrow = T, rel_heights = c(2 / 3, 1 / 3))
  return(p)
}

ClassDimPlot3D <- function(srt, group.by = "orig.ident", reduction = "StandardUMAP3d", split.by = NULL, plotsize = 600) {
  library(plotly)
  for (i in group.by) {
    if (!i %in% colnames(srt@meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (!is.factor(srt[[i, drop = TRUE]])) {
      srt[[i, drop = TRUE]] <- factor(srt[[i, drop = TRUE]], levels = unique(as.character(srt[[i, drop = TRUE]])))
    }
  }
  if (is.null(reduction)) {
    reduction <- Seurat:::DefaultDimReduc(srt)
  }
  if (!reduction %in% Reductions(srt)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }
  if (ncol(Embeddings(srt, reduction = reduction)) < 3) {
    stop("Reduction must be in three dimensions or higher.")
  }

  df0 <- cbind(Embeddings(srt, reduction = reduction), srt@meta.data)
  if (!is.factor(df0[[group.by]])) {
    df0[[group.by]] <- factor(df0[[group.by]], levels = unique(as.character(df0[[group.by]])))
  }
  df0[["group.by"]] <- df0[[group.by]]

  if (!is.null(split.by)) {
    if (!is.factor(df0[[split.by]])) {
      df0[[split.by]] <- factor(df0[[split.by]], levels = unique(as.character(df0[[split.by]])))
    }
    for (i in levels(df0[[split.by]])) {
      df0[[paste0(reduction, "_1_", i)]] <- ifelse(df0[[split.by]] == i, df0[[paste0(reduction, "_1")]], NA)
      df0[[paste0(reduction, "_2_", i)]] <- ifelse(df0[[split.by]] == i, df0[[paste0(reduction, "_2")]], NA)
      df0[[paste0(reduction, "_3_", i)]] <- ifelse(df0[[split.by]] == i, df0[[paste0(reduction, "_3")]], NA)
    }
  } else {
    split.by <- "All_cells"
    df0[[split.by]] <- factor("All_cells")
  }
  df0[[paste0(reduction, "_1_All_cells")]] <- df0[[paste0(reduction, "_1")]]
  df0[[paste0(reduction, "_2_All_cells")]] <- df0[[paste0(reduction, "_2")]]
  df0[[paste0(reduction, "_3_All_cells")]] <- df0[[paste0(reduction, "_3")]]

  p <- plot_ly(data = df0, width = plotsize * 1.5, height = plotsize)
  i <- 0
  sp <- ifelse(i == 0, "All_cells", levels(df0[[split.by]])[i])
  p <- p %>% add_trace(
    x = df0[[paste0(reduction, "_1_", sp)]],
    y = df0[[paste0(reduction, "_2_", sp)]],
    z = df0[[paste0(reduction, "_3_", sp)]],
    colors = palette_zh(df0[["group.by"]]),
    color = df0[["group.by"]],
    text = paste0("Cell:", rownames(df0), "\ngroup.by:", df0[["group.by"]]),
    type = "scatter3d",
    mode = "markers",
    showlegend = ifelse(i == 0, TRUE, FALSE),
    visible = ifelse(i == 0, TRUE, FALSE),
    marker = list(
      size = 1.5
    )
  )

  split_option <- list()
  for (i in 0:nlevels(df0[[split.by]])) {
    sp <- ifelse(i == 0, "All_cells", levels(df0[[split.by]])[i])
    ncells <- ifelse(i == 0, nrow(df0), table(df0[[split.by]])[sp])
    split_option[[i + 1]] <- list(
      method = "update",
      args = list(
        list(
          x = list(df0[[paste0(reduction, "_1_", sp)]]),
          y = list(df0[[paste0(reduction, "_2_", sp)]]),
          z = list(df0[[paste0(reduction, "_3_", sp)]]),
          marker = list(
            size = 1.5,
            color = palette_zh(df0[["group.by"]])[df0[["group.by"]]]
          )
        ),
        list(title = list(
          text = paste0(sp, " (nCells:", ncells, ")"),
          font = list(size = 16, color = "black"),
          y = 0.95
        ))
      ),
      label = sp
    )
  }
  p <- p %>% plotly::layout(
    title = list(
      text = paste0("Total", " (nCells:", nrow(df0), ")"),
      font = list(size = 16, color = "black"),
      y = 0.95
    ),
    font = list(size = 12, color = "black"),
    showlegend = T,
    legend = list(
      itemsizing = "constant",
      y = 0.5,
      x = 1,
      xanchor = "left"
    ),
    scene = list(
      xaxis = list(title = paste0(reduction, "_1"), range = c(min(df0[[paste0(reduction, "_1")]]), max(df0[[paste0(reduction, "_1")]]))),
      yaxis = list(title = paste0(reduction, "_2"), range = c(min(df0[[paste0(reduction, "_2")]]), max(df0[[paste0(reduction, "_2")]]))),
      zaxis = list(title = paste0(reduction, "_3"), range = c(min(df0[[paste0(reduction, "_3")]]), max(df0[[paste0(reduction, "_3")]]))),
      aspectratio = list(x = 1, y = 1, z = 1)
    ),
    updatemenus = list(
      list(
        y = 0.5,
        buttons = split_option
      )
    ),
    autosize = FALSE
  )

  return(p)
}
