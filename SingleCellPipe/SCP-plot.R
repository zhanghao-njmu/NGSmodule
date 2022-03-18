theme_scp <- function(aspect.ratio = 1, ...) {
  library(ggplot2)
  args1 <- list(
    aspect.ratio = aspect.ratio,
    text = element_text(size = 14, color = "black"),
    title = element_text(size = 14, colour = "black"),
    axis.line = element_blank(),
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 14, colour = "black"),
    strip.text = element_text(size = 12, colour = "black", hjust = 0.5),
    strip.background = element_rect(fill = "transparent", linetype = 0),
    legend.text = element_text(size = 12, colour = "black"),
    legend.title = element_text(size = 14, colour = "black", hjust = 0),
    legend.key = element_rect(fill = "transparent", color = "transparent"),
    legend.background = element_blank(),
    plot.subtitle = element_text(size = 12, hjust = 0, margin = margin(b = 2)),
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

palette_scp <- function(x, palette = "Paired", type = "auto", matched = FALSE, reverse = FALSE, NA_keep = FALSE, NA_color = "grey90") {
  library(RColorBrewer)
  library(ggsci)
  library(Redmonder)
  library(rcartocolor)
  library(nord)
  library(viridis)
  brewer.pal.info <- RColorBrewer::brewer.pal.info
  ggsci_db <- ggsci:::ggsci_db
  redmonder.pal.info <- Redmonder::redmonder.pal.info
  metacartocolors <- rcartocolor::metacartocolors
  rownames(metacartocolors) <- metacartocolors$Name
  nord_palettes <- nord::nord_palettes
  viridis_names <- c("magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo")
  viridis_palettes <- lapply(setNames(viridis_names, viridis_names), function(x) viridis::viridis(100, option = x))

  palette_list <- list()
  all_colors <- c(rownames(brewer.pal.info), names(ggsci_db), rownames(redmonder.pal.info), rownames(metacartocolors), names(nord_palettes), names(viridis_palettes))
  for (pal in all_colors) {
    if (!pal %in% all_colors) {
      stop(paste0("Invalid pal Must be one of ", paste0(all_colors, collapse = ",")))
    }
    if (pal %in% rownames(brewer.pal.info)) {
      pal_n <- brewer.pal.info[pal, "maxcolors"]
      pal_category <- brewer.pal.info[pal, "category"]
      if (pal_category == "div") {
        pal_color <- rev(brewer.pal(name = pal, n = pal_n))
      } else {
        if (pal == "Paired") {
          pal_color <- brewer.pal(12, "Paired")[c(1:4, 7, 8, 5, 6, 9, 10, 11, 12)]
        } else {
          pal_color <- brewer.pal(name = pal, n = pal_n)
        }
      }
    } else if (pal %in% names(ggsci_db)) {
      if (pal %in% c("d3", "uchicago", "material")) {
        for (pal in names(ggsci_db[[pal]])) {
          pal_color <- ggsci_db[[pal]][[pal]]
          pal_n <- length(pal_color)
          palette_list[[paste0(pal, "-", pal)]][["pal_color"]] <- pal_color
          palette_list[[paste0(pal, "-", pal)]][["pal_n"]] <- pal_n
        }
        next
      } else {
        pal_color <- ggsci_db[[pal]][[1]]
        pal_n <- length(pal_color)
      }
    } else if (pal %in% rownames(redmonder.pal.info)) {
      pal_n <- redmonder.pal.info[pal, "maxcolors"]
      pal_category <- redmonder.pal.info[pal, "category"]
      if (pal_category == "div") {
        pal_color <- rev(redmonder.pal(name = pal, n = pal_n))
      } else {
        pal_color <- redmonder.pal(name = pal, n = pal_n)
      }
    } else if (pal %in% rownames(metacartocolors)) {
      pal_n <- metacartocolors[pal, "Max_n"]
      pal_color <- carto_pal(name = pal, n = pal_n)
    } else if (pal %in% names(nord_palettes)) {
      pal_color <- nord_palettes[[pal]]
      pal_n <- length(pal_color)
    } else if (pal %in% names(viridis_palettes)) {
      pal_color <- viridis_palettes[[palette]]
      pal_n <- length(pal_color)
    }
    palette_list[[pal]][["pal_color"]] <- pal_color
    palette_list[[pal]][["pal_n"]] <- pal_n
  }
  # save(palette_list, file = "~/Git/SCP/data/palette_list.rda", version = 2)

  if (!palette %in% names(palette_list)) {
    stop(paste0("Invalid palette Must be one of ", paste0(names(palette_list), collapse = ",")))
  }
  pal_color <- palette_list[[palette]][["pal_color"]]
  pal_n <- palette_list[[palette]][["pal_n"]]

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
    if (any(is.na(x))) {
      color <- c(color, setNames(NA_color, "NA"))
    }
    if (isTRUE(matched)) {
      color <- color[x]
      color[is.na(color)] <- NA_color
    }
  } else if (type == "continuous") {
    if (!is.numeric(x)) {
      stop("'x' must be type of numeric when use continuous color palettes.")
    }
    if (all(is.na(x))) {
      values <- as.factor(rep(0, 100))
    } else if (length(unique(na.omit(x))) == 1) {
      values <- as.factor(rep(unique(na.omit(x)), 100))
    } else {
      values <- cut(x, breaks = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = 101), include.lowest = TRUE)
    }
    n <- nlevels(values)
    color <- ifelse(rep(n, n) <= pal_n,
      pal_color[1:n],
      colorRampPalette(pal_color)(n)
    )
    names(color) <- levels(values)
    if (any(is.na(x))) {
      color <- c(color, setNames(NA_color, "NA"))
    }
    if (isTRUE(matched)) {
      color <- color[as.character(values)]
      color[is.na(color)] <- NA_color
    }
  }

  if (isTRUE(reverse)) {
    color <- rev(color)
  }
  if (!isTRUE(NA_keep)) {
    color <- color[names(color) != "NA"]
  }
  return(color)
}

PanelRatio <- function(p, width, height, units = "px", dpi = 300) {
  library(ggplot2)
  library(gtable)
  library(grid)
  gGrob <- ggplotGrob(p)
  tmpfile <- tempfile(pattern = "png")
  png(tmpfile, width = width, height = height, units = units, res = ifelse(units == "px", NA, dpi))
  plot(p)
  legendSize <- as.numeric(convertWidth(grobWidth(gtable_filter(gGrob, pattern = "guide-box")), unitTo = "inches"))
  panelSize <- as.numeric(convertWidth(grobWidth(gtable_filter(gGrob, pattern = "background")), unitTo = "inches"))
  ratio <- panelSize / (panelSize + legendSize)
  dev.off()
  unlink(tmpfile)
  return(ratio)
}

ClassDimPlot <- function(srt, group.by = "orig.ident", split.by = NULL, reduction = NULL,
                         palette = "Paired", NA_color = "grey90", pt.size = NULL, pt.alpha = 1,
                         label = FALSE, label_insitu = FALSE, label.size = 4,
                         cells.highlight = NULL, cols.highlight = "black", sizes.highlight = 1,
                         legend.position = "right", combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE) {
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
  if (!is.null(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt))) {
      stop("No cells in 'cells.highlight' found in srt.")
    }
    if (!all(cells.highlight %in% colnames(srt))) {
      warning("Some cells in 'cells.highlight' not found in srt.")
    }
  }

  dat_meta <- FetchData(object = srt, vars = c(group.by, split.by))
  dat_dim <- Embeddings(srt, reduction = reduction)
  reduction_key <- Key(srt[[reduction]])
  dat_use <- cbind(dat_dim, dat_meta)
  if (is.null(pt.size)) {
    pt.size <- min(2000 / nrow(dat_use), 0.5)
  }

  plist <- list()
  for (g in group.by) {
    pal_color <- palette_scp(levels(dat_use[[g]]), palette = palette)
    for (s in levels(dat_use[[split.by]])) {
      dat <- dat_use
      cells_mask <- dat[[split.by]] != s
      dat[[g]][cells_mask] <- NA
      labels_tb <- table(dat[[g]])
      labels_tb <- labels_tb[labels_tb != 0]
      dat[, "x"] <- dat[, paste0(paste0(reduction_key, "1"))]
      dat[, "y"] <- dat[, paste0(paste0(reduction_key, "2"))]
      dat[, "g"] <- dat[, g]
      dat[, "GroupBy"] <- g
      dat[, "Split"] <- s
      dat <- dat[order(dat[, "g"], decreasing = FALSE, na.last = FALSE), ]
      naindex <- which(is.na(dat[, "g"]))
      naindex <- ifelse(length(naindex) > 0, max(naindex), 1)
      dat <- dat[c(1:naindex, sample((naindex + 1):nrow(dat))), ]
      subtitle <- paste0(s, " (nCell:", sum(!is.na(dat[[g]])), ")")
      p <- ggplot(dat, aes(x = x, y = y)) +
        geom_point(aes(color = g), size = pt.size, alpha = pt.alpha) +
        labs(subtitle = subtitle, x = paste0(paste0(reduction_key, "1")), y = paste0(paste0(reduction_key, "2"))) +
        scale_x_continuous(limits = c(min(dat_dim[, paste0(paste0(reduction_key, "1"))]), max(dat_dim[, paste0(paste0(reduction_key, "1"))]))) +
        scale_y_continuous(limits = c(min(dat_dim[, paste0(paste0(reduction_key, "2"))]), max(dat_dim[, paste0(paste0(reduction_key, "2"))]))) +
        facet_grid(Split ~ GroupBy) +
        theme_scp(aspect.ratio = 1, legend.position = legend.position)
      if (!is.null(cells.highlight)) {
        p$data[, "cells.highlight"] <- (!is.na(p$data[, g])) & rownames(p$data) %in% cells.highlight
        cell_df <- subset(p$data, cells.highlight == TRUE)
        cell_df[, "x"] <- cell_df[, paste0(paste0(reduction_key, "1"))]
        cell_df[, "y"] <- cell_df[, paste0(paste0(reduction_key, "2"))]
        cell_df[, "g"] <- cell_df[, g]
        if (nrow(cell_df) > 0) {
          # point_size <- p$layers[[1]]$aes_params$size
          p <- p + geom_point(
            data = cell_df,
            aes(x = x, y = y, fill = g),
            shape = 21, color = cols.highlight, size = sizes.highlight, show.legend = FALSE
          )
        }
      }
      p$data[, "g"] <- p$data[, g]
      p$data[, "x"] <- p$data[, paste0(paste0(reduction_key, "1"))]
      p$data[, "y"] <- p$data[, paste0(paste0(reduction_key, "2"))]
      if (isTRUE(label)) {
        label_df <- p$data %>%
          dplyr::group_by(g) %>%
          dplyr::summarize(x = median(x), y = median(y)) %>%
          as.data.frame()
        colnames(label_df)[1] <- "label"
        label_df <- label_df[!is.na(label_df[, "label"]), ]
        label_df[, "rank"] <- 1:nrow(label_df)
        if (isTRUE(label_insitu)) {
          p <- p + geom_shadowtext(
            data = label_df, aes(x = x, y = y, label = label),
            color = "white", bg.colour = "black", size = label.size, inherit.aes = FALSE
          ) + scale_color_manual(
            name = "Groups:",
            values = pal_color[names(labels_tb)],
            labels = paste0(names(labels_tb), "(", labels_tb, ")"),
            na.value = NA_color
          ) + scale_fill_manual(
            name = "Groups:",
            values = pal_color[names(labels_tb)],
            labels = paste0(names(labels_tb), "(", labels_tb, ")"),
            na.value = NA_color
          )
        } else {
          p <- p + geom_shadowtext(
            data = label_df, aes(x = x, y = y, label = rank),
            color = "white", bg.colour = "black", size = label.size, inherit.aes = FALSE
          ) + scale_color_manual(
            name = "Groups:",
            values = pal_color[names(labels_tb)],
            labels = paste0(1:length(labels_tb), ": ", names(labels_tb), "(", labels_tb, ")"),
            na.value = NA_color
          ) + scale_fill_manual(
            name = "Groups:",
            values = pal_color[names(labels_tb)],
            labels = paste0(1:length(labels_tb), ": ", names(labels_tb), "(", labels_tb, ")"),
            na.value = NA_color
          )
        }
      } else {
        p <- p + scale_color_manual(
          name = "Groups:",
          values = pal_color[names(labels_tb)],
          labels = paste0(names(labels_tb), "(", labels_tb, ")"),
          na.value = NA_color
        ) + scale_fill_manual(
          name = "Groups:",
          values = pal_color[names(labels_tb)],
          labels = paste0(names(labels_tb), "(", labels_tb, ")"),
          na.value = NA_color
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
    if (length(plist) > 1) {
      plot <- cowplot::plot_grid(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow, align = "hv", axis = "tblr")
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}

ExpDimPlot <- function(srt, features = NULL, reduction = NULL, split.by = NULL,
                       palette = "Spectral", NA_color = "grey90", pt.size = NULL, pt.alpha = 1, keep.scale = NULL,
                       cells.highlight = NULL, cols.highlight = "black", sizes.highlight = 1,
                       combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE) {
  library(shadowtext)
  if (is.null(features)) {
    stop("'features' must be provided.")
  }
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
  if (!is.null(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt))) {
      stop("No cells in 'cells.highlight' found in srt.")
    }
    if (!all(cells.highlight %in% colnames(srt))) {
      warning("Some cells in 'cells.highlight' not found in srt.")
    }
  }

  nafeatures <- features[!features %in% c(rownames(srt), colnames(srt@meta.data))]
  if (length(nafeatures) > 0) {
    warning(paste0(nafeatures, collapse = ","), " are not in the features of srt.")
    features <- features[features != nafeatures]
  }
  dat_exp <- FetchData(object = srt, vars = features)
  if (!all(sapply(dat_exp, is.numeric))) {
    stop("'features' must be type of numeric variable.")
  }
  dat_sp <- FetchData(object = srt, vars = split.by)
  dat_dim <- Embeddings(srt, reduction = reduction)
  reduction_key <- Key(srt[[reduction]])
  dat_use <- cbind(dat_dim, dat_sp)
  if (is.null(pt.size)) {
    pt.size <- min(2000 / nrow(dat_use), 0.5)
  }

  plist <- list()
  for (f in features) {
    for (s in levels(dat_sp[[split.by]])) {
      dat <- cbind(dat_use, dat_exp[, f, drop = FALSE])
      dat[[f]][dat[, f] == 0] <- NA
      dat[, "x"] <- dat[, paste0(paste0(reduction_key, "1"))]
      dat[, "y"] <- dat[, paste0(paste0(reduction_key, "2"))]
      dat[, "f"] <- dat[, f]
      dat[, "Feature"] <- f
      dat[, "Split"] <- s
      cells_keep <- dat[[split.by]] == s
      dat <- dat[cells_keep, ]
      dat <- dat[order(dat[, "f"], decreasing = FALSE, na.last = FALSE), ]
      colors <- palette_scp(dat[, f], type = "continuous", palette = palette, matched = FALSE)
      if (all(is.na(dat[, f]))) {
        colors_value <- rep(0, 100)
      } else {
        if (is.null(keep.scale)) {
          colors_value <- seq(min(dat[, f], na.rm = TRUE), quantile(dat[, f], 0.99, na.rm = T) + 0.001, length.out = 100)
        } else {
          if (keep.scale == "feature") {
            colors_value <- seq(min(dat_exp[, f], na.rm = TRUE), quantile(dat_exp[, f], 0.99, na.rm = T) + 0.001, length.out = 100)
          }
          if (keep.scale == "all") {
            colors_value <- seq(min(dat_exp[, features], na.rm = T), quantile(dat_exp[, features], 0.99, na.rm = T) + 0.001, length.out = 100)
          }
        }
      }
      dat[which(dat[, "f"] > max(colors_value)), "f"] <- max(colors_value)
      subtitle <- paste0(s, " (nPos:", sum(dat$f > 0, na.rm = T), ", ", round(sum(dat$f > 0, na.rm = T) / nrow(dat) * 100, 2), "%)")
      p <- ggplot(dat, aes(x = x, y = y)) +
        geom_point(aes(color = f), size = pt.size, alpha = pt.alpha) +
        labs(subtitle = subtitle, x = paste0(paste0(reduction_key, "1")), y = paste0(paste0(reduction_key, "2"))) +
        scale_x_continuous(limits = c(min(dat_use[, paste0(paste0(reduction_key, "1"))]), max(dat_use[, paste0(paste0(reduction_key, "1"))]))) +
        scale_y_continuous(limits = c(min(dat_use[, paste0(paste0(reduction_key, "2"))]), max(dat_use[, paste0(paste0(reduction_key, "2"))]))) +
        facet_grid(Split ~ Feature) +
        theme_scp(aspect.ratio = 1)

      if (!is.null(cells.highlight)) {
        p$data[, "cells.highlight"] <- rownames(p$data) %in% cells.highlight
        cell_df <- subset(p$data, cells.highlight == TRUE)
        if (nrow(cell_df) > 0) {
          # point_size <- p$layers[[1]]$aes_params$size
          p <- p + geom_point(
            data = cell_df,
            aes(x = x, y = y, fill = f),
            shape = 21, color = cols.highlight, size = sizes.highlight, show.legend = FALSE
          )
        }
      }
      if (all(is.na(dat[, f]))) {
        p <- p + scale_colour_gradient(
          name = "Value", na.value = NA_color,
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
        ) + scale_fill_gradient(
          name = "Value", na.value = NA_color,
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
        )
      } else {
        p <- p + scale_colour_gradientn(
          name = "Value", colours = colors, values = scales::rescale(colors_value), na.value = NA_color,
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
        ) + scale_fill_gradientn(
          name = "Value", colours = colors, values = scales::rescale(colors_value), na.value = NA_color,
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
        )
      }
      plist[[paste0(s, ":", f)]] <- p
    }
  }
  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- cowplot::plot_grid(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow, align = "hv", axis = "tblr")
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}

ClassDimPlot3D <- function(srt, group.by = "orig.ident", reduction = NULL, split.by = group.by,
                           palette = "Paired", NA_color = "grey90", pt.size = 1.5,
                           cells.highlight = NULL, cols.highlight = "black", shape.highlight = "circle-open", sizes.highlight = 2,
                           width = NULL, height = NULL) {
  library(plotly)
  NA_color <- gplots::col2hex(NA_color)
  cols.highlight <- gplots::col2hex(cols.highlight)

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
  if (ncol(Embeddings(srt, reduction = reduction)) < 3) {
    stop("Reduction must be in three dimensions or higher.")
  }
  if (!is.null(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt))) {
      stop("No cells in 'cells.highlight' found in srt.")
    }
    if (!all(cells.highlight %in% colnames(srt))) {
      warning("Some cells in 'cells.highlight' not found in srt.")
    }
  }

  df0 <- cbind(Embeddings(srt, reduction = reduction), srt@meta.data)
  reduction_key <- Key(srt[[reduction]])
  if (!is.factor(df0[[group.by]])) {
    df0[[group.by]] <- factor(df0[[group.by]], levels = unique(as.character(df0[[group.by]])))
  }
  df0[["group.by"]] <- df0[[group.by]]
  if (!is.null(cells.highlight)) {
    cells.highlight <- cells.highlight[cells.highlight %in% rownames(df0)]
    df0_highlight <- df0[cells.highlight, ]
    df0_highlight[["group.by"]] <- "highlight"
  }
  if (any(is.na(df0[[group.by]]))) {
    n <- as.character(df0[[group.by]])
    n[is.na(n)] <- "NA"
    df0[[group.by]] <- factor(n, levels = c(levels(df0[[group.by]]), "NA"))
  }

  df0[["color"]] <- df0[[group.by]]
  colors <- palette_scp(df0[["group.by"]], palette = palette, NA_color = NA_color)

  df0[[paste0(reduction_key, "1", "All_cells")]] <- df0[[paste0(reduction_key, "1")]]
  df0[[paste0(reduction_key, "2", "All_cells")]] <- df0[[paste0(reduction_key, "2")]]
  df0[[paste0(reduction_key, "3", "All_cells")]] <- df0[[paste0(reduction_key, "3")]]
  for (i in levels(df0[[split.by]])) {
    df0[[paste0(reduction_key, "1", i)]] <- ifelse(df0[[split.by]] == i, df0[[paste0(reduction_key, "1")]], NA)
    df0[[paste0(reduction_key, "2", i)]] <- ifelse(df0[[split.by]] == i, df0[[paste0(reduction_key, "2")]], NA)
    df0[[paste0(reduction_key, "3", i)]] <- ifelse(df0[[split.by]] == i, df0[[paste0(reduction_key, "3")]], NA)
  }
  if (!is.null(cells.highlight)) {
    cells.highlight <- cells.highlight[cells.highlight %in% rownames(df1)]
    df0_highlight <- df0[cells.highlight, ]
    for (i in levels(df0_highlight[[split.by]])) {
      df0_highlight[[paste0(reduction_key, "1", i)]] <- ifelse(df0_highlight[[split.by]] == i, df0_highlight[[paste0(reduction_key, "1")]], NA)
      df0_highlight[[paste0(reduction_key, "2", i)]] <- ifelse(df0_highlight[[split.by]] == i, df0_highlight[[paste0(reduction_key, "2")]], NA)
      df0_highlight[[paste0(reduction_key, "3", i)]] <- ifelse(df0_highlight[[split.by]] == i, df0_highlight[[paste0(reduction_key, "3")]], NA)
    }
  }

  df_legend <- df0[which(!duplicated(df0[["color"]])), ]
  p <- plot_ly(data = df0, width = width, height = height)
  p <- p %>% add_trace(
    x = rep(0, nrow(df_legend)),
    y = rep(0, nrow(df_legend)),
    z = rep(0, nrow(df_legend)),
    type = "scatter3d",
    mode = "lines",
    color = df_legend[["color"]],
    colors = colors,
    line = list(size = 0),
    showlegend = TRUE,
    visible = TRUE
  )
  p <- p %>% add_trace(
    data = df0,
    x = df0[[paste0(reduction_key, "1", "All_cells")]],
    y = df0[[paste0(reduction_key, "2", "All_cells")]],
    z = df0[[paste0(reduction_key, "3", "All_cells")]],
    text = paste0(
      "Cell:", rownames(df0),
      "\ngroup.by:", df0[["group.by"]],
      "\nsplit.by:", df0[[split.by]],
      "\ncolor:", df0[["color"]]
    ),
    type = "scatter3d",
    mode = "markers",
    marker = list(color = colors[df0[["color"]]], size = pt.size),
    name = "All_cells",
    showlegend = FALSE,
    visible = TRUE
  )
  if (!is.null(cells.highlight)) {
    p <- p %>% add_trace(
      x = df0_highlight[[paste0(reduction_key, "1", "All_cells")]],
      y = df0_highlight[[paste0(reduction_key, "2", "All_cells")]],
      z = df0_highlight[[paste0(reduction_key, "3", "All_cells")]],
      text = paste0(
        "Cell:", rownames(df0_highlight),
        "\ngroup.by:", df0_highlight[["group.by"]],
        "\nsplit.by:", df0_highlight[[split.by]],
        "\ncolor:", df0_highlight[["color"]]
      ),
      type = "scatter3d",
      mode = "markers",
      marker = list(size = sizes.highlight, color = cols.highlight, symbol = shape.highlight),
      name = "highlight",
      showlegend = TRUE,
      visible = TRUE
    )
  }

  split_option <- list()
  for (i in 0:nlevels(df0[[split.by]])) {
    sp <- ifelse(i == 0, "All_cells", levels(df0[[split.by]])[i])
    if (i != 0 & sp == "All_cells") {
      next
    }
    ncells <- ifelse(i == 0, nrow(df0), table(df0[[split.by]])[sp])
    x <- c(
      rep(list(rep(0, nrow(df_legend))), nrow(df_legend)),
      list(df0[[paste0(reduction_key, "1", sp)]])
    )
    y <- c(
      rep(list(rep(0, nrow(df_legend))), nrow(df_legend)),
      list(df0[[paste0(reduction_key, "2", sp)]])
    )
    z <- c(
      rep(list(rep(0, nrow(df_legend))), nrow(df_legend)),
      list(df0[[paste0(reduction_key, "3", sp)]])
    )
    name <- c(levels(df_legend[["color"]]), sp)
    if (!is.null(cells.highlight)) {
      x <- c(x, list(df0_highlight[[paste0(reduction_key, "1", sp)]]))
      y <- c(y, list(df0_highlight[[paste0(reduction_key, "2", sp)]]))
      z <- c(z, list(df0_highlight[[paste0(reduction_key, "3", sp)]]))
      name <- c(name, "highlight")
      split_option[[i + 1]] <- list(
        method = "update",
        args = list(list(
          x = x,
          y = y,
          z = z,
          name = name,
          visible = TRUE
        ), list(title = list(
          text = paste0(sp, " (nCells:", ncells, ")"),
          font = list(size = 16, color = "black"),
          y = 0.95
        ))),
        label = sp
      )
    } else {
      split_option[[i + 1]] <- list(
        method = "update",
        args = list(list(
          x = x,
          y = y,
          z = z,
          name = name,
          visible = TRUE
        ), list(title = list(
          text = paste0(sp, " (nCells:", ncells, ")"),
          font = list(size = 16, color = "black"),
          y = 0.95
        ))),
        label = sp
      )
    }
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
      xanchor = "left",
      alpha = 1
    ),
    scene = list(
      xaxis = list(title = paste0(reduction_key, "1"), range = c(min(df0[[paste0(reduction_key, "1")]]), max(df0[[paste0(reduction_key, "1")]]))),
      yaxis = list(title = paste0(reduction_key, "2"), range = c(min(df0[[paste0(reduction_key, "2")]]), max(df0[[paste0(reduction_key, "2")]]))),
      zaxis = list(title = paste0(reduction_key, "3"), range = c(min(df0[[paste0(reduction_key, "3")]]), max(df0[[paste0(reduction_key, "3")]]))),
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

ExpDimPlot3D <- function(srt, features = NULL, reduction = NULL, split.by = NULL,
                         pt.size = 1.5,
                         cells.highlight = NULL, cols.highlight = "black", shape.highlight = "circle-open", sizes.highlight = 2,
                         width = NULL, height = NULL) {
  library(plotly)
  cols.highlight <- gplots::col2hex(cols.highlight)

  if (is.null(split.by)) {
    split.by <- "All_cells"
    srt[[split.by]] <- factor("All_cells")
  }
  for (i in split.by) {
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
  if (!is.null(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt))) {
      stop("No cells in 'cells.highlight' found in srt.")
    }
    if (!all(cells.highlight %in% colnames(srt))) {
      warning("Some cells in 'cells.highlight' not found in srt.")
    }
  }
  nafeatures <- features[!features %in% c(rownames(srt), colnames(srt@meta.data))]
  if (length(nafeatures) > 0) {
    warning(paste0(nafeatures, collapse = ","), " are not in the features of srt.")
    features <- features[features != nafeatures]
  }

  dat_exp <- FetchData(object = srt, vars = features)
  if (!all(sapply(dat_exp, is.numeric))) {
    stop("'features' must be type of numeric variable.")
  }
  dat_sp <- FetchData(object = srt, vars = split.by)
  dat_dim <- Embeddings(srt, reduction = reduction)
  reduction_key <- Key(srt[[reduction]])
  df1 <- cbind(dat_exp, dat_dim, dat_sp)

  df1[[paste0(reduction_key, "1", "All_cells")]] <- df1[[paste0(reduction_key, "1")]]
  df1[[paste0(reduction_key, "2", "All_cells")]] <- df1[[paste0(reduction_key, "2")]]
  df1[[paste0(reduction_key, "3", "All_cells")]] <- df1[[paste0(reduction_key, "3")]]
  for (i in levels(df1[[split.by]])) {
    df1[[paste0(reduction_key, "1", i)]] <- ifelse(df1[[split.by]] == i, df1[[paste0(reduction_key, "1")]], NA)
    df1[[paste0(reduction_key, "2", i)]] <- ifelse(df1[[split.by]] == i, df1[[paste0(reduction_key, "2")]], NA)
    df1[[paste0(reduction_key, "3", i)]] <- ifelse(df1[[split.by]] == i, df1[[paste0(reduction_key, "3")]], NA)
  }
  if (!is.null(cells.highlight)) {
    cells.highlight <- cells.highlight[cells.highlight %in% rownames(df1)]
    df1_highlight <- df1[cells.highlight, ]
    for (i in levels(df1_highlight[[split.by]])) {
      df1_highlight[[paste0(reduction_key, "1", i)]] <- ifelse(df1_highlight[[split.by]] == i, df1_highlight[[paste0(reduction_key, "1")]], NA)
      df1_highlight[[paste0(reduction_key, "2", i)]] <- ifelse(df1_highlight[[split.by]] == i, df1_highlight[[paste0(reduction_key, "2")]], NA)
      df1_highlight[[paste0(reduction_key, "3", i)]] <- ifelse(df1_highlight[[split.by]] == i, df1_highlight[[paste0(reduction_key, "3")]], NA)
    }
  }

  p <- plot_ly(data = df1, width = width, height = height)
  p <- p %>% add_trace(
    data = df1,
    x = df1[[paste0(reduction_key, "1", "All_cells")]],
    y = df1[[paste0(reduction_key, "2", "All_cells")]],
    z = df1[[paste0(reduction_key, "3", "All_cells")]],
    text = paste0(
      "Cell:", rownames(df1),
      "\nExp:", round(df1[[features[1]]], 3),
      "\nsplit.by:", df1[[split.by]]
    ),
    type = "scatter3d",
    mode = "markers",
    marker = list(
      color = df1[[features[1]]],
      colorbar = list(title = list(text = features[1], font = list(color = "black", size = 14)), len = 0.5),
      size = pt.size,
      showscale = T
    ),
    name = "All_cells",
    showlegend = TRUE,
    visible = TRUE
  )
  if (!is.null(cells.highlight)) {
    p <- p %>% add_trace(
      x = df1_highlight[[paste0(reduction_key, "1", "All_cells")]],
      y = df1_highlight[[paste0(reduction_key, "2", "All_cells")]],
      z = df1_highlight[[paste0(reduction_key, "3", "All_cells")]],
      text = paste0(
        "Cell:", rownames(df1_highlight),
        "\nExp:", round(df1_highlight[[features[1]]], 3),
        "\nsplit.by:", df1_highlight[[split.by]]
      ),
      type = "scatter3d",
      mode = "markers",
      marker = list(size = sizes.highlight, color = cols.highlight, symbol = shape.highlight),
      name = "highlight",
      showlegend = TRUE,
      visible = TRUE
    )
  }

  split_option <- list()
  genes_option <- list()
  for (i in 0:nlevels(df1[[split.by]])) {
    sp <- ifelse(i == 0, "All_cells", levels(df1[[split.by]])[i])
    ncells <- ifelse(i == 0, nrow(df1), table(df1[[split.by]])[sp])
    if (i != 0 & sp == "All_cells") {
      next
    }
    x <- list(df1[[paste0(reduction_key, "1", sp)]])
    y <- list(df1[[paste0(reduction_key, "2", sp)]])
    z <- list(df1[[paste0(reduction_key, "3", sp)]])
    name <- sp
    if (!is.null(cells.highlight)) {
      x <- c(x, list(df1_highlight[[paste0(reduction_key, "1", sp)]]))
      y <- c(y, list(df1_highlight[[paste0(reduction_key, "2", sp)]]))
      z <- c(z, list(df1_highlight[[paste0(reduction_key, "3", sp)]]))
      name <- c(sp, "highlight")
    }
    split_option[[i + 1]] <- list(
      method = "update",
      args = list(list(
        x = x,
        y = y,
        z = z,
        name = name,
        visible = TRUE
      ), list(title = list(
        text = paste0(sp, " (nCells:", ncells, ")"),
        font = list(size = 16, color = "black"),
        y = 0.95
      ))),
      label = sp
    )
  }
  for (j in 1:length(features)) {
    marker <- list(
      color = df1[[features[j]]],
      colorbar = list(title = list(text = features[j], font = list(color = "black", size = 14)), len = 0.5),
      size = pt.size,
      showscale = T
    )

    if (!is.null(cells.highlight)) {
      marker <- list(marker, list(size = sizes.highlight, color = cols.highlight, symbol = shape.highlight))
    }
    genes_option[[j]] <- list(
      method = "update",
      args = list(list(
        text = list(paste0(
          "Cell:", rownames(df1),
          "\nExp:", round(df1[[features[j]]], 3),
          "\nsplit.by:", df1[[split.by]]
        )),
        marker = marker
      )),
      label = features[j]
    )
  }
  p <- p %>% plotly::layout(
    title = list(
      text = paste0("All_cells", " (nCells:", nrow(df1), ")"),
      font = list(size = 16, color = "black"),
      y = 0.95
    ),
    showlegend = T,
    legend = list(
      itemsizing = "constant",
      y = -0.2,
      x = 0.5,
      xanchor = "center"
    ),
    scene = list(
      xaxis = list(title = paste0(reduction_key, "1"), range = c(min(df1[[paste0(reduction_key, "1")]]), max(df1[[paste0(reduction_key, "1")]]))),
      yaxis = list(title = paste0(reduction_key, "2"), range = c(min(df1[[paste0(reduction_key, "2")]]), max(df1[[paste0(reduction_key, "2")]]))),
      zaxis = list(title = paste0(reduction_key, "3"), range = c(min(df1[[paste0(reduction_key, "3")]]), max(df1[[paste0(reduction_key, "3")]]))),
      aspectratio = list(x = 1, y = 1, z = 1)
    ),
    updatemenus = list(
      list(
        y = 0.67,
        buttons = split_option
      ),
      list(
        y = 0.33,
        buttons = genes_option
      )
    ),
    autosize = FALSE
  )

  return(p)
}

ExpDotPlot <- function(srt, genes = NULL, gene_groups = NULL, columns = NULL, exp_method = c("zscore", "logfc", "mean"),
                       heatmap_palette = "YlOrRd", group_palette = "Paired", column_palette = "Paired",
                       cluster_rows = TRUE, cluster_columns = FALSE, grid_size = 0.4, assay = "RNA") {
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(dplyr)

  exp_method <- match.arg(exp_method)

  if (is.null(genes)) {
    genes <- VariableFeatures(srt)
  }
  if (is.null(columns)) {
    stop("'columns' must be provided.")
  }
  columns <- columns[columns %in% colnames(srt@meta.data)]
  if (length(columns) == 0) {
    stop("Stop plot! 'columns' is invalid!")
  }
  columns_length <- lapply(srt[[columns, drop = FALSE]], function(x) length(unique(x))) %>% unlist()
  if (any(columns_length < 2) & exp_method == "zscore") {
    stop(paste0("'columns' ", columns[columns_length < 2], " has only one group."))
  }
  if (!is.null(gene_groups)) {
    if (length(gene_groups) != length(genes)) {
      stop("length(gene_groups)!=length(genes)")
    }
  }

  index <- genes %in% rownames(srt)
  genes <- genes[index]
  gene_groups <- gene_groups[index]
  if (length(genes) > 500) {
    stop("Too many genes, suggest reducing the number of genes.")
  }
  if (is.null(assay)) {
    assay <- DefaultAssay(srt)
  }

  dat <- FetchData(object = srt, vars = c(columns, unique(genes)), slot = "data")
  dat <- dat[, c(columns, genes)]
  dotHT_list <- list()
  for (i in columns) {
    dat_exp <- aggregate(formula(paste0(".~", i)), dat[, c(i, genes)], FUN = function(x) {
      mean(expm1(x))
    }) %>% t()
    colnames(dat_exp) <- dat_exp[i, ]
    dat_exp <- dat_exp[-which(rownames(dat_exp) == i), ]
    dat_exp <- apply(dat_exp, c(1, 2), as.numeric)
    if (exp_method == "zscore") {
      dat_exp <- t(scale(t(dat_exp)))
    } else if (exp_method == "logfc") {
      dat_exp <- log2(dat_exp + 1) - log2(apply(dat_exp, 1, mean) + 1)
    } else if (exp_method == "mean") {
      dat_exp <- dat_exp
    }
    dat_exp[is.na(dat_exp)] <- 0

    color_palette <- palette_scp(x = 1:100, palette = heatmap_palette)
    colors <- colorRamp2(seq(quantile(dat_exp, 0.01), quantile(dat_exp, 0.99), length = 100), color_palette)
    legend_color <- Legend(
      col_fun = colors,
      labels_gp = gpar(fontsize = 12),
      title = exp_method, title_position = "topleft",
      title_gp = gpar(fontsize = 12, fontfamily = "sans"),
      type = "grid",
      border = TRUE,
      background = "transparent"
    )
    legend_point <- Legend(
      labels = paste0(seq(20, 100, length.out = 5), "%"), labels_gp = gpar(fontsize = 12),
      title = "Percent", title_position = "topleft",
      title_gp = gpar(fontsize = 12, fontfamily = "sans"),
      type = "points", pch = 21,
      size = unit(pi * grid_size^2 * seq(0.2, 1, length.out = 5), "cm"),
      grid_height = unit(grid_size, "cm"),
      grid_width = unit(grid_size, "cm"),
      legend_gp = gpar(fill = "grey90"), border = FALSE,
      background = "transparent"
    )

    dat_pec <- aggregate(formula(paste0(".~", i)), dat, FUN = function(x) {
      sum(x > 0) / length(x)
    }) %>% t()
    colnames(dat_pec) <- dat_pec[i, ]
    dat_pec <- dat_pec[-which(rownames(dat_pec) == i), ]
    dat_pec <- apply(dat_pec, c(1, 2), as.numeric)
    dat_pec <- dat_pec[rownames(dat_exp), ]
    assign(paste0(i, ".dat_pec"), dat_pec)

    if (!is.null(gene_groups)) {
      if (!is.factor(gene_groups)) {
        gene_groups <- factor(gene_groups, levels = unique(gene_groups))
      }
      row_color <- palette_scp(levels(gene_groups), palette = group_palette, matched = T)
      row_split <- gene_groups
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
    df_top_annotation <- HeatmapAnnotation(
      Cell = anno_simple(
        x = colnames(dat_exp),
        col = palette_scp(colnames(dat_exp), palette = column_palette),
        gp = gpar(col = "black"),
      ),
      height = unit(0.5, "cm"), which = "column", show_annotation_name = FALSE
    )
    use_raster <- length(genes) > 2000
    funbody <- paste0("
      p <- ", i, ".dat_pec", "[i, j];
      grid.rect(x, y,
                width = w, height = h,
                gp = gpar(col = 'white', lwd = 1, fill = 'white')
      );
      grid.rect(x, y,
                width = w, height = h,
                gp = gpar(col = fill, lwd = 1, fill = alpha(fill, 0.5))
      );
      grid.circle(x, y,
                  r = h * p / 2,
                  gp = gpar(col = 'black', lwd = 1, fill = fill)
      );
    ")
    funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
    eval(parse(text = paste("cell_fun <- function(j, i, x, y, w, h, fill) {", funbody, "}", sep = "")))

    dotHT_list[[i]] <- Heatmap(dat_exp,
      col = colors,
      show_row_names = TRUE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 10),
      column_names_side = "top",
      column_names_rot = 90,
      cluster_columns = cluster_columns,
      cluster_rows = cluster_rows,
      row_split = row_split,
      column_split = NULL,
      cluster_row_slices = FALSE,
      cluster_column_slices = FALSE,
      row_title_rot = 0,
      cell_fun = cell_fun,
      left_annotation = df_left_annotation,
      top_annotation = df_top_annotation,
      width = unit(ncol(dat_exp) * grid_size, "cm"),
      height = unit(length(genes) * grid_size, "cm"),
      show_heatmap_legend = FALSE,
      use_raster = use_raster
    )
  }
  ht_list <- NULL
  for (ht in dotHT_list) {
    ht_list <- ht_list + ht
  }
  grob <- grid.grabExpr({
    ComplexHeatmap::draw(ht_list,
      heatmap_legend_list = list(legend_color, legend_point),
      padding = unit(c(1, max(c(nchar(genes) * 0.13, 1)), max(c(nchar(as.character(srt[[columns, drop = T]])) * 0.13, 1)), 1), "cm")
    ) # bottom, left, top and right
  })
  p <- as_ggplot(grob)

  return(p)
}

ExpHeatmapPlot <- function(srt, genes = NULL, gene_groups = NULL, cell_size = 100, columns = "orig.ident",
                           heatmap_palette = "Spectral", group_palette = "Paired", column_palette = "Paired", assay = "RNA") {
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
  if (!is.null(gene_groups)) {
    if (length(gene_groups) != length(genes)) {
      stop("length(gene_groups)!=length(genes)")
    }
  }
  index <- genes %in% rownames(srt)
  genes <- genes[index]
  gene_groups <- gene_groups[index]
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

  color_palette <- palette_scp(x = 1:100, palette = heatmap_palette)
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
      size <- ifelse(length(cells) > cell_size, cell_size, length(cells))
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
        gp = gpar(fill = palette_scp(levels(cell), palette = group_palette)),
        show_name = FALSE
      ),
      width = unit(0.5, "cm"), which = "row", border = TRUE
    )
    df_top_annotation <- HeatmapAnnotation(
      Cell = anno_block(
        gp = gpar(fill = palette_scp(levels(cell), palette = column_palette)),
        show_name = FALSE
      ),
      height = unit(0.5, "cm"), which = "column", border = TRUE
    )
    use_raster <- nrow(mat) > 2000
    dotHT_list[[i]] <- Heatmap(
      matrix = mat,
      col = colors,
      row_split = gene_groups,
      row_title = table(gene_groups)[table(gene_groups) != 0],
      row_title_rot = 0,
      column_split = cell,
      column_title_rot = 90,
      cluster_row_slices = FALSE,
      cluster_column_slices = FALSE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      left_annotation = df_left_annotation,
      top_annotation = df_top_annotation,
      show_heatmap_legend = FALSE,
      use_raster = use_raster
    )
  }
  ht_list <- NULL
  for (ht in dotHT_list) {
    ht_list <- ht_list + ht
  }
  grob <- grid.grabExpr({
    ComplexHeatmap::draw(ht_list,
      heatmap_legend_list = list(legend_color),
      padding = unit(c(1, max(c(nchar(genes) * 0.13, 1)), max(c(nchar(as.character(srt[[columns, drop = T]])) * 0.13, 1)), 1), "cm")
    ) # bottom, left, top and right
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

SummaryPlot <- function(srt, groups = "orig.ident", clusters = "Standardpcaclusters",
                        features = c("POU5F1", "percent.mito"),
                        reduction = NULL,
                        Class_palette = "Paired", Exp_palette = "Spectral") {
  if (is.null(clusters) | is.null(groups)) {
    stop("'groups' and 'clusters' must be provided.")
  }
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
    features <- c("nCount_RNA", "nFeature_RNA")
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
    palette = Class_palette,
    combine = FALSE
  )

  p2 <- ClassDimPlot(srt,
    reduction = reduction,
    group.by = clusters,
    palette = Class_palette,
    label = T,
    combine = FALSE
  )
  p_a <- cowplot::plot_grid(plotlist = c(p1, p2), align = "hv", axis = "tblr", nrow = 2, byrow = T)

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
    combine = FALSE
  )
  p_c <- cowplot::plot_grid(plotlist = p4_list, align = "hv", axis = "tblr", nrow = 1, byrow = T)

  p <- cowplot::plot_grid(p_ab, p_c, align = "hv", axis = "tblr", nrow = 2, byrow = T, rel_heights = c(2 / 3, 1 / 3))
  return(p)
}
