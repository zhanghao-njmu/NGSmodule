suppressWarnings(suppressPackageStartupMessages(invisible(lapply(
  c(
    "BiocParallel", "edgeR", "DESeq2", "dplyr", "stringr", "scales", "RColorBrewer",
    "ggpubr", "ggsci", "ggforce", "reshape2", "VennDiagram", "gridExtra",
    "gplots", "openxlsx", "ggalluvial", "ggfittext", "ComplexHeatmap", "circlize", "nord"
  ), require,
  character.only = TRUE
))))

##### plot function #####
theme_use <- pal_nejm(palette = "default")(8)

count_bar <- function(lib_size, label) {
  df <- data.frame(x = names(lib_size), y = lib_size)
  p <- ggplot(data = df, aes(x = x, y = y)) +
    geom_bar(aes(fill = y), color = "black", stat = "identity") +
    geom_text(aes(label = round(y, digits = 2)), vjust = -0.5, size = 4) +
    ylim(0, ceiling(max(df[, "y"])) * 1.05) +
    scale_fill_viridis_c(guide = FALSE) +
    labs(title = paste0("Library size (", label, ")"), x = "", y = "Number of reads/fragments\n(Millions)") +
    theme_pubr(border = T) +
    theme(
      aspect.ratio = 5 / length(unique(df[, "x"])),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major = element_line()
    )
  return(p)
}

var_mean_plot <- function(count_matrix, group, label) {
  cpm_matrix <- as.data.frame(count_matrix * (colSums(count_matrix) / mean(colSums(count_matrix))))

  if (ncol(cpm_matrix) == 2) {
    cpm_mean_group <- as.data.frame(apply(cpm_matrix, 1, mean), stringsAsFactors = F)
    colnames(cpm_mean_group) <- "Merge"
    cpm_mean_group[, "id"] <- rownames(cpm_mean_group)

    cpm_var_group <- as.data.frame(apply(cpm_matrix, 1, var), stringsAsFactors = F)
    colnames(cpm_var_group) <- "Merge"
    cpm_var_group[, "id"] <- rownames(cpm_var_group)
    disp <- estimateCommonDisp(count_matrix, group = rep("Merge", 2))
  } else {
    cpm_mean_group <- t(apply(cpm_matrix, 1, function(x) {
      tapply(as.numeric(x), group, mean)
    }))
    cpm_var_group <- t(apply(cpm_matrix, 1, function(x) {
      tapply(as.numeric(x), group, var)
    }))
    disp <- estimateCommonDisp(count_matrix, group = group)
  }
  cpm_mean_group <- suppressMessages(reshape2::melt(cpm_mean_group))
  colnames(cpm_mean_group) <- c("id", "group", "mean")
  cpm_var_group <- suppressMessages(reshape2::melt(cpm_var_group))
  colnames(cpm_var_group) <- c("id", "group", "var")


  df <- merge(x = cpm_mean_group, y = cpm_var_group)
  df <- df[which(df$mean > 0), ]
  df$group <- sapply(df$group, function(x) {
    switch(as.character(x),
      "Merge" = "Merge",
      "Group1" = str_split(label, "/")[[1]][1],
      "Group2" = str_split(label, "/")[[1]][2]
    )
  })


  p <- ggplot(data = df, aes(x = log10(mean), y = log10(var))) +
    stat_density2d(geom = "tile", aes(fill = ..density..^0.25, alpha = 1), contour = FALSE, show.legend = F) +
    stat_density2d(geom = "tile", aes(fill = ..density..^0.25, alpha = ifelse(..density..^0.25 < 0.4, 0, 1)), contour = FALSE, show.legend = F) +
    geom_point(aes(color = group), shape = 16, size = 0.5, alpha = 0.4) +
    geom_line(data = df, aes(x = log10(mean), y = log10(mean + disp * mean^2)), color = "red", alpha = 0.8, size = 1) +
    geom_abline(intercept = 0, slope = 1, color = "black", alpha = 0.8, size = 1) +
    annotate(
      geom = "text", x = 0, y = Inf, hjust = 0, vjust = 3,
      label = "Black: Poisson distribution", color = "black"
    ) +
    annotate(
      geom = "text", x = 0, y = Inf, hjust = 0, vjust = 5,
      label = "Red: Negative binomial distribution", color = "red"
    ) +
    scale_color_viridis_d(name = "") +
    scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256)) +
    labs(title = paste0("Var-Mean scatter (", label, ")"), x = "log10(mean(CPM))", y = "log10(var(CPM))") +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    # facet_grid(.~group)+
    theme_pubr(border = T) +
    theme(aspect.ratio = 0.6)
  return(p)
}

myClip <- function(x, min, max) {
  ifelse(x <= min, min, ifelse(x >= max, max, x))
}

MDplot <- function(count_matrix, group, up, down, nonsig, drop, label) {
  cpm_matrix <- as.data.frame(cpm(count_matrix))
  cpm_rowmean <- rowMeans(cpm_matrix)
  cpm_mean_group <- t(apply(cpm_matrix, 1, function(x) {
    tapply(as.numeric(x), group, mean)
  }))
  log2fc <- log2(cpm_mean_group[, "Group1"] / cpm_mean_group[, "Group2"])
  df <- data.frame(cpm = cpm_rowmean, log2fc = log2fc)
  df[drop, "group"] <- paste0("Excluded(", length(drop), ")")
  df[up, "group"] <- paste0("Up-regulated(", length(up), ")")
  df[down, "group"] <- paste0("Down-regulated(", length(down), ")")
  df[, "group"] <- factor(df[, "group"], levels = c(
    paste0("Excluded(", length(drop), ")"),
    paste0("Up-regulated(", length(up), ")"),
    paste0("Down-regulated(", length(down), ")")
  ))
  colors <- c("grey", theme_use[1], theme_use[2])
  names(colors) <- c(
    paste0("Excluded(", length(drop), ")"),
    paste0("Up-regulated(", length(up), ")"),
    paste0("Down-regulated(", length(down), ")")
  )
  df2 <- df[which(!is.na(df$group)), ]
  p <- ggplot(data = df, aes(x = cpm, y = myClip(log2fc, min = -2, max = 2))) +
    geom_point(aes(color = log2fc), alpha = 0.2, size = 0.5) +
    geom_point(
      data = df2, aes(x = cpm, y = myClip(log2fc, min = -2, max = 2), fill = group),
      color = "transparent", shape = 21, alpha = 0.5, size = 1, inherit.aes = F
    ) +
    geom_hline(yintercept = 0, alpha = 0.5, size = 2, color = "red") +
    geom_hline(yintercept = c(1, -1), alpha = 1, linetype = 2, size = 1) +
    labs(title = paste0("MD plot (", label, ")"), x = "Average CPM", y = "log2FoldChange") +
    scale_colour_viridis_c(guide = F) +
    scale_fill_manual(values = colors, name = "") +
    scale_y_continuous(breaks = seq(-2, 2, 1)) +
    scale_x_continuous(
      trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    guides(fill = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_pubr(border = T) +
    theme(
      aspect.ratio = 0.6,
      panel.grid.major = element_line()
    )
  return(p)
}

pvalue_distibution <- function(pvalue, padj, label) {
  colors <- c("p-value" = theme_use[2], "p.adjust" = theme_use[1])
  p <- ggplot(data = data.frame(pvalue = pvalue, padj = padj)) +
    geom_histogram(aes(x = pvalue, fill = "p-value"), binwidth = 0.01, color = "black", alpha = 0.6, boundary = 0) +
    geom_histogram(aes(x = padj, fill = "p.adjust"), binwidth = 0.01, color = "black", alpha = 0.6, boundary = 0) +
    geom_vline(xintercept = 0.05, linetype = 2, size = 1) +
    scale_fill_manual(values = colors) +
    labs(title = paste0("Probability distribution (", label, ")"), x = "Probability", y = "Frequency") +
    facet_zoom(xlim = c(0, 0.05), zoom.size = 1) +
    theme_pubr(border = T) +
    theme(
      aspect.ratio = 0.6,
      plot.margin = margin(t = 20, r = 100, b = 20, l = 100),
      panel.grid.major = element_line(),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent")
    )
  return(p)
}

filterNumRej_plot <- function(dds_res, max_padj, count_matrix, label) {
  filterNumRej <- dds_res@metadata$filterNumRej
  filterNumRej[, "group"] <- "Actural"
  fit <- data.frame(dds_res@metadata$lo.fit)
  fit[, "group"] <- "Fitted"
  colnames(filterNumRej) <- colnames(fit) <- c("x", "y", "group")
  df <- rbind(filterNumRej, fit)

  remain_dds <- na.omit(dds_res)
  filterThreshold <- dds_res@metadata$filterThreshold
  labels <- paste0(
    "Quantile: ", gsub(x = names(filterThreshold), perl = T, replacement = "", pattern = "(?<=\\.\\d\\d).*(?=%)"), "\n",
    "baseMean: ", round(filterThreshold, 2), "\n",
    "average CPM: ", round(min(rowMeans(cpm(count_matrix))[rownames(remain_dds)]), 2), "\n",
    "average Counts: ", round(min(rowMeans(count_matrix)[rownames(remain_dds)]), 2), "\n",
    "Number of Rejectons: ", filterNumRej[names(filterThreshold), "y"]
  )

  linetype <- c("Actural" = 1, "Fitted" = 2)
  color <- c("Actural" = "black", "Fitted" = "red")
  p <- ggplot(data = df, aes(x = x, y = y, color = group)) +
    geom_segment(
      x = filterNumRej[names(filterThreshold), "x"], y = -Inf,
      xend = filterNumRej[names(filterThreshold), "x"], yend = filterNumRej[names(filterThreshold), "y"],
      linetype = 3, inherit.aes = F
    ) +
    geom_segment(
      x = -Inf, y = filterNumRej[names(filterThreshold), "y"],
      xend = filterNumRej[names(filterThreshold), "x"], yend = filterNumRej[names(filterThreshold), "y"],
      linetype = 3, inherit.aes = F
    ) +
    annotate("text",
      x = filterNumRej[names(filterThreshold), "x"] + 0.01, y = filterNumRej[names(filterThreshold), "y"] / 10,
      label = labels, vjust = 0, hjust = 0
    ) +
    geom_point(alpha = 0.8) +
    geom_line(aes(linetype = group)) +
    labs(title = paste0("Independent Filtering (", label, ")"), x = "Quantiles of normalized counts used for filering", y = paste0("Number of rejections (padj<", max_padj, ")")) +
    scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, max(filterNumRej[, "y"]))) +
    scale_color_manual(values = color) +
    scale_linetype_manual(values = linetype) +
    theme_pubr(border = T, legend = c(0.85, 0.85)) +
    theme(
      aspect.ratio = 0.8,
      panel.grid.major = element_line(),
      legend.title = element_blank(),
      legend.text = element_text(colour = "black", size = 12),
      legend.background = element_rect(fill = "transparent")
    )

  return(p)
}

volcano_plot <- function(count_matrix, res, col_padj, col_log2fc, col_DEgroup, DEgroup_level, max_padj, min_fc, label) {
  cpm_matrix <- as.data.frame(cpm(count_matrix))
  cpm_cov <- apply(cpm_matrix, 1, function(x) {
    sd(x) / mean(x)
  })
  df <- res
  rownames(df) <- res[, 1]
  colnames(df)[colnames(df) == col_padj] <- "padj"
  colnames(df)[colnames(df) == col_log2fc] <- "log2fc"
  colnames(df)[colnames(df) == col_DEgroup] <- "DEgroup"
  df <- df[, c("padj", "log2fc", "DEgroup")]
  df[, "cov"] <- cpm_cov[rownames(df)]

  df[, "DEgroup"] <- factor(df[, "DEgroup"], levels = DEgroup_level)
  df_table <- table(df[, "DEgroup"])
  df_label <- paste0(names(df_table), ": ", df_table)
  maxFC <- df %>%
    filter(.data[["log2fc"]] != Inf & .data[["log2fc"]] != (-Inf)) %>%
    pull("log2fc") %>%
    abs() %>%
    max() %>%
    ceiling()

  y_max <- ceiling(max(min(max(-log10(df$padj)), 8), 4))
  br <- sort(unique(c(seq(0, y_max, 1), -log10(max_padj))))
  br_label <- br
  br_label[which(br_label == (-log10(max_padj)))] <- paste0("padj=", max_padj)
  col_y <- ifelse(br == (-log10(max_padj)), "red", "black")

  p <- ggplot(data = df, aes(x = log2fc, y = -log10(padj), group = DEgroup)) +
    geom_point(aes(fill = DEgroup, alpha = DEgroup, size = cov), shape = 21, position = "jitter") +
    geom_vline(xintercept = log2(min_fc), linetype = "longdash", color = theme_use[1], size = 1, show.legend = F) +
    geom_vline(xintercept = -log2(min_fc), linetype = "longdash", color = theme_use[2], size = 1, show.legend = F) +
    geom_hline(yintercept = -log10(max_padj), linetype = "longdash", color = "black", size = 1, show.legend = F) +
    geom_rug(aes(color = DEgroup), show.legend = F) +
    scale_x_continuous(breaks = seq(-maxFC, maxFC, 1), limits = c(-maxFC, maxFC)) +
    scale_y_continuous(breaks = br, limits = c(0, y_max), labels = br_label) +
    scale_alpha_manual(breaks = DEgroup_level, values = c(1, 1, 0.5), guide = F) +
    scale_color_manual(breaks = DEgroup_level, values = c(theme_use[1], theme_use[2], "grey80")) +
    scale_fill_manual(
      breaks = DEgroup_level, values = c(theme_use[1], theme_use[2], "grey80"),
      labels = df_label
    ) +
    scale_size_continuous(range = c(0.2, 4)) +
    labs(x = "log2fc", y = paste0("-log10(padj)"), title = paste0("Volcano plot (", label, ")")) +
    guides(
      fill = guide_legend(title = paste("DGEs:", sum(df_table[1:2])), override.aes = list(size = 3), order = 1),
      size = guide_legend(title = "Coefficient of variation", override.aes = list(color = "black", shape = 16))
    ) +
    theme_classic2() +
    theme(
      aspect.ratio = 0.9,
      panel.border = element_rect(color = "black", size = 1, fill = "transparent"),
      panel.grid.major = element_line(size = 0.5, colour = "grey80", linetype = 2),
      line = element_line(size = 1),
      axis.line = element_blank(),
      axis.text = element_text(colour = "black", size = 12),
      axis.title = element_text(colour = "black", size = 14),
      axis.text.y = element_text(colour = col_y),
      legend.title = element_text(colour = "black", size = 12),
      legend.text = element_text(colour = "black", size = 12)
    )
  return(p)
}

heatmap <- function(LogNormCounts, genes, res, res_log2fc_col, res_padj_col, row_split, column_split, label,
                    cluster_row_slices = FALSE, cluster_column_slices = FALSE,
                    row_title = " ", column_title = " ",
                    anno_genebiotype = TRUE, anno_zoom_box = FALSE,
                    return_ggplot = TRUE,
                    ...) {
  df <- as.matrix(LogNormCounts[genes, ])
  df <- t(scale(t(df)))
  q01 <- quantile(unlist(df), 0.01)
  q99 <- quantile(unlist(df), 0.99)
  lim <- max(abs(c(q01, q99)))

  row_baseExp <- rowMeans(LogNormCounts[genes, ])
  row_log2fc <- res[genes, res_log2fc_col]
  row_padj <- (-log10(res[genes, res_padj_col]))
  row_padj[row_padj > 10] <- 10

  row_baseExp_color <- colorRampPalette(brewer.pal(n = 9, name = "Purples")[3:9])(length(row_baseExp))[rank(row_baseExp)]
  row_log2fc_color <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlGn")))(length(row_log2fc))[rank(row_log2fc)]
  row_padj_color <- colorRampPalette(brewer.pal(n = 9, name = "Blues")[3:9])(length(row_padj))[rank(row_padj)]

  color_palette <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(100))
  col_color <- if (length(unique(column_split)) == 2) {
    pal_jama("default")(2)
  } else {
    colorRampPalette(pal_jco("default")(length(unique(column_split))))(length(unique(column_split)))
  }
  col_color <- setNames(object = col_color, levels(column_split)[levels(column_split) %in% unique(column_split)])
  row_color <- if (length(unique(row_split)) == 2) {
    pal_nejm(palette = "default")(8)[1:2]
  } else {
    viridis_pal(option = "E")(length(unique(row_split)))
  }
  row_color <- setNames(object = row_color, levels(row_split)[levels(row_split) %in% unique(row_split)])

  ht_legend <- list(
    title = paste0("Z-score"),
    title_gp = gpar(fontsize = 8, fontfamily = "sans"),
    title_position = "topcenter",
    grid_height = unit(3, "mm"),
    grid_width = unit(3, "mm"),
    border = "black",
    labels_gp = gpar(fontsize = 8),
    legend_direction = "vertical",
    legend_width = unit(3, "cm")
  )


  df_top_annotation <- HeatmapAnnotation(
    foo1 = anno_block(
      gp = gpar(fill = col_color),
      labels = levels(column_split)[levels(column_split) %in% unique(column_split)],
      labels_gp = gpar(col = "white", fontsize = 8)
    ),
    height = unit(0.5, "cm")
  )
  if (length(row_split) > 0) {
    df_left_annotation <- HeatmapAnnotation(
      foo1 = anno_block(
        gp = gpar(fill = row_color),
        labels = levels(row_split)[levels(row_split) %in% unique(row_split)],
        labels_gp = gpar(col = "white", fontsize = 8)
      ),
      width = unit(0.5, "cm"), which = "row"
    )
  } else {
    df_left_annotation <- NULL
  }

  df_right_annotation <- HeatmapAnnotation(
    "baseExp" = anno_barplot(row_baseExp, gp = gpar(fill = row_baseExp_color, col = 0), baseline = 0, border = T, width = unit(1, "cm"), axis_param = list(gp = gpar(fontsize = 7))),
    "log2fc" = anno_barplot(row_log2fc, gp = gpar(fill = row_log2fc_color, col = 0), baseline = 0, border = T, width = unit(1, "cm"), axis_param = list(gp = gpar(fontsize = 7))),
    "log10padj" = anno_barplot(row_padj, gp = gpar(fill = row_padj_color, col = 0), baseline = 0, border = T, width = unit(1, "cm"), axis_param = list(gp = gpar(fontsize = 7))),
    show_annotation_name = FALSE,
    gap = unit(0, "points"), which = "row"
  )

  if (anno_genebiotype) {
    protcode <- res[genes, "gene_biotype"] == "protein_coding"
    lincRNA <- res[genes, "gene_biotype"] == "lincRNA"
    miRNA <- res[genes, "gene_biotype"] == "miRNA"
    pseudogene <- str_detect(string = res[genes, "gene_biotype"], pattern = "pseudogene")
    antisense <- res[genes, "gene_biotype"] == "antisense"
    genebiotype <- HeatmapAnnotation(
      "protein_coding" = anno_simple(protcode + 0, col = c("0" = "white", "1" = nord("polarnight")[1]), width = unit(0.5, "cm")),
      "lincRNA" = anno_simple(lincRNA + 0, col = c("0" = "white", "1" = nord("polarnight")[4]), width = unit(0.5, "cm")),
      "miRNA" = anno_simple(miRNA + 0, col = c("0" = "white", "1" = nord("mountain_forms")[1]), width = unit(0.5, "cm")),
      "pseudogene" = anno_simple(pseudogene + 0, col = c("0" = "white", "1" = nord("mountain_forms")[2]), width = unit(0.5, "cm")),
      "antisense" = anno_simple(antisense + 0, col = c("0" = "white", "1" = nord("mountain_forms")[4]), width = unit(0.5, "cm")),
      show_annotation_name = FALSE,
      gap = unit(0, "points"), which = "row"
    )
    df_right_annotation <- c(genebiotype, df_right_annotation)
  }
  if (anno_zoom_box & !is.null(row_split) & !is.null(column_split)) {
    panel_fun <- function(index, nm) {
      if (all(!index %in% which(rownames(df) %in%
        names(row_split)[row_split == tail(levels(row_split), 1)]))) {
        axis_text_x <- element_blank()
      } else {
        axis_text_x <- element_text(size = 6, angle = 45, hjust = 1, vjust = 1)
      }
      d <- melt(df[index, , drop = FALSE], varnames = c("row_index", "col_index"))
      d[, "group"] <- column_split[d[, "col_index"]]
      g <- ggplot(d, aes(x = group, y = value, fill = group)) +
        geom_boxplot(outlier.size = 0, outlier.alpha = 0, show.legend = FALSE) +
        stat_summary(fun = median, geom = "line", size = 0.8, aes(group = 1), show.legend = FALSE) +
        stat_summary(fun = median, geom = "point", size = 1, shape = 20, color = "white", show.legend = FALSE) +
        labs(x = "", y = "Z-score") +
        scale_fill_manual(values = col_color, guide = FALSE) +
        theme_pubr(border = T) +
        theme(
          aspect.ratio = 5 / length(unique(column_split)),
          axis.title = element_text(size = 8),
          axis.text.x = axis_text_x,
          axis.text.y = element_text(size = 6),
          panel.grid.major = element_line()
        )
      g <- grid.grabExpr(print(g))
      pushViewport(viewport())
      grid.rect()
      grid.draw(g)
      popViewport()
    }
    zoom_box <- HeatmapAnnotation(
      "anno_zoom_box" = anno_zoom(
        align_to = row_split, which = "row", panel_fun = panel_fun,
        size = unit(4, "cm"), gap = unit(1, "cm"), width = unit(5, "cm")
      ),
      show_annotation_name = FALSE, gap = unit(0, "points"), which = "row"
    )
    df_right_annotation <- c(df_right_annotation, zoom_box)
  }

  ht <- ComplexHeatmap::Heatmap(
    matrix = df,
    col = colorRamp2(seq(-lim, lim, length = 100), color_palette),

    row_split = row_split,
    column_split = column_split,
    cluster_row_slices = cluster_row_slices,
    cluster_column_slices = cluster_column_slices,

    row_title = row_title,
    column_title = column_title,

    show_row_names = F,
    column_names_gp = gpar(rot = 60, fontsize = 8),

    top_annotation = df_top_annotation,
    left_annotation = df_left_annotation,
    right_annotation = df_right_annotation,

    border = T,
    heatmap_legend_param = ht_legend,
    ...
  )

  draw_ht <- function() {
    ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
    annotation_titles <- c(
      protein_coding = "Protein-coding",
      lincRNA = "lincRNA",
      miRNA = "miRNA",
      pseudogene = "pseudogene",
      antisense = "antisense",
      baseExp = "baseExp",
      log2fc = "log2(fc)",
      log10padj = "-log10(padj)"
    )
    for (an in names(annotation_titles)) {
      if (an %in% names(df_right_annotation@anno_list)) {
        decorate_annotation(an,
          {
            grid.text(annotation_titles[an], y = 1.05, rot = 90, just = "left", gp = gpar(fontsize = 8))
            grid.rect(gp = gpar(fill = NA, col = "black"))
          },
          slice = 1
        )
        if (length(unique(row_split)) == 2) {
          decorate_annotation(an,
            {
              grid.rect(gp = gpar(fill = NA, col = "black"))
            },
            slice = 2
          )
        }
      }
    }
    if (anno_genebiotype) {
      decorate_annotation("protein_coding", {
        grid.lines(x = unit(c(0, 0), "npc"), y = unit.c(unit(1, "npc"), unit(1, "npc") + unit(20, "mm")))
      })
      decorate_annotation("antisense", {
        grid.lines(x = unit(c(1, 1), "npc"), y = unit.c(unit(1, "npc"), unit(1, "npc") + unit(20, "mm")))
      })
    }
    decorate_annotation(
      annotation = "baseExp",
      {
        grid.lines(x = unit(c(0, 0), "native"), y = c(0, 1), gp = gpar(col = "black", lty = 2))
      },
      slice = 1
    )
    decorate_annotation(
      annotation = "log2fc",
      {
        grid.lines(x = unit(c(0, 0), "native"), y = c(0, 1), gp = gpar(col = "black", lty = 2))
      },
      slice = 1
    )
    if (length(unique(row_split)) == 2) {
      decorate_annotation(
        annotation = "baseExp",
        {
          grid.lines(x = unit(c(0, 0), "native"), y = c(0, 1), gp = gpar(col = "black", lty = 2))
        },
        slice = 2
      )
      decorate_annotation(
        annotation = "log2fc",
        {
          grid.lines(x = unit(c(0, 0), "native"), y = c(0, 1), gp = gpar(col = "black", lty = 2))
        },
        slice = 2
      )
    }
  }
  if (return_ggplot) {
    grob <- grid.grabExpr({
      draw_ht()
    })
    grob <- addGrob(
      gTree = grob,
      child = textGrob(
        label = paste0("DGEs Heatmap\n(", label, ")"),
        x = unit(0, "npc"), y = unit(0.98, "npc"), just = "left"
      )
    )
    p <- as_ggplot(grob)
    if (!anno_zoom_box) {
      add_width <- anno_genebiotype * 5
      p <- p + theme(
        aspect.ratio = 40 / (ncol(df) + 36 + add_width),
        plot.margin = margin(t = 15, b = 15)
      )
    }
  } else {
    p <- draw_ht()
  }
  return(p)
}

pie_chart <- function(res, annotation_matrix, label) {
  df <- res
  stat <- table(df[, "DifferentialExpression"])
  df[, "group"] <- paste0(df[, "DifferentialExpression"], "(", stat[df[, "DifferentialExpression"]], ")")
  df <- df[with(df, order(group)), ]

  df <- df %>%
    group_by(gene_biotype, group) %>%
    mutate(count = n())
  df <- df %>%
    group_by(group) %>%
    mutate(perc = count / n())
  df1 <- unique(df[, c("gene_biotype", "group", "count", "perc")])

  group <- rev(sort(unique(annotation_matrix$gene_biotype)))
  color <- pal_igv(palette = c("default"), alpha = 1)(length(group))
  names(color) <- group

  cp <- coord_polar(theta = "y")
  cp$is_free <- function() TRUE

  p <- ggplot(df1, aes(x = 2, y = count, fill = gene_biotype)) +
    geom_bar(color = "black", stat = "identity") +
    geom_text(aes(label = ifelse(perc > 0.15, paste0(round(perc * 100, 1), "%\n", "(", count, ")"), "")),
      color = "white", size = 3.5, stat = "identity", position = position_stack(vjust = 0.5)
    ) +
    labs(title = paste0("Gene biotype (", label, ")")) +
    xlim(0.5, 2.5) +
    scale_fill_manual(values = color) +
    cp +
    facet_wrap(~group, scales = "free") +
    guides(fill = guide_legend(ncol = 4, byrow = F)) +
    theme_pubr(border = T) +
    theme(
      aspect.ratio = 1,
      legend.key.width = unit(4, "mm"),
      legend.key.height = unit(4, "mm"),
      strip.text = element_text(size = 12),
      panel.grid.major = element_line(),
      panel.grid.minor = element_line(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_blank()
    )

  return(p)
}

gene_test_plot <- function(test, drop_low, drop_outlier = NULL, annotation_matrix, label) {
  df <- data.frame(
    gene = c(test, drop_low, drop_outlier),
    class = c(
      rep("Tested", length(test)),
      rep("Excluded.LowExp", length(drop_low)),
      rep("Excluded.Outlier", length(drop_outlier))
    ), stringsAsFactors = F
  )
  df[, "class"] <- factor(df[, "class"], levels = c("Excluded.Outlier", "Excluded.LowExp", "Tested"))
  df <- merge(x = df, by.x = "gene", y = annotation_matrix, by.y = "row.names")

  df <- df %>%
    group_by(gene_biotype, class) %>%
    mutate(count = n())
  df <- df %>%
    group_by(class) %>%
    mutate(perc = count / n(), total_count = n())
  df <- unique(df[, c("gene_biotype", "class", "count", "perc", "total_count")])

  group <- rev(sort(unique(annotation_matrix$gene_biotype)))
  color <- pal_igv(palette = c("default"), alpha = 1)(length(group))
  names(color) <- group

  p <- ggplot(df, aes(x = class, y = count, fill = gene_biotype)) +
    geom_bar(color = "black", stat = "identity") +
    geom_text(aes(label = ifelse(perc > 0.15, paste0(round(perc * 100, 1), "%\n", "(", count, ")"), "")),
      color = "white", size = 3.5, stat = "identity", position = position_stack(vjust = 0.5)
    ) +
    geom_text(
      data = unique(df[, c("class", "total_count")]),
      aes(x = class, y = total_count, label = total_count),
      vjust = -0.5, angle = -90, size = 5, inherit.aes = FALSE
    ) +
    labs(title = paste0("Gene-filtering stat (", label, ")"), x = "", y = "Count") +
    scale_fill_manual(values = color) +
    scale_y_continuous(expand = expansion(mult = 0.1, add = 0.6)) +
    coord_flip() +
    guides(fill = guide_legend(ncol = 4, byrow = F)) +
    theme_pubr(border = T) +
    theme(
      aspect.ratio = 0.5, plot.margin = margin(r = 100),
      legend.key.width = unit(4, "mm"),
      legend.key.height = unit(4, "mm"),
      strip.text = element_text(size = 12),
      panel.grid.major = element_line(),
      panel.grid.minor = element_line(),
      legend.title = element_blank()
    )

  return(p)
}

scale_factor_compare <- function(edgeR_sf, DESeq2_sf, lib_size) {
  df <- data.frame(edgeR_ScalingFactor = edgeR_sf, DESeq2_ScalingFactor = DESeq2_sf, lib_Size = lib_size)
  df[, "edgeR_Strength"] <- df[, "edgeR_ScalingFactor"] / df[, "lib_Size"] * mean(df[, "lib_Size"])
  df[, "DESeq2_Strength"] <- df[, "DESeq2_ScalingFactor"] / df[, "lib_Size"] * mean(df[, "lib_Size"])
  df[, "Sample"] <- rownames(df)
  df <- reshape2::melt(df, id.vars = "Sample")
  df[, "type"] <- str_extract(string = df[, "variable"], pattern = "(?<=_).*")
  df[, "type"] <- factor(df[, "type"], levels = c("Size", "ScalingFactor", "Strength"))
  df[, "tool"] <- str_extract(string = df[, "variable"], pattern = ".*(?=_)")

  df1 <- df[which(df[, "variable"] == "lib_Size"), ]
  p1 <- ggplot(data = df1, aes(x = Sample, y = value)) +
    geom_bar(aes(fill = value), color = "black", stat = "identity") +
    geom_text(aes(label = round(value, digits = 2)), vjust = -0.5, size = 2.5) +
    ylim(0, ceiling(max(df1[, "value"])) * 1.05) +
    scale_fill_viridis_c(guide = FALSE) +
    labs(title = paste0("Library size"), x = "", y = "Number of reads/fragments\n(Millions)") +
    theme_pubr(border = T) +
    theme(
      aspect.ratio = 6 / length(unique(df1[, "Sample"])),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 10),
      panel.grid.major = element_line()
    )

  df2 <- df[which(df[, "variable"] != "lib_Size"), ]
  p2 <- ggplot(data = df2, aes(x = Sample, y = value, fill = tool)) +
    geom_bar(color = "black", stat = "identity", position = position_dodge()) +
    geom_hline(yintercept = 1, linetype = 2, color = "red") +
    geom_text(aes(label = round(value, digits = 2)), position = position_dodge(width = 1), vjust = -0.5, size = 2.5) +
    labs(title = "Scaling factors", x = "", y = "") +
    scale_fill_manual(name = "", values = nord("mountain_forms")[c(1, 3)]) +
    scale_y_continuous(breaks = seq(round(min(df2$value), digits = 1), round(max(df2$value), digits = 1), 0.1)) +
    coord_cartesian(ylim = c(round(min(df2$value), digits = 1) - 0.1, round(max(df2$value), digits = 1) + 0.1)) +
    facet_wrap(type ~ ., scales = "free", nrow = 2) +
    theme_pubr(border = T, legend = "top") +
    theme(
      panel.grid.major = element_line(),
      strip.background = element_rect(fill = "transparent"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 10)
    )
  p <- ggarrange(p1, p2, widths = c(2, 3))
  return(p)
}

cpm_gene_curve <- function(count_matrix, cpm_seq, label, minCPM, minCPM_label) {
  cpm_matrix <- as.data.frame(cpm(count_matrix))
  df <- apply(cpm_matrix, 2, function(y) {
    sapply(cpm_seq, FUN = function(x) {
      sum(y > x)
    })
  })
  df <- as.data.frame(df)
  df[, "CPM"] <- cpm_seq
  df <- reshape2::melt(df, id.vars = "CPM")
  p <- ggplot(data = df, aes(x = CPM, y = value, color = variable)) +
    geom_point(shape = 21, alpha = 0.8) +
    geom_path(alpha = 0.8) +
    scale_color_viridis_d() +
    geom_vline(
      data = data.frame(minCPM = minCPM, minCPM_label = minCPM_label),
      aes(xintercept = minCPM, linetype = minCPM_label), alpha = 0.8
    ) +
    scale_linetype_manual(values = c("edgeR" = 3, "DESeq2" = 2)) +
    guides(linetype = guide_legend(order = 1)) +
    labs(title = paste0("CPM-Gene curve (", label, ")"), x = "CPM", y = "Number of genes >CPM cut-off") +
    facet_zoom(xlim = c(0, ceiling(max(minCPM)) + 1), zoom.size = 1) +
    theme_pubr(border = T, legend = "top") +
    theme(
      aspect.ratio = 0.6,
      plot.margin = margin(t = 20, r = 100, b = 20, l = 100),
      panel.grid.major = element_line(),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent")
    )

  return(p)
}

dges_compare <- function(data, col_log2fc = c("log2FoldChange.x", "log2FoldChange.y"), col_label = "label", comparsion = "", min_fc = 2, tool = "edgeR") {
  groups <- sort(unique(data[, col_label]))

  colors <- rep("purple4", length(groups))
  names(colors) <- groups
  color_spec <- c("#3288BD", "#66C2A5", "#ABDDA4", "grey80", "#FDAE61", "#F46D43", "#D53E4F")
  colors[grep(x = groups, pattern = "Down-regulated.*Down-regulated")] <- color_spec[1]
  colors[grep(x = groups, pattern = "Down-regulated.*Non-significant")] <- color_spec[2]
  colors[grep(x = groups, pattern = "Non-significant.*Down-regulated")] <- color_spec[3]
  colors[grep(x = groups, pattern = "Non-significant.*Non-significant")] <- color_spec[4]
  colors[grep(x = groups, pattern = "Non-significant.*Up-regulated")] <- color_spec[5]
  colors[grep(x = groups, pattern = "Up-regulated.*Non-significant")] <- color_spec[6]
  colors[grep(x = groups, pattern = "Up-regulated.*Up-regulated")] <- color_spec[7]

  colors[grep(x = groups, pattern = "Up-regulated.*Down-regulated")] <- "black"
  colors[grep(x = groups, pattern = "Down-regulated.*Up-regulated")] <- "black"
  colors[grep(x = groups, pattern = "Excluded")] <- "purple4"

  alpha <- rep(1, length(groups))
  names(alpha) <- groups
  alpha[grep(x = groups, pattern = "Down-regulated.*Down-regulated")] <- 1
  alpha[grep(x = groups, pattern = "Down-regulated.*Non-significant")] <- 0.6
  alpha[grep(x = groups, pattern = "Non-significant.*Down-regulated")] <- 0.6
  alpha[grep(x = groups, pattern = "Non-significant.*Non-significant")] <- 0.2
  alpha[grep(x = groups, pattern = "Non-significant.*Up-regulated")] <- 0.6
  alpha[grep(x = groups, pattern = "Up-regulated.*Non-significant")] <- 0.6
  alpha[grep(x = groups, pattern = "Up-regulated.*Up-regulated")] <- 1

  alpha[grep(x = groups, pattern = "Up-regulated.*Down-regulated")] <- 1
  alpha[grep(x = groups, pattern = "Down-regulated.*Up-regulated")] <- 1

  x_max <- ceiling(max(abs(data[, col_log2fc[1]]), na.rm = TRUE))
  y_max <- ceiling(max(abs(data[, col_log2fc[2]]), na.rm = TRUE))

  data1 <- data[-grep(x = data[, col_label], pattern = "Excluded"), ]
  p1 <- ggplot(
    data = data1,
    aes(x = data1[, col_log2fc[1]], y = data1[, col_log2fc[2]], fill = data1[, col_label], alpha = data1[, col_label])
  ) +
    geom_point(shape = 21) +
    geom_rug(aes(color = data1[, col_label], alpha = data1[, col_label])) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = c(-log2(min_fc), log2(min_fc)), linetype = 2) +
    geom_hline(yintercept = c(-log2(min_fc), log2(min_fc)), linetype = 2) +
    labs(title = paste0(tool, ": ", paste0(comparsion, collapse = " vs ")), x = paste0(comparsion[1], "\nlog2FC"), y = paste0(comparsion[2], "\nlog2FC")) +
    scale_x_continuous(breaks = seq(-x_max, x_max, 1), limits = c(-x_max, x_max)) +
    scale_y_continuous(breaks = seq(-y_max, y_max, 1), limits = c(-y_max, y_max)) +
    scale_fill_manual(name = "", values = colors) +
    scale_color_manual(name = "", values = colors, guide = F) +
    scale_alpha_manual(values = alpha, guide = F) +
    guides(fill = guide_legend(override.aes = list(size = 4))) +
    theme_pubr(border = T, legend = "right") +
    theme(
      aspect.ratio = 1,
      panel.grid.major = element_line()
    )

  size <- table(data[, col_label])
  data2 <- data
  data2[, "point_size"] <- size[data2[, col_label]]
  data2 <- unique(data2[, c("DifferentialExpression.x", "DifferentialExpression.y", col_label, "point_size")])
  data2[, "DifferentialExpression.x"] <- factor(data2[, "DifferentialExpression.x"],
    levels = sort(unique(data2[, "DifferentialExpression.x"]))
  )
  data2[, "DifferentialExpression.y"] <- factor(data2[, "DifferentialExpression.y"],
    levels = sort(unique(data2[, "DifferentialExpression.y"]))
  )
  p2 <- ggplot(data = data2, aes(x = data2[, "DifferentialExpression.x"], y = data2[, "DifferentialExpression.y"])) +
    geom_point(aes(size = as.numeric(log2(data2[, "point_size"])), fill = data2[, "label"]),
      color = "black", shape = 21, stroke = 1
    ) +
    geom_text(aes(label = as.character(data2[, "point_size"])), color = "white", fontface = 2) +
    labs(title = paste0(tool, ": ", paste0(comparsion, collapse = " vs ")), x = comparsion[1], y = comparsion[2]) +
    scale_size_continuous(range = c(8, 16), guide = F) +
    scale_fill_manual(name = "", values = colors, guide = F) +
    theme_pubr(border = T) +
    theme(
      aspect.ratio = 1,
      panel.grid.major = element_line(),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    )

  return(list(p1, p2))
}

alluvial_plot <- function(res_list, annotation_matrix, tool, stratum_col) {
  compare <- names(res_list)[str_detect(string = names(res_list), paste0(tool, "_res"))]
  compare_label <- gsub(x = compare, pattern = paste0(tool, "_res@"), replacement = "")
  compare_df <- lapply(compare, function(name) {
    df <- res_list[[name]]
    name <- gsub(x = name, pattern = paste0(tool, "_res@"), replacement = "")
    all_gene <- rownames(counts_matrix)
    df <- df[, c("RowID", "DifferentialExpression")]
    drop <- setdiff(x = all_gene, y = df[, "RowID"])
    df <- rbind(df, data.frame(RowID = drop, DifferentialExpression = "Excluded"), stringsAsFactors = F)
    df <- cbind(df, name, stringsAsFactors = F)
    return(df)
  })
  compare_df <- bind_rows(compare_df)
  compare_df <- compare_df %>%
    group_by(RowID, DifferentialExpression) %>%
    mutate(count = n())
  DGEs <- unique(compare_df[compare_df$DifferentialExpression %in% c("Up-regulated", "Down-regulated"), "RowID", ])
  compare_df <- compare_df[which(compare_df$RowID %in% DGEs$RowID), ]
  compare_df <- merge(x = compare_df, y = annotation_matrix, by.x = "RowID", by.y = "row.names", all.x = T)

  compare_df$name <- gsub(x = compare_df$name, pattern = "_vs_", replacement = "/")
  compare_df$name <- factor(compare_df$name, levels = gsub(x = compare_label, pattern = "_vs_", replacement = "/"))
  compare_df$DifferentialExpression <- factor(compare_df$DifferentialExpression,
    levels = (c("Up-regulated", "Down-regulated", "Non-significant", "Excluded"))
  )
  compare_df$unit <- 1

  p <- ggplot(
    compare_df,
    aes(
      x = name, y = unit, stratum = compare_df[, stratum_col], alluvium = RowID,
      fill = compare_df[, stratum_col], label = compare_df[, stratum_col]
    )
  ) +
    labs(title = paste0("DGEs alluvial plot (", tool, ")"), x = "", y = "Count") +
    scale_x_discrete(expand = c(0.1, 0.1)) +
    # scale_fill_viridis_d()+
    scale_fill_brewer(type = "div", palette = "Spectral") +
    geom_stratum(width = 0.4, alpha = 0.6, color = "black") +
    geom_flow(width = 0.4, alpha = 0.6, color = "black", aes.flow = "forward") +
    ggfittext::geom_fit_text(stat = "stratum", width = 0.4, min.size = 4.5) +
    theme_pubr() +
    theme(
      aspect.ratio = 2 / length(compare),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      legend.title = element_blank()
    )

  return(p)
}
