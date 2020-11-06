library(plyr)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(ggpubr)
library(factoextra)
library(gtools)
library(cowplot)
library(RColorBrewer)
library(future)
library(CellChat)
library(ggplot2)
library(ggalluvial)

srt.integrated <- readRDS("/data/lab/WangMengJie/scRNAseq_data_human_testis/analysis/Results/All_cells/rds_generated_by_run/resolution-1/yo_1-yo_7-yo_11-yo_13-yo_14-yo_25.regress_umi_cellcycle.mod.rds")
srt_use <- srt.integrated
DefaultAssay(srt_use) <- "RNA"
new.cluster.ids <- c("Peritubular Myoid & Leydig",
                     "Peritubular Myoid & Leydig",
                     "Endothelial cells",
                     "Germ cells",
                     "Sertoli",
                     "Sertoli",
                     "Endothelial cells",
                     "Germ cells",
                     "Peritubular Myoid & Leydig",
                     "Germ cells",
                     "Peritubular Myoid & Leydig",
                     "Macrophage",
                     "Sertoli",
                     "Germ cells",
                     "Endothelial cells",
                     "Germ cells",
                     "Endothelial cells",
                     "Germ cells",
                     "Sertoli",
                     "Smooth muscle",
                     "Germ cells",
                     "Endothelial cells",
                     "Uncertain1",
                     "Uncertain2")
names(new.cluster.ids) <- levels(srt_use)
srt_use <- RenameIdents(srt_use, new.cluster.ids)
DimPlot(srt_use,group.by = "ident",label = TRUE)
















# LR calculate ------------------------------------------------------------
LR_score <- function(norm_data, ligand = c("AMH"), receptor = c("AMHR2", "BMPR1A"), LR_vec = c(1, -1),
                     cluster1 = germ, cluster2 = sertoli, clusters_name = NULL,
                     do_plot = FALSE) {
  require(dplyr)
  require(ggplot2)
  require(ggpubr)
  require(ggsci)
  require(ggridges)
  require(cowplot)

  if (LR_vec[1] <= 0 | LR_vec[2] >= 0) {
    stop("LR_vec[1] should be a positive value for ligand and LR_vec[2] be a negative value for receptor.",
      call. = FALSE
    )
  }

  df <- norm_data[c(ligand, receptor), c(cluster1, cluster2)]
  exp <- apply(df, 1, mean)
  exp_coef <- exp / mean(exp)
  weight <- c(LR_vec[1] / length(ligand) / exp_coef[ligand], LR_vec[2] / length(receptor) / exp_coef[receptor])
  weight[is.infinite(weight)] <- 0

  LRpreference <- apply(df, 2, function(x) {
    res <- sum(x * weight) / length(x)
    return(res)
  })

  comb <- expand.grid(ligand, receptor)
  autocrine <- apply(df, 2, function(x) {
    res <- 0
    for (i in 1:nrow(comb)) {
      l <- as.character(comb[i, 1])
      r <- as.character(comb[i, 2])
      lexp <- x[l] * weight[l]
      rexp <- (-x[r] * weight[r])
      gm <- exp(mean(log(c(lexp, rexp))))
      res <- res + gm
    }
    return(res)
  })

  if (isTRUE(do_plot)) {
    if (is.null(clusters_name)) {
      clusters_name <- c("cluster1", "cluster2")
    }

    df_plot <- as.data.frame(t(df))
    df_plot[, "cluster"] <- c(rep(clusters_name[1], length(cluster1)), rep(clusters_name[2], length(cluster2)))
    df_plot[, "LRpreference"] <- LRpreference
    df_plot[, "autocrine"] <- autocrine


    p1 <- reshape2::melt(df_plot, measure.vars = c(ligand, receptor), variable.name = "LR", value.name = "exp") %>%
      ggplot(aes(x = cluster, y = exp)) +
      geom_violin(aes(fill = LR)) +
      geom_boxplot(aes(color = LR), position = position_dodge(0.9), width = 0.2, fill = "grey", alpha = 0.5, show.legend = FALSE) +
      stat_summary(aes(group = LR),
        fun = median, geom = "point",
        size = 2, shape = 21, color = "black", fill = "white", position = position_dodge(0.9)
      ) +
      scale_fill_nejm() +
      scale_color_manual(values = rep("black", length(c(ligand, receptor)))) +
      labs(title = "Normalized expression", subtitle = "1st:ligand; other:receptor") +
      theme_classic(base_size = 13) +
      theme(aspect.ratio = 0.8, axis.title.x = element_blank())

    p2 <- reshape2::melt(df_plot, measure.vars = c(ligand, receptor), variable.name = "LR", value.name = "exp") %>%
      ggplot(aes(x = cluster, y = exp, fill = LR)) +
      stat_summary(fun = mean, geom = "bar", color = "black", position = position_dodge(0.9)) +
      stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", position = position_dodge(0.9), width = 0.5) +
      scale_fill_nejm() +
      scale_color_manual(values = rep("black", length(c(ligand, receptor)))) +
      labs(title = "Normalized expression", subtitle = "1st:ligand; other:receptor") +
      theme_classic(base_size = 13) +
      theme(aspect.ratio = 0.8, axis.title.x = element_blank())

    p3 <- ggplot(data = df_plot, aes(x = LRpreference, y = cluster, fill = cluster)) +
      geom_vline(xintercept = 0, size = 1, color = "grey50") +
      geom_density_ridges(
        jittered_points = TRUE,
        position = position_points_jitter(width = 0.05, height = 0, yoffset = -0.05),
        point_shape = "|", point_size = 3, point_alpha = 1, alpha = 0.7, show.legend = FALSE
      ) +
      scale_fill_npg() +
      scale_x_continuous(limits = c(-max(abs(df_plot$LRpreference)), max(abs(df_plot$LRpreference)))) +
      labs(title = "LR preference distribution", subtitle = "ligand>0; receptor<0", x = "LR preference") +
      theme_classic(base_size = 13) +
      theme(aspect.ratio = 0.8, axis.title.y = element_blank())

    df_plot_summary <- df_plot %>%
      group_by(cluster) %>%
      summarise(new = list(mean_se(LRpreference))) %>%
      unnest(new) %>%
      mutate(Type = ifelse(y > 0, "ligand", "receptor"))

    p4 <- ggplot(df_plot_summary, aes(x = cluster, y = y, fill = Type)) +
      geom_col(color = "black") +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.5) +
      geom_hline(yintercept = 0, size = 1, color = "grey50") +
      stat_compare_means(
        data = df_plot, mapping = aes(x = cluster, y = LRpreference, label = paste0("p = ", ..p.format..)),
        method = "t.test", label.y = 1.05 * max(df_plot_summary$ymax), inherit.aes = FALSE
      ) +
      scale_fill_manual(values = setNames(
        nm = c("ligand", "receptor"),
        object = c("red3", "royalblue4")
      )) +
      labs(title = "LR preference stat", subtitle = "ligand>0; receptor<0", y = "LR preference") +
      theme_classic(base_size = 13) +
      theme(aspect.ratio = 0.8, axis.title.x = element_blank())

    p5 <- ggplot(data = df_plot, aes(x = autocrine, y = df_plot$cluster, fill = df_plot$cluster)) +
      geom_vline(xintercept = 0, size = 1, color = "grey50") +
      geom_density_ridges(
        jittered_points = TRUE,
        position = position_points_jitter(width = 0.05, height = 0, yoffset = -0.05),
        point_shape = "|", point_size = 3, point_alpha = 1, alpha = 0.7, show.legend = FALSE
      ) +
      scale_fill_npg() +
      scale_x_continuous(limits = c(-max(abs(df_plot$autocrine)), max(abs(df_plot$autocrine)))) +
      labs(title = "Autocrine score distribution", subtitle = "high-value indicated autocrine signals", x = "Autocrine score") +
      theme_classic(base_size = 13) +
      theme(aspect.ratio = 0.8, axis.title.y = element_blank())

    df_plot_summary2 <- df_plot %>%
      group_by(cluster) %>%
      summarise(new = list(mean_se(autocrine))) %>%
      unnest(new)

    p6 <- ggplot(df_plot_summary2, aes(x = df_plot$cluster, y = y, fill = y)) +
      geom_col(color = "black") +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.5) +
      stat_compare_means(
        data = df_plot, mapping = aes(x = df_plot$cluster, y = autocrine, label = paste0("p = ", ..p.format..)),
        method = "t.test", label.y = 1.05 * max(df_plot_summary2$ymax), inherit.aes = FALSE
      ) +
      scale_fill_material("blue-grey", guide = FALSE) +
      labs(title = "LR autocrine score stat", subtitle = "high-value indicated autocrine signals", y = "autocrine score") +
      theme_classic(base_size = 13) +
      theme(aspect.ratio = 0.8, axis.title.x = element_blank())

    plot <- plot_grid(p1, p2, p3, p4, p5, p6, align = "hv", axis = "tblr", nrow = 3, ncol = 2, labels = "AUTO")
    print(plot)
  }
  return(data.frame(
    row.names = colnames(df),
    "preference" = LRpreference,
    "autocrine" = autocrine,
    stringsAsFactors = FALSE
  ))
}

norm_data <- srt_use@assays$RNA@data
germ <- WhichCells(object = srt_use, idents = "Germ cells")
sertoli <- WhichCells(object = srt_use, idents = "Sertoli")
pair <- readRDS("/data/lab/WangMengJie/scRNAseq_data_human_testis/analysis/Results/All_cells/rds_generated_by_run/resolution-1/Allages_L_Sertoli_R_Germ_cells_rough_RNAslot.rds")

LRpreference <- LR_score(
  norm_data = norm_data, ligand = c("AMH"), receptor = c("AMHR2", "BMPR1A"),LR_vec = c(1, -1),
  cluster1 = germ, cluster2 = sertoli, clusters_name = c("germ","sertoli"),
  do_plot = TRUE
)
LRpreference <- LR_score(
  norm_data = norm_data, ligand = c("BMP2"), receptor = c("BMPR1A", "ACVR2B"),LR_vec = c(1, -1),
  cluster1 = germ, cluster2 = sertoli, clusters_name = c("germ","sertoli"),
  do_plot = TRUE
)
LRpreference <- LR_score(
  norm_data = norm_data, ligand = c("WNT2"), receptor = c("FZD1", "LRP6"),LR_vec = c(1, -1),
  cluster1 = germ, cluster2 = sertoli, clusters_name = c("germ","sertoli"),
  do_plot = TRUE
)

all_pair <- unlist(pair, recursive = FALSE) %>% unique()
df <- bplapply(all_pair, function(LR) {
  ligand <- LR[1]
  receptor <- LR[2:length(LR)]
  df_score <- LR_score(
    norm_data = norm_data, ligand = ligand, receptor = receptor,LR_vec = c(1, -1),
    cluster1 = germ, cluster2 = sertoli, clusters_name = c("germ","sertoli"),
    do_plot = FALSE
  )
  return(df_score)
}, BPPARAM = MulticoreParam())
names(df) <- sapply(all_pair, function(x)paste0(x,collapse = ";"))


i <- 500
LR <- all_pair[[i]]
df_score <- LR_score(
  norm_data = norm_data, ligand = LR[1], receptor = LR[2:length(LR)], LR_vec = c(1, -1),
  cluster1 = germ, cluster2 = sertoli, clusters_name = c("germ", "sertoli"),
  do_plot = TRUE
)






























