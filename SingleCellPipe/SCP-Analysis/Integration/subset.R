
# parameters: global settings ---------------------------------------------
work_dir <- "/data/lab/LiLaiHua/scRNA-seq/Gonadal_ridge/analysis_zh/"
cellranger_dir <- "/data/lab/LiLaiHua/scRNA-seq/Gonadal_ridge/cellranger/"
threads <- 80
dataset_id <- c("LLH_5D", "LLH_7D_2")
dataset_name <- c("5d", "7d")
compare <- c("5d", "7d")

exogenous_gene <- "GFP"

# parameters: cell filtering ----------------------------------------------
minGene <- 1000
maxGene <- 10000
minUMI <- 4000
maxUMI <- 40000
maxMitoPercent <- 20

# parameters: integration -------------------------------------------------
HVF_source <- "separate"
nHVF <- 3000
anchor_dims <- 1:30
integrate_dims <- 1:30

# parameters: clustering --------------------------------------------------
maxPC <- 100
resolution <- 0.8


# parameters: subset ------------------------------------------------------
integration_method <- "Standard"
subset_cluster <- setNames(
  object = list(
    c(4, 5, 6, 9),
    c(3, 14, 10),
    c(1, 13),
    c(2, 7, 8, 11, 12, 15, 16, 17, 18, 19)
  ),
  nm = c("Germ", "Somatic1", "Somatic2", "Somatic3")
)


########################### Start the workflow ############################
setwd(work_dir)

srt_list <- readRDS(paste0("HVF_",HVF_source, "/", paste0("srt_list_", integration_method),".rds"))
srt_use <- srt_list[[paste0(compare, collapse = "-")]]

srt_subset_list <- list()
for (i in 1:length(subset_cluster)) {
  cat(names(subset_cluster)[i],":",subset_cluster[[i]], "\n")
  srt_subset <- subset(srt_use,seurat_clusters %in% subset_cluster[[i]])
  srt_subset <- Standard_integrate(sc_list = SplitObject(srt_subset,split.by = "orig.ident"),
                                   nHVF = nHVF,anchor_dims = anchor_dims,integrate_dims = integrate_dims,
                                   maxPC = maxPC,resolution = resolution,HVF_source = HVF_source)
  srt_subset_list[[paste0(names(subset_cluster)[i],":",paste0(subset_cluster[[i]],collapse = ","))]] <- srt_subset
}

if (!dir.exists(paste0("HVF_",HVF_source,"/subset"))) {
  dir.create(paste0("HVF_",HVF_source,"/subset"), recursive = T)
}
saveRDS(srt_subset_list,file = paste0("HVF_",HVF_source,"/subset","/srt_subset_list.rds"))


