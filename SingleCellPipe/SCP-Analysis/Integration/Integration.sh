#!/usr/bin/env bash

#######################################################################################
trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

Rscript &>/dev/null
[ $? -eq 127 ] && {
  color_echo "red" "Cannot find the command Rscript.\n"
  exit 1
}

R_packages=("sctransform" "Seurat" "SeuratWrappers" "SeuratDisk" "intrinsicDimension" "scater" "Matrix" "BiocParallel" "future" "reticulate" "harmony" "simspec" "plyr" "dplyr" "RColorBrewer" "scales" "gtools" "ggsci" "ggpubr" "ggplot2" "ggtree" "cowplot" "reshape2" "stringr" "velocyto.R" "scDblFinder" "biomaRt" "rvest" "xml2")
for package in "${R_packages[@]}"; do
  Rscript -e "installed.packages()" | awk '{print $1}' | grep $package &>/dev/null
  if [ $? -ne 0 ]; then
    color_echo "red" "Cannot find the R package $package.\n"
    if [[ $package == "refGenome" ]]; then
      color_echo "red" "Please install it in the R environment using \"remotes::install_version('$package')\" "
    else
      color_echo "red" "Please install it in the R environment using \"BiocManager::install('$package')\" "
    fi
    exit 1
  fi
done

echo -e "########################## Quantification Parameters ###########################\n"
echo -e "*** base parameters ***"
echo -e "    Rscript path: $(which Rscript)\n    datasets: ${datasets}\n    species: ${species}\n    exogenous_genes: ${exogenous_genes}"
echo -e "*** cell-filtering parameters ***"
echo -e "    cell_calling_methodNum: ${cell_calling_methodNum}\n    mito_threshold: ${mito_threshold}"
echo -e "*** integration parameters ***"
echo -e "    HVF_source: ${HVF_source}\n    nHVF: ${nHVF}\n    anchor_dims: ${anchor_dims}\n    integrate_dims: ${integrate_dims}"
echo -e "*** clustering parameters ***"
echo -e "    maxPC: ${maxPC}\n    resolution: ${resolution}\n    reduction: ${reduction}\n"
echo -e "################################################################################\n"

echo -e "****************** Start Integration ******************\n"
SECONDS=0

echo -e "Integrating the data....\n"
mkdir -p $maindir/NGSmodule_SCP_analysis/Integration
cd $maindir/NGSmodule_SCP_analysis/Integration

Rscript $1 $1 $maindir/NGSmodule_SCP_analysis/Integration "${work_dir}" "${threads}" "${datasets}" \
  "${species}" "${exogenous_genes}" "${cell_calling_methodNum}" "${mito_threshold}" "${HVF_source}" \
  "${nHVF}" "${anchor_dims}" "${integrate_dims}" "${maxPC}" "${resolution}" \
  "${reduction}" 2>&1 |tee Integration.log 
echo -e "Integration completed.\n"

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "****************** Quantification Done ******************\n"
