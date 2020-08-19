#!/usr/bin/env bash

#######################################################################################
trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

Rscript &>/dev/null
[ $? -eq 127 ] && {
  color_echo "red" "Cannot find the command Rscript.\n"
  exit 1
}
R_packages=("BiocParallel" "edgeR" "DESeq2" "stringr" "scales" "RColorBrewer" "ggpubr" "ggsci" "ggforce" "reshape2" "VennDiagram" "gridExtra" "gplots" "dplyr" "openxlsx" "ggalluvial" "ggfittext" "ComplexHeatmap" "circlize" "nord" "ggupset")
for package in "${R_packages[@]}"; do
  Rscript -e "installed.packages()" | awk '{print $1}' | grep $package &>/dev/null
  [ $? -ne 0 ] && {
    color_echo "red" "Cannot find the R package $package.\n"
    exit 1
  }
done

echo -e "####################### DifferentialExpression Parameters ######################\n"
echo -e "  Rscript path: $(which Rscript)"
echo -e "  max_padj: $max_padj\n  min_fc: $min_fc\n  min_count:$min_count\n  group_compare:$group_compare\n  DGEs_multi_compare:$DGEs_multi_compare\n"
echo -e "################################################################################\n"

echo -e "****************** Start DifferentialExpression analysis ******************\n"
SECONDS=0

mkdir -p $maindir/NGSmodule_analysis/DifferentialExpression/DGEs_data
mkdir -p $maindir/NGSmodule_analysis/DifferentialExpression/DGEs_plot
mkdir -p $maindir/NGSmodule_analysis/DifferentialExpression/DGEs_rds
cd $maindir/NGSmodule_analysis/DifferentialExpression

Rscript $1 $maindir $Aligner $SampleInfoFile $group_compare $max_padj $min_fc $min_count $DGEs_multi_compare $1

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "*********************** DifferentialExpression Done ***********************\n"
