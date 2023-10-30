#!/usr/bin/env bash

#######################################################################################
trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

Rscript &>/dev/null
[ $? -eq 127 ] && {
    color_echo "red" "Cannot find the command Rscript.\n"
    exit 1
}

R_packages=("dplyr" "stringr" "ggplot2" "ggsci" "ggtree" "RColorBrewer" "cowplot" "aplot" "ggplotify" "edgeR" "sva" "limma" "patchwork" "ggrepel" "Rtsne" "plotly" "plot3D" "grid" "ggforce" "aplot" "factoextra" "ComplexHeatmap" "circlize")

all_installed=$(Rscript -e "installed.packages()" | awk '{print $1}')
Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager',repos='https://cran.r-project.org/')"
Rscript -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes',repos='https://cran.r-project.org/')"
for package in "${R_packages[@]}"; do
  if ! echo "$all_installed" | grep -q "$package"; then
    color_echo "yellow" "Install the R package: '$package'.\n"
    Rscript -e "BiocManager::install('$package')"
  fi
done

echo -e "########################## BatchCorrection Parameters ##########################\n"
echo -e "  Rscript path: $(which Rscript)"
echo -e "  Aligner: ${Aligner}\n"
echo -e "  SampleInfoFile: ${SampleInfoFile}\n"
echo -e "################################################################################\n"

echo -e "****************** Start BatchCorrection ******************\n"
SECONDS=0

mkdir -p $maindir/NGSmodule_analysis/Quantification/BatchCorrection
cd $maindir/NGSmodule_analysis/Quantification/BatchCorrection

Rscript $1 $maindir $Aligner $SampleInfoFile $1

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "****************** BatchCorrection Done *******************\n"
