#!/usr/bin/env bash


#######################################################################################
trap 'j=`ps aux | grep -P "$work_dir" |grep -P "Rscript"| awk '"'"'{print $2}'"'"'`;kill -9 $j;kill -9 $(jobs -p);echo -e "\nKilling all background processes......\nExiting the script......\n";exit 1' SIGINT

$Rscript &>/dev/null;[ $? -eq 127 ] && { echo -e "Cannot find the command Rscript.\n";exit 1; }
R_packages=("limma" "edgeR" "data.table" "gplots" "stringr" "ComplexHeatmap" "ggsci" "ggpubr" "RColorBrewer" "circlize" "ggrepel" "GGally" "factoextra" "nord")
for package in ${R_packages[@]};do
  $Rscript -e "installed.packages()" |awk '{print $1}' |grep $package &>/dev/null;[ $? -ne 0 ] && { echo -e "Cannot find the R package $package.\n";exit 1; }
done

echo -e "######################## postQuantificationQC Parameters #######################\n"
echo -e "  Aligner: ${Aligner}\n"
echo -e "  SampleInfoFile: ${SampleInfoFile}\n"
echo -e "################################################################################\n"


echo -e "****************** Start postQuantificationQC ******************\n"
SECONDS=0

mkdir -p $maindir/NGSmodule_analysis/Quantification/postQuantificationQC 
cd $maindir/NGSmodule_analysis/Quantification/postQuantificationQC

$Rscript $1 $maindir $Aligner $SampleInfoFile


ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "****************** postQuantificationQC Done *******************\n"





