#!/usr/bin/env bash

#######################################################################################
trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

Rscript &>/dev/null
[ $? -eq 127 ] && {
  color_echo "red" "Cannot find the command Rscript.\n"
  exit 1
}

R_packages=("biomaRt" "Biostrings")
for package in "${R_packages[@]}"; do
  Rscript -e "installed.packages()" | awk '{print $1}' | grep $package &>/dev/null
  if [ $? -ne 0 ]; then
    color_echo "red" "Cannot find the R package $package.\n"
    color_echo "red" "Please install it in the R environment using \"BiocManager::install('$package')\" "
    exit 1
  fi
done

echo -e "########################## Prefetch_Sequence Parameters ###########################\n"
echo -e "  Rscript path: $(which Rscript)"
echo -e "  Species: ${species}"
echo -e "  Sequence_type: ${sequence_type}"
echo -e "################################################################################\n"

echo -e "****************** Start prefetching sequence ******************\n"
SECONDS=0

Rscript $2 $work_dir $gtf $Aligner $Species $Source &>Annotation.R.log
echo -e "Integrated quantification matrix: $maindir/NGSmodule_analysis/Quantification/Quantification.${Aligner}.*.tab\n"

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "****************** Quantification Done ******************\n"
