#!/usr/bin/env bash

#######################################################################################
trap 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM EXIT

#featureCounts &>/dev/null;[ $? -eq 127 ] && { echo -e "Cannot find the command featureCounts.\n";exit 1; }
Rscript &>/dev/null
[ $? -eq 127 ] && {
  color_echo "red" "Cannot find the command Rscript.\n"
  exit 1
}

R_packages=("Rsubread" "edgeR" "Rsamtools" "refGenome" "AnnotationDbi" "org.Hs.eg.db" "org.Mm.eg.db" "org.Mmu.eg.db" "org.Dm.eg.db")
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

if [[ ! -f $gtf ]]; then
  color_echo "red" "ERROR! Cannot find the gtf file: $gtf\nPlease check the Alignment Paramaters in your ConfigFile.\n"
  exit 1
fi

echo -e "########################## Quantification Parameters ###########################\n"
echo -e "  featurecounts_threads: ${threads_featurecounts}\n  Strand_Specific: ${strandspecific} (0=unstranded,1=stranded,2=reversely stranded)\n"
echo -e "  GTF_File: ${gtf}\n "
echo -e "################################################################################\n"

echo -e "****************** Start Quantification ******************\n"
SECONDS=0

for sample in "${arr[@]}"; do
  read -u1000
  {
    echo "+++++ $sample +++++"
    dir=$work_dir/$sample
    bam=$dir/$Aligner/${sample}.${Aligner}.bam
    if [[ ! -f $bam ]]; then
      echo -e "ERROR: Bam file:$bam do not exist. Please check the file.\n"
      exit 1
    fi

    mkdir -p $dir/$Aligner/Quantification
    cd $dir/$Aligner/Quantification
    Rscript $1 $threads_featurecounts $gtf $strandspecific $bam ${sample}.${Aligner} &>Quantification.R.log

    echo "Completed: $sample" >>$TMPFILE
    color_echo "green" "***** Completed:$(cat $TMPFILE | grep "Completed" | wc -l) | Interrupted:$(cat $TMPFILE | grep "Interrupted" | wc -l) | Total:$total_task *****"

    echo >&1000
  } &
  ((bar++))
  processbar $bar $total_task
done
wait

echo -e "\nIntegrating and annotating the matrix....\n"
mkdir -p $maindir/NGSmodule_analysis/Quantification
cd $maindir/NGSmodule_analysis/Quantification

species_arr=('human' 'mouse' 'rhesus' 'fly')
if [[ "${species_arr[*]}" =~ ${Species} ]]; then
  Species_anno=$Species
else
  Species_anno=""
  echo -e "No additional annotation for species: $Species\n"
fi

Rscript $2 $work_dir $gtf $Aligner $Species_anno $Database &>Annotation.R.log
echo -e "Integrated quantification matrix: $maindir/NGSmodule_analysis/Quantification/Quantification.${Aligner}.*.tab\n"

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "****************** Quantification Done ******************\n"
