#!/usr/bin/env bash

#######################################################################################
trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

#featureCounts &>/dev/null;[ $? -eq 127 ] && { echo -e "Cannot find the command featureCounts.\n";exit 1; }
Rscript &>/dev/null
[ $? -eq 127 ] && {
  color_echo "red" "Cannot find the command Rscript.\n"
  exit 1
}

R_packages=("Rsubread" "edgeR" "Rsamtools" "data.table" "refGenome" "AnnotationDbi" "org.Hs.eg.db" "org.Mm.eg.db" "org.Mmu.eg.db" "org.Dm.eg.db")

all_installed=$(Rscript -e "installed.packages()" | awk '{print $1}')
Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager',repos='https://cran.r-project.org/')"
Rscript -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes',repos='https://cran.r-project.org/')"
for package in "${R_packages[@]}"; do
  if ! echo "$all_installed" | grep -q "$package"; then
    color_echo "yellow" "Install the R package: '$package'.\n"
    if [[ $package == "refGenome" ]]; then
      Rscript -e "remotes::install_version('refGenome',repos='https://cran.r-project.org/')"
    else
      Rscript -e "BiocManager::install('$package')"
    fi
  fi
done

if [[ ! -f $gtf ]]; then
  color_echo "red" "ERROR! Cannot find the gtf file: $gtf\nPlease check the Alignment Paramaters in your ConfigFile.\n"
  exit 1
fi

echo -e "########################## Quantification Parameters ###########################\n"
echo -e "  Rscript path: $(which Rscript)"
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
    
    if [[ ${Aligner} == "kallisto" ]]; then
      if [[ -f $dir/Alignment-$Aligner/abundance.tsv ]]; then
        color_echo "green" "kallisto quant output (abundance.tsv) detected."
      else
        echo -e "ERROR: kallisto quant output do not exist. Please re-run Alignment.\n"
        exit 1
      fi
    else
      bam=$dir/Alignment-$Aligner/${sample}.${Aligner}.bam
      if [[ ! -f $bam ]]; then
        echo -e "ERROR: Bam file:$bam do not exist. Please re-run Alignment.\n"
        exit 1
      fi

      mkdir -p $dir/Alignment-$Aligner/Quantification
      rm -f $dir/Alignment-$Aligner/Quantification/*.tmp
      
      cd $dir/Alignment-$Aligner/Quantification
      Rscript $1 $threads_featurecounts $gtf $strandspecific $bam ${sample}.${Aligner} &>Quantification.R.log
    fi
    
    echo "Completed: $sample" >>$tmpfile
    color_echo "green" "***** Completed:$(cat "$tmpfile" | grep "Completed" | uniq | wc -l) | Interrupted:$(cat "$tmpfile" | grep "Interrupted" | uniq | wc -l) | Total:$total_task *****"

    echo >&1000
  } &
  ((bar++))
  processbar $bar $total_task
done
wait

echo -e "\nIntegrating and annotating the matrix....\n"
mkdir -p $maindir/NGSmodule_analysis/Quantification
cd $maindir/NGSmodule_analysis/Quantification

Rscript $2 $work_dir $gtf $Aligner $Species $Source &>Annotation.R.log
echo -e "Integrated quantification matrix: $maindir/NGSmodule_analysis/Quantification/Quantification.${Aligner}.*.tab\n"

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "****************** Quantification Done ******************\n"
