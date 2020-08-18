#!/usr/bin/env bash

#######################################################################################
trap_add 'trap - SIGTERM && kill -- -$$;kill $(jobs -pr);kill 0' SIGINT SIGTERM

bam_stat.py &>/dev/null
[ $? -eq 127 ] && {
  color_echo "red" "Cannot find the package RSeQC. User can install RSeQC by 'conda install -c bioconda rseqc'.\n"
  exit 1
}
preseq &>/dev/null
[ $? -eq 127 ] && {
  color_echo "red" "Cannot find the package preseq. User can install preseq by 'conda install -c bioconda preseq'.\n"
  exit 1
}
goleft &>/dev/null
[ $? -eq 127 ] && {
  color_echo "red" "Cannot find the command goleft. User can install goleft by 'conda install -c bioconda goleft'.\n"
  exit 1
}
mosdepth &>/dev/null
[ $? -eq 127 ] && {
  color_echo "red" "Cannot find the command mosdepth. User can install mosdepth by 'conda install -c bioconda mosdepth'.\n"
  exit 1
}

Rscript &>/dev/null
[ $? -eq 127 ] && {
  color_echo "red" "Cannot find the command Rscript.\n"
  exit 1
}

R_packages=("dupRadar" "parallel")
for package in "${R_packages[@]}"; do
  Rscript -e "installed.packages()" | awk '{print $1}' | grep $package &>/dev/null
  [ $? -ne 0 ] && {
    color_echo "red" "Cannot find the R package $package.\n"
    exit 1
  }
done

if [[ ! -f $gtf ]]; then
  color_echo "red" "ERROR! Cannot find the gtf file: $gtf\nPlease check the paramaters in your ConfigFile.\n"
  exit 1
fi

echo -e "########################## postAlignmentQC Parameters ##########################\n"
echo -e "  Rscript path: $(which Rscript)"
echo -e "  SequenceType: ${SequenceType}\n  Aligner: ${Aligner}\n"
echo -e "  GTF_File: ${gtf}\n "
echo -e "################################################################################\n"

echo -e "****************** Start postAlignmentQC ******************\n"
SECONDS=0

genes_bed=${gtf%/*}/genes.bed
if [[ ! -f $genes_bed ]]; then
  gtfToGenePred $gtf tmp.genePhred && genePredToBed tmp.genePhred $genes_bed && rm tmp.genePhred
fi

for sample in "${arr[@]}"; do
  read -u1000
  {
    echo "+++++ $sample +++++"
    Layout=${Layout_dict[$sample]}

    dir=$work_dir/$sample
    if [[ "$SequenceType" == "BSdna" ]] && [[ "$Aligner" =~ bismark_* ]]; then
      bam=$(ls $dir/$Aligner/*.bam)
    else
      bam=$dir/$Aligner/${sample}.${Aligner}.bam
    fi
    if [[ ! -f $bam ]]; then
      color_echo "red" "ERROR! Bam file:$bam do not exist. Please check the file.\n"
      exit 1
    fi

    mkdir -p $dir/$Aligner/postAlignmentQC/RSeQC
    cd $dir/$Aligner/postAlignmentQC/RSeQC
    bam_stat.py -i $bam >${sample}.${Aligner}.bam_stat.txt 2>bam_stat.log
    infer_experiment.py -r $genes_bed -i $bam >${sample}.${Aligner}.infer_experiment.txt 2>infer_experiment.log
    inner_distance.py -r $genes_bed -i $bam -o ${sample}.${Aligner} &>inner_distance.log
    read_distribution.py -r $genes_bed -i $bam >${sample}.${Aligner}.read_distribution.txt 2>read_distribution.log
    read_duplication.py -i $bam -o ${sample}.${Aligner} &>read_duplication.log
    read_GC.py -i $bam -o ${sample}.${Aligner} &>read_GC.log
    if [[ $SequenceType == "rna" ]]; then
      geneBody_coverage.py -r $genes_bed -i $bam -o ${sample}.${Aligner} &>geneBody_coverage.log
      junction_annotation.py -r $genes_bed -i $bam -o ${sample}.${Aligner} &>${sample}.${Aligner}.log
      junction_saturation.py -r $genes_bed -i $bam -o ${sample}.${Aligner} &>junction_saturation.log
    fi

    if [[ "$SequenceType" == "BSdna" ]] && [[ "$Aligner" =~ bismark_* ]]; then
      echo "+++++ Waiting for the background processes..........+++++"
    else
      mkdir -p $dir/$Aligner/postAlignmentQC/Preseq
      cd $dir/$Aligner/postAlignmentQC/Preseq
      preseq lc_extrap -B $bam -o ${sample}.${Aligner}.txt &>preseq_lc_extrap.log
      mkdir -p $dir/$Aligner/postAlignmentQC/goleft
      cd $dir/$Aligner/postAlignmentQC/goleft
      goleft indexcov --directory ./ $bam &>goleft_indexcov.log
      mkdir -p $dir/$Aligner/postAlignmentQC/mosdepth
      cd $dir/$Aligner/postAlignmentQC/mosdepth
      mosdepth -t $threads -n --fast-mode ${sample}.${Aligner} $bam &>mosdepth.log
      mkdir -p $dir/$Aligner/postAlignmentQC/dupRadar
      cd $dir/$Aligner/postAlignmentQC/dupRadar
      Rscript $1 $bam $gtf $strandspecific $Layout $threads_featurecounts $dir/$Aligner/postAlignmentQC/dupRadar ${sample}.${Aligner}
    fi

    ##  mkdir -p $dir/$Aligner/Qualimap; cd $dir/$Aligner/Qualimap #### too slow!!!
    ##  unset DISPLAY
    ##  if [[ $Layout == "SE" && $SequenceType == "rna" ]];then
    ##    qualimap rnaseq -bam $bam -gtf $gtf -outdir ${sample}.${Aligner} -outformat HTML --java-mem-size=10G  &
    ##  elif [[ $Layout == "PE" && $SequenceType == "rna" ]];then
    ##    qualimap rnaseq -bam $bam -gtf $gtf --paired -outdir ${sample}.${Aligner} -outformat HTML --java-mem-size=10G &
    ##  fi
    ##  wait
    ##  qualimap bamqc -bam $bam -gff $gtf -nt $threads -nr 100000 -nw 300 -outdir ${sample}.${Aligner} -outformat HTML java_options="-Djava.awt.headless=true -Xmx$JAVA_MEM_SIZE -XX:MaxPermSize=10G"
    #
    ##  mkdir -p $dir/$Aligner/deepTools; cd $dir/$Aligner/deepTools
    ##  bamPEFragmentSize -p $threads -b $bam  --table ${sample}.${Aligner}.table.txt --outRawFragmentLengths ${sample}.${Aligner}.RawFragmentLengths.txt
    ##  estimateReadFiltering -p $threads -b $bam >${sample}.${Aligner}.estimateRead.txt
    ##  plotCoverage -p $threads -b $bam --outRawCounts ${sample}.${Aligner}.CoverageRawCounts.txt
    ##  plotFingerprint -p $threads -b $bam --outRawCounts ${sample}.${Aligner}.FingerprintRawCounts.txt  --outQualityMetrics ${sample}.${Aligner}.FingerprintQualityMetrics.txt

    wait

    echo "Completed: $sample" >>$tmpfile
    color_echo "green" "***** Completed:$(cat "$tmpfile" | grep "Completed" | uniq | wc -l) | Interrupted:$(cat "$tmpfile" | grep "Interrupted" | uniq | wc -l) | Total:$total_task *****"

    echo >&1000
  } &
  ((bar++))
  processbar $bar $total_task
done
wait

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "****************** postAlignmentQC Done ******************\n"
