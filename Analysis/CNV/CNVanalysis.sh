#!/usr/bin/bash
trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

#######################################################################################

Rscript &>/dev/null
[ $? -eq 127 ] && {
  color_echo "red" "Cannot find the command Rscript.\n"
  exit 1
}
R_packages=("HMMcopy" "DNAcopy" "scales" "vcfR" "stringr" "dplyr" "ggpubr")
for package in "${R_packages[@]}"; do
  Rscript -e "installed.packages()" | awk '{print $1}' | grep $package &>/dev/null
  [ $? -ne 0 ] && {
    color_echo "red" "Cannot find the R package $package.\n"
    exit 1
  }
done

readCounter --help &>/dev/null
[ $? -ne 0 ] && {
  color_echo "red" "Cannot find the command readCounter. User can install it from 'https://github.com/shahcompbio/hmmcopy_utils'.\n"
  exit 1
}

freec --help &>/dev/null
[ $? -ne 0 ] && {
  color_echo "red" "Cannot find the tool Control-FREEC. User can install it with the command 'conda install -c bioconda control-freec'.\n"
  exit 1
}

bcftools --version &>/dev/null
[ $? -ne 0 ] && {
  color_echo "red" "Cannot find the tool bcftools. User can install it with the command 'conda install -c bioconda bcftools'.\n"
  exit 1
}

$GATK3 --help &>/dev/null
[ $? -ne 0 ] && {
  color_echo "red" "Cannot find the tool GATK3.\n"
  exit 1
}

configureStrelkaGermlineWorkflow.py --version &>/dev/null
[ $? -ne 0 ] && {
  color_echo "red" "Cannot find the tool Strelka2.\n"
  exit 1
}

GC_bin="$iGenomes_Dir/$Species/$Source/$Build/Sequence/GemIndex/Mappability/${Kmer}mer/windows/win$Window/genome.win${Window}.gc.wig"
if [[ ! -f $GC_bin ]]; then
  color_echo "red" "ERROR! Cannot find the wig file containing GC content per bin: ${GC_bin}\n"
  exit 1
fi

Map_bin="$iGenomes_Dir/$Species/$Source/$Build/Sequence/GemIndex/Mappability/${Kmer}mer/windows/win$Window/genome.win$Window.${Kmer}mer.gem.wig"
if [[ ! -f $Map_bin ]]; then
  color_echo "red" "ERROR! Cannot find the wig file containing average mappability per bin: ${Map_bin}\n"
  exit 1
fi

echo -e "############################# CNVanalysis Parameters #############################\n"
echo -e "  Window: ${Window}\n  Kmer: ${Kmer}\n  PloidyAssumed: ${PloidyAssumed}"
echo -e "  GATK3: \"${GATK3}\"\n"
echo -e "################################################################################\n"

echo -e "****************** Start CNVanalysis ******************\n"
SECONDS=0

for sample in "${arr[@]}"; do
  read -u1000
  {
    dir=${work_dir}/${sample}
    mkdir -p $dir/$Aligner
    cd $dir/$Aligner

    force=${force_complete}
    status="uncompleted"
    attempt=0

    echo "===== $sample ====="

    while [[ $status == "uncompleted" ]] && (("$attempt" <= 1)); do
      ((attempt++))
      if [[ $attempt != 1 ]]; then
        echo -e "+++++ ${sample}: Number of attempts: $attempt +++++"
      fi

      logfiles=("HMMcopyStatus.log" "GATK3Status.log" "Strelka2Status.log")
      globalcheck_logfile "$dir/$Aligner" logfiles[@] "$force" "$error_pattern" "$complete_pattern" "$sample"

      ##### HMMcopy #####
      check_logfile "$sample" "HMMcopy" "$dir/$Aligner/CNV/HMMcopy/HMMcopyStatus.log" "$error_pattern" "$complete_pattern" "precheck"
      if [[ $? == 1 ]]; then
        rm -rf $dir/$Aligner/CNV/HMMcopy
        mkdir -p $dir/$Aligner/CNV/HMMcopy
        cd $dir/$Aligner/CNV/HMMcopy
        readCounter -w $Window ${dir}/${Aligner}/${sample}.${Aligner}.dedup.bam >${sample}.${Aligner}.w$Window.wig 2>>HMMcopyStatus.log
        cmd="Rscript $1 ${sample}.${Aligner}.w$Window.wig $GC_bin $Map_bin $PloidyAssumed ${sample}.${Aligner}.HMMcopy"
        echo $cmd >>HMMcopyStatus.log
        eval $cmd &>>HMMcopyStatus.log

        check_logfile "$sample" "HMMcopy" "$dir/$Aligner/CNV/HMMcopy/HMMcopyStatus.log" "$error_pattern" "$complete_pattern" "postcheck"
        if [[ $? == 1 ]]; then
          continue
        fi
      fi

      #####
      ##### BaseqCNV #####
      #mkdir -p $dir/baseqCNV
      #cd $dir/baseqCNV
      ###baseqCNV align -t $threads -1 $fq1 -2 $fq2 -g hg19 -c $baseqCNV_config -o ${sample}.hg19.bam
      #baseqCNV bincount -g hg19 -i $dir/${sample}.hg19.bwamem.rmdup.bam -c $baseqCNV_config -o ${sample}.bincounts
      #baseqCNV normalize -g hg19 -i ${sample}.bincounts -c $baseqCNV_config -o ${sample}.bincounts_norm.txt
      #baseqCNV cbs -i ${sample}.bincounts_norm.txt -o ${sample}.cbs.txt
      #baseqCNV plotgenome -i ${sample}.bincounts_norm.txt -c ${sample}.cbs.txt -o ${sample}.baseqCNV

      #####
      ##### Control-FREEC #####
      #mkdir -p $dir/Control_FREEC
      #cd $dir/Control_FREEC
      #sed "s#PE_sample.bam#$dir/${sample}.hg19.bowtie2.rmdup.bam#" $FREEC_config >FREEC_config.ini
      #sed -i "s#thread#$thread#" $FREEC_config
      #freec -conf FREEC_config.ini
      #makeGraph.R ${sample}.hg19.bowtie2.rmdup.bam_ratio.txt 2 ${sample}.FREEC

      #####
      ##### nQuire estimate ploidy #####
      #nQuire create -f 0 -b ${sample}.bowtie2.hg19.filter.rmdup.bam -x -o ${sample}.nQuire
      #nQuire denoise ${sample}.nQuire.bin -o ${sample}.nQuire.denoise
      #nQuire histo ${sample}.nQuire.denoise.bin > ${sample}.nQuire.histo.txt
      #nQuire lrdmodel -t 88 ${sample}.nQuire.denoise.bin
      #nQuire modeltest ${sample}.nQuire.denoise.bin
      #nQuire histotest ${sample}.nQuire.denoise.bin

      #####
      ##### ploidyNGS estimate ploidy ##### Need large memory
      #mkdir -p $dir/ploidyNGS
      #cd $dir/ploidyNGS
      #ploidyNGS.py --out ${sample}.ploidyNGS --bam $dir/${sample}.hg19.bwa.rmdup.bam

      ##### call SNV
      ### bcftools
      #bcftools mpileup --threads $threads -d 500 -Ou -f $genome ${sample}.bowtie2.hg19.filter.rmdup.bam | bcftools call --threads $threads -mv -Oz | bcftools view -i '%QUAL>=20' -Oz -o ${sample}.calls.vcf

      ## GATK3 #####
      check_logfile "$sample" "GATK3" "$dir/$Aligner/SNV/GATK3/GATK3Status.log" "$error_pattern" "$complete_pattern" "precheck"
      if [[ $? == 1 ]]; then
        rm -rf $dir/$Aligner/SNV/GATK3
        mkdir -p $dir/$Aligner/SNV/GATK3
        cd $dir/$Aligner/SNV/GATK3
        eval "$GATK3 -T HaplotypeCaller -nct $threads -R $genome -I ${dir}/${Aligner}/${sample}.${Aligner}.dedup.bam -o $dir/$Aligner/SNV/GATK3/${sample}.${Aligner}.GATK3.vcf.gz" &>>$dir/$Aligner/SNV/GATK3/GATK3Status.log
        if [[ $? == 0 ]]; then
          bcftools view $dir/$Aligner/SNV/GATK3/${sample}.${Aligner}.GATK3.vcf.gz | bcftools filter -i 'TYPE="snp" && MIN(FORMAT/DP)>=5 && QUAL>=20' -Oz -o $dir/$Aligner/SNV/GATK3/${sample}.${Aligner}.GATK3.filter.vcf.gz &>>$dir/$Aligner/SNV/GATK3/GATK3Status.log
          Rscript $2 $dir/$Aligner/SNV/GATK3/${sample}.${Aligner}.GATK3.filter.vcf.gz $dir/$Aligner/CNV/HMMcopy/${sample}.${Aligner}.HMMcopy ${sample}.${Aligner}.GATK3 &>>$dir/$Aligner/SNV/GATK3/GATK3Status.log
        fi

        check_logfile "$sample" "GATK3" "$dir/$Aligner/SNV/GATK3/GATK3Status.log" "$error_pattern" "$complete_pattern" "postcheck"
        if [[ $? == 1 ]]; then
          continue
        fi
      fi

      ### Strelka2 #####
      check_logfile "$sample" "Strelka2" "$dir/$Aligner/SNV/Strelka2/Strelka2Status.log" "$error_pattern" "$complete_pattern" "precheck"
      if [[ $? == 1 ]]; then
        rm -rf $dir/$Aligner/SNV/Strelka2
        mkdir -p $dir/$Aligner/SNV/Strelka2
        cd $dir/$Aligner/SNV/Strelka2

        configureStrelkaGermlineWorkflow.py \
        --bam ${dir}/${Aligner}/${sample}.${Aligner}.dedup.bam \
        --referenceFasta $genome \
        --runDir $dir/$Aligner/SNV/Strelka2 &>>$dir/$Aligner/SNV/Strelka2/Strelka2Status.log
        $dir/$Aligner/SNV/Strelka2/runWorkflow.py -m local -j $threads &>>$dir/$Aligner/SNV/Strelka2/Strelka2Status.log
        if [[ $? == 0 ]]; then
          bcftools view $dir/$Aligner/SNV/Strelka2/results/variants/variants.vcf.gz | bcftools filter -i 'TYPE="snp" && MIN(FORMAT/DP)>=5 && QUAL>=20' -Oz -o $dir/$Aligner/SNV/Strelka2/${sample}.${Aligner}.Strelka2.filter.vcf.gz &>>$dir/$Aligner/SNV/Strelka2/Strelka2Status.log
          Rscript $2 $dir/$Aligner/SNV/Strelka2/${sample}.${Aligner}.Strelka2.filter.vcf.gz $dir/$Aligner/CNV/HMMcopy/${sample}.${Aligner}.HMMcopy ${sample}.${Aligner}.Strelka2 &>>$dir/$Aligner/SNV/Strelka2/Strelka2Status.log

          if [[ $? != 0 ]]; then
            continue
          else
            echo -e "NGSmodule finished the job[Strelka2]" &>>$dir/$Aligner/SNV/Strelka2/Strelka2Status.log
          fi
        fi

      fi

      status="completed"
      color_echo "blue" "+++++ ${sample}: CNVanalysis completed +++++"
    done

    if [[ "$status" == "completed" ]]; then
      echo "Completed: $sample" >>"$tmpfile"
    else
      echo "Interrupted: $sample" >>"$tmpfile"
      color_echo "red" "ERROR! ${sample} interrupted! Please check the processing log and your raw fastq file."
    fi

    color_echo "green" "***** Completed:$(cat "$tmpfile" | grep "Completed" | uniq | wc -l) | Interrupted:$(cat "$tmpfile" | grep "Interrupted" | uniq | wc -l) | Total:$total_task *****"

    echo >&1000
  } &
  ((bar++))
  processbar $bar $total_task
done
wait

ninterrupted=$(cat "$tmpfile" | grep "Interrupted" | uniq | wc -l)
if [[ $ninterrupted != 0 ]]; then
  cat "$tmpfile" | grep "Interrupted" | uniq >$maindir/CNVanalysis.Interrupted.txt
  color_echo "red" "\n\n################################################################################"
  color_echo "red" "    $ninterrupted of $total_task tasks interrupted."
  color_echo "red" "    Please check the samples in $maindir/CNVanalysis.Interrupted.txt"
  color_echo "red" "################################################################################\n\n"
fi

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "****************** CNVanalysis Done ******************\n"
