#!/usr/bin/env bash

#######################################################################################
##### PreAlignment-QC&Trim.sh include four QC or trimming steps:
##### 1. Do QC for the untrimmed reads with FastQC
##### 2. Trim and filter the reads with fastp
##### 3. Align trimmed reads on multi-genomes to detect the contaminantion with fastq-screen
##### 4. Remove the rRNA with SortmeRNA

trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

fastqc --version &>/dev/null
[ $? -ne 0 ] && {
  color_echo "red" "Cannot find the command fastqc.\n"
  exit 1
}
fastp --version &>/dev/null
[ $? -ne 0 ] && {
  color_echo "red" "Cannot find the command fastp.\n"
  exit 1
}
fastq_screen --version &>/dev/null
[ $? -ne 0 ] && {
  color_echo "red" "Cannot find the command fastq_screen.\n"
  exit 1
}
reformat.sh --version &>/dev/null
[ $? -ne 0 ] && {
  color_echo "red" "Cannot find the command reformat.sh which is a tool in BBmap.\n"
  exit 1
}
sortmerna --version &>/dev/null
[ $? -ne 0 ] && {
  color_echo "red" "Cannot find the command sortmerna.\n"
  exit 1
}
pigz --version &>/dev/null
[ $? -ne 0 ] && {
  color_echo "red" "Cannot find the command pigz.\n"
  exit 1
}
v=$(sortmerna --version 2>/dev/null | grep -oP "(?<=SortMeRNA version ).*" | cut -d. -f1)
((v < 4)) || [[ $v == "" ]] && {
  color_echo "red" "sortmerna vertsion need to be updated to 4.x.x(https://github.com/biocore/sortmerna).\n"
  exit 1
}

force_complete_option=("TRUE" "FALSE")
if [[ " ${force_complete_option[*]} " != *" $force_complete "* ]]; then
  color_echo "red" "ERROR! force_complete must be TRUE or FALSE.\nPlease check theParamaters in your ConfigFile.\n"
  exit 1
fi

if [[ ! -f $SortmeRNA_ref ]]; then
  color_echo "red" "ERROR! Cannot find the SortmeRNA_ref file: ${SortmeRNA_ref}\n"
  exit 1
fi

if [[ ! -f $FastqScreen_config ]]; then
  color_echo "red" "ERROR! Cannot find the FastqScreen_config file: ${FastqScreen_config}\n"
  exit 1
fi

echo -e "########################### preAlignmentQC Parameters ##########################\n"
echo -e "Fastp\n  trim_front1: ${trim_front1}\n  trim_tail1: ${trim_tail1}\n  trim_front2: ${trim_front2}\n  trim_tail2: ${trim_tail2}\n  qualified_quality_phred: ${qualified_quality_phred}\n  unqualified_percent_limit: ${unqualified_percent_limit}\n  read_cutting: ${read_cutting}\n  cut_window_size: ${cut_window_size}\n  cut_mean_quality: ${cut_mean_quality}\n  length_required: ${length_required}\n"
echo -e "FastqScreen\n  FastqScreen_config: ${FastqScreen_config}\n"
echo -e "SortmeRNA\n  SortmeRNA_ref: ${SortmeRNA_ref}\n"
echo -e "################################################################################\n\n\n"

echo -e "****************** Start preAlignmentQC ******************\n"
SECONDS=0

for sample in "${arr[@]}"; do
  read -u1000
  {
    dir=${work_dir}/${sample}
    cd ${dir}
    mkdir -p ${dir}/PreAlignmentQC
    echo "+++++ ${sample} +++++"
    Layout=${Layout_dict[${sample}]}
    force=${force_complete}

    status="uncompleted"
    while [[ $status == "uncompleted" ]]; do

      if [[ $Layout == "SE" ]]; then

        if [[ ! -f $dir/fq.log ]] || [[ $force == "TRUE" ]]; then
          if (($(ls ${dir}/run*_${sample}.fq.gz | wc -l) == 1)); then
            cp -fa ${dir}/run1_${sample}.fq.gz ${dir}/${sample}.fq.gz
            echo -e "Fastq files for ${sample} is ready.\n====== ${sample}.fq.gz ======\n${dir}/run1_${sample}.fq.gz" >$dir/fq.log
          else
            runs=$(ls -lL ${dir}/run*_${sample}.fq.gz)
            while [[ ! -f ${dir}/${sample}.fq.gz ]] || [[ ! $(echo ${runs[*]} | awk 'BEGIN{sum=0}{sum+=$5}END{print sum}') == $(ls -lL ${dir}/${sample}.fq.gz | awk '{print$5}') ]]; do
              echo "${runs[*]}" | awk '{print$9}' | xargs cat >"${dir}"/"${sample}".fq.gz
            done
            echo -e "Fastq files for ${sample} is ready.\n====== ${sample}.fq.gz ======\n${runs[*]}" >$dir/fq.log
          fi
        fi
        fq1=${dir}/${sample}.fq.gz

        check_logfile $sample "FastQC" $dir/PreAlignmentQC/fastqc/fastqc.log "$error_pattern" "$complete_pattern"
        if [[ $? == 1 ]] || [[ $force == "TRUE" ]]; then
          mkdir -p $dir/PreAlignmentQC/fastqc
          fastqc -o $dir/PreAlignmentQC/fastqc -t $threads ${fq1} &>$dir/PreAlignmentQC/fastqc/fastqc.log
        fi
        check_logfile $sample "FastQC" $dir/PreAlignmentQC/fastqc/fastqc.log "$error_pattern" "$complete_pattern"
        if [[ $? == 1 ]]; then
          echo "Interrupted: $sample" >>"$tmpfile"
          break
        fi

        check_logfile $sample "Fastp" $dir/PreAlignmentQC/fastp/fastp.log "$error_pattern" "$complete_pattern"
        if [[ $? == 1 ]] || [[ $force == "TRUE" ]]; then
          mkdir -p $dir/PreAlignmentQC/fastp
          fastp --thread $threads_fastp --trim_front1 $trim_front1 --trim_tail1 $trim_tail1 \
          --qualified_quality_phred $qualified_quality_phred --unqualified_percent_limit $unqualified_percent_limit \
          $read_cutting --cut_window_size $cut_window_size --cut_mean_quality $cut_mean_quality \
          --low_complexity_filter --trim_poly_x --trim_poly_g --overrepresentation_analysis \
          --length_required $length_required \
          --in1 ${fq1} \
          --out1 ${sample}.fq \
          -j $dir/PreAlignmentQC/fastp/${sample}.fastp.json \
          -h $dir/PreAlignmentQC/fastp/${sample}.fastp.html 2>$dir/PreAlignmentQC/fastp/fastp.log
        fi
        check_logfile $sample "Fastp" $dir/PreAlignmentQC/fastp/fastp.log "$error_pattern" "$complete_pattern"
        if [[ $? == 1 ]]; then
          echo "Interrupted: $sample" >>"$tmpfile"
          break
        fi

        fq1=$dir/${sample}.fq

        if [[ -f $fq1 ]]; then

          check_logfile $sample "FastQ_Screen" $dir/PreAlignmentQC/fastq_screen/fastq_screen.log "$error_pattern" "$complete_pattern"
          if [[ $? == 1 ]] || [[ $force == "TRUE" ]]; then
            mkdir -p $dir/PreAlignmentQC/fastq_screen
            fastq_screen --force --Aligner bowtie2 $FastqScreen_mode --conf $FastqScreen_config --threads $threads $fq1 \
            --outdir $dir/PreAlignmentQC/fastq_screen 2>$dir/PreAlignmentQC/fastq_screen/fastq_screen.log
          fi
          check_logfile $sample "FastQ_Screen" $dir/PreAlignmentQC/fastq_screen/fastq_screen.log "$error_pattern" "$complete_pattern"
          if [[ $? == 1 ]]; then
            echo "Interrupted: $sample" >>"$tmpfile"
            break
          fi

          if [[ $SequenceType == "rna" ]]; then
            check_logfile $sample "SortMeRNA" $dir/PreAlignmentQC/sortmerna/sortmerna.process.log "$error_pattern" "$complete_pattern"
            if [[ $? == 1 ]] || [[ $force == "TRUE" ]]; then
              rm -rf $dir/PreAlignmentQC/sortmerna_tmp
              mkdir -p $dir/PreAlignmentQC/sortmerna_tmp
              mkdir -p $dir/PreAlignmentQC/sortmerna
              sortmerna --ref ${SortmeRNA_ref} \
              --reads ${sample}.fq \
              --threads $threads \
              --workdir $dir/PreAlignmentQC/sortmerna_tmp \
              --fastx \
              --num_alignments 1 \
              --aligned aligned \
              --other other \
              -v &>$dir/PreAlignmentQC/sortmerna/sortmerna.process.log
            fi
            check_logfile $sample "SortMeRNA" $dir/PreAlignmentQC/sortmerna/sortmerna.process.log "$error_pattern" "$complete_pattern"
            if [[ $? == 1 ]]; then
              echo "Interrupted: $sample" >>"$tmpfile"
              break
            fi
            mv other.fq $dir/${sample}_trim.fq
            rm -rf $fq1 aligned.fq $dir/PreAlignmentQC/sortmerna_tmp
            mv aligned.log $dir/PreAlignmentQC/sortmerna/sortmerna.log
          else
            mv $fq1 $dir/${sample}_trim.fq
            color_echo "yellow" "+++++ ${sample}: SequenceType='rna'. SortMeRNA skipped. +++++"
          fi

        fi

        if [[ -f $dir/${sample}_trim.fq ]]; then
          pigz -p $threads -f $dir/${sample}_trim.fq
          status="completed"
          color_echo "blue" "+++++ ${sample}: Processing complete +++++"
        elif [[ -f $dir/${sample}_trim.fq.gz ]]; then
          status="completed"
          color_echo "blue" "+++++ ${sample}: Processing complete +++++"
        else
          status="uncompleted"
          force="TRUE"
          color_echo "yellow" "+++++ ${sample}: cannot generate fq.gz files. Start to re-analysis. +++++"
        fi

      elif [[ $Layout == "PE" ]]; then

        if [[ ! -f $dir/fq.log ]] || [[ $force == "TRUE" ]]; then
          if (($(ls ${dir}/run*_"${sample}"_1.fq.gz | wc -l) == 1)); then
            cp -fa "${dir}"/run1_"${sample}"_1.fq.gz "${dir}"/"${sample}"_1.fq.gz
            cp -fa "${dir}"/run1_"${sample}"_2.fq.gz "${dir}"/"${sample}"_2.fq.gz
            echo -e "Fastq files for ${sample} is ready.\n====== ${sample}_1.fq.gz ======\n${dir}/run1_${sample}_1.fq.gz\n====== ${sample}_2.fq.gz ======\n${dir}/run1_${sample}_2.fq.gz" >$dir/fq.log
          else
            runs1=$(ls -lL "${dir}"/run*_"${sample}"_1.fq.gz)
            runs2=$(ls -lL "${dir}"/run*_"${sample}"_2.fq.gz)
            while [[ ! -f ${dir}/${sample}_1.fq.gz ]] || [[ ! $(echo ${runs1[*]} | awk 'BEGIN{sum=0}{sum+=$5}END{print sum}') == $(ls -lL ${dir}/${sample}_1.fq.gz | awk '{print$5}') ]]; do
              echo "${runs1[*]}" | awk '{print$9}' | xargs cat >"${dir}"/"${sample}"_1.fq.gz
            done
            while [[ ! -f ${dir}/${sample}_2.fq.gz ]] || [[ ! $(echo ${runs2[*]} | awk 'BEGIN{sum=0}{sum+=$5}END{print sum}') == $(ls -lL ${dir}/${sample}_2.fq.gz | awk '{print$5}') ]]; do
              echo "${runs2[*]}" | awk '{print$9}' | xargs cat >"${dir}"/"${sample}"_2.fq.gz
            done
            echo -e "Fastq files for ${sample} is ready.\n====== ${sample}_1.fq.gz ======\n${runs1[*]}\n====== ${sample}_2.fq.gz ======\n${runs2[*]}" >$dir/fq.log
          fi
        fi
        fq1=${dir}/${sample}_1.fq.gz
        fq2=${dir}/${sample}_2.fq.gz

        ##To verify that reads appear to be correctly paired
        check_logfile $sample "reformat" $dir/reformat_vpair.log "$error_pattern" "$complete_pattern"
        if [[ $? == 1 ]] || [[ $force == "TRUE" ]]; then
          reformat.sh in1=$fq1 in2=$fq2 vpair allowidenticalnames=t 2>$dir/reformat_vpair.log
        fi

        check_logfile $sample "reformat" $dir/reformat_vpair.log "$error_pattern" "$complete_pattern"
        if [[ $? == 1 ]]; then
          fq1_nlines=$(zcat $fq1 | wc -l)
          fq2_nlines=$(zcat $fq2 | wc -l)
          if [[ $fq1_nlines == "$fq2_nlines" ]]; then
            echo -e "fq1_nlines:$fq1_nlines\nfq2_nlines:$fq2_nlines\nNames appear to be correctly paired(custom)" >>$dir/reformat_vpair.log
          else
            color_echo "red" "ERROR! R1 and R2 for ${sample} have different numbers of reads."
            echo -e "fq1_nlines:$fq1_nlines\nfq2_nlines:$fq2_nlines\n" >>$dir/reformat_vpair.log
            echo -e "ERROR! R1 and R2 for ${sample} have different numbers of reads." >>$dir/reformat_vpair.log
            echo "Interrupted: $sample" >>"$tmpfile"
            break
          fi
        fi

        check_logfile $sample "FastQC" $dir/PreAlignmentQC/fastqc/fastqc.log "$error_pattern" "$complete_pattern"
        if [[ $? == 1 ]] || [[ $force == "TRUE" ]]; then
          mkdir -p $dir/PreAlignmentQC/fastqc
          fastqc -o $dir/PreAlignmentQC/fastqc -t $threads ${fq1} ${fq2} &>$dir/PreAlignmentQC/fastqc/fastqc.log
        fi
        check_logfile $sample "FastQC" $dir/PreAlignmentQC/fastqc/fastqc.log "$error_pattern" "$complete_pattern"
        if [[ $? == 1 ]]; then
          echo "Interrupted: $sample" >>"$tmpfile"
          break
        fi

        check_logfile $sample "Fastp" $dir/PreAlignmentQC/fastp/fastp.log "$error_pattern" "$complete_pattern"
        if [[ $? == 1 ]] || [[ $force == "TRUE" ]]; then
          mkdir -p $dir/PreAlignmentQC/fastp
          fastp --thread $threads_fastp --trim_front1 $trim_front1 --trim_tail1 $trim_tail1 --trim_front2 $trim_front2 --trim_tail2 $trim_tail2 \
          --qualified_quality_phred $qualified_quality_phred --unqualified_percent_limit $unqualified_percent_limit \
          $read_cutting --cut_window_size $cut_window_size --cut_mean_quality $cut_mean_quality \
          --low_complexity_filter --trim_poly_x --trim_poly_g --overrepresentation_analysis \
          --length_required $length_required --detect_adapter_for_pe --correction \
          --in1 ${fq1} --in2 ${fq2} \
          --out1 ${sample}_1.fq --out2 ${sample}_2.fq \
          -j $dir/PreAlignmentQC/fastp/${sample}.fastp.json \
          -h $dir/PreAlignmentQC/fastp/${sample}.fastp.html 2>$dir/PreAlignmentQC/fastp/fastp.log
        fi
        check_logfile $sample "Fastp" $dir/PreAlignmentQC/fastp/fastp.log "$error_pattern" "$complete_pattern"
        if [[ $? == 1 ]]; then
          echo "Interrupted: $sample" >>"$tmpfile"
          break
        fi

        fq1=$dir/${sample}_1.fq
        fq2=$dir/${sample}_2.fq

        if [[ -f $fq1 ]] && [[ -f $fq2 ]]; then

          check_logfile $sample "FastQ_Screen" $dir/PreAlignmentQC/fastq_screen/fastq_screen.log "$error_pattern" "$complete_pattern"
          if [[ $? == 1 ]] || [[ $force == "TRUE" ]]; then
            mkdir -p $dir/PreAlignmentQC/fastq_screen
            fastq_screen --force --Aligner bowtie2 $FastqScreen_mode --conf $FastqScreen_config --threads $threads $fq1 $fq2 \
            --outdir $dir/PreAlignmentQC/fastq_screen 2>$dir/PreAlignmentQC/fastq_screen/fastq_screen.log
          fi
          check_logfile $sample "FastQ_Screen" $dir/PreAlignmentQC/fastq_screen/fastq_screen.log "$error_pattern" "$complete_pattern"
          if [[ $? == 1 ]]; then
            echo "Interrupted: $sample" >>"$tmpfile"
            break
          fi

          if [[ $SequenceType == "rna" ]]; then
            check_logfile $sample "SortMeRNA" $dir/PreAlignmentQC/sortmerna/sortmerna.process.log "$error_pattern" "$complete_pattern"
            if [[ $? == 1 ]] || [[ $force == "TRUE" ]]; then
              rm -rf $dir/PreAlignmentQC/sortmerna_tmp
              mkdir -p $dir/PreAlignmentQC/sortmerna_tmp
              mkdir -p $dir/PreAlignmentQC/sortmerna
              reformat.sh in1=$fq1 in2=$fq2 out=$dir/${sample}.fq overwrite=true 2>$dir/PreAlignmentQC/sortmerna/reformat_merge.log
              sortmerna --ref ${SortmeRNA_ref} \
              --reads ${sample}.fq --paired_in \
              --threads $threads \
              --workdir $dir/PreAlignmentQC/sortmerna_tmp \
              --fastx \
              --num_alignments 1 \
              --aligned aligned \
              --other other \
              -v &>$dir/PreAlignmentQC/sortmerna/sortmerna.process.log
            fi
            check_logfile $sample "SortMeRNA" $dir/PreAlignmentQC/sortmerna/sortmerna.process.log "$error_pattern" "$complete_pattern"
            if [[ $? == 1 ]]; then
              echo "Interrupted: $sample" >>"$tmpfile"
              break
            fi
            reformat.sh in=other.fq out1=$dir/${sample}_1_trim.fq out2=$dir/${sample}_2_trim.fq overwrite=true 2>$dir/PreAlignmentQC/sortmerna/reformat_split.log
            rm -rf $fq1 $fq2 aligned.fq other.fq $dir/PreAlignmentQC/sortmerna_tmp
            mv aligned.log $dir/PreAlignmentQC/sortmerna/sortmerna.log
          else
            mv $fq1 $dir/${sample}_1_trim.fq
            mv $fq2 $dir/${sample}_2_trim.fq
            color_echo "yellow" "+++++ ${sample}: SequenceType='rna'. SortMeRNA skipped. +++++"
          fi

        fi

        if [[ -f ${sample}_1_trim.fq ]] && [[ -f ${sample}_2_trim.fq ]]; then
          pigz -p $threads -f ${sample}_1_trim.fq ${sample}_2_trim.fq
          status="completed"
          color_echo "blue" "+++++ ${sample}: Processing complete +++++"
        elif [[ -f ${sample}_1_trim.fq.gz ]] && [[ -f ${sample}_2_trim.fq.gz ]]; then
          status="completed"
          color_echo "blue" "+++++ ${sample}: Processing complete +++++"
        else
          status="uncompleted"
          force="TRUE"
          color_echo "yellow" "+++++ ${sample}: cannot generate fq.gz files. Start to re-analysis. +++++"
        fi

      fi

    done
    
    if [[ "$status" == "completed" ]];then
      echo "Completed: $sample" >>"$tmpfile"
    fi

    color_echo "green" "***** Completed:$(cat "$tmpfile" | grep "Completed" | uniq | wc -l) | Interrupted:$(cat "$tmpfile" | grep "Interrupted" | uniq | wc -l) | Total:$total_task *****"

    echo >&1000
  } &
  ((bar++))
  processbar $bar $total_task
done
wait

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "****************** preAlignmentQC Done ******************\n"
