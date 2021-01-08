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
reformat.sh &>/dev/null
[ $? -ne 0 ] && {
  color_echo "red" "Cannot find the command reformat.sh.\n"
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
    dir="${work_dir}/${sample}"
    cd "${dir}"
    mkdir -p "${dir}"/PreAlignmentQC

    Layout=${Layout_dict[${sample}]}
    force=${force_complete}
    status="uncompleted"
    attempt=0

    echo "+++++ ${sample} +++++"

    while [[ $status == "uncompleted" ]] && (("$attempt" <= 1)); do
      ((attempt++))
      if [[ $attempt != 1 ]]; then
        echo -e "+++++ ${sample}: Number of attempts: $attempt +++++"
      fi

      ### check existed logs
      logfiles=("fqCheck.log" "fastqc.log" "fastp.log" "fastq_screen.log" "sortmerna.process.log")
      globalcheck_logfile "$dir" logfiles[@] "$force" "$error_pattern" "$complete_pattern" "$sample"
      if [[ -f "$dir"/PreAlignmentQC/fastp/fastp.log ]]; then
        if [[ ! -f "$dir"/PreAlignmentQC/fastq_screen/fastq_screen.log ]]; then
          color_echo "yellow" "+++++ ${sample}: fastq_screen is uncompleted. Will force to run Fastp. +++++"
          rm -f "$dir"/PreAlignmentQC/fastp/fastp.log
        fi
        if [[ $SequenceType == "rna" ]] && ([[ ! -f "$dir"/PreAlignmentQC/sortmerna/sortmerna.result.log ]] || [[ ! -s "$dir"/PreAlignmentQC/sortmerna/sortmerna.result.log ]]); then
          color_echo "yellow" "+++++ ${sample}: SortMeRNA is uncompleted. Will force to run Fastp. +++++"
          rm -f "$dir"/PreAlignmentQC/fastp/fastp.log "$dir"/PreAlignmentQC/sortmerna/sortmerna.process.log
        fi
      fi

      if [[ $Layout == "SE" ]]; then

        if (($(ls "${dir}"/run*_"${sample}".fq.gz | wc -l) == 1)); then
          cp -fa "${dir}"/run1_"${sample}".fq.gz "${dir}"/"${sample}".fq.gz
          echo -e "Fastq files for ${sample} is ready.\n====== ${sample}.fq.gz ======\n${dir}/run1_${sample}.fq.gz" >"$dir"/fqPrepare.log
        else
          runs=$(ls -lL "${dir}"/run*_"${sample}".fq.gz)
          if [[ ! -f ${dir}/${sample}.fq.gz ]] || [[ ! $(echo "${runs[*]}" | awk 'BEGIN{sum=0}{sum+=$5}END{print sum}') == $(ls -lL "${dir}"/"${sample}".fq.gz | awk '{print$5}') ]]; then
            rm -f "$dir"/fqPrepare.log
            echo "${runs[*]}" | awk '{print$9}' | xargs cat >"${dir}"/"${sample}".fq.gz
          fi
          echo -e "Fastq files for ${sample} is ready.\n====== ${sample}.fq.gz ======\n${runs[*]}" >"$dir"/fqPrepare.log
        fi

        fq1=${dir}/${sample}.fq.gz

        check_logfile "$sample" "FastqCheck(raw)" "$dir"/fqCheck.log "$error_pattern" "$complete_pattern" "precheck"
        if [[ $? == 1 ]]; then
          fqCheck_SE "$sample" "$fq1" "$dir"/fqCheck.log
          check_logfile "$sample" "FastqCheck(raw)" "$dir"/fqCheck.log "$error_pattern" "$complete_pattern" "postcheck" $?
          if [[ $? == 1 ]]; then
            force="TRUE"
            continue
          fi
        fi

        # FastQC
        check_logfile "$sample" "FastQC" "$dir"/PreAlignmentQC/fastqc/fastqc.log "$error_pattern" "$complete_pattern" "precheck"
        if [[ $? == 1 ]]; then
          mkdir -p "$dir"/PreAlignmentQC/fastqc
          fastqc -o "$dir"/PreAlignmentQC/fastqc -t "$threads" "${fq1}" &>"$dir"/PreAlignmentQC/fastqc/fastqc.log

          check_logfile "$sample" "FastQC" "$dir"/PreAlignmentQC/fastqc/fastqc.log "$error_pattern" "$complete_pattern" "postcheck" $?
          if [[ $? == 1 ]]; then
            force="TRUE"
            continue
          fi
        fi

        # Fastp
        check_logfile "$sample" "Fastp" "$dir"/PreAlignmentQC/fastp/fastp.log "$error_pattern" "$complete_pattern" "precheck"
        if [[ $? == 1 ]]; then
          mkdir -p "$dir"/PreAlignmentQC/fastp
          fastp --thread "$threads_fastp" --trim_front1 "$trim_front1" --trim_tail1 "$trim_tail1" \
            --qualified_quality_phred "$qualified_quality_phred" --unqualified_percent_limit "$unqualified_percent_limit" \
            "$read_cutting" --cut_window_size "$cut_window_size" --cut_mean_quality "$cut_mean_quality" \
            --low_complexity_filter --trim_poly_x --trim_poly_g --overrepresentation_analysis \
            --length_required "$length_required" \
            --in1 "${fq1}" \
            --out1 "${sample}".fq \
            -j "$dir"/PreAlignmentQC/fastp/"${sample}".fastp.json \
            -h "$dir"/PreAlignmentQC/fastp/"${sample}".fastp.html 2>"$dir"/PreAlignmentQC/fastp/fastp.log

          check_logfile "$sample" "Fastp" "$dir"/PreAlignmentQC/fastp/fastp.log "$error_pattern" "$complete_pattern" "postcheck" $?
          if [[ $? == 1 ]]; then
            force="TRUE"
            continue
          fi
        fi

        fq1=$dir/${sample}.fq

        if [[ -f $fq1 ]]; then
          #FastQ_Screen
          check_logfile "$sample" "FastQ_Screen" "$dir"/PreAlignmentQC/fastq_screen/fastq_screen.log "$error_pattern" "$complete_pattern" "precheck"
          if [[ $? == 1 ]]; then
            mkdir -p "$dir"/PreAlignmentQC/fastq_screen
            fastq_screen --force --Aligner bowtie2 "$FastqScreen_mode" --conf "$FastqScreen_config" --threads "$threads" "$fq1" \
              --outdir "$dir"/PreAlignmentQC/fastq_screen 2>"$dir"/PreAlignmentQC/fastq_screen/fastq_screen.log

            check_logfile "$sample" "FastQ_Screen" "$dir"/PreAlignmentQC/fastq_screen/fastq_screen.log "$error_pattern" "$complete_pattern" "postcheck" $?
            if [[ $? == 1 ]]; then
              force="TRUE"
              continue
            fi
          fi

          #SortMeRNA
          if [[ $SequenceType == "rna" ]]; then
            check_logfile "$sample" "SortMeRNA" "$dir"/PreAlignmentQC/sortmerna/sortmerna.process.log "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
              rm -rf "$dir"/PreAlignmentQC/sortmerna_tmp
              mkdir -p "$dir"/PreAlignmentQC/sortmerna_tmp
              mkdir -p "$dir"/PreAlignmentQC/sortmerna
              sortmerna --ref "${SortmeRNA_ref}" \
                --reads "${fq1}" \
                --threads "$threads" \
                --workdir "$dir"/PreAlignmentQC/sortmerna_tmp \
                --fastx \
                --num_alignments 1 \
                --aligned aligned \
                --other other \
                -v &>"$dir"/PreAlignmentQC/sortmerna/sortmerna.process.log
              check_logfile "$sample" "SortMeRNA" "$dir"/PreAlignmentQC/sortmerna/sortmerna.process.log "$error_pattern" "$complete_pattern" "postcheck" $?
              if [[ $? == 1 ]]; then
                force="TRUE"
                continue
              fi

              mv other.fq "$dir"/"${sample}"_trim.fq
              rm -rf "${sample}".fq aligned.fq "$dir"/PreAlignmentQC/sortmerna_tmp
              mv aligned.log "$dir"/PreAlignmentQC/sortmerna/sortmerna.result.log

            fi

          else
            mv "$fq1" "$dir"/"${sample}"_trim.fq
            color_echo "blue" "+++++ ${sample}: SequenceType is not 'rna'. SortMeRNA skipped. +++++"
          fi

        elif [[ -f "${sample}"_trim.fq ]]; then
          color_echo "blue" "+++++ ${sample}: ${sample}_trim.fq existed. +++++"
        elif [[ -f "${sample}"_trim.fq.gz ]]; then
          color_echo "blue" "+++++ ${sample}: ${sample}_trim.fq.gz existed. +++++"
        else
          color_echo "yellow" "+++++ ${sample}: Cannot find ${sample}.fq or ${sample}_trim.fq or ${sample}_trim.fq.gz . Start a complete preAlignmentQC.+++++"
          force="TRUE"
          continue
        fi

        if [[ -f $dir/${sample}_trim.fq ]]; then
          pigz -p "$threads" -f "$dir"/"${sample}"_trim.fq
        fi

        if [[ -f $dir/${sample}_trim.fq.gz ]]; then
          fqCheck_SE "$sample" "$dir/${sample}_trim.fq.gz" "$dir"/fqCheck.log
          check_logfile "$sample" "FastqCheck(trimmed)" "$dir"/fqCheck.log "$error_pattern" "$complete_pattern" "postcheck" $?
          if [[ $? == 1 ]]; then
            force="TRUE"
            continue
          else
            status="completed"
            color_echo "blue" "+++++ ${sample}: preAlignmentQC completed +++++"
          fi
        else
          force="TRUE"
          color_echo "yellow" "[INFO] ${sample}: trim.fq.gz files not found. Force to do a complete preAlignmentQC. +++++"
        fi

      elif [[ $Layout == "PE" ]]; then

        if (($(ls "${dir}"/run*_"${sample}"_1.fq.gz | wc -l) == 1)); then
          cp -fa "${dir}"/run1_"${sample}"_1.fq.gz "${dir}"/"${sample}"_1.fq.gz
          cp -fa "${dir}"/run1_"${sample}"_2.fq.gz "${dir}"/"${sample}"_2.fq.gz
          echo -e "Fastq files for ${sample} is ready.\n====== ${sample}_1.fq.gz ======\n${dir}/run1_${sample}_1.fq.gz\n====== ${sample}_2.fq.gz ======\n${dir}/run1_${sample}_2.fq.gz" >"$dir"/fqPrepare.log
        else
          runs1=$(ls -lL "${dir}"/run*_"${sample}"_1.fq.gz)
          runs2=$(ls -lL "${dir}"/run*_"${sample}"_2.fq.gz)
          if [[ ! -f ${dir}/${sample}_1.fq.gz ]] || [[ ! $(echo "${runs1[*]}" | awk 'BEGIN{sum=0}{sum+=$5}END{print sum}') == $(ls -lL "${dir}"/"${sample}"_1.fq.gz | awk '{print$5}') ]]; then
            rm -f "$dir"/fqPrepare.log
            echo "${runs1[*]}" | awk '{print$9}' | xargs cat >"${dir}"/"${sample}"_1.fq.gz
          fi
          if [[ ! -f ${dir}/${sample}_2.fq.gz ]] || [[ ! $(echo "${runs2[*]}" | awk 'BEGIN{sum=0}{sum+=$5}END{print sum}') == $(ls -lL "${dir}"/"${sample}"_2.fq.gz | awk '{print$5}') ]]; then
            rm -f "$dir"/fqPrepare.log
            echo "${runs2[*]}" | awk '{print$9}' | xargs cat >"${dir}"/"${sample}"_2.fq.gz
          fi
          echo -e "Fastq files for ${sample} is ready.\n====== ${sample}_1.fq.gz ======\n${runs1[*]}\n====== ${sample}_2.fq.gz ======\n${runs2[*]}" >"$dir"/fqPrepare.log
        fi

        fq1=${dir}/${sample}_1.fq.gz
        fq2=${dir}/${sample}_2.fq.gz

        ##To verify that reads appear to be correctly paired
        check_logfile "$sample" "FastqCheck(raw)" "$dir"/fqCheck.log "$error_pattern" "$complete_pattern" "precheck"
        if [[ $? == 1 ]]; then
          fqCheck_PE "$sample" "$fq1" "$fq2" "$dir"/fqCheck.log
          check_logfile "$sample" "FastqCheck(raw)" "$dir"/fqCheck.log "$error_pattern" "$complete_pattern" "postcheck" $?
          if [[ $? == 1 ]]; then
            force="TRUE"
            continue
          fi
        fi

        # FastQC
        check_logfile "$sample" "FastQC" "$dir"/PreAlignmentQC/fastqc/fastqc.log "$error_pattern" "$complete_pattern" "precheck"
        if [[ $? == 1 ]]; then
          mkdir -p "$dir"/PreAlignmentQC/fastqc
          fastqc -o "$dir"/PreAlignmentQC/fastqc -t "$threads" "${fq1}" "${fq2}" &>"$dir"/PreAlignmentQC/fastqc/fastqc.log

          check_logfile "$sample" "FastQC" "$dir"/PreAlignmentQC/fastqc/fastqc.log "$error_pattern" "$complete_pattern" "postcheck" $?
          if [[ $? == 1 ]]; then
            force="TRUE"
            continue
          fi
        fi

        # Fastp
        check_logfile "$sample" "Fastp" "$dir"/PreAlignmentQC/fastp/fastp.log "$error_pattern" "$complete_pattern" "precheck"
        if [[ $? == 1 ]]; then
          mkdir -p "$dir"/PreAlignmentQC/fastp
          fastp --thread "$threads_fastp" --trim_front1 "$trim_front1" --trim_tail1 "$trim_tail1" --trim_front2 "$trim_front2" --trim_tail2 "$trim_tail2" \
            --qualified_quality_phred "$qualified_quality_phred" --unqualified_percent_limit "$unqualified_percent_limit" \
            "$read_cutting" --cut_window_size "$cut_window_size" --cut_mean_quality "$cut_mean_quality" \
            --low_complexity_filter --trim_poly_x --trim_poly_g --overrepresentation_analysis \
            --length_required "$length_required" --detect_adapter_for_pe --correction \
            --in1 "${fq1}" --in2 "${fq2}" \
            --out1 "${sample}"_1.fq --out2 "${sample}"_2.fq \
            -j "$dir"/PreAlignmentQC/fastp/"${sample}".fastp.json \
            -h "$dir"/PreAlignmentQC/fastp/"${sample}".fastp.html 2>"$dir"/PreAlignmentQC/fastp/fastp.log

          check_logfile "$sample" "Fastp" "$dir"/PreAlignmentQC/fastp/fastp.log "$error_pattern" "$complete_pattern" "postcheck" $?
          if [[ $? == 1 ]]; then
            force="TRUE"
            continue
          fi
        fi

        fq1=$dir/${sample}_1.fq
        fq2=$dir/${sample}_2.fq

        if [[ -f $fq1 ]] && [[ -f $fq2 ]]; then
          #FastQ_Screen
          check_logfile "$sample" "FastQ_Screen" "$dir"/PreAlignmentQC/fastq_screen/fastq_screen.log "$error_pattern" "$complete_pattern" "precheck"
          if [[ $? == 1 ]]; then
            mkdir -p "$dir"/PreAlignmentQC/fastq_screen
            fastq_screen --force --Aligner bowtie2 "$FastqScreen_mode" --conf "$FastqScreen_config" --threads "$threads" "$fq1" "$fq2" \
              --outdir "$dir"/PreAlignmentQC/fastq_screen 2>"$dir"/PreAlignmentQC/fastq_screen/fastq_screen.log

            check_logfile "$sample" "FastQ_Screen" "$dir"/PreAlignmentQC/fastq_screen/fastq_screen.log "$error_pattern" "$complete_pattern" "postcheck" $?
            if [[ $? == 1 ]]; then
              force="TRUE"
              continue
            fi
          fi

          #SortMeRNA
          if [[ $SequenceType == "rna" ]]; then
            check_logfile "$sample" "SortMeRNA" "$dir"/PreAlignmentQC/sortmerna/sortmerna.process.log "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
              rm -rf "$dir"/PreAlignmentQC/sortmerna_tmp
              mkdir -p "$dir"/PreAlignmentQC/sortmerna_tmp
              mkdir -p "$dir"/PreAlignmentQC/sortmerna
              reformat.sh in1="$fq1" in2="$fq2" out="$dir"/"${sample}".fq overwrite=true 2>"$dir"/PreAlignmentQC/sortmerna/reformat_merge.log
              sortmerna --ref "${SortmeRNA_ref}" \
                --reads "${sample}".fq --paired_in \
                --threads "$threads" \
                --workdir "$dir"/PreAlignmentQC/sortmerna_tmp \
                --fastx \
                --num_alignments 1 \
                --aligned aligned \
                --other other \
                -v &>"$dir"/PreAlignmentQC/sortmerna/sortmerna.process.log
              check_logfile "$sample" "SortMeRNA" "$dir"/PreAlignmentQC/sortmerna/sortmerna.process.log "$error_pattern" "$complete_pattern" "postcheck" $?
              if [[ $? == 1 ]]; then
                force="TRUE"
                continue
              fi

              reformat.sh in=other.fq out1="$dir"/"${sample}"_1_trim.fq out2="$dir"/"${sample}"_2_trim.fq overwrite=true 2>"$dir"/PreAlignmentQC/sortmerna/reformat_split.log
              rm -rf "$fq1" "$fq2" "${sample}".fq "${sample}".fq.gz aligned.fq other.fq "$dir"/PreAlignmentQC/sortmerna_tmp
              mv aligned.log "$dir"/PreAlignmentQC/sortmerna/sortmerna.result.log

            fi

          else
            mv "$fq1" "$dir"/"${sample}"_1_trim.fq
            mv "$fq2" "$dir"/"${sample}"_2_trim.fq
            color_echo "blue" "+++++ ${sample}: SequenceType is not 'rna'. SortMeRNA skipped. +++++"
          fi
        elif [[ -f "${sample}"_1_trim.fq ]] && [[ -f "${sample}"_2_trim.fq ]]; then
          color_echo "blue" "+++++ ${sample}: ${sample}_[1,2]_trim.fq existed. +++++"
        elif [[ -f "${sample}"_1_trim.fq.gz ]] && [[ -f "${sample}"_2_trim.fq.gz ]]; then
          color_echo "blue" "+++++ ${sample}: ${sample}_[1,2]_trim.fq.gz existed. +++++"
        else
          color_echo "yellow" "+++++ ${sample}: Cannot find ${sample}_[1,2].fq or ${sample}_[1,2]_trim.fq or ${sample}_[1,2]_trim.fq.gz . Start a complete preAlignmentQC.+++++"
          force="TRUE"
          continue
        fi

        if [[ -f $dir/${sample}_1_trim.fq ]] && [[ -f $dir/${sample}_2_trim.fq ]]; then
          pigz -p "$threads" -f "$dir"/"${sample}"_1_trim.fq "$dir"/"${sample}"_2_trim.fq
        fi

        if [[ -f $dir/${sample}_1_trim.fq.gz ]] && [[ -f $dir/${sample}_2_trim.fq.gz ]]; then
          fqCheck_PE "$sample" "$dir/${sample}_1_trim.fq.gz" "$dir/${sample}_2_trim.fq.gz" "$dir"/fqCheck.log
          check_logfile "$sample" "FastqCheck(trimmed)" "$dir"/fqCheck.log "$error_pattern" "$complete_pattern" "postcheck" $?
          if [[ $? == 1 ]]; then
            force="TRUE"
            continue
          else
            status="completed"
            color_echo "blue" "+++++ ${sample}: preAlignmentQC completed +++++"
          fi
        else
          force="TRUE"
          color_echo "yellow" "+++++ ${sample}: trim.fq.gz files not found. Force to do a complete preAlignmentQC. +++++"
          continue
        fi

      else
        color_echo "yellow" "+++++ ${sample}: Cannot determine the layout of sequencing data! +++++"
        attempt=2
        continue
      fi

    done

    if [[ "$status" == "completed" ]]; then
      echo "Completed: $sample" >>"$tmpfile"
    else
      echo "Interrupted: $sample" >>"$tmpfile"
      color_echo "red" "ERROR! ${sample} interrupted! Please check the processing log and your raw fastq files."
    fi

    color_echo "green" "***** Completed:$(cat "$tmpfile" | grep "Completed" | uniq | wc -l) | Interrupted:$(cat "$tmpfile" | grep "Interrupted" | uniq | wc -l) | Total:$total_task *****"

    echo >&1000
  } &
  ((bar++))
  processbar "$bar" "$total_task"
done
wait

ninterrupted=$(cat "$tmpfile" | grep "Interrupted" | uniq | wc -l)
if [[ $ninterrupted != 0 ]]; then
  cat "$tmpfile" | grep "Interrupted" | uniq >$maindir/preAlignmentQC.Interrupted.txt
  color_echo "red" "\n\n################################################################################"
  color_echo "red" "    $ninterrupted of $total_task tasks interrupted."
  color_echo "red" "    Please check the samples in $maindir/preAlignmentQC.Interrupted.txt"
  color_echo "red" "################################################################################\n\n"
fi

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "\n****************** preAlignmentQC Finished ******************\n"
