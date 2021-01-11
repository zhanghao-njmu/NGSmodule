#!/usr/bin/env bash
trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

pigz --version &>/dev/null
[ $? -ne 0 ] && {
    color_echo "red" "Cannot find the command pigz.\n"
    exit 1
}
fastqc --version &>/dev/null
[ $? -ne 0 ] && {
    color_echo "red" "Cannot find the command fastqc.\n"
    exit 1
}
fastq_screen --version &>/dev/null
[ $? -ne 0 ] && {
    color_echo "red" "Cannot find the command fastq_screen.\n"
    exit 1
}
cellranger &>/dev/null
[ $? -eq 127 ] && {
    echo -e "Cannot find the command cellranger.\n"
    exit 1
}
velocyto &>/dev/null
[ $? -eq 127 ] && {
    echo -e "Cannot find the command velocyto.\n"
    exit 1
}
dropest &>/dev/null
[ $? -eq 127 ] && {
    echo -e "Cannot find the command dropest.\n"
    exit 1
}
dropReport.Rsc &>/dev/null
[ $? -eq 127 ] && {
    echo -e "Cannot find the command dropReport.Rsc.\n"
    exit 1
}
Rscript &>/dev/null
[ $? -eq 127 ] && {
    echo -e "Cannot find the command Rscript.\n"
    exit 1
}
force_complete_option=("TRUE" "FALSE")
if [[ " ${force_complete_option[*]} " != *" $force_complete "* ]]; then
    color_echo "red" "ERROR! force_complete must be TRUE or FALSE.\nPlease check theParamaters in your SCPConfigFile.\n"
    exit 1
fi

R_packages=("DropletUtils" "dropestr" "ggplot2" "ggupset" "ggsci" "cowplot" "dplyr" "reshape2" "SingleCellExperiment" "grid" "png" "gridExtra" "scales" "inflection")
for package in "${R_packages[@]}"; do
    Rscript -e "installed.packages()" | awk '{print $1}' | grep "$package" &>/dev/null
    [ $? -ne 0 ] && {
        color_echo "red" "Cannot find the R package $package.\n"
        exit 1
    }
done

if [[ ! -f $FastqScreen_config ]]; then
    color_echo "red" "ERROR! Cannot find the FastqScreen_config file: ${FastqScreen_config}\nPlease check the path in your SCPConfigFile.\n"
    exit 1
elif [[ ! -d $cellranger_ref ]]; then
    color_echo "red" "ERROR! Cannot find the cellranger_ref directory: $cellranger_ref\nPlease check the path in your SCPConfigFile.\n"
    exit 1
elif [[ ! -f $gene_gtf ]]; then
    color_echo "red" "ERROR! Cannot find the gene_gtf file: $gene_gtf\nPlease check the path in your SCPConfigFile.\n"
    exit 1
elif [[ ! -f $rmsk_gtf ]]; then
    color_echo "red" "ERROR! Cannot find the rmsk_gtf file: $rmsk_gtf\nPlease check the path in your SCPConfigFile.\n"
    exit 1
elif [[ ! -f $dropEst_config ]]; then
    color_echo "red" "ERROR! Cannot find the dropEst_config file: $dropEst_config\nPlease check the path in your SCPConfigFile.\n"
    exit 1
fi

echo -e "########################### RunCellranger Parameters ###########################\n"
echo -e "  FastqScreen\n  FastqScreen_config: ${FastqScreen_config}\n"
echo -e "  cellranger_ref: ${cellranger_ref}"
echo -e "  gene_gtf: ${gene_gtf}\n  rmsk_gtf: ${rmsk_gtf}\n  dropEst_config: ${dropEst_config}\n"
echo -e "################################################################################\n\n"

echo -e "****************** Start RunCellranger ******************\n"
SECONDS=0
SCP_path=$1

for sample in "${arr[@]}"; do
    read -u1000
    {
        dir=$work_dir/$sample
        cd "$dir"

        force=${force_complete}
        status="uncompleted"
        attempt=0

        echo "+++++ ${sample} +++++"

        while [[ $status == "uncompleted" ]] && (("$attempt" <= 1)); do
            ((attempt++))
            if [[ $attempt != 1 ]]; then
                echo -e "+++++ ${sample}: Number of attempts: $attempt +++++"
            fi

            ### clear existed logs
            logfiles=("fqCheck.log" "fastqc.log" "fastq_screen.log" "cellranger.log" "velocyto.log" "dropEst.log" "CellCalling.log")
            globalcheck_logfile "$dir" logfiles[@] "$force" "$error_pattern" "$complete_pattern" "$sample"

            if (($(ls "${dir}"/run*_"${sample}"_S1_L001_R1_001.fastq.gz | wc -l) == 1)); then
                cp -fa "${dir}"/run1_"${sample}"_S1_L001_R1_001.fastq.gz "${dir}"/"${sample}"_S1_L001_R1_001.fastq.gz
                cp -fa "${dir}"/run1_"${sample}"_S1_L001_R2_001.fastq.gz "${dir}"/"${sample}"_S1_L001_R2_001.fastq.gz
                echo -e "Fastq files for ${sample} is ready.\n====== ${sample}_S1_L001_R1_001.fastq.gz ======\n${dir}/run1_${sample}_S1_L001_R1_001.fastq.gz\n====== ${sample}_S1_L001_R2_001.fastq.gz ======\n${dir}/run1_${sample}_S1_L001_R2_001.fastq.gz" >"$dir"/fqPrepare.log
            else
                runs1=$(ls -lL "${dir}"/run*_"${sample}"_S1_L001_R1_001.fastq.gz)
                runs2=$(ls -lL "${dir}"/run*_"${sample}"_S1_L001_R2_001.fastq.gz)
                if [[ ! -f ${dir}/${sample}_S1_L001_R1_001.fastq.gz ]] || [[ ! $(echo "${runs1[*]}" | awk 'BEGIN{sum=0}{sum+=$5}END{print sum}') == $(ls -lL "${dir}"/"${sample}"_S1_L001_R1_001.fastq.gz | awk '{print$5}') ]]; then
                    rm -f "$dir"/fqPrepare.log
                    echo "${runs1[*]}" | awk '{print$9}' | xargs cat >"${dir}"/"${sample}"_S1_L001_R1_001.fastq.gz
                fi
                if [[ ! -f ${dir}/${sample}_S1_L001_R2_001.fastq.gz ]] || [[ ! $(echo "${runs2[*]}" | awk 'BEGIN{sum=0}{sum+=$5}END{print sum}') == $(ls -lL "${dir}"/"${sample}"_S1_L001_R2_001.fastq.gz | awk '{print$5}') ]]; then
                    rm -f "$dir"/fqPrepare.log
                    echo "${runs2[*]}" | awk '{print$9}' | xargs cat >"${dir}"/"${sample}"_S1_L001_R2_001.fastq.gz
                fi
                echo -e "Fastq files for ${sample} is ready.\n====== ${sample}_S1_L001_R1_001.fastq.gz ======\n${runs1[*]}\n====== ${sample}_S1_L001_R2_001.fastq.gz ======\n${runs2[*]}" >"$dir"/fqPrepare.log
            fi

            fq1=${dir}/${sample}_S1_L001_R1_001.fastq.gz
            fq2=${dir}/${sample}_S1_L001_R2_001.fastq.gz

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

                check_logfile "$sample" "FastQC" "$dir"/PreAlignmentQC/fastqc/fastqc.log "$error_pattern" "$complete_pattern" "postcheck"
                if [[ $? == 1 ]]; then
                    force="TRUE"
                    continue
                fi
            fi

            #FastQ_Screen
            check_logfile "$sample" "FastQ_Screen" "$dir"/PreAlignmentQC/fastq_screen/fastq_screen.log "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
                mkdir -p "$dir"/PreAlignmentQC/fastq_screen
                fastq_screen --force --Aligner bowtie2 "$FastqScreen_mode" --conf "$FastqScreen_config" --threads "$threads" "$fq2" \
                    --outdir "$dir"/PreAlignmentQC/fastq_screen 2>"$dir"/PreAlignmentQC/fastq_screen/fastq_screen.log

                check_logfile "$sample" "FastQ_Screen" "$dir"/PreAlignmentQC/fastq_screen/fastq_screen.log "$error_pattern" "$complete_pattern" "postcheck"
                if [[ $? == 1 ]]; then
                    force="TRUE"
                    continue
                fi
            fi

            # cellranger
            check_logfile "$sample" "cellranger" "$dir"/Alignment-Cellranger/cellranger.log "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
                rm -rf "$dir"/Alignment-Cellranger
                mkdir -p "$dir"/Alignment-Cellranger
                cd "$dir"/Alignment-Cellranger
                # sample_run=($(find "$dir" -type l | grep -P "$SufixPattern" | grep -oP "(?<=$dir/).*(?=_S\d+_L\d+)" | sort | uniq))
                # sample_run=$(printf ",%s" "${sample_run[@]}")
                # sample_run=${sample_run:1}
                # echo -e "$sample: $sample_run"

                cellranger count --id "${sample}" \
                    --fastqs "${dir}" \
                    --sample "${sample}" \
                    --include-introns \
                    --disable-ui \
                    --localcores "$threads" \
                    --localmem "$memory" \
                    --transcriptome "$cellranger_ref" &>"$dir"/Alignment-Cellranger/cellranger.log

                check_logfile "$sample" "cellranger" "$dir"/Alignment-Cellranger/cellranger.log "$error_pattern" "$complete_pattern" "postcheck"
                if [[ $? == 1 ]]; then
                    continue
                fi
            fi

            # velocyto
            # check_logfile "$sample" "velocyto" "$dir"/Alignment-Cellranger/"$sample"/velocyto/velocyto.log "$error_pattern" "$complete_pattern" "precheck"
            # if [[ $? == 1 ]]; then
            #     rm -rf "$dir"/Alignment-Cellranger/"$sample"/velocyto
            #     mkdir -p "$dir"/Alignment-Cellranger/"$sample"/velocyto
            #     cd "$dir"/Alignment-Cellranger/"$sample"/velocyto
            #     #velocyto run10x -m "$rmsk_gtf" --samtools-threads "$threads" "$dir"/Alignment-Cellranger/"$sample" "$gene_gtf" &>"$dir"/Alignment-Cellranger/"$sample"/velocyto/velocyto.log
            #     cp -f "$dir"/Alignment-Cellranger/"$sample"/outs/raw_feature_bc_matrix/barcodes.tsv.gz "$dir"/Alignment-Cellranger/"$sample"/velocyto/barcodes.tsv.gz
            #     gunzip -f "$dir"/Alignment-Cellranger/"$sample"/velocyto/barcodes.tsv.gz
            #     velocyto run -b "$dir"/Alignment-Cellranger/"$sample"/velocyto/barcodes.tsv -o "$dir"/Alignment-Cellranger/"$sample"/velocyto -m "$rmsk_gtf" --samtools-threads "$threads" "$dir"/Alignment-Cellranger/"$sample"/outs/possorted_genome_bam.bam "$gene_gtf" &>"$dir"/Alignment-Cellranger/"$sample"/velocyto/velocyto.log
            #     rm -f "$dir"/Alignment-Cellranger/"$sample"/velocyto/barcodes.tsv

            #     check_logfile "$sample" "velocyto" "$dir"/Alignment-Cellranger/"$sample"/velocyto/velocyto.log "$error_pattern" "$complete_pattern" "postcheck"
            #     if [[ $? == 1 ]]; then
            #         continue
            #     fi
            # fi

            # dropEst
            check_logfile "$sample" "dropEst" "$dir"/Alignment-Cellranger/"$sample"/dropEst/dropEst.log "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
                rm -rf "$dir"/Alignment-Cellranger/"$sample"/dropEst
                mkdir -p "$dir"/Alignment-Cellranger/"$sample"/dropEst
                cd "$dir"/Alignment-Cellranger/"$sample"/dropEst
                dropest -V -f -g "$gene_gtf" -c "$dropEst_config" "$dir"/Alignment-Cellranger/"$sample"/outs/possorted_genome_bam.bam &>"$dir"/Alignment-Cellranger/"$sample"/dropEst/dropEst.log

                check_logfile "$sample" "dropEst" "$dir"/Alignment-Cellranger/"$sample"/dropEst/dropEst.log "$error_pattern" "$complete_pattern" "postcheck"
                if [[ $? == 1 ]]; then
                    continue
                fi
            fi

            check_logfile "$sample" "dropReport" "$dir"/Alignment-Cellranger/"$sample"/dropEst/dropReport.log "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
                mkdir -p "$dir"/Alignment-Cellranger/"$sample"/dropEst
                cd "$dir"/Alignment-Cellranger/"$sample"/dropEst
                dropReport.Rsc "$dir"/Alignment-Cellranger/"$sample"/dropEst/cell.counts.rds &>"$dir"/Alignment-Cellranger/"$sample"/dropEst/dropReport.log

                check_logfile "$sample" "dropReport" "$dir"/Alignment-Cellranger/"$sample"/dropEst/dropReport.log "$error_pattern" "$complete_pattern" "postcheck"
                if [[ $? == 1 ]]; then
                    continue
                fi
            fi

            # Cell-calling
            check_logfile "$sample" "CellCalling" "$dir"/Alignment-Cellranger/"$sample"/CellCalling/CellCalling.log "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
                mkdir -p "$dir"/Alignment-Cellranger/"$sample"/CellCalling
                cd "$dir"/Alignment-Cellranger/"$sample"/CellCalling
                script=$SCP_path/SCP-GeneralSteps/RunCellranger/CellCalling.R
                Rscript "$script" "$dir"/Alignment-Cellranger "$sample" "$threads" "$EmptyThreshold" "$CellLabel" &>"$dir"/Alignment-Cellranger/"$sample"/CellCalling/CellCalling.log

                check_logfile "$sample" "CellCalling" "$dir"/Alignment-Cellranger/"$sample"/CellCalling/CellCalling.log "$error_pattern" "$complete_pattern" "postcheck"
                if [[ $? == 1 ]]; then
                    continue
                fi
            fi

            status="completed"
        done

        if [[ "$status" == "completed" ]]; then
            echo "Completed: $sample" >>"$tmpfile"
        else
            echo "Interrupted: $sample" >>"$tmpfile"
            color_echo "red" "ERROR! ${sample} interrupted! Please check the processing log."
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
echo -e "****************** RunCellranger Done ******************\n"
