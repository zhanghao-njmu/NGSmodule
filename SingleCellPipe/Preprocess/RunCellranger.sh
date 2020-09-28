#!/usr/bin/env bash
trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

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

R_packages=("DropletUtils" "dropestr" "ggplot2" "ggupset" "ggsci" "cowplot" "dplyr" "reshape2" "SingleCellExperiment" "grid" "png" "gridExtra" "scales")
for package in "${R_packages[@]}"; do
    Rscript -e "installed.packages()" | awk '{print $1}' | grep "$package" &>/dev/null
    [ $? -ne 0 ] && {
        color_echo "red" "Cannot find the R package $package.\n"
        exit 1
    }
done

if [[ ! -d $cellranger_ref ]]; then
    color_echo "red" "ERROR! Cannot find the cellranger_ref directory: $cellranger_ref\nPlease check the path.\n"
    exit 1
elif [[ ! -f $gene_gtf ]]; then
    color_echo "red" "ERROR! Cannot find the gene_gtf file: $gtf\nPlease check the path.\n"
    exit 1
elif [[ ! -f $rmsk_gtf ]]; then
    color_echo "red" "ERROR! Cannot find the rmsk_gtf file: $rmsk_gtf\nPlease check the path.\n"
    exit 1
elif [[ ! -f $dropEst_config ]]; then
    color_echo "red" "ERROR! Cannot find the dropEst_config file: $dropEst_config\nPlease check the path.\n"
    exit 1
fi

echo -e "########################### RunCellranger Parameters ###########################\n"
echo -e "  cellranger_ref: ${cellranger_ref}"
echo -e "  gene_gtf: ${gene_gtf}\n  rmsk_gtf: ${rmsk_gtf}\n  dropEst_config: ${dropEst_config}\n"
echo -e "################################################################################\n"

echo -e "****************** Start Alignment ******************\n"
SECONDS=0

for sample in "${arr[@]}"; do
    read -u1000
    {
        dir=$work_dir/$sample
        cd "$dir"

        force=${force_complete}
        status="uncompleted"
        attempt=0

        while [[ $status == "uncompleted" ]] && (("$attempt" <= 1)); do
            ((attempt++))
            if [[ $attempt != 1 ]]; then
                echo -e "+++++ ${sample}: Number of attempts: $attempt +++++"
            fi

            ### clear existed logs
            existlogs=()
            while IFS='' read -r line; do
                existlogs+=("$line")
            done < <(find "$dir" -name "fastqc.log" -o -name "cellranger.log" -o -name "velocyto.log" -o -name "dropEst.log" -o -name "CellCalling.log")
            if ((${#existlogs} >= 1)); then
                for existlog in "${existlogs[@]}"; do
                    if [[ $force == "TRUE" ]] || [[ $(grep -iP "${error_pattern}" "${existlog}") ]] || [[ ! $(grep -iP "${complete_pattern}" "${existlog}") ]]; then
                        color_echo "yellow" "Warning! ${sample}: Detected problems in logfile: ${existlog}."
                        rm -f ${existlog}
                    fi
                done
            fi

            # fastqc
            check_logfile "$sample" "FastQC" "$dir"/PreAlignmentQC/fastqc/fastqc.log "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
                mkdir -p "$dir"/PreAlignmentQC/fastqc
                files=($(find "$dir" -type l | grep -P "$SufixPattern"))
                fastqc -o "$dir"/PreAlignmentQC/fastqc -t "$threads" $(printf " %s" "${files[@]}") &>"$dir"/PreAlignmentQC/fastqc/fastqc.log

                check_logfile "$sample" "FastQC" "$dir/PreAlignmentQC/fastqc/fastqc.log" "$error_pattern" "$complete_pattern" "postcheck"
                if [[ $? == 1 ]]; then
                    continue
                fi
            fi

            # cellranger
            check_logfile "$sample" "cellranger" "$dir"/Alignment/Cellranger/cellranger.log "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
                rm -rf "$dir"/Alignment/Cellranger
                mkdir -p "$dir"/Alignment/Cellranger
                cd "$dir"/Alignment/Cellranger
                sample_run=($(find "$dir" -type l | grep -P "$SufixPattern" | grep -oP "(?<=$dir/).*(?=_S\d+_L\d+)" | sort | uniq))
                sample_run=$(printf ",%s" "${sample_run[@]}")
                sample_run=${sample_run:1}
                echo -e "$sample: $sample_run"

                cellranger count --id "${sample}" \
                --fastqs "${dir}" \
                --sample "${sample_run}" \
                --localcores "$threads" \
                --localmem "$memory" \
                --transcriptome "$cellranger_ref" &>"$dir"/Alignment/Cellranger/cellranger.log

                check_logfile "$sample" "cellranger" "$dir"/Alignment/Cellranger/cellranger.log "$error_pattern" "$complete_pattern" "postcheck"
                if [[ $? == 1 ]]; then
                    continue
                fi
            fi

            # velocyto
            check_logfile "$sample" "velocyto" "$dir"/Alignment/Cellranger/"$sample"/velocyto/velocyto.log "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
                mkdir -p "$dir"/Alignment/Cellranger/"$sample"/velocyto
                cd "$dir"/Alignment/Cellranger/"$sample"/velocyto
                velocyto run10x -m "$rmsk_gtf" --samtools-threads "$threads" "$dir"/Alignment/Cellranger/"$sample" "$gene_gtf" &>"$dir"/Alignment/Cellranger/"$sample"/velocyto/velocyto.log

                check_logfile "$sample" "velocyto" "$dir"/Alignment/Cellranger/"$sample"/velocyto/velocyto.log "$error_pattern" "$complete_pattern" "postcheck"
                if [[ $? == 1 ]]; then
                    continue
                fi
            fi

            # dropEst
            check_logfile "$sample" "dropEst" "$dir"/Alignment/Cellranger/"$sample"/dropEst/dropEst.log "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
                mkdir -p "$dir"/Alignment/Cellranger/"$sample"/dropEst
                cd "$dir"/Alignment/Cellranger/"$sample"/dropEst
                dropest -f -g "$gene_gtf" -c "$dropEst_config" "$dir"/Alignment/Cellranger/"$sample"/outs/possorted_genome_bam.bam &>"$dir"/Alignment/Cellranger/"$sample"/dropEst/dropEst.log

                check_logfile "$sample" "dropEst" "$dir"/Alignment/Cellranger/"$sample"/dropEst/dropEst.log "$error_pattern" "$complete_pattern" "postcheck"
                if [[ $? == 1 ]]; then
                    continue
                fi
            fi

            check_logfile "$sample" "dropReport" "$dir"/Alignment/Cellranger/"$sample"/dropEst/dropReport.log "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
                mkdir -p "$dir"/Alignment/Cellranger/"$sample"/dropEst
                cd "$dir"/Alignment/Cellranger/"$sample"/dropEst
                dropReport.Rsc "$dir"/Alignment/Cellranger/"$sample"/dropEst/cell.counts.rds &>"$dir"/Alignment/Cellranger/"$sample"/dropEst/dropReport.log

                check_logfile "$sample" "dropReport" "$dir"/Alignment/Cellranger/"$sample"/dropEst/dropReport.log "$error_pattern" "$complete_pattern" "postcheck"
                if [[ $? == 1 ]]; then
                    continue
                fi
            fi

            # Cell-calling
            check_logfile "$sample" "CellCalling" "$dir"/Alignment/Cellranger/"$sample"/CellCalling/CellCalling.log "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
                mkdir -p "$dir"/Alignment/Cellranger/"$sample"/CellCalling
                cd "$dir"/Alignment/Cellranger/"$sample"/CellCalling
                Rscript "$1" "$dir"/Alignment/Cellranger "$sample" "$threads" &>"$dir"/Alignment/Cellranger/"$sample"/CellCalling/CellCalling.log

                check_logfile "$sample" "CellCalling" "$dir"/Alignment/Cellranger/"$sample"/CellCalling/CellCalling.log "$error_pattern" "$complete_pattern" "postcheck"
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
echo -e "****************** Alignment Done ******************\n"
