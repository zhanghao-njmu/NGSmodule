#!/usr/bin/env bash
trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

if [[ $mode == "rna" ]]; then
    cellranger &>/dev/null
    [ $? -eq 127 ] && {
        echo -e "Cannot find the command cellranger.\n"
        exit 1
    }
elif [[ $mode == "atac" ]]; then
    cellranger-atac &>/dev/null
    [ $? -eq 127 ] && {
        echo -e "Cannot find the command cellranger-atac.\n"
        exit 1
    }
elif [[ $mode == "arc" ]]; then
    cellranger-arc &>/dev/null
    [ $? -eq 127 ] && {
        echo -e "Cannot find the command cellranger-arc.\n"
        exit 1
    }
fi

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

R_packages=("SCP" "Seurat" "SeuratWrappers" "future" "Signac" "rtracklayer" "Rsamtools" "GenomicRanges")
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
fi

genome="$cellranger_ref/fasta/genome.fa"
if [[ -f $gene_gtf.gz ]]; then
    gunzip -f $gene_gtf.gz
fi
if [[ ! -f $gene_gtf ]]; then
    color_echo "red" "ERROR! Cannot find the gene_gtf file: $gene_gtf\nPlease check the path in your SCPConfigFile.\n"
    exit 1
fi

if [[ $mode == "rna" ]]; then
    if [[ ! -f $rmsk_gtf ]]; then
        color_echo "red" "ERROR! Cannot find the rmsk_gtf file: $rmsk_gtf\nPlease check the path in your SCPConfigFile.\n"
        exit 1
    fi
fi

echo -e "########################### RunCellranger Parameters ###########################\n"
echo -e "  FastqScreen_config: ${FastqScreen_config}\n"
echo -e "  cellranger_ref: ${cellranger_ref}"
echo -e "  gene_gtf: ${gene_gtf}\n  rmsk_gtf: ${rmsk_gtf}\n"
echo -e "################################################################################\n\n"

echo -e "****************** Start RunCellranger ******************\n"
SECONDS=0
SCP_path=$1

for sample in "${arr[@]}"; do
    read -u1000
    {
        dir="$work_dir/$sample"
        cd "$dir"

        force=${force_complete}
        status="uncompleted"
        attempt=0

        echo "+++++ ${sample} +++++"

        while [[ $status == "uncompleted" ]] && (("$attempt" <= $retry)); do
            ((attempt++))
            if [[ $attempt != 1 ]]; then
                echo -e "+++++ ${sample}: Number of retries: $attempt +++++"
            fi

            ### clear existed logs
            logfiles=("fqCheck.log" "fastqc.log" "fastq_screen.log" "cellranger.log" "velocyto.log" "dropEst.log" "CellCalling.log")
            globalcheck_logfile "$dir" logfiles[@] "$force" "$error_pattern" "$complete_pattern" "$sample"

            if [[ $mode == "rna" ]]; then
                if (($(ls "${dir}"/run*_"${sample}"_S1_L001_R1_001.fastq.gz | wc -l) == 1)); then
                    cp -fa "${dir}"/run1_"${sample}"_S1_L001_R1_001.fastq.gz "${dir}"/"${sample}"_S1_L001_R1_001.fastq.gz
                    cp -fa "${dir}"/run1_"${sample}"_S1_L001_R2_001.fastq.gz "${dir}"/"${sample}"_S1_L001_R2_001.fastq.gz
                    echo -e "Fastq files for ${sample} is ready.\n====== ${sample}_S1_L001_R1_001.fastq.gz ======\n${dir}/run1_${sample}_S1_L001_R1_001.fastq.gz\n====== ${sample}_S1_L001_R2_001.fastq.gz ======\n${dir}/run1_${sample}_S1_L001_R2_001.fastq.gz" >"$dir"/fqPrepare.log
                else
                    runs1=$(ls -lL "${dir}"/run*_"${sample}"_S1_L001_R1_001.fastq.gz)
                    runs2=$(ls -lL "${dir}"/run*_"${sample}"_S1_L001_R2_001.fastq.gz)
                    if [[ ! -f "${dir}"/"${sample}_S1_L001_R1_001.fastq.gz" ]] || [[ ! $(echo "${runs1[*]}" | awk 'BEGIN{sum=0}{sum+=$5}END{print sum}') == $(ls -lL "${dir}"/"${sample}"_S1_L001_R1_001.fastq.gz | awk '{print$5}') ]]; then
                        rm -f "$dir"/fqPrepare.log
                        echo "${runs1[*]}" | awk '{print$9}' | xargs cat >"${dir}"/"${sample}"_S1_L001_R1_001.fastq.gz
                    fi
                    if [[ ! -f "${dir}"/"${sample}_S1_L001_R2_001.fastq.gz" ]] || [[ ! $(echo "${runs2[*]}" | awk 'BEGIN{sum=0}{sum+=$5}END{print sum}') == $(ls -lL "${dir}"/"${sample}"_S1_L001_R2_001.fastq.gz | awk '{print$5}') ]]; then
                        rm -f "$dir"/fqPrepare.log
                        echo "${runs2[*]}" | awk '{print$9}' | xargs cat >"${dir}"/"${sample}"_S1_L001_R2_001.fastq.gz
                    fi
                    echo -e "Fastq files for ${sample} is ready.\n====== ${sample}_S1_L001_R1_001.fastq.gz ======\n${runs1[*]}\n====== ${sample}_S1_L001_R2_001.fastq.gz ======\n${runs2[*]}" >"$dir"/fqPrepare.log
                fi
                fq1=${dir}/${sample}_S1_L001_R1_001.fastq.gz
                fq2=${dir}/${sample}_S1_L001_R2_001.fastq.gz
            elif [[ $mode == "atac" ]]; then
                if (($(ls "${dir}"/run*_"${sample}"_S1_L001_R1_001.fastq.gz | wc -l) == 1)); then
                    cp -fa "${dir}"/run1_"${sample}"_S1_L001_R1_001.fastq.gz "${dir}"/"${sample}"_S1_L001_R1_001.fastq.gz
                    cp -fa "${dir}"/run1_"${sample}"_S1_L001_R2_001.fastq.gz "${dir}"/"${sample}"_S1_L001_R2_001.fastq.gz
                    cp -fa "${dir}"/run1_"${sample}"_S1_L001_R3_001.fastq.gz "${dir}"/"${sample}"_S1_L001_R3_001.fastq.gz
                    cp -fa "${dir}"/run1_"${sample}"_S1_L001_I1_001.fastq.gz "${dir}"/"${sample}"_S1_L001_I1_001.fastq.gz
                    echo -e "Fastq files for ${sample} is ready." >"$dir"/fqPrepare.log
                    echo -e "====== ${sample}_S1_L001_R1_001.fastq.gz ======\n${dir}/run1_${sample}_S1_L001_R1_001.fastq.gz\n====== ${sample}_S1_L001_R2_001.fastq.gz ======\n${dir}/run1_${sample}_S1_L001_R2_001.fastq.gz" >>"$dir"/fqPrepare.log
                    echo -e "====== ${sample}_S1_L001_R3_001.fastq.gz ======\n${dir}/run1_${sample}_S1_L001_R3_001.fastq.gz\n====== ${sample}_S1_L001_I1_001.fastq.gz ======\n${dir}/run1_${sample}_S1_L001_I1_001.fastq.gz" >>"$dir"/fqPrepare.log
                else
                    runs1=$(ls -lL "${dir}"/run*_"${sample}"_S1_L001_R1_001.fastq.gz)
                    runs2=$(ls -lL "${dir}"/run*_"${sample}"_S1_L001_R2_001.fastq.gz)
                    runs3=$(ls -lL "${dir}"/run*_"${sample}"_S1_L001_R3_001.fastq.gz)
                    runs4=$(ls -lL "${dir}"/run*_"${sample}"_S1_L001_I1_001.fastq.gz)
                    if [[ ! -f "${dir}"/"${sample}_S1_L001_R1_001.fastq.gz" ]] || [[ ! $(echo "${runs1[*]}" | awk 'BEGIN{sum=0}{sum+=$5}END{print sum}') == $(ls -lL "${dir}"/"${sample}"_S1_L001_R1_001.fastq.gz | awk '{print$5}') ]]; then
                        rm -f "$dir"/fqPrepare.log
                        echo "${runs1[*]}" | awk '{print$9}' | xargs cat >"${dir}"/"${sample}"_S1_L001_R1_001.fastq.gz
                    fi
                    if [[ ! -f "${dir}"/"${sample}_S1_L001_R2_001.fastq.gz" ]] || [[ ! $(echo "${runs2[*]}" | awk 'BEGIN{sum=0}{sum+=$5}END{print sum}') == $(ls -lL "${dir}"/"${sample}"_S1_L001_R2_001.fastq.gz | awk '{print$5}') ]]; then
                        rm -f "$dir"/fqPrepare.log
                        echo "${runs2[*]}" | awk '{print$9}' | xargs cat >"${dir}"/"${sample}"_S1_L001_R2_001.fastq.gz
                    fi
                    if [[ ! -f "${dir}"/"${sample}_S1_L001_R3_001.fastq.gz" ]] || [[ ! $(echo "${runs3[*]}" | awk 'BEGIN{sum=0}{sum+=$5}END{print sum}') == $(ls -lL "${dir}"/"${sample}"_S1_L001_R3_001.fastq.gz | awk '{print$5}') ]]; then
                        rm -f "$dir"/fqPrepare.log
                        echo "${runs3[*]}" | awk '{print$9}' | xargs cat >"${dir}"/"${sample}"_S1_L001_R3_001.fastq.gz
                    fi
                    if [[ ! -f "${dir}"/"${sample}_S1_L001_I1_001.fastq.gz" ]] || [[ ! $(echo "${runs4[*]}" | awk 'BEGIN{sum=0}{sum+=$5}END{print sum}') == $(ls -lL "${dir}"/"${sample}"_S1_L001_I1_001.fastq.gz | awk '{print$5}') ]]; then
                        rm -f "$dir"/fqPrepare.log
                        echo "${runs4[*]}" | awk '{print$9}' | xargs cat >"${dir}"/"${sample}"_S1_L001_R3_001.fastq.gz
                    fi
                    echo -e "Fastq files for ${sample} is ready." >"$dir"/fqPrepare.log
                    echo -e "====== ${sample}_S1_L001_R1_001.fastq.gz ======\n${runs1[*]}\n====== ${sample}_S1_L001_R2_001.fastq.gz ======\n${runs2[*]}" >>"$dir"/fqPrepare.log
                    echo -e "====== ${sample}_S1_L001_R3_001.fastq.gz ======\n${runs3[*]}\n====== ${sample}_S1_L001_I1_001.fastq.gz ======\n${runs4[*]}" >>"$dir"/fqPrepare.log
                fi
                fq1="${dir}"/"${sample}_S1_L001_R1_001.fastq.gz"
                fq2="${dir}"/"${sample}_S1_L001_R3_001.fastq.gz"
            fi

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
                rm -rf "$dir"/PreAlignmentQC/fastqc
                mkdir -p "$dir"/PreAlignmentQC/fastqc
                fastqc -o "$dir"/PreAlignmentQC/fastqc -t "$threads" "${fq1}" "${fq2}" &>"$dir"/PreAlignmentQC/fastqc/fastqc.log

                check_logfile "$sample" "FastQC" "$dir"/PreAlignmentQC/fastqc/fastqc.log "$error_pattern" "$complete_pattern" "postcheck" $?
                if [[ $? == 1 ]]; then
                    force="TRUE"
                    continue
                fi
            fi

            #FastQ_Screen
            check_logfile "$sample" "FastQ_Screen" "$dir"/PreAlignmentQC/fastq_screen/fastq_screen.log "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
                rm -rf "$dir"/PreAlignmentQC/fastq_screen
                mkdir -p "$dir"/PreAlignmentQC/fastq_screen
                fastq_screen --force --Aligner bowtie2 "$FastqScreen_mode" --conf "$FastqScreen_config" --threads "$threads" "${fq1}" "${fq2}" \
                    --outdir "$dir"/PreAlignmentQC/fastq_screen 2>"$dir"/PreAlignmentQC/fastq_screen/fastq_screen.log

                check_logfile "$sample" "FastQ_Screen" "$dir"/PreAlignmentQC/fastq_screen/fastq_screen.log "$error_pattern" "$complete_pattern" "postcheck" $?
                if [[ $? == 1 ]]; then
                    force="TRUE"
                    continue
                fi
            fi

            if [[ $mode == "rna" ]]; then
                # cellranger
                check_logfile "$sample" "cellranger" "$dir"/Alignment-Cellranger/cellranger.log "$error_pattern" "$complete_pattern" "precheck"
                if [[ $? == 1 ]]; then
                    rm -rf "$dir"/Alignment-Cellranger
                    mkdir -p "$dir"/Alignment-Cellranger
                    cd "$dir"/Alignment-Cellranger

                    cellranger count \
                        --id "${sample}" \
                        --sample "${sample}" \
                        --fastqs "${dir}" \
                        --transcriptome "$cellranger_ref" \
                        --disable-ui \
                        --localcores "$threads" \
                        --localmem "$memory" \
                        ${cellranger_param} &>"$dir"/Alignment-Cellranger/cellranger.log

                    rm -rf "$dir"/Alignment-Cellranger/SC_RNA_COUNTER_CS

                    check_logfile "$sample" "cellranger" "$dir"/Alignment-Cellranger/cellranger.log "$error_pattern" "$complete_pattern" "postcheck" $?
                    if [[ $? == 1 ]]; then
                        continue
                    fi
                    mv "$dir"/Alignment-Cellranger/"$sample"/outs/web_summary.html "$dir"/Alignment-Cellranger/"$sample"/outs/"$sample".cellranger.html
                    mv "$dir"/Alignment-Cellranger/"$sample"/outs/cloupe.cloupe "$dir"/Alignment-Cellranger/"$sample"/outs/"$sample".cloupe
                fi

                # dropEst
                # check_logfile "$sample" "dropEst" "$dir"/Alignment-Cellranger/"$sample"/dropEst/dropEst.log "$error_pattern" "$complete_pattern" "precheck"
                # if [[ $? == 1 ]]; then
                #     rm -rf "$dir"/Alignment-Cellranger/"$sample"/dropEst
                #     mkdir -p "$dir"/Alignment-Cellranger/"$sample"/dropEst
                #     cd "$dir"/Alignment-Cellranger/"$sample"/dropEst
                #     dropest -V -f -g "$gene_gtf" -c "$dropEst_config" "$dir"/Alignment-Cellranger/"$sample"/outs/possorted_genome_bam.bam &>"$dir"/Alignment-Cellranger/"$sample"/dropEst/dropEst.log

                #     check_logfile "$sample" "dropEst" "$dir"/Alignment-Cellranger/"$sample"/dropEst/dropEst.log "$error_pattern" "$complete_pattern" "postcheck" $?
                #     if [[ $? == 1 ]]; then
                #         continue
                #     fi
                # fi

                # check_logfile "$sample" "dropReport" "$dir"/Alignment-Cellranger/"$sample"/dropEst/dropReport.log "$error_pattern" "$complete_pattern" "precheck"
                # if [[ $? == 1 ]]; then
                #     mkdir -p "$dir"/Alignment-Cellranger/"$sample"/dropEst
                #     cd "$dir"/Alignment-Cellranger/"$sample"/dropEst
                #     dropReport.Rsc "$dir"/Alignment-Cellranger/"$sample"/dropEst/cell.counts.rds &>"$dir"/Alignment-Cellranger/"$sample"/dropEst/dropReport.log

                #     check_logfile "$sample" "dropReport" "$dir"/Alignment-Cellranger/"$sample"/dropEst/dropReport.log "$error_pattern" "$complete_pattern" "postcheck" $?
                #     if [[ $? == 1 ]]; then
                #         continue
                #     fi
                # fi

                # Cell-calling
                # check_logfile "$sample" "CellCalling" "$dir"/Alignment-Cellranger/"$sample"/CellCalling/CellCalling.log "$error_pattern" "$complete_pattern" "precheck"
                # if [[ $? == 1 ]]; then
                #     if [[ -f "$dir"/Alignment-Cellranger/"$sample"/velocyto/velocyto.log ]]; then
                #         rm -f "$dir"/Alignment-Cellranger/"$sample"/velocyto/velocyto.log
                #     fi
                #     rm -rf "$dir"/Alignment-Cellranger/"$sample"/CellCalling
                #     mkdir -p "$dir"/Alignment-Cellranger/"$sample"/CellCalling
                #     cd "$dir"/Alignment-Cellranger/"$sample"/CellCalling
                #     script=$SCP_path/SCP-GeneralSteps/RunCellranger/CellCalling.R
                #     Rscript "$script" "$dir"/Alignment-Cellranger "$sample" "$threads" "$EmptyThreshold" "$CellLabel" &>"$dir"/Alignment-Cellranger/"$sample"/CellCalling/CellCalling.log

                #     check_logfile "$sample" "CellCalling" "$dir"/Alignment-Cellranger/"$sample"/CellCalling/CellCalling.log "$error_pattern" "$complete_pattern" "postcheck"
                #     if [[ $? == 1 ]]; then
                #         continue
                #     fi
                # fi

                # velocyto
                check_logfile "$sample" "velocyto" "$dir"/Alignment-Cellranger/"$sample"/velocyto/velocyto.log "$error_pattern" "$complete_pattern" "precheck"
                if [[ $? == 1 ]]; then
                    rm -rf "$dir"/Alignment-Cellranger/"$sample"/velocyto
                    mkdir -p "$dir"/Alignment-Cellranger/"$sample"/velocyto
                    cd "$dir"/Alignment-Cellranger/"$sample"/velocyto

                    BAM="$dir"/Alignment-Cellranger/"$sample"/outs/cellsorted_possorted_genome_bam.bam
                    if [ -f ${BAM} ]; then
                        samtools quickcheck -v ${BAM} &>/dev/null
                        if [[ $? != 0 ]]; then
                            color_echo "yellow" "[INFO] $sample: BAM file checked failed. Remove the cellsorted_possorted_genome_bam.bam"
                            rm -f ${BAM}
                            continue
                        fi
                    fi
                    if (($(find "$dir"/Alignment-Cellranger/"$sample"/outs/ -name "*bam.tmp.*.bam" | wc -l) > 1)); then
                        rm -f "$dir"/Alignment-Cellranger/"$sample"/outs/*bam.tmp.*.bam
                    fi

                    # cp "$dir"/Alignment-Cellranger/"$sample"/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ./
                    # gunzip ./barcodes.tsv.gz
                    ### samtools sort -t CB -O BAM -o "$dir"/Alignment-Cellranger/"$sample"/outs/cellsorted_possorted_genome_bam.bam "$dir"/Alignment-Cellranger/"$sample"/outs/possorted_genome_bam.bam
                    #velocyto run -b "$dir"/Alignment-Cellranger/"$sample"/CellCalling/barcodes.tsv -e $sample -o "$dir"/Alignment-Cellranger/"$sample"/velocyto -m "$rmsk_gtf" --samtools-threads "$threads" "$dir"/Alignment-Cellranger/"$sample"/outs/possorted_genome_bam.bam "$gene_gtf" &>"$dir"/Alignment-Cellranger/"$sample"/velocyto/velocyto.log
                    velocyto run10x -m "$rmsk_gtf" --samtools-threads "$threads" "$dir"/Alignment-Cellranger/"$sample" "$gene_gtf" &>"$dir"/Alignment-Cellranger/"$sample"/velocyto/velocyto.log

                    check_logfile "$sample" "velocyto" "$dir"/Alignment-Cellranger/"$sample"/velocyto/velocyto.log "$error_pattern" "$complete_pattern" "postcheck" $?
                    if [[ $? == 1 ]]; then
                        continue
                    else
                        rm -f "$dir"/Alignment-Cellranger/"$sample"/outs/cellsorted_possorted_genome_bam.bam
                    fi
                fi
            elif [[ $mode == "atac" ]]; then
                # cellranger-atac
                check_logfile "$sample" "cellranger-atac" "$dir"/Alignment-Cellranger-atac/cellranger-atac.log "$error_pattern" "$complete_pattern" "precheck"
                if [[ $? == 1 ]]; then
                    rm -rf "$dir"/Alignment-Cellranger-atac
                    mkdir -p "$dir"/Alignment-Cellranger-atac
                    cd "$dir"/Alignment-Cellranger-atac

                    cellranger-atac count \
                        --id "${sample}" \
                        --sample "${sample}" \
                        --fastqs "${dir}" \
                        --reference "$cellranger_ref" \
                        --disable-ui \
                        --localcores "$threads" \
                        --localmem "$memory" \
                        ${cellranger_param} &>"$dir"/Alignment-Cellranger-atac/"${sample}"/cellranger-atac.log

                    rm -rf "$dir"/Alignment-Cellranger-atac/SC_ATAC_COUNTER_CS

                    check_logfile "$sample" "cellranger-atac" "$dir"/Alignment-Cellranger-atac/cellranger-atac.log "$error_pattern" "$complete_pattern" "postcheck" $?
                    if [[ $? == 1 ]]; then
                        continue
                    fi
                    mv "$dir"/Alignment-Cellranger-atac/"$sample"/outs/web_summary.html "$dir"/Alignment-Cellranger-atac/"$sample"/outs/"$sample".cellranger-atac.html
                    mv "$dir"/Alignment-Cellranger-atac/"$sample"/outs/cloupe.cloupe "$dir"/Alignment-Cellranger-atac/"$sample"/outs/"$sample".cloupe
                fi
            elif [[ $mode == "arc" ]]; then
                # cellranger-arc
                check_logfile "$sample" "cellranger-arc" "$dir"/Alignment-Cellranger-arc/cellranger-arc.log "$error_pattern" "$complete_pattern" "precheck"
                if [[ $? == 1 ]]; then
                    rm -rf "$dir"/Alignment-Cellranger-arc
                    mkdir -p "$dir"/Alignment-Cellranger-arc
                    cd "$dir"/Alignment-Cellranger-arc

                    cellranger-arc count \
                        --id "${sample}" \
                        --sample "${sample}" \
                        --libraries=/home/jdoe/runs/libraries.csv \
                        --reference "$cellranger_ref" \
                        --disable-ui \
                        --localcores "$threads" \
                        --localmem "$memory" \
                        ${cellranger_param} &>"$dir"/Alignment-Cellranger-arc/cellranger-arc.log

                    check_logfile "$sample" "cellranger-arc" "$dir"/Alignment-Cellranger-arc/cellranger-arc.log "$error_pattern" "$complete_pattern" "postcheck" $?
                    if [[ $? == 1 ]]; then
                        continue
                    fi
                    mv "$dir"/Alignment-Cellranger-arc/"$sample"/outs/web_summary.html "$dir"/Alignment-Cellranger-arc/"$sample"/outs/"$sample".cellranger-arc.html
                    mv "$dir"/Alignment-Cellranger-arc/"$sample"/outs/cloupe.cloupe "$dir"/Alignment-Cellranger-arc/"$sample"/outs/"$sample".cloupe
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

logfiles=("preparation.log")
globalcheck_logfile "$maindir/NGSmodule_SCP_analysis" logfiles[@] "$force_complete" "$error_pattern" "$complete_pattern" "All"
check_logfile "All" "Preparation" $maindir/NGSmodule_SCP_analysis/Preparation/preparation.log "$error_pattern" "$complete_pattern" "precheck"
if [[ $? == 1 ]]; then
    rm -rf $maindir/NGSmodule_SCP_analysis/Preparation
    mkdir -p $maindir/NGSmodule_SCP_analysis/Preparation
    cd $maindir/NGSmodule_SCP_analysis/Preparation
    script=$SCP_path/SCP-GeneralSteps/RunCellranger/Preparation.R
    Rscript $script $work_dir $mode $genome $gene_gtf &>"$maindir/NGSmodule_SCP_analysis/Preparation/preparation.log"

    check_logfile "All" "Preparation" $maindir/NGSmodule_SCP_analysis/Preparation/preparation.log "$error_pattern" "$complete_pattern" "postcheck" $?
    if [[ $? == 1 ]]; then
        color_echo "red" "ERROR! Preparation failed! Please check the $maindir/NGSmodule_SCP_analysis/Preparation/preparation.log."
    fi
fi

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "****************** RunCellranger Done ******************\n"
