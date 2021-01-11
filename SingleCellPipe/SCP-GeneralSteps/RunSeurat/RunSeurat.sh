#!/usr/bin/env bash
trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

pigz --version &>/dev/null
[ $? -ne 0 ] && {
    color_echo "red" "Cannot find the command pigz.\n"
    exit 1
}
velocyto &>/dev/null
[ $? -eq 127 ] && {
    echo -e "Cannot find the command velocyto.\n"
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

echo -e "########################### RunSeurat Parameters ##############################\n"
echo -e "*** base parameters ***"
echo -e "    Rscript_path: $(which Rscript)\n     Rscript_threads: ${Rscript_threads}\n    species: ${species}\n    exogenous_genes: ${exogenous_genes}"
echo -e "*** cell-filtering parameters ***"
echo -e "    cell_calling_methodNum: ${cell_calling_methodNum}\n    mito_threshold: ${mito_threshold}\n    gene_threshold: ${gene_threshold}\n    UMI_threshold: ${UMI_threshold}"
echo -e "*** Seurat parameters ***"
echo -e "    normalization_method: ${normalization_method}\n    nHVF: ${nHVF}\n    maxPC: ${maxPC}\n    resolution: ${resolution}\n    reduction: ${reduction}"
echo -e "################################################################################\n"

echo -e "****************** Start RunSeurat ******************\n"
SECONDS=0

for sample in "${arr[@]}"; do
    read -u1000
    {
        dir=$work_dir/$sample/Alignment/Cellranger/$sample/Seurat
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
            logfiles=("RunSeuratStatus.log" )
            globalcheck_logfile "$dir" logfiles[@] "$force" "$error_pattern" "$complete_pattern" "$sample"

            fq1=${dir}/${sample}_S1_L001_R1_001.fastq.gz
            fq2=${dir}/${sample}_S1_L001_R2_001.fastq.gz

            ##To verify that reads appear to be correctly paired
            check_logfile "$sample" "RunSeurat" "$dir"/RunSeurat "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
                fqCheck_PE "$sample" "$fq1" "$fq2" "$dir"/fqCheck.log
                check_logfile "$sample" "RunSeurat" "$dir"/fqCheck.log "$error_pattern" "$complete_pattern" "postcheck" $?
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
            check_logfile "$sample" "cellranger" "$dir"/Alignment/Cellranger/cellranger.log "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
                rm -rf "$dir"/Alignment/Cellranger
                mkdir -p "$dir"/Alignment/Cellranger
                cd "$dir"/Alignment/Cellranger
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
                    --transcriptome "$cellranger_ref" &>"$dir"/Alignment/Cellranger/cellranger.log

                check_logfile "$sample" "cellranger" "$dir"/Alignment/Cellranger/cellranger.log "$error_pattern" "$complete_pattern" "postcheck"
                if [[ $? == 1 ]]; then
                    continue
                fi
            fi

            # velocyto
            # check_logfile "$sample" "velocyto" "$dir"/Alignment/Cellranger/"$sample"/velocyto/velocyto.log "$error_pattern" "$complete_pattern" "precheck"
            # if [[ $? == 1 ]]; then
            #     rm -rf "$dir"/Alignment/Cellranger/"$sample"/velocyto
            #     mkdir -p "$dir"/Alignment/Cellranger/"$sample"/velocyto
            #     cd "$dir"/Alignment/Cellranger/"$sample"/velocyto
            #     #velocyto run10x -m "$rmsk_gtf" --samtools-threads "$threads" "$dir"/Alignment/Cellranger/"$sample" "$gene_gtf" &>"$dir"/Alignment/Cellranger/"$sample"/velocyto/velocyto.log
            #     cp -f "$dir"/Alignment/Cellranger/"$sample"/outs/raw_feature_bc_matrix/barcodes.tsv.gz "$dir"/Alignment/Cellranger/"$sample"/velocyto/barcodes.tsv.gz
            #     gunzip -f "$dir"/Alignment/Cellranger/"$sample"/velocyto/barcodes.tsv.gz
            #     velocyto run -b "$dir"/Alignment/Cellranger/"$sample"/velocyto/barcodes.tsv -o "$dir"/Alignment/Cellranger/"$sample"/velocyto -m "$rmsk_gtf" --samtools-threads "$threads" "$dir"/Alignment/Cellranger/"$sample"/outs/possorted_genome_bam.bam "$gene_gtf" &>"$dir"/Alignment/Cellranger/"$sample"/velocyto/velocyto.log
            #     rm -f "$dir"/Alignment/Cellranger/"$sample"/velocyto/barcodes.tsv

            #     check_logfile "$sample" "velocyto" "$dir"/Alignment/Cellranger/"$sample"/velocyto/velocyto.log "$error_pattern" "$complete_pattern" "postcheck"
            #     if [[ $? == 1 ]]; then
            #         continue
            #     fi
            # fi

            # dropEst
            check_logfile "$sample" "dropEst" "$dir"/Alignment/Cellranger/"$sample"/dropEst/dropEst.log "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
                rm -rf "$dir"/Alignment/Cellranger/"$sample"/dropEst
                mkdir -p "$dir"/Alignment/Cellranger/"$sample"/dropEst
                cd "$dir"/Alignment/Cellranger/"$sample"/dropEst
                dropest -V -f -g "$gene_gtf" -c "$dropEst_config" "$dir"/Alignment/Cellranger/"$sample"/outs/possorted_genome_bam.bam &>"$dir"/Alignment/Cellranger/"$sample"/dropEst/dropEst.log

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
                Rscript "$1" "$dir"/Alignment/Cellranger "$sample" "$threads" "$EmptyThreshold" "$CellLabel" &>"$dir"/Alignment/Cellranger/"$sample"/CellCalling/CellCalling.log

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
