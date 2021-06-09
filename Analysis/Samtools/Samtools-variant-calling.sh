#!/usr/bin/env bash

#######################################################################################
trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

samtools &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command samtools.\n"
  exit 1
}

bcftools &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command bcftools.\n"
  exit 1
}

tabix &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command tabix.\n"
  exit 1
}

echo -e "############################# Samtools Parameters #############################\n"
echo -e "  Genome: ${genome} \n"
echo -e "################################################################################\n"

echo -e "****************** Samtools-variant-calling ******************\n"
SECONDS=0

for sample in "${arr[@]}"; do
    read -u1000
    {
        dir=${work_dir}/${sample}
        if [ ! -d "$dir/Alignment-$Aligner" ]; then
            color_echo "red" "Cannot find alignment result. Please run 'Alignment' first.\n"
            exit 1
        fi

        case ${Deduplication} in
        FALSE)
            BAM="$dir/Alignment-$Aligner/${sample}.${Aligner}.markdup.bam"
            ;;
        TRUE)
            BAM="$dir/Alignment-$Aligner/${sample}.${Aligner}.dedup.bam"
            ;;
        *)
            BAM="$dir/Alignment-$Aligner/${sample}.${Aligner}.bam"
            ;;
        esac
        prefix=${BAM%%.bam}
        prefix=${prefix##*/}

        dir_result="$dir/Alignment-$Aligner/Samtools-variant-calling"
        mkdir -p $dir_result
        cd $dir_result

        force=${force_complete}
        status="uncompleted"
        attempt=0

        echo "+++++ ${sample} +++++"

        while [[ $status == "uncompleted" ]] && (("$attempt" <= 0)); do
            ((attempt++))
            if [[ $attempt != 1 ]]; then
                echo -e "+++++ ${sample}: Number of attempts: $attempt +++++"
            fi

            logfiles=("Samtools.log")
            globalcheck_logfile "$dir_result" logfiles[@] "$force" "$error_pattern" "$complete_pattern" "$sample"

            ##### Samtools #####
            check_logfile "$sample" "Samtools" "$dir_result/Samtools/Samtools.log" "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
                rm -rf $dir_result/Samtools
                mkdir -p $dir_result/Samtools
                cd $dir_result/Samtools
                bcftools mpileup --threads $threads -Ou -f $genome $BAM | bcftools call -vmO z -o ${prefix}.Samtools.vcf.gz &>$dir_result/Samtools/Samtools.log
                tabix -p vcf ${prefix}.Samtools.vcf.gz &>$dir_result/Samtools/Samtools.log
                bcftools stats -F $genome -s - ${prefix}.Samtools.vcf.gz > ${prefix}.Samtools.vcf.gz.stats 2>$dir_result/Samtools/Samtools.log

                check_logfile "$sample" "Samtools" "$dir_result/Samtools/Samtools.log" "$error_pattern" "$complete_pattern" "postcheck" $?
                if [[ $? == 1 ]]; then
                    continue
                fi
            fi

            ##### VariantFiltration #####
            check_logfile "$sample" "VariantFiltration" "$dir_result/VariantFiltration/VariantFiltration.log" "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
                rm -rf $dir_result/VariantFiltration
                mkdir -p $dir_result/VariantFiltration
                cd $dir_result/VariantFiltration
                bcftools view --threads $threads --types snps --output-type z --output-file ${prefix}.Samtools.snps.vcf.gz $dir_result/Samtools/${prefix}.Samtools.vcf.gz &>>$dir_result/VariantFiltration/VariantFiltration.log
                bcftools view --threads $threads --types indels  --output-type z --output-file ${prefix}.Samtools.indels.vcf.gz $dir_result/Samtools/${prefix}.Samtools.vcf.gz &>>$dir_result/VariantFiltration/VariantFiltration.log
                bcftools filter --threads $threads -e 'MQ < 40.0' --output-type z --output ${prefix}.Samtools.snps.filter.vcf.gz ${prefix}.Samtools.snps.vcf.gz &>>$dir_result/VariantFiltration/VariantFiltration.log
                bcftools filter --threads $threads -e 'MQ < 20.0' --output-type z --output ${prefix}.Samtools.indels.filter.vcf.gz ${prefix}.Samtools.indels.vcf.gz &>>$dir_result/VariantFiltration/VariantFiltration.log
                bcftools concat --threads $threads --output-type z --output ${prefix}.Samtools.filter.vcf.gz ${prefix}.Samtools.snps.filter.vcf.gz ${prefix}.Samtools.indels.filter.vcf.gz &>>$dir_result/VariantFiltration/VariantFiltration.log
                
                check_logfile "$sample" "VariantFiltration" "$dir_result/VariantFiltration/VariantFiltration.log" "$error_pattern" "$complete_pattern" "postcheck" $?
                if [[ $? == 1 ]]; then
                    continue
                fi
            fi

            status="completed"
            color_echo "blue" "+++++ ${sample}: Samtools-variant-calling completed +++++"
        done

        if [[ "$status" == "completed" ]]; then
            echo "Completed: $sample" >>"$tmpfile"
        else
            echo "Interrupted: $sample" >>"$tmpfile"
            color_echo "red" "ERROR! ${sample} interrupted! Please check the processing log and your raw bam file."
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
    cat "$tmpfile" | grep "Interrupted" | uniq >$maindir/Samtools-variant-calling.Interrupted.txt
    color_echo "red" "\n\n################################################################################"
    color_echo "red" "    $ninterrupted of $total_task tasks interrupted."
    color_echo "red" "    Please check the samples in $maindir/Samtools-variant-calling.Interrupted.txt"
    color_echo "red" "################################################################################\n\n"
fi

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "****************** Samtools-variant-calling Done ******************\n"
