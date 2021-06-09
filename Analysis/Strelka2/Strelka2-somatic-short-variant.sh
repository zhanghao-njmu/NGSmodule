#!/usr/bin/env bash

#######################################################################################
trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

configureStrelkaSomaticWorkflow.py --version &>/dev/null
[ $? -ne 0 ] && {
    color_echo "red" "Cannot find the command 'configureStrelkaSomaticWorkflow.py' from Strelka2.\n"
    exit 1
}

echo -e "############################# Strelka2 Parameters #############################\n"
echo -e "  Strelka2: $(which configureStrelkaSomaticWorkflow.py)\n"
echo -e "  Genome: ${genome} \n"
echo -e "  Exome: ${exome} \n"
echo -e "  Targeted: ${targeted} \n"
echo -e "################################################################################\n"

if [[ $exome == "TRUE" ]]; then
    par_exome="--exome"
else
    par_exome=""
fi

if [[ $targeted == "TRUE" ]]; then
    par_targeted="--targeted"
else
    par_targeted=""
fi

echo -e "****************** Start Strelka2-somatic-short-variant ******************\n"
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

        dir_result="$dir/Alignment-$Aligner/Strelka2-somatic-short-variant"
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

            logfiles=("Strelka2.log")
            globalcheck_logfile "$dir_result" logfiles[@] "$force" "$error_pattern" "$complete_pattern" "$sample"

            ### Strelka2 #####
            check_logfile "$sample" "Strelka2" "$dir_result/Strelka2/Strelka2.log" "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
                rm -rf $dir_result/Strelka2
                mkdir -p $dir_result/Strelka2
                cd $dir_result/Strelka2
                configureStrelkaSomaticWorkflow.py \
                    --tumorBam $BAM \
                    --referenceFasta $genome \
                    $par_exome $par_targeted \
                    --runDir $dir_result/Strelka2 &>>$dir_result/Strelka2/Strelka2.log
                $dir_result/Strelka2/runWorkflow.py -m local -j $threads &>>$dir_result/Strelka2/Strelka2.log
                status=$?
                cp $dir_result/Strelka2/results/variants/variants.vcf.gz $dir_result/Strelka2/${prefix}.Strelka2.vcf.gz
                cp $dir_result/Strelka2/results/variants/variants.vcf.gz.tbi $dir_result/Strelka2/${prefix}.Strelka2.vcf.gz.tbi

                check_logfile "$sample" "Strelka2" "$dir_result/Strelka2/Strelka2.log" "$error_pattern" "$complete_pattern" "postcheck" $status
                if [[ $? == 1 ]]; then
                    continue
                fi
            fi

            status="completed"
            color_echo "blue" "+++++ ${sample}: Strelka2-somatic-short-variant completed +++++"
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
    cat "$tmpfile" | grep "Interrupted" | uniq >$maindir/Strelka2-somatic-short-variant.Interrupted.txt
    color_echo "red" "\n\n################################################################################"
    color_echo "red" "    $ninterrupted of $total_task tasks interrupted."
    color_echo "red" "    Please check the samples in $maindir/Strelka2-somatic-short-variant.Interrupted.txt"
    color_echo "red" "################################################################################\n\n"
fi

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "****************** Strelka2-somatic-short-variant Done ******************\n"
