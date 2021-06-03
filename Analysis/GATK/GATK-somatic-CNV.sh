#!/usr/bin/env bash

#######################################################################################
trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

$GATK3 --help &>/dev/null
[ $? -ne 0 ] && {
    color_echo "red" "Cannot find the tool GATK3.\n"
    exit 1
}

if [[ ${Species} == "Homo_sapiens" ]]; then
    case ${Build} in
    hg19)
        Hapmap_snps="$iGenomes_Dir/$Species/GATK/hg19/hapmap_3.3.hg19.sites.vcf.gz"
        Omni_snps="$iGenomes_Dir/$Species/GATK/hg19/1000G_omni2.5.hg19.sites.vcf.gz"
        KG_indels="$iGenomes_Dir/$Species/GATK/hg19/1000G_phase1.indels.hg19.sites.vcf.gz"
        KG_snps="$iGenomes_Dir/$Species/GATK/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz"
        Mills_indels="$iGenomes_Dir/$Species/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
        dbSNP_snps="$iGenomes_Dir/$Species/GATK/hg19/dbsnp_138.hg19.vcf.gz"

        known_indels+=($KG_indels $Mills_indels)
        known_snps+=($dbSNP_snps)
        ;;
    hg38)
        Hapmap_snps="$iGenomes_Dir/$Species/GATK/hg38/hapmap_3.3.hg38.vcf.gz"
        Omni_snps="$iGenomes_Dir/$Species/GATK/hg38/1000G_omni2.5.hg38.vcf.gz"
        KG_snps="$iGenomes_Dir/$Species/GATK/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
        Mills_indels="$iGenomes_Dir/$Species/GATK/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
        dbSNP_snps="$iGenomes_Dir/$Species/GATK/hg38/dbsnp_146.hg38.vcf.gz"

        known_indels+=($Mills_indels)
        known_snps+=($dbSNP_snps)
        ;;
    GRCh37)
        Hapmap_snps="$iGenomes_Dir/$Species/GATK/GRCh37/Annotation/GATKBundle/hapmap_3.3.b37.vcf.gz"
        Omni_snps="$iGenomes_Dir/$Species/GATK/GRCh37/Annotation/GATKBundle/1000G_omni2.5.b37.vcf.gz"
        KG_snps="$iGenomes_Dir/$Species/GATK/GRCh37/Annotation/GATKBundle/1000G_phase1.snps.high_confidence.b37.vcf.gz"
        KG_indels="$iGenomes_Dir/$Species/GATK/GRCh37/Annotation/GATKBundle/1000G_phase1.indels.b37.vcf.gz"
        Mills_indels="$iGenomes_Dir/$Species/GATK/GRCh37/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
        dbSNP_snps="$iGenomes_Dir/$Species/GATK/GRCh37/Annotation/GATKBundle/dbsnp_138.b37.vcf.gz"

        known_indels+=($KG_indels $Mills_indels)
        known_snps+=($dbSNP_snps)
        ;;
    GRCh38)
        Hapmap_snps="$iGenomes_Dir/$Species/GATK/GRCh38/Annotation/GATKBundle/hapmap_3.3.hg38.vcf.gz"
        Omni_snps="$iGenomes_Dir/$Species/GATK/GRCh38/Annotation/GATKBundle/1000G_omni2.5.hg38.vcf.gz"
        KG_snps="$iGenomes_Dir/$Species/GATK/GRCh38/Annotation/GATKBundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
        Mills_indels="$iGenomes_Dir/$Species/GATK/GRCh38/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
        dbSNP_snps="$iGenomes_Dir/$Species/GATK/GRCh38/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz"

        known_indels+=($Mills_indels)
        known_snps+=($dbSNP_snps)
        ;;
    *)
        color_echo "yellow" "No support vcf for the ${Build}.\n"
        ;;
    esac

    resource_snps+=("hapmap,known=false,training=true,truth=true,prior=15.0 $Hapmap_snps")
    resource_snps+=("omini,known=false,training=true,truth=false,prior=12.0 $Omni_snps")
    resource_snps+=("1000G,known=false,training=true,truth=false,prior=10.0 $KG_snps")
    resource_snps+=("dbsnp,known=true,training=false,truth=false,prior=2.0 $dbSNP_snps")
    resource_indels+=("mills,known=true,training=true,truth=true,prior=12.0 $Mills_indels")
    resource_indels+=("dbsnp,known=true,training=false,truth=false,prior=2.0 $dbSNP_snps")
fi
if [[ ${Species} == "Mus_musculus" ]] && [[ ${Source} == "Ensembl" ]] && [[ ${Build} == "GRCm38" ]]; then
    dbSNP_indels="$iGenomes_Dir/$Species/Ensembl/GRCm38/MouseGenomeProject/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz"
    dbSNP_snps="$iGenomes_Dir/$Species/Ensembl/GRCm38/MouseGenomeProject/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"

    known_indels+=($dbSNP_indels)
    known_snps+=($dbSNP_snps)
    resource_snps+=("dbsnp,known=true,training=false,truth=false,prior=2.0 $dbSNP_snps")
    resource_indels+=("dbsnp,known=true,training=false,truth=false,prior=2.0 $dbSNP_snps")
fi

echo -e "############################# GATK Parameters #############################\n"
echo -e "  GATK3: \"${GATK3}\"\n"
echo -e "  Genome: ${genome} \n"
echo -e "  known_indels: ${known_indels[*]} \n"
echo -e "  known_snps: ${known_snps[*]} \n"
echo -e "  resource_snps: ${resource_snps[*]} \n"
echo -e "  resource_indels: ${resource_indels[*]} \n"
echo -e "################################################################################\n"

echo -e "****************** Start GATK-germline-short-variant ******************\n"
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

        dir_result="$dir/Alignment-$Aligner/GATK-somatic-short-variant"
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

            logfiles=("Realigner.log" "BQSR.log" "HaplotypeCaller.log")
            globalcheck_logfile "$dir_result" logfiles[@] "$force" "$error_pattern" "$complete_pattern" "$sample"

            ##### Realigner #####
            # check_logfile "$sample" "Realigner" "$dir_result/Realigner/Realigner.log" "$error_pattern" "$complete_pattern" "precheck"
            # if [[ $? == 1 ]]; then
            #     rm -rf $dir_result/Realigner
            #     mkdir -p $dir_result/Realigner
            #     cd $dir_result/Realigner
            #     par_known_indels=$(printf -- " -known '%s'" "${known_indels[@]}")
            #     eval "$GATK3 -T RealignerTargetCreator -nt $threads -R $genome -I $BAM $par_known_indels -o ${prefix}.IndelRealigner.intervals" &>>$dir_result/Realigner/Realigner.log
            #     eval "$GATK3 -T IndelRealigner -R $genome -I $BAM $par_known_indels -o ${prefix}.realign.bam --targetIntervals ${prefix}.IndelRealigner.intervals" &>>$dir_result/Realigner/Realigner.log

            #     check_logfile "$sample" "Realigner" "$dir_result/Realigner/Realigner.log" "$error_pattern" "$complete_pattern" "postcheck"
            #     if [[ $? == 1 ]]; then
            #         continue
            #     fi
            # fi

            ##### BQSR #####
            check_logfile "$sample" "BQSR" "$dir_result/BQSR/BQSR.log" "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
                rm -rf $dir_result/BQSR
                mkdir -p $dir_result/BQSR
                cd $dir_result/BQSR
                par_known_indels=$(printf -- " --knownSites '%s'" "${known_indels[@]}")
                par_known_snps=$(printf -- " --knownSites '%s'" "${known_snps[@]}")
                eval "$GATK3 -T BaseRecalibrator -nct $threads -R $genome -I $BAM $par_known_indels $par_known_snps -o ${prefix}.BQSR.table" &>>$dir_result/BQSR/BQSR.log
                eval "$GATK3 -T PrintReads -nct $threads -R $genome -I $BAM -BQSR ${prefix}.BQSR.table -o ${prefix}.BQSR.bam" &>>$dir_result/BQSR/BQSR.log

                check_logfile "$sample" "BQSR" "$dir_result/BQSR/BQSR.log" "$error_pattern" "$complete_pattern" "postcheck"
                if [[ $? == 1 ]]; then
                    continue
                fi
            fi

            ##### HaplotypeCaller #####
            check_logfile "$sample" "HaplotypeCaller" "$dir_result/HaplotypeCaller/HaplotypeCaller.log" "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
                rm -rf $dir_result/HaplotypeCaller
                mkdir -p $dir_result/HaplotypeCaller
                cd $dir_result/HaplotypeCaller
                eval "$GATK3 -T HaplotypeCaller --emitRefConfidence GVCF -nct $threads -R $genome -I ${dir_result}/BQSR/${prefix}.BQSR.bam -variant_index_type LINEAR -variant_index_parameter 128000 -o ${prefix}.gvcf.gz" &>>$dir_result/HaplotypeCaller/HaplotypeCaller.log
                eval "$GATK3 -T GenotypeGVCFs -nct $threads -R $genome --variant ${prefix}.gvcf.gz -o ${prefix}.vcf.gz" &>>$dir_result/HaplotypeCaller/HaplotypeCaller.log

                check_logfile "$sample" "HaplotypeCaller" "$dir_result/HaplotypeCaller/HaplotypeCaller.log" "$error_pattern" "$complete_pattern" "postcheck"
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
                eval "$GATK3 -T SelectVariants -R $genome -V ${prefix}.Mutect2.vcf.gz -selectType SNP -o ${prefix}.Mutect2.snps.vcf.gz" &>>$dir_result/VariantFiltration/VariantFiltration.log
                eval "$GATK3 -T SelectVariants -R $genome -V ${prefix}.Mutect2.vcf.gz -selectType INDEL -o ${prefix}.Mutect2.indels.vcf.gz" &>>$dir_result/VariantFiltration/VariantFiltration.log
                eval "$GATK3 -T VariantFiltration -R $genome -V ${prefix}.Mutect2.snps.vcf.gz --filterExpression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'" &>>$dir_result/VariantFiltration/VariantFiltration.log
                eval "$GATK3 -T VariantFiltration -R $genome -V ${prefix}.Mutect2.indels.vcf.gz --filterExpression 'QD < 2.0 || ReadPosRankSum < -20.0 || InbreedingCoeff < -0.8 || FS > 200.0 || SOR > 10.0'" &>>$dir_result/VariantFiltration/VariantFiltration.log
                eval "$GATK3 -T CombineVariants -V ${prefix}.Mutect2.snps.vcf.gz -V ${prefix}.Mutect2.indels.vcf.gz -o ${prefix}.Mutect2.VariantFiltration.vcf.gz" &>>$dir_result/VariantFiltration/VariantFiltration.log

                check_logfile "$sample" "VariantFiltration" "$dir_result/VariantFiltration/VariantFiltration.log" "$error_pattern" "$complete_pattern" "postcheck"
                if [[ $? == 1 ]]; then
                    continue
                fi
            fi
            
            ##### VQSR #####
            # check_logfile "$sample" "VQSR" "$dir_result/VQSR/VQSR.log" "$error_pattern" "$complete_pattern" "precheck"
            # if [[ $? == 1 ]]; then
            #     rm -rf $dir_result/VQSR
            #     mkdir -p $dir_result/VQSR
            #     cd $dir_result/VQSR
            #     par_resource_snps=$(printf -- "-resource:%s " "${resource_snps[@]}")
            #     par_resource_indels=$(printf -- "-resource:%s " "${resource_indels[@]}")
            #     eval "$GATK3 -T VariantRecalibrator -nct $threads -R $genome -input $dir_result/Mutect2/${prefix}.Mutect2.vcf.gz $par_resource_snps -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -mode SNP -recalFile ${prefix}.Mutect2.snps.recal -tranchesFile ${prefix}.Mutect2.snps.tranches -rscriptFile ${prefix}.Mutect2.snps.plots.R" &>>$dir_result/VQSR/VQSR.log
            #     eval "$GATK3 -T ApplyRecalibration -nct $threads -R $genome -input $dir_result/Mutect2/${prefix}.Mutect2.vcf.gz --ts_filter_level 99.0 -recalFile ${prefix}.Mutect2.snps.recal -tranchesFile ${prefix}.Mutect2.snps.tranches -mode SNP -o ${prefix}.Mutect2.snps.VQSR.vcf.gz" &>>$dir_result/VQSR/VQSR.log
            #     eval "$GATK3 -T VariantRecalibrator -nct $threads -R $genome -input ${prefix}.Mutect2.snps.VQSR.vcf.gz $par_resource_indels -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -mode INDEL -recalFile ${prefix}.Mutect2.snps.indels.recal -tranchesFile ${prefix}.Mutect2.snps.indels.tranches -rscriptFile ${prefix}.Mutect2.snps.indels.plots.R" &>>$dir_result/VQSR/VQSR.log
            #     eval "$GATK3 -T ApplyRecalibration -nct $threads -R $genome -input ${prefix}.Mutect2.snps.VQSR.vcf.gz --ts_filter_level 99.0 -recalFile ${prefix}.Mutect2.snps.indels.recal -tranchesFile ${prefix}.Mutect2.snps.indels.tranches -mode INDEL -o ${prefix}.Mutect2.snps.indels.VQSR.vcf.gz" &>>$dir_result/VQSR/VQSR.log
            #     rm -f ${prefix}.Mutect2.snps.VQSR.vcf.gz
            #     mv ${prefix}.Mutect2.snps.indels.VQSR.vcf.gz ${prefix}.Mutect2.VQSR.vcf.gz

            #     check_logfile "$sample" "VQSR" "$dir_result/VQSR/VQSR.log" "$error_pattern" "$complete_pattern" "postcheck"
            #     if [[ $? == 1 ]]; then
            #         continue
            #     fi
            # fi

            status="completed"
            color_echo "blue" "+++++ ${sample}: GATK-germline-short-variant completed +++++"
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
    cat "$tmpfile" | grep "Interrupted" | uniq >$maindir/GATK-germline-short-variant.Interrupted.txt
    color_echo "red" "\n\n################################################################################"
    color_echo "red" "    $ninterrupted of $total_task tasks interrupted."
    color_echo "red" "    Please check the samples in $maindir/GATK-germline-short-variant.Interrupted.txt"
    color_echo "red" "################################################################################\n\n"
fi

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "****************** GATK-germline-short-variant Done ******************\n"
