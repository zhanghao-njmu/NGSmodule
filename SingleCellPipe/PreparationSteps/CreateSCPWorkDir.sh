#!/usr/bin/env bash

echo -e "****************** Start PrepareWorkDir ******************\n"
#######################################################################################
if [[ -d $work_dir ]]; then
    tmp=($(date +"%Y%m%d%H%M%S"))
    mv "${work_dir}" "${work_dir}"/../bk_"${tmp}"_work
    mkdir "${work_dir}"
fi

PE_pattern="(^${RunIDPattern}${R1_SufixPattern}$)|(^${RunIDPattern}${R2_SufixPattern}$)"
PE_RunID=($(find $rawdata_dir -type f | grep -P $PE_pattern | sort | sed "s/.*\///g" | perl -pe "s/(${R1_SufixPattern})|(${R2_SufixPattern})//g" | sort | uniq))

if [[ ${#PE_RunID} == 0 ]]; then
    color_echo "red" "Error! Cannot find any file matched the pattern!\nPlease check the RunIDPattern and SufixPattern in the ConfigFile!\n"
    echo $PE_RunID
    echo ${#PE_RunID}
    echo $PE_pattern
    echo $((find $rawdata_dir -type f | grep -P $PE_pattern | sort | sed "s/.*\///g"))
    exit 1
fi

for RunID in "${PE_RunID[@]}"; do
    R1_pattern="^${RunID}${R1_SufixPattern}$"
    R1_arr=($(find $rawdata_dir -type f | grep -P $R1_pattern | sort))
    R2_pattern="^${RunID}${R2_SufixPattern}$"
    R2_arr=($(find $rawdata_dir -type f | grep -P $R2_pattern | sort))

    if ((${#R1_arr} == 0)); then
        color_echo "red" "Error! RunID: $RunID have no R1 fastq file!"
        exit 1
    fi
    if ((${#R1_arr} > 1)); then
        color_echo "red" "Error! RunID: $RunID have more than one R1 fastq file: ${R1_arr[*]}"
        exit 1
    fi
    if ((${#R2_arr} == 0)); then
        color_echo "red" "Error! RunID: $RunID have no R2 fastq file!"
        exit 1
    fi
    if ((${#R2_arr} > 1)); then
        color_echo "red" "Error! RunID: $RunID have more than one R2 fastq file: ${R2_arr[*]}"
        exit 1
    fi

    if [[ "${#Sample_dict[@]}" != 0 ]]; then
        use_run="FALSE"
        if [[ ${Sample_dict[$RunID]} ]]; then
            SampleID=${Sample_dict[$RunID]}
            if [[ $SampleID == "" ]]; then
                SampleID=$RunID
                color_echo "yellow" "Warning! Cannot find the SampleID for RunID: $RunID. Use '$RunID' as its SampleID."
            fi
            use_run="TRUE"
            break
        fi
    else
        color_echo "red" "Error! Cannot find the RunID-SampleID matching from the SampleInfoFile."
        exit 1
    fi

    if [[ $use_run == "FALSE" ]]; then
        color_echo "yellow" "Warning! SampleInfoFile have no RunID-SampleID matching information for the RunID: $RunID "
        continue
    else

        echo "RunID: ${RunID}  SampleID: ${SampleID}"
        mkdir -p "${work_dir}"/"${SampleID}"

        R1_raw=${R1_arr[1]}
        R2_raw=${R2_arr[1]}
        R1_new=run1_${SampleID}_S1_L001_R1_001.fastq.gz
        R2_new=run1_${SampleID}_S1_L001_R2_001.fastq.gz

        if [[ ! -f ${work_dir}/$SampleID/${R1_new} ]] && [[ ! -f ${work_dir}/$SampleID/${R2_new} ]]; then
            ln -s "$R1_raw" "${work_dir}"/"$SampleID"/"${R1_new}"
            ln -s "$R2_raw" "${work_dir}"/"$SampleID"/"${R2_new}"
        else
            num1=$(ls ${work_dir}/$SampleID/run*_${SampleID}_S1_L001_R1_001.fastq.gz | wc -l)
            R1_new=run$(($num1 + 1))_${SampleID}_S1_L001_R1_001.fastq.gz
            num2=$(ls ${work_dir}/$SampleID/run*_${SampleID}_S1_L001_R2_001.fastq.gz | wc -l)
            R2_new=run$(($num2 + 1))_${SampleID}_S1_L001_R2_001.fastq.gz

            if (($num1==$num2));then
                ln -s "$R1_raw" "${work_dir}"/"$SampleID"/"${R1_new}"
                ln -s "$R2_raw" "${work_dir}"/"$SampleID"/"${R2_new}"
            else
                color_echo "red" "Error! RunID:$RunID have different R1/R2 file in the dir: ${work_dir}/$SampleID "
                exit 1
            fi
        fi

    fi
done

echo -e "\n****************** CreateWorkDir Finished ******************\n"
