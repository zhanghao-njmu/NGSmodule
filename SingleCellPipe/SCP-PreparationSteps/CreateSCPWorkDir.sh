#!/usr/bin/env bash

echo -e "****************** Start PrepareWorkDir ******************\n"
#######################################################################################
if [[ -d $work_dir ]]; then
    tmp=($(date +"%Y%m%d%H%M%S"))
    mv "${work_dir}" "${work_dir}"/../bk_"${tmp}"_work
    mkdir "${work_dir}"
fi

R2_SufixPattern=$(echo ${R1_SufixPattern} | perl -pe "s/$R1_to_R2/g")
R3_SufixPattern=$(echo ${R1_SufixPattern} | perl -pe "s/$R1_to_R3/g")
I1_SufixPattern=$(echo ${R1_SufixPattern} | perl -pe "s/$R1_to_I1/g")
reads_pattern="(/${RunIDPattern}${R1_SufixPattern}$)|(/${RunIDPattern}${R2_SufixPattern}$)|(/${RunIDPattern}${R3_SufixPattern}$)|(/${RunIDPattern}${I1_SufixPattern}$)"
color_echo "green" ">>> grep_pattern=$reads_pattern"
reads_RunID=()
while IFS='' read -r line; do reads_RunID+=("$line"); done < <(find "$rawdata_dir" -type f | grep -P "$reads_pattern" | sort | sed "s/.*\///g" | perl -pe "s/(${R1_SufixPattern})|(${R2_SufixPattern})|(${R3_SufixPattern})|(${I1_SufixPattern})//g" | sort | uniq)

if [[ ${#reads_RunID} == 0 ]]; then
    color_echo "red" "Error! Cannot find any file matched the pattern!\nPlease check the RunIDPattern and SufixPattern in the ConfigFile!\n"
    exit 1
fi

for RunID in "${reads_RunID[@]}"; do
    R1_pattern="/${RunID}${R1_SufixPattern}$"
    R1_arr=()
    while IFS='' read -r line; do R1_arr+=("$line"); done < <(find "$rawdata_dir" -type f | grep -P "$R1_pattern" | sort)

    R2_pattern="/${RunID}${R2_SufixPattern}$"
    R2_arr=()
    while IFS='' read -r line; do R2_arr+=("$line"); done < <(find "$rawdata_dir" -type f | grep -P "$R2_pattern" | sort)

    R3_pattern="/${RunID}${R2_SufixPattern}$"
    R3_arr=()
    while IFS='' read -r line; do R3_arr+=("$line"); done < <(find "$rawdata_dir" -type f | grep -P "$R3_pattern" | sort)

    I1_pattern="/${RunID}${R2_SufixPattern}$"
    I1_arr=()
    while IFS='' read -r line; do I1_arr+=("$line"); done < <(find "$rawdata_dir" -type f | grep -P "$I1_pattern" | sort)

    if ((${#R1_arr[*]} == 0)); then
        color_echo "red" "Error! RunID: $RunID has no R1 fastq file!"
        exit 1
    fi
    if ((${#R2_arr[*]} == 0)); then
        color_echo "red" "Error! RunID: $RunID has no R2 fastq file!"
        exit 1
    fi

    if ((${#R1_arr[*]} > 1)) || ((${#R2_arr[*]} > 1)); then
        color_echo "yellow" "Warning! RunID: $RunID has more than one R1(${#R1_arr[*]}) or R2(${#R2_arr[*]}) fastq file."
    fi
    if ((${#R1_arr[*]} != ${#R2_arr[*]})); then
        color_echo "red" "Error! RunID: $RunID has diffent number of R1,R2 fastq file!"
        exit 1
    fi

    if [[ $mode == "atac" ]]; then
        if ((${#R3_arr[*]} == 0)); then
            color_echo "red" "Error! RunID: $RunID has no R3 fastq file!"
            exit 1
        fi
        if ((${#I1_arr[*]} == 0)); then
            color_echo "red" "Error! RunID: $RunID has no I1 fastq file!"
            exit 1
        fi
        if ((${#R1_arr[*]} != ${#R3_arr[*]})); then
            color_echo "red" "Error! RunID: $RunID has diffent number of R1,R3 fastq file!"
            exit 1
        fi
        if ((${#R1_arr[*]} != ${#I1_arr[*]})); then
            color_echo "red" "Error! RunID: $RunID has diffent number of R1,I1 fastq file!"
            exit 1
        fi
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
        fi
    else
        color_echo "red" "Error! Cannot find the RunID-SampleID matching from the SampleInfoFile."
        exit 1
    fi

    if [[ $use_run == "FALSE" ]]; then
        color_echo "yellow" "Warning! SampleInfoFile has no RunID-SampleID matching information for the RunID: $RunID "
        continue
    else

        color_echo "green" "RunID: ${RunID}  SampleID: ${SampleID}"
        mkdir -p "${work_dir}"/"${SampleID}"
        for run in "${R1_arr[@]}"; do
            R1_raw=$run
            R2_raw=$(echo ${run} | perl -pe "s/$R1_to_R2/g")
            R3_raw=$(echo ${run} | perl -pe "s/$R1_to_R3/g")
            I1_raw=$(echo ${run} | perl -pe "s/$R1_to_I1/g")

            if [[ $mode == "rna" ]]; then
                echo -e "R1:${R1_raw##*/}  R2:${R2_raw##*/}"
                if [[ ! -f $R2_raw ]]; then
                    color_echo "red" "Error! Cannot find the R1 corresponding R2 file: $R2_raw"
                    exit 1
                fi
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

                    if (($num1 == $num2)); then
                        ln -s "$R1_raw" "${work_dir}"/"$SampleID"/"${R1_new}"
                        ln -s "$R2_raw" "${work_dir}"/"$SampleID"/"${R2_new}"
                    else
                        color_echo "red" "Error! RunID:$RunID has different R1/R2 file in the dir: ${work_dir}/$SampleID "
                        exit 1
                    fi
                fi
            elif [[ $mode == "atac" ]]; then
                echo -e "R1:${R1_raw##*/}  R2:${R2_raw##*/}  R3:${R3_raw##*/}  I1:${I1_raw##*/}"
                if [[ ! -f $R2_raw ]]; then
                    color_echo "red" "Error! Cannot find the R1 corresponding R2 file: $R2_raw"
                    exit 1
                fi
                if [[ ! -f $R2_raw ]]; then
                    color_echo "red" "Error! Cannot find the R1 corresponding R2 file: $R2_raw"
                    exit 1
                fi
                if [[ ! -f $R3_raw ]]; then
                    color_echo "red" "Error! Cannot find the R1 corresponding R3 file: $R3_raw"
                    exit 1
                fi
                if [[ ! -f $I1_raw ]]; then
                    color_echo "red" "Error! Cannot find the I1 corresponding I1 file: $I1_raw"
                    exit 1
                fi

                R1_new=run1_${SampleID}_S1_L001_R1_001.fastq.gz
                R2_new=run1_${SampleID}_S1_L001_R2_001.fastq.gz
                R3_new=run1_${SampleID}_S1_L001_R3_001.fastq.gz
                I1_new=run1_${SampleID}_S1_L001_I1_001.fastq.gz

                if [[ ! -f ${work_dir}/$SampleID/${R1_new} ]] && [[ ! -f ${work_dir}/$SampleID/${R2_new} ]] && [[ ! -f ${work_dir}/$SampleID/${R3_new} ]] && [[ ! -f ${work_dir}/$SampleID/${I1_new} ]]; then
                    ln -s "$R1_raw" "${work_dir}"/"$SampleID"/"${R1_new}"
                    ln -s "$R2_raw" "${work_dir}"/"$SampleID"/"${R2_new}"
                    ln -s "$R3_raw" "${work_dir}"/"$SampleID"/"${R3_new}"
                    ln -s "$I1_raw" "${work_dir}"/"$SampleID"/"${I1_new}"
                else
                    num1=$(ls ${work_dir}/$SampleID/run*_${SampleID}_S1_L001_R1_001.fastq.gz | wc -l)
                    R1_new=run$(($num1 + 1))_${SampleID}_S1_L001_R1_001.fastq.gz
                    num2=$(ls ${work_dir}/$SampleID/run*_${SampleID}_S1_L001_R2_001.fastq.gz | wc -l)
                    R2_new=run$(($num2 + 1))_${SampleID}_S1_L001_R2_001.fastq.gz
                    num3=$(ls ${work_dir}/$SampleID/run*_${SampleID}_S1_L001_R3_001.fastq.gz | wc -l)
                    R3_new=run$(($num1 + 1))_${SampleID}_S1_L001_R3_001.fastq.gz
                    num4=$(ls ${work_dir}/$SampleID/run*_${SampleID}_S1_L001_I1_001.fastq.gz | wc -l)
                    I1_new=run$(($num2 + 1))_${SampleID}_S1_L001_I1_001.fastq.gz

                    if (($num1 == $num2 && $num1 == $num3 && $num1 == $num4)); then
                        ln -s "$R1_raw" "${work_dir}"/"$SampleID"/"${R1_new}"
                        ln -s "$R2_raw" "${work_dir}"/"$SampleID"/"${R2_new}"
                        ln -s "$R3_raw" "${work_dir}"/"$SampleID"/"${R3_new}"
                        ln -s "$I1_raw" "${work_dir}"/"$SampleID"/"${I1_new}"
                    else
                        color_echo "red" "Error! RunID:$RunID has different R1/R2/R3/I1 file in the dir: ${work_dir}/$SampleID "
                        exit 1
                    fi
                fi
            fi


        done

    fi
done

echo -e "\n****************** CreateWorkDir Finished ******************\n"
