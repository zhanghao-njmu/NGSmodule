#!/usr/bin/env bash

echo -e "****************** Start PrepareWorkDir ******************\n"
#######################################################################################
if [[ -d $work_dir ]]; then
    tmp=($(date +"%Y%m%d%H%M%S"))
    mv "${work_dir}" "${work_dir}"/../bk_"${tmp}"_work
    mkdir "${work_dir}"
fi

grep_pattern="(${LibraryIdPattern}${SufixPattern})"
color_echo "green" "Grep pattern: $grep_pattern \n"
arr=($(find $rawdata_dir -type f | grep -P $grep_pattern))
if [[ ${#arr} == 0 ]]; then
    color_echo "red" "Error! Cannot find any file matched the pattern!\nPlease check the LibraryIdPattern and SufixPattern in the ConfigFile!\n"
    exit 1
fi

for file in "${arr[@]}"; do
    file_sim=${file##*/}
    Sufix=($(echo "${file_sim}" | grep -oP "$SufixPattern"))

    if [[ "${#Sample_dict[@]}" != 0 ]] && [[ "${#Layout_dict[@]}" != 0 ]]; then
        use_run="FALSE"
        if [[ $Sufix ]] && [[ ${Sample_dict[${file_sim%%$Sufix}]} ]]; then
            LibraryId=${file_sim%%$Sufix}
            SampleID=${Sample_dict[$LibraryId]}
            Layout=${Layout_dict[$SampleID]}
            if [[ $SampleID == "" ]]; then
                SampleID=$LibraryId
                color_echo "yellow" "Warning! Cannot find the SampleID for LibraryId: $LibraryId. Use '$LibraryId' as its SampleID."
            fi
            use_run="TRUE"
        fi
    else
        color_echo "red" "Error! Cannot find the SampleID or Layout information. Please check the SampleInfoFile."
        exit 1
    fi

    if [[ $use_run == "FALSE" ]]; then
        continue
    else

        if [[ $Sufix == $R1_Sufix ]]; then
            fq=run1_${SampleID}_1.fq.gz
            fq_Layout="PE"
        elif [[ $Sufix == $R2_Sufix ]]; then
            fq=run1_${SampleID}_2.fq.gz
            fq_Layout="PE"
        elif [[ $Sufix == $SE_Sufix ]]; then
            fq=run1_${SampleID}.fq.gz
            fq_Layout="SE"
        fi

        if [[ $Layout != $fq_Layout ]]; then
            color_echo "red" "Error caused by the file:$file\nThe layout of the fastq file:$fq_Layout is conflict with the Layout information of the SampleInfoFile:$Layout"
            exit 1
        fi

        echo "File: ${file_sim}  RunId: ${RunId}  SampleID: ${SampleID}"
        mkdir -p ${work_dir}/${SampleID}
        if [[ ! -f ${work_dir}/$SampleID/${fq} ]]; then
            ln -s $file ${work_dir}/$SampleID/${fq}
        else
            num=$(ls ${work_dir}/$SampleID/run*_${SampleID}${fq##run1_${SampleID}} | wc -l)
            ln -s $file ${work_dir}/$SampleID/run$(($num + 1))_${fq##run1_}
        fi

    fi

done

echo -e "\n****************** CreateWorkDir Finished ******************\n"
