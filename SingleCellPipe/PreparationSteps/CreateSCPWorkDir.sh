#!/usr/bin/env bash

echo -e "****************** Start PrepareWorkDir ******************\n"
#######################################################################################
if [[ -d $work_dir ]]; then
    tmp=($(date +"%Y%m%d%H%M%S"))
    mv "${work_dir}" "${work_dir}"/../bk_"${tmp}"_work
    mkdir "${work_dir}"
fi

grep_pattern="(^${LibraryIdPattern}${SufixPattern}$)"
color_echo "green" "Grep pattern: $grep_pattern \n"
arr=($(find $rawdata_dir -type f | grep -P $grep_pattern))
if [[ ${#arr} == 0 ]]; then
    color_echo "red" "Error! Cannot find any file matched the pattern!\nPlease check the LibraryIdPattern and SufixPattern in the ConfigFile!\n"
    exit 1
fi

for file in "${arr[@]}"; do
    file_sim=${file##*/}
    Sufix=($(echo "${file_sim}" | grep -oP "$SufixPattern"))

    if [[ "${#Sample_dict[@]}" != 0 ]]; then
        use_run="FALSE"
        if [[ $Sufix ]] && [[ ${Sample_dict[${file_sim%%$Sufix}]} ]]; then
            LibraryId=${file_sim%%$Sufix}
            SampleID=${Sample_dict[$LibraryId]}
            if [[ $SampleID == "" ]]; then
                SampleID=$LibraryId
                color_echo "yellow" "Warning! Cannot find the SampleID for LibraryId: $LibraryId. Use '$LibraryId' as its SampleID."
            fi
            use_run="TRUE"
        fi
    else
        color_echo "red" "Error! Cannot find the SampleID information!"
        exit 1
    fi

    if [[ $use_run == "FALSE" ]]; then
        color_echo "yellow" "Warning! Cannot find the LibraryId information for the file: $file "
        continue
    else
        SampleID=$(echo $SampleID | xargs)
        echo "File: ${file_sim}  LibraryId: ${LibraryId}  SampleID: ${SampleID}"
        mkdir -p ${work_dir}/${SampleID}
        ln -fs $file ${work_dir}/$SampleID/${SampleID}${Sufix}
    fi

done

echo -e "\n****************** CreateWorkDir Finished ******************\n"
