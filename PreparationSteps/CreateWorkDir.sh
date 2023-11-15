#!/usr/bin/env bash

echo -e "****************** Start PrepareWorkDir ******************\n"
#######################################################################################
if [[ -d $work_dir ]]; then
  tmp=($(date +"%Y%m%d%H%M%S"))
  mv "${work_dir}" "${work_dir}"/../bk_"${tmp}"_work
  mkdir "${work_dir}"
fi

grep_pattern="(${RunIDPattern}${SE_SufixPattern}$)|(${RunIDPattern}${R1_SufixPattern}$)|(${RunIDPattern}${R2_SufixPattern}$)"

color_echo "green" ">>> Please make sure your file SufixPattern matched with the following pattern:\n"
color_echo "green" "    SE_SufixPattern=$SE_SufixPattern"
color_echo "green" "    R1_SufixPattern=$R1_SufixPattern"
color_echo "green" "    R2_SufixPattern=$R2_SufixPattern\n"

arr=($(find $rawdata_dir -type f | grep -P $grep_pattern | sort))
if [[ ${#arr} == 0 ]]; then
  color_echo "red" "Error! Cannot find any file matched the pattern!\nPlease check the RunIDPattern and SufixPattern in the ConfigFile!\n"
  exit 1
fi

for file in "${arr[@]}"; do
  file_sim=${file##*/}
  SE_Sufix=($(echo "${file_sim}" | grep -oP "$SE_SufixPattern"))
  R1_Sufix=($(echo "${file_sim}" | grep -oP "$R1_SufixPattern"))
  R2_Sufix=($(echo "${file_sim}" | grep -oP "$R2_SufixPattern"))

  if [[ "${#Sample_dict[@]}" != 0 ]] && [[ "${#Layout_dict[@]}" != 0 ]]; then
    use_run="FALSE"
    for map_Sufix in $SE_Sufix $R1_Sufix $R2_Sufix; do
      id=$(echo "${file_sim%%$map_Sufix}" | grep -oP "$RunIDPattern")
      if [[ $map_Sufix ]] && [[ ${Sample_dict[${id}]} ]]; then
        RunID=$id 
        SampleID=${Sample_dict[$RunID]}
        Layout=${Layout_dict[$SampleID]}
        if [[ $SampleID == "" ]]; then
          SampleID=$RunID
          color_echo "yellow" "Warning! Cannot find the SampleID for RunID: $RunID. Use '$RunID' as its SampleID."
        fi
        use_run="TRUE"
      fi
    done
  else
    color_echo "red" "Error! Cannot find the RunID-SampleID information(${#Sample_dict[@]}) or Layout information(${#Layout_dict[@]}) from the SampleInfoFile."
    exit 1
  fi

  if [[ $use_run == "FALSE" ]]; then
    color_echo "yellow" "Warning! SampleInfoFile have no RunID-SampleID matching information for the file: $file "
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

    echo "File: ${file_sim}  RunID: ${RunID}  SampleID: ${SampleID}"
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
