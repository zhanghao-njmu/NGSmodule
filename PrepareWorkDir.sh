#!/usr/bin/env bash

echo -e "****************** Start PrepareWorkDir ******************\n"
#######################################################################################
if [[ -d $work_dir ]];then
  tmp=(`date +"%Y%m%d%H%M%S"`)
  mv ${work_dir} ${work_dir}/../bk_${tmp}_work
  mkdir ${work_dir}
fi

arr=(` find $rawdata_dir -type f  | grep -P "(${RunIdPattern}${SE_SufixPattern})|(${RunIdPattern}${R1_SufixPattern})|(${RunIdPattern}${R2_SufixPattern})" `)
if [[ ${#arr} == 0 ]];then
  echo -e "Error! Cannot find the rawdata!\nPlease check the RunIdPattern and SufixPattern in the ConfigFile!\n"
  exit 1
fi

for file in ${arr[@]};do
  file_sim=${file##*/}
  SE_Sufix=(`echo "${file_sim}" | grep -oP "$SE_SufixPattern"`)
  R1_Sufix=(`echo "${file_sim}" | grep -oP "$R1_SufixPattern"`)
  R2_Sufix=(`echo "${file_sim}" | grep -oP "$R2_SufixPattern"`)
  
  if [[ "${#Sample_dict[@]}" != 0 ]] && [[ "${#Layout_dict[@]}" != 0 ]];then
    for map_Sufix in $SE_Sufix $R1_Sufix $R2_Sufix;do
      if [[ $map_Sufix ]] && [[ ${Sample_dict[${file_sim%%$map_Sufix}]} ]];then
        Sufix=$map_Sufix
        RunId=${file_sim%%$map_Sufix}
        SampleID=${Sample_dict[$RunId]}
        Layout=${Layout_dict[$SampleID]}
        if [[ $SampleID == "" ]];then
          SampleID=$RunId
          echo -e "Warning! Cannot find the SampleID for RunId: $RunId. Use '$RunId' as its SampleID."
        fi
        continue;
      fi
    done
  else
    echo -e "Error! Cannot find the SampleID or Layout information. Please check the SampleInfoFile."
    exit 1
  fi
  
  if [[ $Sufix == $SE_Sufix ]];then
    fq=run1_${SampleID}.fq.gz
    fq_Layout="SE"
  elif [[ $Sufix == $R1_Sufix ]];then
    fq=run1_${SampleID}_1.fq.gz
    fq_Layout="PE"
  elif [[ $Sufix == $R2_Sufix ]];then
    fq=run1_${SampleID}_2.fq.gz
    fq_Layout="PE"
  fi

  if [[ $Layout != $fq_Layout ]];then
    echo -e "Error caused by the file:$file\nThe layout of the fastq file:$fq_Layout is conflict with the Layout information of the SampleInfoFile:$Layout"
    exit 1
  fi

  echo "File: ${file_sim}  RunId: ${RunId}  SampleID: ${SampleID}"
  mkdir -p ${work_dir}/${SampleID}
  if [[ ! -f ${work_dir}/$SampleID/${fq} ]];then
    ln -s $file ${work_dir}/$SampleID/${fq}
  else
    num=$(ls ${work_dir}/$SampleID/run*_${SampleID}${fq##run1_${SampleID}} |wc -l)
    ln -s $file ${work_dir}/$SampleID/run$(($num+1))_${fq##run1_}
  fi

done
echo -e ""
echo -e "\n****************** PrepareWorkDir Done ******************\n"
