#!/usr/bin/env bash

echo -e "****************** Start PrepareWorkDir ******************\n"
#######################################################################################
if [[ -d $work_dir ]];then
  tmp=(`date +"%Y%m%d%H%M%S"`)
  mv ${work_dir} ${work_dir}/../bk_${tmp}_work
  mkdir ${work_dir}
fi

arr=(` find $rawdata_dir -type f  | grep -P "${SampleIdPattern}${SampleSufixPattern}" `)
if [[ ${#arr} == 0 ]];then
  echo -e "Cannot find the rawdata.\nPlease check the SampleIdPattern and SampleSufixPattern in the ConfigFile!\n"
fi

for file in ${arr[@]};do
  file_sim=${file##*/}
  SampleSufix=(`echo "${file_sim}" | grep -oP "$SampleSufixPattern"`)
  sampleid=${file_sim%%$SampleSufix}

  if [[ "${#Sample_dict[@]}" != 0 ]];then
    samplename=${Sample_dict[$sampleid]}
    layout=${Layout_dict[$samplename]}    
    if [[ $samplename == "" ]];then
      samplename=$sampleid
      echo -e "Warning! Cannot find the SampleName for SampleID: $sampleid. Please check the SampleInfoFile."
    fi
  else
    echo -e "Error! Cannot find the layout info. Please check the SampleInfoFile."
    exit 1
  fi
  
  fqU=$(echo "${file_sim}"  |grep -P "(.fastq.gz)|(.fq.gz)" | grep -Pv "(_R\d.fastq.gz)|(_R\d.fq.gz)|(_trim.fq.gz)")
  fq1=$(echo "${file_sim}"  |grep -P "(_1.fastq.gz)|(_R1.fastq.gz)|(_1.fq.gz)|(_R1.fq.gz)" | grep -Pv "_trim.fq.gz")
  fq2=$(echo "${file_sim}"  |grep -P "(_2.fastq.gz)|(_R2.fastq.gz)|(_2.fq.gz)|(_R2.fq.gz)" | grep -Pv "_trim.fq.gz")
  
  if [[ $fqU ]]  && [[ "${layout}" == "SE" ]];then
    fq=${samplename}.fq.gz
  elif [[ $fq1 ]] && [[ "${layout}" == "PE" ]];then
    fq=${samplename}_1.fq.gz
  elif [[ $fq2 ]] && [[ "${layout}" == "PE" ]];;then
    fq=${samplename}_2.fq.gz
  else
    echo -e "Error! Cannot get the layout information of the file: ${file_sim}"
  fi

  echo "File: ${file_sim}  SampleID: ${sampleid}  SampleName: ${samplename}"
  mkdir -p ${work_dir}/${samplename}
  ln -s $file ${work_dir}/$samplename/${fq}

done
echo -e ""
echo -e "\n****************** PrepareWorkDir Done ******************\n"
