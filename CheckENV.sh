#!/usr/bin/env bash

#######################################################################################
trap_add 'conda deactivate' SIGINT SIGTERM EXIT

conda &>/dev/null
[ $? -eq 127 ] && {
  color_echo "red" "Cannot find the conda. Please install conda and add conda to your PATH environment variable.\n"
  exit 1
}
eval "$(conda shell.bash hook)"

declare -A env_packages
env_packages["NGSmodule-PrefetchData"]="-c conda-forge -c bioconda -c dranew samtools pyfaidx gem2 genmap ucsc-wigtobigwig hmmcopy_utils picard awscli bowtie bowtie2 hisat2 star bismark sra-tools"
env_packages["NGSmodule-SCP"]="-c conda-forge -c bioconda pigz fastqc fastq-screen"
env_packages["NGSmodule-preAlignmentQC"]="-c conda-forge -c bioconda fastqc fastp fastq-screen bbmap sortmerna pigz"
env_packages["NGSmodule-Alignment"]="-c conda-forge -c bioconda bwa bowtie hisat2 star bismark samtools sambamba"
env_packages["NGSmodule-postAlignmentQC"]="-c conda-forge -c bioconda rseqc preseq goleft mosdepth"

env_name=$1

if [[ ! $(conda env list | grep $env_name) ]]; then
  color_echo "green" ">>> Create the environment: $env_name\n"
  conda create -y -q --name $env_name ${env_packages[$env_name]}
  if [[ $? != 0 ]];then
    color_echo "red" ">>> Failed to create the environment: $env_name\n"
    exit 1
  fi
else
  color_echo "green" ">>> Found the environment: $env_name\n"
fi

conda activate $env_name
