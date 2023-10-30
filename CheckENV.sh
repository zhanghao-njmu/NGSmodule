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
env_packages["NGSmodule-preAlignmentQC"]="-c conda-forge -c bioconda fastqc fastp fastq-screen bbmap sortmerna pigz"
env_packages["NGSmodule-Alignment"]="-c conda-forge -c bioconda bwa bowtie hisat2 star bismark kallisto==0.48.0 samtools sambamba"
env_packages["NGSmodule-postAlignmentQC"]="-c conda-forge -c bioconda rseqc preseq goleft mosdepth"
env_packages["NGSmodule-Analysis"]="-c conda-forge -c bioconda r-base"

env_packages["NGSmodule-SCP"]="-c conda-forge -c bioconda fastqc fastq-screen velocyto.py r-base r-seurat r-signac r-future r-hdf5r r-seuratdisk r-velocyto.r bioconductor-rtracklayer bioconductor-rsamtools bioconductor-genomicranges bioconductor-rhdf5 bioconductor-hdf5array"

env_name=$1

# 获取系统环境中的Rscript路径  
system_rscript=$(whereis Rscript)  

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

# 激活conda环境  
conda activate $env_name

# 获取conda环境中的Rscript路径  
conda_rscript=$(whereis Rscript)  
  
# 将路径字符串转换为数组  
IFS=' ' read -r -a system_rscript_arr <<< "$system_rscript"  
IFS=' ' read -r -a conda_rscript_arr <<< "$conda_rscript"  
  
# 遍历conda环境中的Rscript路径数组  
for path in "${conda_rscript_arr[@]}"; do  
    # 检查路径是否在系统环境的Rscript路径数组中出现  
    if [[ ! "${system_rscript_arr[@]}" =~ "${path}" ]]; then  
        # 如果没有出现，那么这是conda环境中的Rscript路径  
        alias Rscript="$path"
    fi  
done