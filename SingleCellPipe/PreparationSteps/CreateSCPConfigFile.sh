#!/usr/bin/env bash

if [[ $1 != "" ]];then
    ConfigFile=$1
    if [[ -f $ConfigFile ]];then
        tmp=(`date +"%Y%m%d%H%M%S"`)
        mv ${ConfigFile} bk_${tmp}_${ConfigFile}
    fi
else
    tmp=(`date +"%Y%m%d%H%M%S"`)
    ConfigFile=temp_${tmp}.config
fi

cat <<- EOF >$ConfigFile
############# Global Paramaters ###########################################################################
maindir="$(pwd)"        ## Absolute path of your project directory.
rawdata_dir="$(pwd)/rawdata/"
total_threads=$(grep 'processor' /proc/cpuinfo | sort -u | wc -l)                             ## Total threads.
total_memory=$(free -g | awk 'NR==2 {print $7}')                             ## Total memory.
ntask_per_run="ALL"
SampleInfoFile="$(pwd)/temp_Sample_info.csv" ## Absolute path of a .csv SampleInfoFile.
SampleGrepPattern=""                          ## Optional. Perl-compatible regexps used for matching the SampleID under the NGSmodule_work directory.
force_complete="FALSE"

############# CreateWorkDir Paramaters ###################################################################
### raw_fastq_file_name= LibraryIdPattern  + SufixPattern
### Example: 1_P31_20200911NB_S244_L004_R2_001.fastq.gz
### LibraryIdPattern=".*"
### SufixPattern="_2020.*\.fastq\.gz"
LibraryIdPattern=".*"                      ## This pattern must could be matched with the RunId in the SampleInfoFile after excluding Sufix.
SufixPattern="_S\d+_L\d+.*\.((fastq)|(fq))\.gz"         ## Must start with and end with a fixed string pattern. 


############# Cellranger Paramaters #######################################################################
### cellranger ###
cellranger_ref="/data/database/CellRanger/custom/GFP_and_refdata-gex-GRCh38-2020-A"

### velocyto and dropEst ###
dropEst_config="/home/zhanghao/Program/NGS/SingleCell/dropEst/configs/10x_v3.xml"
gene_gtf="\$cellranger_ref/genes/genes.gtf"
rmsk_gtf="\$cellranger_ref/genes/rmsk.gtf"



EOF

echo -e "ConfigFile: $ConfigFile\n"
echo -e "\n****************** CreateSCPConfigFile Finished ******************\n"
