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
total_memory=$(free -g | awk 'NR==2 {print $7}')                             ## Total memory(gigabytes).
ntask_per_run="ALL"
SampleInfoFile="$(pwd)/temp_Sample_info.csv" ## Absolute path of a .csv SampleInfoFile.
SampleGrepPattern=""                          ## Optional. Perl-compatible regexps used for matching the SampleID under the NGSmodule_work directory.
force_complete="FALSE"

############# CreateWorkDir Paramaters ###################################################################
### raw_fastq_file_name= LibraryIdPattern  + SufixPattern
### Example: ESC_20200911NB_S244_L004_R1_001.fastq.gz
### LibraryIdPattern=".*"
### SufixPattern="_2020.*\.fastq\.gz"
RunIDPattern=".*"                                                          ## This pattern must could be matched with the LibraryId after excluding Sufix.    
R1_SufixPattern="_2020.*_S\d+_L\d+_R1_\d+\.((fastq)|(fq))\.gz"             ## Must start with and end with a fixed string. 
R2_SufixPattern="_2020.*_S\d+_L\d+_R2_\d+\.((fastq)|(fq))\.gz"             ## Must start with and end with a fixed string. 
R1_to_R2="_R1_/_R2_"

############# Cellranger Paramaters #######################################################################
### FastqScreen ###
FastqScreen_config="/archive/reference/FastQ_Screen/FastQ_Screen_Genomes/fastq_screen.conf"

### cellranger ###
cellranger_ref="/archive/reference/CellRanger/refdata-gex-GRCh38-2020-A_addGFP"

### velocyto and dropEst ###
dropEst_config="/home/zhanghao/Program/NGS/SingleCell/dropEst/configs/10x_v3.xml"
gene_gtf="\$cellranger_ref/genes/genes.gtf"
rmsk_gtf="\$cellranger_ref/genes/hg38_rmsk.gtf"

### Cell-Calling ###
EmptyThreshold="AUTO"                               # an integer number or 'AUTO'
CellLabel="NULL"                                    # a gene name or 'NULL'

############# Intergration #######################################################################
### base parameters ###
Rscript_threads=120
datasets="ESC,iMeLC,PGC;Testis1,Testis2;"
species="Homo_sapiens"                              # Homo_sapiens,Mus_musculus,Macaca_fascicularis,Macaca_mulatta,Drosophila_melanogaster 
exogenous_genes="NULL"                              # a gene name, e.g. 'GFP' or 'NULL'

### cell-filtering ###
cell_calling_methodNum=3
mito_threshold=0.2
gene_threshold=1000
UMI_threshold=3000

### basic ###
normalization_method="logCPM,SCT"       # logCPM,SCT
nHVF=3000
maxPC=100
resolution=0.8
reduction="umap"                        # umap,tsne

### integration ###
HVF_source="separate"                  # global,separate
integration_method="Uncorrected,Seurat,fastMNN,Harmony,Scanorama,BBKNN,CSS,LIGER,scMerge" # Uncorrected,Seurat,fastMNN,Harmony,Scanorama,BBKNN,CSS,LIGER,scMerge


EOF

echo -e "ConfigFile: $ConfigFile\n"
echo -e "\n****************** CreateSCPConfigFile Finished ******************\n"
