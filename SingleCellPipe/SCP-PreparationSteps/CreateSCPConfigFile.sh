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
total_memory=$(free -g | awk 'NR==2 {print ($7/3)}')                             ## Total memory(gigabytes).
ntask_per_run="ALL"
SampleInfoFile="$(pwd)/temp_Sample_info.csv" ## Absolute path of a .csv SampleInfoFile.
SampleGrepPattern=""                          ## Optional. Perl-compatible regexps used for matching the SampleID under the NGSmodule_work directory.
force_complete="FALSE"
retry=0

############# CreateWorkDir Paramaters ###################################################################
### raw_fastq_file_name= RunIDPattern + SufixPattern
### Example: ESC_20200911NB_S244_L004_R1_001.fastq.gz
### RunIDPattern="ESC_.*"
### SufixPattern="_S\d+_L\d+_R1_001\.fastq\.gz"
RunIDPattern=".*"                                                          ## This pattern must could be matched with the LibraryId after excluding Sufix.    
R1_SufixPattern="_S\d+_L\d+_R1_\d+\.((fastq)|(fq))\.gz"             ## Must start with and end with a fixed string. 
R2_SufixPattern="_S\d+_L\d+_R2_\d+\.((fastq)|(fq))\.gz"             ## Must start with and end with a fixed string. 
R1_to_R2="_R1_/_R2_"

############# Cellranger Paramaters #######################################################################
### FastqScreen ###
FastqScreen_config="/data/reference/FastQ_Screen/FastQ_Screen_Genomes/fastq_screen.conf"

### CellRanger ###
cellranger_ref="/data/reference/CellRanger/refdata-gex-GRCh38-and-mm10-2020-A"
include_introns="TRUE"

############# CellCalling Paramaters #######################################################################
# ### dropEst ### (test)
# dropEst_config="/home/zhanghao/Program/NGS/SingleCell/dropEst/configs/10x_v3.xml"

# ### CellCalling ### (test)
# EmptyThreshold="AUTO"                               # an integer number or "AUTO"
# CellLabel="NULL"                                    # a gene name or "NULL"

### Velocyto ###
gene_gtf="\$cellranger_ref/genes/genes.gtf"
rmsk_gtf="\$cellranger_ref/genes/genes_rmsk.gtf"

############# Analysis Paramaters #######################################################################
### CellQC ###
db_method="scDblFinder" ## Doublet-calling methods used. Can be one of scDblFinder, Scrublet, DoubletDetection, scds_cxds, scds_bcds, scds_hybrid.
outlier_cutoff="log10_nCount:both:2.5,log10_nFeature:both:2.5,featurecount_dist:lower:2.5"
outlier_n=1
gene_threshold=1000 ## 1000. Minimum threshold for the cell gene count.
UMI_threshold=3000 ## 3000. Minimum threshold for the cell UMI count.
mito_threshold=20 ## 20. Maximum threshold for the count proportion of mitochondrial genes.
ribo_threshold=50 ## 50. Maximum threshold for the count proportion of ribosomal genes.
species="" ## Leave blank or comma-separated species names, e.g. "Homo_sapiens,Mus_musculus". The first is the species of interest.
species_gene_prefix="" ## Leave blank or comma-separated prefixes, e.g. "GRCh38-,mm10-". The first is the species of interest.
species_percent=95 ## Count proportion thresholds for species of interest.
exogenous_genes="" ## Leave blank or or comma-separated gene symbol.
features_inspect="nCount_RNA,nFeature_RNA,percent.mito,percent.ribo" ## Comma-separated gene symbol or other features.
samples="" ## Leave blank or comma-separated names of samples to be analyzed. Leave blank means analyze all samples.
simple_analysis="TRUE" ## Whether to do a simple analysis.

### StandardSCP ###
normalization_method="logCPM"  ## "logCPM,SCT", comma-separated.
vars_to_regress=NULL  ## "nCount_RNA,nFeature_RNA,percent.mt,percent.ribo" or NULL, comma-separated.
nHVF=3000  ## 3000. Number of high-variable features to use.
liner_reduction="pca"  ## "pca,ica,nmf,mds", comma-separated. Used to compute neighbors and nonlinear dimensionality reduction. The first method will be used in the Integration_SCP step.
liner_reduction_dims=50  ## 50. Number of dimensions to calculate in nonlinear dimensionality reduction.
liner_reduction_dims_use=NULL  ## A number or NULL (automatic). The top N dimensions used to compute neighbors and nonlinear dimensionality reduction.
liner_reduction_distance="cosine"  ## "cosine" or other method in parallelDist::parDist. The distance measure to be used in the MDS or other liner dimensionality reduction.
nonliner_reduction="umap"  ## "umap,tsne,dm", comma-separated. Methods for nonlinear dimensionality reduction.
nonliner_reduction_dims="2,3"  ## "2,3", comma-separated. Dimensional numbers in descending space.
nonliner_reduction_distance="euclidean"  ## "euclidean" or other method in parallelDist::parDist. The distance measure to be used in the DiffusionMap or other nonliner dimensionality reduction.
cluster_algorithm="louvain"  ## One of "louvain","slm" and "leiden". The algorithm used to identify clusters of cells.
cluster_resolution=0.8  ## The resolution used to identify clusters of cells.

### IntegrationSCP ###
HVF_source="separate"  ## "separate,global", comma-separated. Source of high variable featues.
datasets_integration="all;"  ## "all;LLH_A,LLH_B;", comma-separated. Comma separates datasets in one integration, semicolon separates integrations.
datasets_names="all;"  ## "all;LLH_sets;", semicolon-separated. Name for each integration.
integration_method="Uncorrected,Seurat,fastMNN,Harmony,Scanorama,BBKNN,CSS,LIGER"  ## "Uncorrected,Seurat,fastMNN,Harmony,Scanorama,BBKNN,CSS,LIGER", comma-separated. Methods used for integration.



EOF

echo -e "ConfigFile: $ConfigFile\n"
echo -e "\n****************** CreateSCPConfigFile Finished ******************\n"







