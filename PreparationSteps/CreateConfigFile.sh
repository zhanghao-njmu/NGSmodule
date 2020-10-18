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
rawdata_dir="$(pwd)/rawdata/"   ## Absolute path of directory containing the raw fastq.gz data.
SequenceType="rna"                            ## rna,dna,BSdna
total_threads=$(grep 'processor' /proc/cpuinfo | sort -u | wc -l)                             ## Total threads.
ntask_per_run="ALL"                           ## "ALL" or numeric value to specify the number of tasks run simultaneously at the backend.
SampleInfoFile="$(pwd)/temp_Sample_info.csv" ## Absolute path of a .csv SampleInfoFile.
SampleGrepPattern=""                          ## Optional. Perl-compatible regexps used for matching the SampleID under the NGSmodule_work directory.
force_complete="FALSE"                        ## A global option to determine whether to execute a complete process for any MODE .


############# CreateWorkDir Paramaters ###################################################################
### raw_run_file_name= RunIDPattern  + SufixPattern
### Example: Sample1: R19051060_2020_L001_1.fq.gz, R19051060_2020_L001_2.fq.gz; Sample2: R19051061_2020_L002.fq.gz
### RunIDPattern=".*"
### SE_SufixPattern="_2020.*\.fastq\.gz"; R1_SufixPattern="_2020.*_1\.fastq\.gz"; R2_SufixPattern="_2020.*_2\.fastq\.gz"
RunIDPattern=".*"                      ## This pattern must could be matched with the RunId in the SampleInfoFile after excluding Sufix.
SE_SufixPattern="\.((fastq)|(fq))\.gz"          # Must start with and end with a fixed string. 
R1_SufixPattern="_1\.((fastq)|(fq))\.gz"        # Must start with and end with a fixed string.     
R2_SufixPattern="_2\.((fastq)|(fq))\.gz"        # Must start with and end with a fixed string. 


############# preAlignmentQC Paramaters ###################################################################
### Fastp ###
trim_front1=2                  ## trimming how many bases in front for read1. e.g. 1-4 bp for RNAseq and 9-10 bp for WGBS.
trim_tail1=0                   ## trimming how many bases in tail for read1.
trim_front2=2                  ## trimming how many bases in front for read2. e.g. 1-4 bp for RNAseq and 9-10 bp for WGBS. Only valid when Layout=PE.
trim_tail2=0                   ## trimming how many bases in tail for read2.
qualified_quality_phred=20     ## base quaility over this threshold value will be qualified.
unqualified_percent_limit=50   ## how many percents of bases are allowed to be unqualified. Otherwise the reads will be dropped.
read_cutting="--cut_right"     ## one or more options in the following: --cut_front --cut_tail --cut_right. Multiple selection should be seperated by spaces.
cut_window_size=4              ## the window size option shared by cut_front, cut_tail or cut_right. Range: 1~1000.
cut_mean_quality=20            ## the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. Range: 1~36.
length_required=20             ## reads shorter than length_required will be discarded.

### FastqScreen ###
FastqScreen_config="/data/database/FastQ_Screen/FastQ_Screen_Genomes/fastq_screen.conf"

### SortmeRNA ###
SortmeRNA_Dir="/data/database/SortmeRNA"       ## the path of the dir containing the reference sequence.
SortmeRNA_Type="rRNA"                          ## rRNA,Mt_tRNA,Mt_rRNA
SortmeRNA_Species="Homo_sapiens"               ## Homo_sapiens,Mus_musculus,Macaca_fascicularis,Macaca_mulatta,Drosophila_melanogaster 
SortmeRNA_DataVersion="EnsemblGenes98"         ## the version of the sequence 
SortmeRNA_ref_direct=""                        ## Optional. Specify the path of the reference sequence file for SortmeRNA. 


############# Alignment Paramaters ##########################################################################
iGenomes_Dir="/data/database/iGenomes"           ## The iGenomes dir containing the index under a directory tree: {iGenomes_Dir}/{Species}/{Source}/Sequence/{Aligner}
Species="Homo_sapiens"                           ## Homo_sapiens,Mus_musculus,Macaca_fascicularis,Macaca_mulatta,Drosophila_melanogaster 
Source="Ensembl"                                 ## Ensembl,NCBI,UCSC
Build="GRCh38"                                   ## The genome build version.
Aligner="hisat2"                                 ## bwa,bowtie,bowtie2,hisat2,tophat2,star,bismark_bowtie2,bismark_hisat2
Aligner_parament=""                              ## Optional. Specify custom parameters and overwrite the default patameters excluding the index path and threads number.  
Genome_direct=""                                 ## Optional. Specify a genome file path used for the alignment.  
GTF_direct=""                                    ## Optional. Specify a gtf file path used for the alignment.  
Index_direct=""                                  ## Optional. Specify the index path used for the alignment.  


############# Quantification Paramaters ######################################################################
strandspecific=0                                   ## 0(unstranded),1(stranded),2(reversely stranded)


############# DifferentialExpression Paramaters ##############################################################
max_padj=0.05                                      ## e.g. 0.05 or 0.01 or 0.001
min_fc=2                                           ## e.g. 2
min_count=10                                       ## Minimum count required for at least n samples (n is the smallest group sample size).
group_compare="Hom-80S,WT-80S;Hom-Input,WT-Input;" ## Groups are seperated by comma(,). Different comparisons are seperated by semicolon(;). 
DGEs_multi_compare=1                               ## Whether to compare DGEs among different comparisons. 0(not to do),1(do).


############# CNVanalysis Paramaters ##############################################################
Window=1000000
Kmer=130
PloidyAssumed=2

############# SNV Paramaters ##############################################################
GATK3="java -jar -Xmx32g /path/to/GenomeAnalysisTK.jar"

EOF

echo -e "ConfigFile: $ConfigFile\n"
echo -e "\n****************** CreateConfigFile Finished ******************\n"
