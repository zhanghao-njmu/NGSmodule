#!/usr/bin/env bash

cat << EOF >$1
#!/usr/bin/env bash

############# Global Paramaters ###########################################################################
maindir="$(pwd)"        ## Absolute path.
rawdata_dir="$(pwd)/rawdata/"   ## Absolute path of dir containing the raw fastq.gz data.
SequenceType="rna"                            ## rna,dna,BSdna
total_threads=$(grep 'processor' /proc/cpuinfo | sort -u | wc -l)                             ## Total threads for use.
ntask_per_run="ALL"                           ## "ALL" or numeric value to specify the number of tasks run simultaneously at the backend.
SampleInfoFile=""                             ## Absolute path of a .csv file or leave with blank when there is no need to rename the sample.
SampleGrepPattern=""                          ## Optional. Perl-compatible regexps used for matching the SampleName under the work dir.


############# Rscript path ################################################################################
Rscript="/usr/local/bin/Rscript"


############# PrepareWorkDir Paramaters ###################################################################
### raw_file_name= SampleIdPattern  + SampleSufixPattern
### Example: R19051060_BKDL190818861-1a_1.fq.gz and R19051060_BKDL190818861-1a_2.fq.gz
### SampleIdPattern="R.*-1a"
### SampleSufixPattern="_BKDL.*_\d\.fq\.gz"
SampleIdPattern="R.*"                                              ## This argument must be same pattern with the SampleID column in the SampleInfoFile.
SampleSufixPattern="(_BKDL.*_\d\.fastq\.gz|_BKDL.*\.fastq\.gz)"    ## SE must end with fq.gz or .fastq.gz. PE must end with _1.fastq.gz,_1.fq.gz,_R1.fastq.gz,_R1.fq.gz


############# preAlignmentQC Paramaters ###################################################################
### Fastp ###
trim_front1=10                 ## trimming how many bases in front for read1.6-10 bp for RNAseq and WGBS.
trim_tail1=0                   ## trimming how many bases in tail for read1.
trim_front2=10                 ## trimming how many bases in front for read2. Only valid when layout=PE.
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
SortmeRNA_Dir="/data/database/SortmeRNA"       ## SortmeRNA_ref: the dir containing the reference sequence.
SortmeRNA_Type="rRNA"                          ## SortmeRNA_ref: rRNA,Mt_tRNA,Mt_rRNA
SortmeRNA_Species="mouse"                      ## SortmeRNA_ref: human,mouse,machin,rhesus,fly
SortmeRNA_DataVersion="EnsemblGenes98"         ## SortmeRNA_ref: version of sortmerna fetched sequence 

SortmeRNA_ref_direct=""                        ## Optional. Directly specify the path of the SortmeRNA_ref sequence file. 


############# Alignment Paramaters ##########################################################################
iGenomes_Dir="/data/database/iGenomes"         ## The iGenomes dir
Species="mouse"                                ## human,mouse,machin,rhesus,fly
Database="Ensembl"                             ## Ensembl,NCBI,UCSC
Genome_build="GRCm38"                          ## The genome version under the dir GenomeDir/Species_arr[Species]/Database.
Genome_name="genome.fa"                        ## genome.fa,genome_main.fa
Aligner="hisat2"                               ## bwa,bowtie,bowtie2,hisat2,tophat2,star,bismark_bowtie2,bismark_hisat2


############# Quantification Paramaters ######################################################################
strandspecific=0                               ## 0(unstranded),1(stranded),2(reversely stranded)


############# DifferentialExpression Paramaters ##############################################################
max_padj=0.05                                      ## Typically 0.05 or 0.01 or 0.001
min_fc=2                                           ## Typically 2
min_count=10                                       ## Minimum count required for at least n samples (n is the smallest group sample size).

group_compare="Hom-80S,WT-80S;Hom-Input,WT-Input;" ## Groups are seperated by comma(,). Different comparisons are seperated by semicolon(;). 

DGEs_multi_compare=1                               ## Whether to compare DGEs among different comparisons. 0(not to do),1(do).

EOF
echo -e "Task finished \nConfigFile: $1\n"

