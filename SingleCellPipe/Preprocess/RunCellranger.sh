#!/usr/bin/env bash
trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

cellranger &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command cellranger.\n"
  exit 1
}
velocyto &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command velocyto.\n"
  exit 1
}
dropest &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command dropest.\n"
  exit 1
}
dropReport.Rsc &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command dropReport.Rsc.\n"
  exit 1
}
Rscript &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command Rscript.\n"
  exit 1
}

R_packages=("DropletUtils" "dropestr" "ggplot2" "ggupset" "ggsci" "cowplot" "dplyr" "reshape2" "SingleCellExperiment" "grid" "png" "gridExtra" "scales")
for package in "${R_packages[@]}"; do 
  Rscript -e "installed.packages()" | awk '{print $1}' | grep $package &>/dev/null
  [ $? -ne 0 ] && {
    color_echo "red" "Cannot find the R package $package.\n"
    exit 1
  }
done

if [[ ! -d $cellranger_ref ]]; then
  color_echo "red" "ERROR! Cannot find the cellranger_ref directory: $cellranger_ref\nPlease check the path.\n"
  exit 1
elif [[ ! -f $gene_gtf ]]; then
  color_echo "red" "ERROR! Cannot find the gene_gtf file: $gtf\nPlease check the path.\n"
  exit 1
elif [[ ! -f $rmsk_gtf ]]; then
  color_echo "red" "ERROR! Cannot find the rmsk_gtf file: $rmsk_gtf\nPlease check the path.\n"
  exit 1
elif [[ ! -f $dropEst_config ]]; then
  color_echo "red" "ERROR! Cannot find the dropEst_config file: $dropEst_config\nPlease check the path.\n"
  exit 1
fi

echo -e "############################# Alignment Parameters #############################\n"
echo -e "  SequenceType: ${SequenceType}\n  Aligner: ${Aligner}\n"
echo -e "  Genome_File: ${genome}\n  GTF_File: ${gtf}\n  Aligner_Index: ${index}\n"
echo -e "################################################################################\n"

echo -e "****************** Start Alignment ******************\n"
SECONDS=0


arr=($(find $rawdata_dir -name "*.fastq.gz" | sed "s/\/.*\///g"| sed "s/_2020.*//g" | sort | uniq ))

for sample in "${arr[@]}"; do
do
{
cd $cellranger_dir
files=($(find $rawdata_dir -name "${id}_*.fastq.gz" |grep -P "(?<=rawdata\/).*(?=_S\d+)" -o | sort | uniq ))
sample_run=$(printf ",%s" "${files[@]}")
sample_run=${sample_run:1}
echo -e "$id: $sample_run"

# cellranger
cellranger count --id ${id} \
                --fastqs ${rawdata_dir} \
                --sample ${sample_run} \
                --localcores $threads \
                --localmem $((threads*2)) \
                --transcriptome $cellranger_ref

# velocyto                 
velocyto run10x -m $rmsk_gtf --samtools-threads $threads $cellranger_dir/$id $gene_gtf
 
# dropEst
mkdir -p $cellranger_dir/$id/dropEst 
cd $cellranger_dir/$id/dropEst
dropest -f -g $gene_gtf -c $dropEst_config $cellranger_dir/$id/outs/possorted_genome_bam.bam
dropReport.Rsc  $cellranger_dir/$id/dropEst/cell.counts.rds

# droplet-QC
Rscript $1

}&
done

wait 
echo "Done"

