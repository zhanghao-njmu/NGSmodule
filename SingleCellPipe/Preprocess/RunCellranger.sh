#!/usr/bin/env bash
trap "trap - SIGTERM && kill -- -$$" SIGINT SIGTERM


################ START ################
arr=($(find $rawdata_dir -name "*.fastq.gz" | sed "s/\/.*\///g"| sed "s/_2020.*//g" | sort | uniq ))

for id in ${arr[@]}
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

