#!/usr/bin/env bash

#######################################################################################
trap 'j=`ps aux | grep -P "$maindir" |grep -P "(bwa)|(bowtie)|(hisat)|(tophat)|(STAR)|(bismark)|(samtools)|(sambamba)|(picard)"| awk '"'"'{print $2}'"'"'`;kill -9 $j;kill -9 $(jobs -p);echo -e "\nKilling all background processes......\nExiting the script......\n";exit 1' SIGINT


bwa &>/dev/null;[ $? -eq 127 ] && { echo -e "Cannot find the command bwa.\n";exit 1; }
bowtie --version &>/dev/null;[ $? -ne 0 ] && { echo -e "Cannot find the command bowtie.\n";exit 1; }
hisat2 --version &>/dev/null;[ $? -ne 0 ] && { echo -e "Cannot find the command hisat2.\n";exit 1; }
STAR --version  &>/dev/null;[ $? -ne 0 ] && { echo -e "Cannot find the command STAR.\n";exit 1; }
bismark --version &>/dev/null;[ $? -ne 0 ] && { echo -e "Cannot find the command bismark.\n";exit 1; }
samtools --version &>/dev/null;[ $? -ne 0 ] && { echo -e "Cannot find the command samtools.\n";exit 1; }
sambamba --version &>/dev/null;[ $? -ne 0 ] && { echo -e "Cannot find the command sambamba.\n";exit 1; }
picard &>/dev/null;[ $? -eq 127 ] && { echo -e "Cannot find the command picard.\n";exit 1; }
bam &>/dev/null;[ $? -eq 127 ] && { echo -e "Cannot find the command bam. User can install mosdepth by 'conda install -c bioconda bamutil'.\n";exit 1; }

eval "index=\${${aligner}_index}"
echo -e "############################# Alignment Parameters #############################\n"
echo -e "  Sequencing: ${Sequencing}\n  Aligner: ${aligner}\n  iGenomes_Dir: ${iGenomes_Dir}\n  Species: ${Species}(${Species_arr[$Species]})\n  Database: ${Database}\n  Genome_build: ${Genome_build}\n  Genome_name: ${Genome_name}\n"
echo -e "  Genome_File: ${genome}\n  GTF_File: ${gtf}\n  Aligner_Index: ${index}\n"
echo -e "################################################################################\n"

echo -e "****************** Start Alignment ******************\n"
SECONDS=0

if [[ ! -d $work_dir ]];then
  echo -e "Cannot find the workdir: ${work_dir}\n"
  exit 1
fi

for sample in ${arr[@]};do
  read -u1000
  {
  dir=$work_dir/$sample
  mkdir -p $dir/$aligner; cd $dir/$aligner
  
 	echo "+++++ Alignment: $sample +++++"
  if [[ $layout == "SE" ]]; then
    fq1=$dir/${sample}_trim.fq.gz
    if [[ "$aligner" = "bwa" ]];then
      bwa mem -t $threads -M $index ${fq1} | samtools view -@ $threads -Shb - | samtools sort -@ $threads - >${sample}.${aligner}.bam 2>/dev/null
    elif [[ "$aligner" == "bowtie" ]];then
      bowtie -p $threads -1 ${fq1} -l 22 --fullref --chunkmbs 512 --best --strata -m 20 -n 2 --mm $index -S | samtools view -@ $threads -Shb - | samtools sort -@ $threads - >${sample}.${aligner}.bam 2>/dev/null
    elif [[ "$aligner" == "bowtie2" ]];then
      bowtie2 -p $threads -x $index -1 ${fq1} 2>${sample}.${aligner}.log | samtools view -@ $threads -Shb - | samtools sort -@ $threads - >${sample}.${aligner}.bam 2>/dev/null
    elif [[ "$aligner" == "hisat2" ]];then
      hisat2 -p $threads -x $index -U ${fq1} --new-summary 2>${sample}.${aligner}.log | samtools view -@ $threads -Shb - | samtools sort -@ $threads - >${sample}.${aligner}.bam 2>/dev/null
    elif [[ "$aligner" == "tophat2" ]];then
      tophat2 -p $threads --GTF $gtf --output-dir ./ $index ${fq1}
      mv accepted_hits.sam ${sample}.${aligner}.bam
    elif [[ "$aligner" == "star" ]];then
      STAR --runThreadN $threads --genomeDir $index --readFilesIn ${fq1} --genomeLoad LoadAndKeep  --limitBAMsortRAM 10000000000 \
           --outSAMunmapped Within  --outFilterType BySJout  --outSAMattributes NH HI AS NM MD  \
           --outFilterMultimapNmax 20  --outFilterMismatchNmax 999  --outFilterMismatchNoverReadLmax 0.04 \
           --alignIntronMin 20  --alignIntronMax 1000000  --alignMatesGapMax 1000000   \
           --alignSJoverhangMin 8   --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesCommand zcat \
           --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
      mv Aligned.sortedByCoord.out.bam ${sample}.${aligner}.bam  
    elif [[ "$aligner" == "bismark_bowtie2" ]];then
      bismark --bowtie2 --multicore $((($threads)/8)) -p 3 --genome $index ${fq1} --quiet \
              --non_directional --nucleotide_coverage \
              --output_dir $dir/$aligner 2>$dir/$aligner/bismark.log
      for file in ./*_trim*; do mv $file ${file//_trim/};done
    elif [[ "$aligner" == "bismark_hisat2" ]];then
      bismark --hisat2 --multicore $((($threads)/8)) -p 3 --genome $index ${fq1} --quiet \
              --non_directional --nucleotide_coverage \
              --output_dir $dir/$aligner 2>$dir/$aligner/bismark.log
      for file in ./*_trim*; do mv $file ${file//_trim/};done
    fi
    
  elif [[ $layout == "PE" ]]; then
    fq1=$dir/${sample}_1_trim.fq.gz
    fq2=$dir/${sample}_2_trim.fq.gz
    if [[ "$aligner" == "bwa" ]];then
      bwa mem -t $threads -M $index ${fq1} ${fq2} | samtools view -@ $threads -Shb - | samtools sort -@ $threads - >${sample}.${aligner}.bam 2>/dev/null
    elif [[ "$aligner" = "bowtie" ]];then
      bowtie -p $threads -1 ${fq1} -2 ${fq2} -l 22 --fullref --chunkmbs 512 --best --strata -m 20 -n 2 --mm $index -S | samtools view -@ $threads -Shb - | samtools sort -@ $threads - >${sample}.${aligner}.bam 2>/dev/null
    elif [[ "$aligner" == "bowtie2" ]];then
      bowtie2 -p $threads -x $index -1 ${fq1} -2 ${fq2} 2>${sample}.${aligner}.log | samtools view -@ $threads -Shb - | samtools sort -@ $threads - >${sample}.${aligner}.bam 2>/dev/null
    elif [[ "$aligner" == "hisat2" ]];then
      hisat2 -p $threads -x $index -1 ${fq1} -2 ${fq2} --new-summary 2>${sample}.${aligner}.log | samtools view -@ $threads -Shb - | samtools sort -@ $threads - >${sample}.${aligner}.bam 2>/dev/null
    elif [[ "$aligner" == "tophat2" ]];then
      tophat2 -p $threads --GTF $gtf --output-dir ./ $index ${fq1} ${fq2} 
      mv accepted_hits.sam ${sample}.${aligner}.bam
    elif [[ "$aligner" == "star" ]];then
      STAR --runThreadN $threads --genomeDir $index --readFilesIn ${fq1} ${fq2} --genomeLoad LoadAndKeep  --limitBAMsortRAM 10000000000 \
           --outSAMunmapped Within  --outFilterType BySJout  --outSAMattributes NH HI AS NM MD  \
           --outFilterMultimapNmax 20  --outFilterMismatchNmax 999  --outFilterMismatchNoverReadLmax 0.04 \
           --alignIntronMin 20  --alignIntronMax 1000000  --alignMatesGapMax 1000000   \
           --alignSJoverhangMin 8   --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesCommand zcat \
           --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
      mv Aligned.sortedByCoord.out.bam ${sample}.${aligner}.bam  
    elif [[ "$aligner" == "bismark_bowtie2" ]];then
      bismark --bowtie2 --multicore $((($threads)/8)) -p 3 --genome $index -1 ${fq1} -2 ${fq2} --quiet \
              --non_directional --nucleotide_coverage \
              --output_dir $dir/$aligner 2>$dir/$aligner/bismark.log
      for file in ./*_1_trim*; do mv $file ${file//_1_trim/};done
    elif [[ "$aligner" == "bismark_hisat2" ]];then
      bismark --hisat2 --multicore $((($threads)/8)) -p 3 --genome $index -1 ${fq1} -2 ${fq2} --quiet \
              --non_directional --nucleotide_coverage \
              --output_dir $dir/$aligner 2>$dir/$aligner/bismark.log
      for file in ./*_1_trim*; do mv $file ${file//_1_trim/};done
    fi
  fi
  
  echo "+++++ $sample: $aligner done +++++"

 	echo "+++++ Bam processing: $sample +++++"
  if [[ "$Sequencing" == "bsseq" ]] && [[ "$aligner" =~ bismark_* ]];then
    bam=$(ls ./*.bam)
    samtools stats -@ $threads $bam >${bam}.stats
    samtools flagstat -@ $threads $bam >${bam}.flagstat
  else
    samtools index -@ $threads ${sample}.${aligner}.bam
    samtools stats -@ $threads ${sample}.${aligner}.bam >${sample}.${aligner}.bam.stats
    samtools idxstats -@ $threads ${sample}.${aligner}.bam >${sample}.${aligner}.bam.idxstats
    samtools flagstat -@ $threads ${sample}.${aligner}.bam >${sample}.${aligner}.bam.flagstat
  fi
  
  if [[ "$Sequencing" == "wgs" ]];then
    echo "+++++ WGS deduplication: $sample +++++"
    sambamba markdup -r -t $threads ${sample}.${aligner}.bam ${sample}.${aligner}.dedup.bam
    picard AddOrReplaceReadGroups I=${sample}.${aligner}.dedup.bam O=${sample}.${aligner}.dedup.RG.bam RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$sample
    picard FixMateInformation I=${sample}.${aligner}.dedup.RG.bam O=${sample}.${aligner}.dedup.bam ADD_MATE_CIGAR=true;
    rm -f ${sample}.${aligner}.dedup.RG.bam 
    samtools index -@ $threads ${sample}.${aligner}.dedup.bam
  fi
  
  if [[ "$Sequencing" == "rnaseq" ]];then
    echo "+++++ RNAseq Mark Duplicates: $sample +++++"
    bam dedup --force --noPhoneHome --in ${sample}.${aligner}.bam --out ${sample}.${aligner}.markdup.bam --log ${sample}.${aligner}.markdup.log
    mv ${sample}.${aligner}.markdup.bam ${sample}.${aligner}.bam
    samtools index -@ $threads ${sample}.${aligner}.bam
  fi
  
  if [[ "$Sequencing" == "bsseq" ]] && [[ "$aligner" =~ bismark_* ]];then
    echo "+++++ BS-seq deduplication: $sample +++++"
    mkdir -p $dir/$aligner/deduplicate_bismark
    bam=$(ls $dir/$aligner/*.bam)
    deduplicate_bismark --bam $bam --output_dir $dir/$aligner/deduplicate_bismark 2>$dir/$aligner/deduplicate_bismark/deduplicate_bismark.log
    
    echo "+++++ BS-seq methylation extractor: $sample +++++"
    mkdir -p $dir/$aligner/bismark_methylation_extractor
    bam=$(ls $dir/$aligner/deduplicate_bismark/*.deduplicated.bam)
    bismark_methylation_extractor --multicore $((($threads)/2)) --gzip --comprehensive --merge_non_CpG \
                                  --bedGraph --buffer_size 10G \
                                  --cytosine_report --genome_folder $index \
                                  --output $dir/$aligner/bismark_methylation_extractor $bam 2>$dir/$aligner/bismark_methylation_extractor/bismark_methylation_extractor.log
                                  
    echo "+++++ BS-seq html processing report: $sample +++++"
    mkdir -p $dir/$aligner/bismark2report
    alignment_report=$(ls $dir/$aligner/*_[SP]E_report.txt)
    dedup_report=$(ls $dir/$aligner/deduplicate_bismark/*.deduplication_report.txt)
    splitting_report=$(ls $dir/$aligner/bismark_methylation_extractor/*_splitting_report.txt)
    mbias_report=$(ls $dir/$aligner/bismark_methylation_extractor/*M-bias.txt)
    nucleotide_report=$(ls $dir/$aligner/*.nucleotide_stats.txt)
    bismark2report --dir $dir/$aligner/bismark2report \
                   --alignment_report $alignment_report \
                   --dedup_report $dedup_report \
                   --splitting_report $splitting_report \
                   --mbias_report $mbias_report \
                   --nucleotide_report $nucleotide_report
    
  fi
  

  echo >&1000
  }&
  ((bar++))
  processbar $bar $total_task
done
wait

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "****************** Alignment Done ******************\n"
