#!/usr/bin/env bash

#######################################################################################
trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

bwa &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command bwa.\n"
  exit 1
}
bowtie &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command bowtie.\n"
  exit 1
}
hisat2 &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command hisat2.\n"
  exit 1
}
STAR &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command STAR.\n"
  exit 1
}
bismark &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command bismark.\n"
  exit 1
}
samtools &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command samtools.\n"
  exit 1
}
sambamba &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command sambamba.\n"
  exit 1
}
picard &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command picard.\n"
  exit 1
}

if [[ "$SequenceType" == "BSdna" ]] && [[ ! "$Aligner" =~ bismark_* ]]; then
  color_echo "red" "ERROR! Aligner must be bismark_bowtie2 or bismark_hisat2 for the SequenceType 'BSdna'."
  exit 1
fi

if [[ ! "$SequenceType" == "BSdna" ]] && [[ "$Aligner" =~ bismark_* ]]; then
  color_echo "red" "ERROR! SequenceType must be BSdna for the Aligner '$Aligner'."
  exit 1
fi

aligners=("bwa" "bowtie" "bowtie2" "hisat2" "tophat2" "star" "bismark_bowtie2" "bismark_hisat2")
if [[ " ${aligners[*]} " != *" $Aligner "* ]]; then
  color_echo "red" "ERROR! Aligner is wrong.\nPlease check theParamaters in your ConfigFile.\n"
  exit 1
fi

de_option=("TRUE" "FALSE")
if [[ ${Deduplication} == "" ]]; then
  case ${SequenceType} in
  rna)
    Deduplication="FALSE"
    ;;
  dna)
    Deduplication="TRUE"
    ;;
  BSdna)
    Deduplication="TRUE"
    ;;
  *)
    Deduplication="FALSE"
    ;;
  esac
elif [[ " ${de_option[*]} " != *" $Deduplication "* ]]; then
  color_echo "red" "ERROR! Deduplication must be empty or one of 'TRUE' and 'FALSE'.\nPlease check the paramaters in your ConfigFile.\n"
  exit 1
fi

if [[ $Genome_direct == "" ]]; then
  genome="$iGenomes_Dir/$Species/$Source/$Build/Sequence/WholeGenomeFasta/genome.fa"
else
  genome=$Genome_direct
fi

if [[ $GTF_direct == "" ]]; then
  gtf="$iGenomes_Dir/$Species/$Source/$Build/Annotation/Genes/genes.gtf"
else
  gtf=$GTF_direct
fi

if [[ ! -f $genome ]]; then
  color_echo "red" "ERROR! Cannot find the genome file: $genome\nPlease check the Alignment Paramaters in your ConfigFile.\n"
  exit 1
elif [[ ! -f $gtf ]]; then
  color_echo "red" "ERROR! Cannot find the gtf file: $gtf\nPlease check the Alignment Paramaters in your ConfigFile.\n"
  exit 1
fi

bwa_index="$iGenomes_Dir/$Species/$Source/$Build/Sequence/BWAIndex/genome.fa"
bowtie_index="$iGenomes_Dir/$Species/$Source/$Build/Sequence/BowtieIndex/genome"
bowtie2_index="$iGenomes_Dir/$Species/$Source/$Build/Sequence/Bowtie2Index/genome"
hisat2_index="$iGenomes_Dir/$Species/$Source/$Build/Sequence/Hisat2Index/genome"
star_index="$iGenomes_Dir/$Species/$Source/$Build/Sequence/STARIndex/genome"
bismark_bowtie2_index="$iGenomes_Dir/$Species/$Source/$Build/Sequence/BismarkIndex/bowtie2"
bismark_hisat2_index="$iGenomes_Dir/$Species/$Source/$Build/Sequence/BismarkIndex/hisat2"
tophat2_index=$bowtie2_index

if [[ $Index_direct == "" ]]; then
  eval "index=\${${Aligner}_index}"
else
  index=$Index_direct
fi

echo -e "############################# Alignment Parameters #############################\n"
echo -e "  SequenceType: ${SequenceType}\n  Aligner: ${Aligner}\n  Deduplication: ${Deduplication}\n"
echo -e "  Genome_File: ${genome}\n  GTF_File: ${gtf}\n  Aligner_Index: ${index}\n"
echo -e "################################################################################\n"

echo -e "****************** Start Alignment ******************\n"
SECONDS=0

for sample in "${arr[@]}"; do
  read -u1000
  {
    dir="$work_dir/$sample"
    mkdir -p "$dir/Alignment-$Aligner"
    cd "$dir/Alignment-$Aligner"

    Layout=${Layout_dict[${sample}]}
    force=${force_complete}
    status="uncompleted"
    attempt=0

    echo "+++++ ${sample} +++++"

    while [[ $status == "uncompleted" ]] && (("$attempt" <= 1)); do
      ((attempt++))
      if [[ $attempt != 1 ]]; then
        echo -e "+++++ ${sample}: Number of attempts: $attempt +++++"
      fi

      logfiles=("AlignmentStatus.log" "BAMprocessStatus.log")
      globalcheck_logfile "$dir/Alignment-$Aligner" logfiles[@] "$force" "$error_pattern" "$complete_pattern" "$sample"

      check_logfile "$sample" "Alignment" "$dir/Alignment-$Aligner/AlignmentStatus.log" "$error_pattern" "$complete_pattern" "precheck"
      if [[ $? == 1 ]]; then
        rm -rf $dir/Alignment-$Aligner/*
        touch "$dir/Alignment-$Aligner/AlignmentStatus.log"

        if [[ $Layout == "SE" ]]; then
          fq1=$dir/${sample}_trim.fq.gz
          if [[ "$Aligner" = "bwa" ]]; then
            bwa mem -t $threads -M $index ${fq1} |
              samtools view -@ $threads -Shb - |
              samtools sort -@ $threads - >${sample}.${Aligner}.bam 2>>"$dir/Alignment-$Aligner/AlignmentStatus.log"
          elif [[ "$Aligner" == "bowtie" ]]; then
            bowtie -p $threads -1 ${fq1} -l 22 --fullref --chunkmbs 512 --best --strata -m 20 -n 2 --mm $index -S |
              samtools view -@ $threads -Shb - |
              samtools sort -@ $threads - >${sample}.${Aligner}.bam 2>>"$dir/Alignment-$Aligner/AlignmentStatus.log"
          elif [[ "$Aligner" == "bowtie2" ]]; then
            bowtie2 -p $threads -x $index -1 ${fq1} 2>${sample}.${Aligner}.log |
              samtools view -@ $threads -Shb - |
              samtools sort -@ $threads - >${sample}.${Aligner}.bam 2>>"$dir/Alignment-$Aligner/AlignmentStatus.log"
          elif [[ "$Aligner" == "hisat2" ]]; then
            hisat2 -p $threads -x $index -U ${fq1} --new-summary 2>${sample}.${Aligner}.log |
              samtools view -@ $threads -Shb - |
              samtools sort -@ $threads - >${sample}.${Aligner}.bam 2>>"$dir/Alignment-$Aligner/AlignmentStatus.log"
          elif [[ "$Aligner" == "tophat2" ]]; then
            tophat2 -p $threads --GTF $gtf --output-dir ./ $index ${fq1}
            samtools view -@ $threads -Shb accepted_hits.bam |
              samtools sort -@ $threads - >${sample}.${Aligner}.bam 2>>"$dir/Alignment-$Aligner/AlignmentStatus.log"
            rm -f accepted_hits.bam
          elif [[ "$Aligner" == "star" ]]; then
            STAR --runThreadN $threads --genomeDir $index --readFilesIn ${fq1} --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 \
              --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
              --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
              --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
              --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesCommand zcat \
              --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM
            samtools view -@ $threads -Shb Aligned.sortedByCoord.out.bam |
              samtools sort -@ $threads - >${sample}.${Aligner}.bam 2>>"$dir/Alignment-$Aligner/AlignmentStatus.log"
            rm -f Aligned.sortedByCoord.out.bam
          elif [[ "$Aligner" == "bismark_bowtie2" ]]; then
            bismark --bowtie2 --multicore $bismark_threads -p 4 --genome $index ${fq1} --quiet \
              --non_directional --nucleotide_coverage \
              --output_dir $dir/Alignment-$Aligner &>>"$dir/Alignment-$Aligner/AlignmentStatus.log"
            for file in ./*_trim*; do mv $file ${file//_trim/}; done
          elif [[ "$Aligner" == "bismark_hisat2" ]]; then
            bismark --hisat2 --multicore $bismark_threads -p 4 --genome $index ${fq1} --quiet \
              --non_directional --nucleotide_coverage \
              --output_dir $dir/Alignment-$Aligner &>>"$dir/Alignment-$Aligner/AlignmentStatus.log"
            for file in ./*_trim*; do mv $file ${file//_trim/}; done
          fi

        elif [[ $Layout == "PE" ]]; then
          fq1=$dir/${sample}_1_trim.fq.gz
          fq2=$dir/${sample}_2_trim.fq.gz
          if [[ "$Aligner" == "bwa" ]]; then
            bwa mem -t $threads -M $index ${fq1} ${fq2} |
              samtools view -@ $threads -Shb - |
              samtools sort -@ $threads - >${sample}.${Aligner}.bam 2>>"$dir/Alignment-$Aligner/AlignmentStatus.log"
          elif [[ "$Aligner" = "bowtie" ]]; then
            bowtie -p $threads -1 ${fq1} -2 ${fq2} -l 22 --fullref --chunkmbs 512 --best --strata -m 20 -n 2 --mm $index -S |
              samtools view -@ $threads -Shb - |
              samtools sort -@ $threads - >${sample}.${Aligner}.bam 2>>"$dir/Alignment-$Aligner/AlignmentStatus.log"
          elif [[ "$Aligner" == "bowtie2" ]]; then
            bowtie2 -p $threads -x $index -1 ${fq1} -2 ${fq2} 2>${sample}.${Aligner}.log |
              samtools view -@ $threads -Shb - |
              samtools sort -@ $threads - >${sample}.${Aligner}.bam 2>>"$dir/Alignment-$Aligner/AlignmentStatus.log"
          elif [[ "$Aligner" == "hisat2" ]]; then
            hisat2 -p $threads -x $index -1 ${fq1} -2 ${fq2} --new-summary 2>${sample}.${Aligner}.log |
              samtools view -@ $threads -Shb - |
              samtools sort -@ $threads - >${sample}.${Aligner}.bam 2>>"$dir/Alignment-$Aligner/AlignmentStatus.log"
          elif [[ "$Aligner" == "tophat2" ]]; then
            tophat2 -p $threads --GTF $gtf --output-dir ./ $index ${fq1} ${fq2}
            samtools view -@ $threads -Shb accepted_hits.bam |
              samtools sort -@ $threads - >${sample}.${Aligner}.bam 2>>"$dir/Alignment-$Aligner/AlignmentStatus.log"
            rm -f accepted_hits.bam
          elif [[ "$Aligner" == "star" ]]; then
            STAR --runThreadN $threads --genomeDir $index --readFilesIn ${fq1} ${fq2} --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 \
              --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
              --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
              --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
              --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesCommand zcat \
              --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM
            samtools view -@ $threads -Shb Aligned.sortedByCoord.out.bam |
              samtools sort -@ $threads - >${sample}.${Aligner}.bam 2>>"$dir/Alignment-$Aligner/AlignmentStatus.log"
            rm -f Aligned.sortedByCoord.out.bam
          elif [[ "$Aligner" == "bismark_bowtie2" ]]; then
            bismark --bowtie2 --multicore $bismark_threads -p 4 --genome $index -1 ${fq1} -2 ${fq2} --quiet \
              --non_directional --nucleotide_coverage \
              --output_dir $dir/Alignment-$Aligner &>>"$dir/Alignment-$Aligner/AlignmentStatus.log"
            for file in ./*_1_trim*; do mv $file ${file//_1_trim/}; done
          elif [[ "$Aligner" == "bismark_hisat2" ]]; then
            bismark --hisat2 --multicore $bismark_threads -p 4 --genome $index -1 ${fq1} -2 ${fq2} --quiet \
              --non_directional --nucleotide_coverage \
              --output_dir $dir/Alignment-$Aligner &>>"$dir/Alignment-$Aligner/AlignmentStatus.log"
            for file in ./*_1_trim*; do mv $file ${file//_1_trim/}; done
          fi

        else
          color_echo "yellow" "ERROR! ${sample}: Cannot determine the layout of sequencing data!"
          attempt=2
          echo "ERROR! ${sample}: Cannot determine the layout of sequencing data!" >>"$dir/Alignment-$Aligner/AlignmentStatus.log"
          continue
        fi

        check_logfile "$sample" "Alignment" "$dir/Alignment-$Aligner/AlignmentStatus.log" "$error_pattern" "$complete_pattern" "postcheck"
        if [[ $? == 1 ]]; then
          force="TRUE"
          continue
        fi
      fi

      check_logfile "$sample" "BamProcessing" "$dir/Alignment-$Aligner/BamProcessingStatus.log" "$error_pattern" "$complete_pattern" "precheck"
      if [[ $? == 1 ]]; then
        rm -f "$dir/Alignment-$Aligner/BamProcessingStatus.log"
        touch "$dir/Alignment-$Aligner/BamProcessingStatus.log"

        BAM=$(ls $dir/Alignment-$Aligner/*.bam)
        samtools quickcheck -v ${BAM} &>/dev/null
        if [[ $? != 0 ]]; then
          color_echo "yellow" "[INFO] $sample: BAM file checked failed."
          force="TRUE"
          continue
        fi

        echo "+++++ Samtools stat: $sample +++++" | tee -a "$dir/Alignment-$Aligner/BamProcessingStatus.log"
        if [[ "$SequenceType" == "BSdna" ]] && [[ "$Aligner" =~ bismark_* ]]; then
          samtools stats -@ $threads $BAM >${BAM}.stats 2>>"$dir/Alignment-$Aligner/BamProcessingStatus.log"
          samtools flagstat -@ $threads $BAM >${BAM}.flagstat 2>>"$dir/Alignment-$Aligner/BamProcessingStatus.log"
        else
          samtools index -@ $threads $BAM 2>"$dir/Alignment-$Aligner/BamProcessingStatus.log"
          samtools stats -@ $threads $BAM >$BAM.stats 2>>"$dir/Alignment-$Aligner/BamProcessingStatus.log"
          samtools idxstats -@ $threads $BAM >$BAM.idxstats 2>>"$dir/Alignment-$Aligner/BamProcessingStatus.log"
          samtools flagstat -@ $threads $BAM >$BAM.flagstat 2>>"$dir/Alignment-$Aligner/BamProcessingStatus.log"
        fi

        if [[ "$SequenceType" == "dna" ]]; then
          if [[ "$Deduplication" == "TRUE" ]]; then
            echo "+++++ DNA remove duplicates: $sample +++++" | tee -a "$dir/Alignment-$Aligner/BamProcessingStatus.log"
            sambamba markdup -r -t $threads --hash-table-size=3000000 --overflow-list-size=3000000 $BAM ${sample}.${Aligner}.dedup.bam &>>"$dir/Alignment-$Aligner/BamProcessingStatus.log"
            BAM=${sample}.${Aligner}.dedup.bam
          else
            echo "+++++ DNA mark duplicates: $sample +++++" | tee -a "$dir/Alignment-$Aligner/BamProcessingStatus.log"
            sambamba markdup -t $threads --hash-table-size=3000000 --overflow-list-size=3000000 $BAM ${sample}.${Aligner}.markdup.bam &>>"$dir/Alignment-$Aligner/BamProcessingStatus.log"
            BAM=${sample}.${Aligner}.markdup.bam
          fi
            picard AddOrReplaceReadGroups I=$BAM O=${sample}.${Aligner}.RG.bam RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$sample &>>"$dir/Alignment-$Aligner/BamProcessingStatus.log"
            picard FixMateInformation I=${sample}.${Aligner}.RG.bam O=${sample}.${Aligner}.RG.FMI.bam ADD_MATE_CIGAR=true &>>"$dir/Alignment-$Aligner/BamProcessingStatus.log"
            rm -f ${sample}.${Aligner}.RG.bam
            mv ${sample}.${Aligner}.RG.FMI.bam $BAM
            samtools index -@ $threads $BAM 2>>"$dir/Alignment-$Aligner/BamProcessingStatus.log"
        fi

        if [[ "$SequenceType" == "rna" ]]; then
          if [[ "$Deduplication" == "TRUE" ]]; then
            echo "+++++ RNA remove duplicates: $sample +++++" | tee -a "$dir/Alignment-$Aligner/BamProcessingStatus.log"
            sambamba markdup -r -t $threads --hash-table-size=3000000 --overflow-list-size=3000000 $BAM ${sample}.${Aligner}.dedup.bam &>>"$dir/Alignment-$Aligner/BamProcessingStatus.log"
            BAM=${sample}.${Aligner}.dedup.bam
          else
            echo "+++++ RNA mark duplicates: $sample +++++" | tee -a "$dir/Alignment-$Aligner/BamProcessingStatus.log"
            sambamba markdup -t $threads --hash-table-size=3000000 --overflow-list-size=3000000 $BAM ${sample}.${Aligner}.markdup.bam &>>"$dir/Alignment-$Aligner/BamProcessingStatus.log"
            BAM=${sample}.${Aligner}.markdup.bam
          fi
          samtools index -@ $threads $BAM 2>>"$dir/Alignment-$Aligner/BamProcessingStatus.log"
        fi

        if [[ "$SequenceType" == "BSdna" ]] && [[ "$Aligner" =~ bismark_* ]]; then
          if [[ "$Deduplication" == "TRUE" ]]; then
            echo "+++++ BSdna remove duplicates: $sample +++++" | tee -a "$dir/Alignment-$Aligner/BamProcessingStatus.log"
            mkdir -p $dir/Alignment-$Aligner/deduplicate_bismark
            deduplicate_bismark --bam $BAM --output_dir $dir/Alignment-$Aligner/deduplicate_bismark &>>"$dir/Alignment-$Aligner/BamProcessingStatus.log"
            BAM=$(ls $dir/Alignment-$Aligner/deduplicate_bismark/*.deduplicated.bam)
            dedup_report=$(ls $dir/Alignment-$Aligner/deduplicate_bismark/*.deduplication_report.txt)
            samtools quickcheck -v ${BAM} &>/dev/null
            if [[ $? != 0 ]]; then
              color_echo "yellow" "[INFO] $sample: BS-seq deduplicated.bam check failed."
              continue
            fi
          else
            dedup_report='none'
          fi

          echo "+++++ BS-seq methylation extractor: $sample +++++" | tee -a "$dir/Alignment-$Aligner/BamProcessingStatus.log"
          mkdir -p $dir/Alignment-$Aligner/bismark_methylation_extractor
          bismark_methylation_extractor --multicore $bismark_threads --gzip --comprehensive --merge_non_CpG \
            --bedGraph --buffer_size 10G \
            --cytosine_report --genome_folder $index \
            --output $dir/Alignment-$Aligner/bismark_methylation_extractor $BAM &>>"$dir/Alignment-$Aligner/BamProcessingStatus.log"

          echo "+++++ BS-seq html processing report: $sample +++++" | tee -a "$dir/Alignment-$Aligner/BamProcessingStatus.log"
          mkdir -p $dir/Alignment-$Aligner/bismark2report
          alignment_report=$(ls $dir/Alignment-$Aligner/*_[SP]E_report.txt)
          splitting_report=$(ls $dir/Alignment-$Aligner/bismark_methylation_extractor/*_splitting_report.txt)
          mbias_report=$(ls $dir/Alignment-$Aligner/bismark_methylation_extractor/*M-bias.txt)
          nucleotide_report=$(ls $dir/Alignment-$Aligner/*.nucleotide_stats.txt)
          bismark2report --dir $dir/Alignment-$Aligner/bismark2report \
            --alignment_report $alignment_report \
            --dedup_report $dedup_report \
            --splitting_report $splitting_report \
            --mbias_report $mbias_report \
            --nucleotide_report $nucleotide_report &>>"$dir/Alignment-$Aligner/BamProcessingStatus.log"
        fi

        check_logfile "$sample" "BamProcessing" "$dir/Alignment-$Aligner/BamProcessingStatus.log" "$error_pattern" "$complete_pattern" "postcheck"
        if [[ $? == 1 ]]; then
          continue
        fi
      fi

      status="completed"
      color_echo "blue" "+++++ ${sample}: Alignment completed +++++"
    done

    if [[ "$status" == "completed" ]]; then
      echo "Completed: $sample" >>"$tmpfile"
    else
      echo "Interrupted: $sample" >>"$tmpfile"
      color_echo "red" "ERROR! ${sample} interrupted! Please check the processing log and your raw fastq file."
    fi

    color_echo "green" "***** Completed:$(cat "$tmpfile" | grep "Completed" | uniq | wc -l) | Interrupted:$(cat "$tmpfile" | grep "Interrupted" | uniq | wc -l) | Total:$total_task *****"

    echo >&1000
  } &
  ((bar++))
  processbar $bar $total_task
done
wait

ninterrupted=$(cat "$tmpfile" | grep "Interrupted" | uniq | wc -l)
if [[ $ninterrupted != 0 ]]; then
  cat "$tmpfile" | grep "Interrupted" | uniq >$maindir/Alignment.Interrupted.txt
  color_echo "red" "\n\n################################################################################"
  color_echo "red" "    $ninterrupted of $total_task tasks interrupted."
  color_echo "red" "    Please check the samples in $maindir/Alignment.Interrupted.txt"
  color_echo "red" "################################################################################\n\n"
fi

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "****************** Alignment Done ******************\n"
