#!/usr/bin/env bash
trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

############# Paramaters #############################################################
# iGenomes_dir="/reference/iGenomes/"
# Species=("Homo_sapiens" "Mus_musculus" "Macaca_mulatta")
# Sources=("Ensembl" "NCBI" "UCSC")
# kmers=(150 100 50)
# windows=(2000000 1000000 500000 100000)
# total_threads=384
# ntask_per_run="ALL"

######## check command available ##########
samtools &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the tool samtools. User can install it with the command 'conda install -c bioconda samtools'.\n"
  exit 1
}
faidx &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the tool pyfaidx. User can install it with the command 'conda install -c bioconda pyfaidx'.\n"
  exit 1
}
gem-mappability &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command gem-mappability. User can install it from 'https://sourceforge.net/projects/gemlibrary'.\n"
  exit 1
}
genmap &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command genmap. User can install it with the command 'conda install -c bioconda genmap'.\n"
  exit 1
}
wigToBigWig &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command wigToBigWig. User can install it with the command 'conda install -c bioconda ucsc-wigtobigwig'.\n"
  exit 1
}
mapCounter &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command mapCounter. User can install it from 'https://github.com/shahcompbio/hmmcopy_utils'.\n"
  exit 1
}
picard &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command picard.\n"
  exit 1
}

if [[ -d $iGenomes_dir ]];then
  echo -e ">>> iGenomes_dir exist: $iGenomes_dir \n"
else
  color_echo "yellow" ">>> iGenomes_dir does not exist and will be created: $iGenomes_dir \n"
fi

######## Download the iGenomes #####
for s in "${Species[@]}"; do
  for i in "${Sources[@]}"; do
    color_echo "green" "Downloading the iGenomes: $s/$i"
    igenome="s3://ngi-igenomes/igenomes/$s/$i"
    par=""
    if [[ -d $iGenomes_dir/$s/$i ]]; then
      bismark_exist=($(find $iGenomes_dir/$s/$i -name "IndexStatus.log" | grep -oP "(?<=$i/).*/BismarkIndex/(?=bowtie2)"))
      if [[ "${#bismark_exist[@]}" != 0 ]]; then
        par=$(printf -- " --exclude '*%s*'" "${bismark_exist[@]}")
      fi
    fi
    cmd="aws s3 --no-sign-request sync $igenome $iGenomes_dir/$s/$i --exclude '*/genome.fa' --include '*/WholeGenomeFasta/genome.fa' $par"
    #echo "$cmd"
    eval $cmd

    if [[ ! "$(ls -A $iGenomes_dir/$s/$i)" ]]; then
      color_echo "yellow" "iGenomes do not exist: $s/$i"
      rm -rf $iGenomes_dir/$s/$i
    fi
    if [[ -d $iGenomes_dir/$s/$i ]]; then
      index_dir=($(find $iGenomes_dir/$s/$i -name "*Index" -type d))
      if [[ "${#index_dir[@]}" != 0 ]]; then
        for index in "${index_dir[@]}"; do
          if [[ "$(ls -A $index | grep -v 'IndexStatus.log')" ]]; then
            echo "NGSmodule finished the job [Index]" >$index/IndexStatus.log
          elif [[ ! "$index" =~ *BismarkIndex* ]]; then
            rm -rf $index
          fi
        done
      fi
    fi

  done
done

#aws s3 --no-sign-request sync s3://ngi-igenomes/igenomes $iGenomes_dir

######## Start buiding #########
arr=($(find $iGenomes_dir -name "genome.fa" | grep "WholeGenomeFasta"))
total_task=${#arr[@]}

if [[ "$total_task" == 0 ]]; then
  color_echo "red" "ERROR! No sample sub-directory found in the iGenomes_dir:$iGenomes_dir\n"
  exit 1
fi

if [[ "$ntask_per_run" =~ ^[0-9]+$ ]]; then
  ntask_per_run=$ntask_per_run
elif [[ "$ntask_per_run" = "ALL" ]]; then
  if ((total_task > total_threads)); then
    ntask_per_run=$total_threads
  else
    ntask_per_run=$total_task
  fi
else
  color_echo "red" "ERROR! ntask_per_run should be 'ALL' or an interger!\n"
  exit 1
fi

threads=$((total_threads / ntask_per_run))

if ((threads == 0)); then
  threads=1
else
  threads=$threads
fi

if ((threads > 120)); then
  hisat2_threads=120
else
  hisat2_threads=$threads
fi

if ((((threads / ${#kmers})) == 0)); then
  map_threads=1
else
  map_threads=$((threads / ${#kmers}))
fi

###### fifo ######
fifo $ntask_per_run

echo -e "****************** Start buiding the remaining index ******************\n"
SECONDS=0

for genome in "${arr[@]}"; do
  read -u1000
  {
    SequenceDir=${genome%%/WholeGenomeFasta/genome.fa}
    genome_size=$(ls -lL $genome | awk '{print$5}')
    id=${SequenceDir##$iGenomes_dir/}
    id=${id%%Sequence}
    echo "+++++ ID: $id +++++"
    cd $SequenceDir

    arr1=($(find $SequenceDir -mindepth 1 -maxdepth 1 -name "*Index" -type d))
    if [[ "${#arr1[@]}" != 0 ]]; then
      for index in "${arr1[@]}"; do
        ln -fs $genome $index/genome.fa
      done
    fi

    arr2=($(find $SequenceDir -name "genome.fa" | grep -v "WholeGenomeFasta"))
    if [[ "${#arr2[@]}" != 0 ]]; then
      for genomeln in "${arr2[@]}"; do
        if [[ -L $genomeln ]]; then
          #color_echo "green" "Soft link ok: $genomeln"
          ln -fs $genome $genomeln
        else
          genomeln_size=$(ls -lL $genomeln | awk '{print$5}')
          if [[ "$genomeln_size" == "$genome_size" ]]; then
            #color_echo "green" "Soft link ok: $genomeln"
            ln -fs $genome $genomeln
          else
            color_echo "red" "Soft link fail: $genomeln genomeln_size:$genomeln_size genome_size:$genome_size"
          fi
        fi
      done
    fi

    gtf="$SequenceDir/../Annotation/Genes/genes.gtf"
    BWAIndex="$SequenceDir/BWAIndex"
    BowtieIndex="$SequenceDir/BowtieIndex"
    Bowtie2Index="$SequenceDir/Bowtie2Index"
    STARIndex="$SequenceDir/STARIndex"
    BismarkIndex="$SequenceDir/BismarkIndex"
    Hisat2Index="$SequenceDir/Hisat2Index"
    GemIndex="$SequenceDir/GemIndex/"
    GenmapIndex="$SequenceDir/GenmapIndex/"

    force="FALSE"
    status="uncompleted"
    attempt=0

    while [[ $status == "uncompleted" ]] && (("$attempt" <= 1)); do
      ((attempt++))
      if [[ $attempt != 1 ]]; then
        echo -e "+++++ ${id}: Number of attempts: $attempt +++++"
      fi

      ####### clear existed logs #######
      logfiles=("IndexStatus.log" "KmerStatus.log" "WindowStatus.log")
      globalcheck_logfile "$SequenceDir" logfiles[@] "$force" "$error_pattern" "$complete_pattern" "$id"
      find "$SequenceDir" -size 0 -print | while IFS= read -r pathname; do
        basedir=${pathname%/*}
        rm -rf $basedir
      done

      ####### genome.fa index #####
      check_logfile "$id" "genome_index" "${genome%%genome.fa}IndexStatus.log" "$error_pattern" "$complete_pattern" "precheck"
      if [[ $? == 1 ]]; then
        rm -f ${genome}.fai ${genome%%fa}dict
        samtools faidx ${genome} &>${genome%%genome.fa}IndexStatus.log
        picard CreateSequenceDictionary R=$genome &>${genome%%genome.fa}IndexStatus.log
        # faidx --regex "^(chr)*(([1-9][0-9]*)|([X,Y]))$" $genome >$SequenceDir/WholeGenomeFasta/genome_main.fa

        check_logfile "$id" "genome_index" "${genome%%genome.fa}IndexStatus.log" "$error_pattern" "$complete_pattern" "postcheck"
        if [[ $? == 1 ]]; then
          continue
        fi
      fi

      ###### BWA index #####
      check_logfile "$id" "BWA_index" "$BWAIndex/IndexStatus.log" "$error_pattern" "$complete_pattern" "precheck"
      if [[ $? == 1 ]]; then
        rm -rf $BWAIndex
        mkdir -p $BWAIndex
        ln -fs $genome $BWAIndex/genome.fa
        bwa index $BWAIndex/genome.fa &>$BWAIndex/IndexStatus.log

        check_logfile "$id" "BWA_index" "$BWAIndex/IndexStatus.log" "$error_pattern" "$complete_pattern" "postcheck"
        if [[ $? == 1 ]]; then
          continue
        fi
      fi

      ###### Bowtie index #####
      check_logfile "$id" "Bowtie_index" "$BowtieIndex/IndexStatus.log" "$error_pattern" "$complete_pattern" "precheck"
      if [[ $? == 1 ]]; then
        rm -rf $BowtieIndex
        mkdir -p $BowtieIndex
        ln -fs $genome $BowtieIndex/genome.fa
        bowtie-build --quiet --threads $threads $BowtieIndex/genome.fa $BowtieIndex/genome &>$BowtieIndex/IndexStatus.log

        check_logfile "$id" "Bowtie_index" "$BowtieIndex/IndexStatus.log" "$error_pattern" "$complete_pattern" "postcheck"
        if [[ $? == 1 ]]; then
          continue
        fi
      fi

      ###### Bowtie2 index #####
      check_logfile "$id" "Bowtie2_index" "$Bowtie2Index/IndexStatus.log" "$error_pattern" "$complete_pattern" "precheck"
      if [[ $? == 1 ]]; then
        rm -rf $Bowtie2Index
        mkdir -p $Bowtie2Index
        ln -fs $genome $Bowtie2Index/genome.fa
        bowtie2-build --quiet --threads $threads $Bowtie2Index/genome.fa $Bowtie2Index/genome &>$Bowtie2Index/IndexStatus.log

        check_logfile "$id" "Bowtie2_index" "$Bowtie2Index/IndexStatus.log" "$error_pattern" "$complete_pattern" "postcheck"
        if [[ $? == 1 ]]; then
          continue
        fi
      fi

      ###### Hisat2 index #####  Segmentation fault when thread is too large (>120?)
      check_logfile "$id" "Hisat2_index" "$Hisat2Index/IndexStatus.log" "$error_pattern" "$complete_pattern" "precheck"
      if [[ $? == 1 ]]; then
        rm -rf $Hisat2Index
        mkdir -p $Hisat2Index
        ln -fs $genome $Hisat2Index/genome.fa
        hisat2-build --quiet -p $hisat2_threads $Hisat2Index/genome.fa $Hisat2Index/genome &>$Hisat2Index/IndexStatus.log

        check_logfile "$id" "Hisat2_index" "$Hisat2Index/IndexStatus.log" "$error_pattern" "$complete_pattern" "postcheck"
        if [[ $? == 1 ]]; then
          continue
        fi
      fi

      ###### STAR index #####
      check_logfile "$id" "STAR_index" "$STARIndex/IndexStatus.log" "$error_pattern" "$complete_pattern" "precheck"
      if [[ $? == 1 ]]; then
        rm -rf $STARIndex
        mkdir -p $STARIndex
        ln -fs $genome $STARIndex/genome.fa
        if [[ -f $gtf ]]; then
          STAR --runMode genomeGenerate --runThreadN $threads \
          --genomeDir $STARIndex \
          --genomeFastaFiles $STARIndex/genome.fa \
          --sjdbGTFfile $gtf \
          --sjdbOverhang 100 &>$STARIndex/IndexStatus.log
        else
          STAR --runMode genomeGenerate --runThreadN $threads \
          --genomeDir $STARIndex \
          --genomeFastaFiles $STARIndex/genome.fa &>$STARIndex/IndexStatus.log
        fi

        check_logfile "$id" "STAR_index" "$STARIndex/IndexStatus.log" "$error_pattern" "$complete_pattern" "postcheck"
        if [[ $? == 1 ]]; then
          continue
        fi
      fi

      ###### rebuild bismark index from iGenome ######
      if [[ -f $BismarkIndex/genome.fa ]] && [[ -d $BismarkIndex/Bisulfite_Genome ]] && [[ -f $BismarkIndex/IndexStatus.log ]]; then
        rm -rf $BismarkIndex/bowtie2
        mkdir -p $BismarkIndex/bowtie2
        mv $BismarkIndex/genome.fa $BismarkIndex/bowtie2/
        mv $BismarkIndex/Bisulfite_Genome $BismarkIndex/bowtie2/
        mv $BismarkIndex/IndexStatus.log $BismarkIndex/bowtie2/
      else
        find "$BismarkIndex" -maxdepth 1 -type f -o -type l -print | xargs -i rm -f {}
      fi

      ###### bismark index #####
      check_logfile "$id" "bismarkBowtie2_index" "$BismarkIndex/bowtie2/IndexStatus.log" "$error_pattern" "$complete_pattern" "precheck"
      if [[ $? == 1 ]]; then
        rm -rf $BismarkIndex/bowtie2
        mkdir -p $BismarkIndex/bowtie2
        ln -fs $genome $BismarkIndex/bowtie2/genome.fa
        bismark_genome_preparation --genomic_composition --bowtie2 --parallel $threads $BismarkIndex/bowtie2 &>$BismarkIndex/bowtie2/IndexStatus.log

        check_logfile "$id" "bismarkBowtie2_index" "$BismarkIndex/bowtie2/IndexStatus.log" "$error_pattern" "$complete_pattern" "postcheck"
        if [[ $? == 1 ]]; then
          continue
        fi
      fi

      check_logfile "$id" "bismarkHisat2_index" "$BismarkIndex/hisat2/IndexStatus.log" "$error_pattern" "$complete_pattern" "precheck"
      if [[ $? == 1 ]]; then
        rm -rf $BismarkIndex/hisat2
        mkdir -p $BismarkIndex/hisat2
        ln -fs $genome $BismarkIndex/hisat2/genome.fa
        bismark_genome_preparation --genomic_composition --hisat2 --parallel $threads $BismarkIndex/hisat2 &>$BismarkIndex/hisat2/IndexStatus.log

        check_logfile "$id" "bismarkHisat2_index" "$BismarkIndex/hisat2/IndexStatus.log" "$error_pattern" "$complete_pattern" "postcheck"
        if [[ $? == 1 ]]; then
          continue
        fi
      fi

      #### Gem #####
      check_logfile "$id" "Gem_index" "$GemIndex/IndexStatus.log" "$error_pattern" "$complete_pattern" "precheck"
      if [[ $? == 1 ]]; then
        rm -rf $GemIndex
        mkdir -p $GemIndex
        gem-indexer --threads $threads -i $genome -o $GemIndex/genome &>$GemIndex/IndexStatus.log

        check_logfile "$id" "Gem_index" "$GemIndex/IndexStatus.log" "$error_pattern" "$complete_pattern" "postcheck"
        if [[ $? == 1 ]]; then
          continue
        fi
      fi

      for kmer in "${kmers[@]}"; do
        {
          check_logfile "$id" "Gem_Kmer_$kmer" "$GemIndex/Mappability/${kmer}mer/KmerStatus.log" "$error_pattern" "$complete_pattern" "precheck"
          if [[ $? == 1 ]]; then
            rm -rf $GemIndex/Mappability/${kmer}mer
            mkdir -p $GemIndex/Mappability/${kmer}mer
            cd $GemIndex/Mappability/${kmer}mer
            gem-mappability -T $map_threads -I $GemIndex/genome.gem -l ${kmer} -o genome.${kmer}mer.gem &>KmerStatus.log
            gem-2-wig -I $GemIndex/genome.gem -i genome.${kmer}mer.gem.mappability -o genome.${kmer}mer.gem &>>KmerStatus.log
            wigToBigWig genome.${kmer}mer.gem.wig genome.${kmer}mer.gem.sizes genome.${kmer}mer.gem.bigwig &>>KmerStatus.log
            rm -f genome.${kmer}mer.gem.mappability

            check_logfile "$id" "Gem_Kmer_$kmer" "$GemIndex/Mappability/${kmer}mer/KmerStatus.log" "$error_pattern" "$complete_pattern" "postcheck"
            if [[ $? == 1 ]]; then
              continue 2
            fi
          fi

          for window in "${windows[@]}"; do
            check_logfile "$id" "Gem_Kmer_${kmer}_Windows_${window}" "$GemIndex/Mappability/${kmer}mer/windows/$window/WindowStatus.log" "$error_pattern" "$complete_pattern" "precheck"
            if [[ $? == 1 ]]; then
              rm -rf $GemIndex/Mappability/${kmer}mer/windows/win$window
              mkdir -p $GemIndex/Mappability/${kmer}mer/windows/win$window
              cd $GemIndex/Mappability/${kmer}mer/windows/win$window
              gcCounter -w $window --forgiving $genome >genome.win${window}.gc.wig &>WindowStatus.log
              mapCounter -w $window $GemIndex/Mappability/${kmer}mer/genome.${kmer}mer.gem.bigwig >genome.win${window}.${kmer}mer.gem.wig &>>WindowStatus.log

              check_logfile "$id" "Gem_Kmer_${kmer}_Windows_${window}" "$GemIndex/windows/$window/WindowStatus.log" "$error_pattern" "$complete_pattern" "postcheck"
              if [[ $? == 1 ]]; then
                continue 3
              fi
            fi
          done

        } &
      done
      wait
      echo -e "\033[32mComplete Gem mappability building.\033[0m"

      # ##### Genmap #####
      # if [[ ! -d $GenmapIndex ]] || [[ ! "$(ls -A $GenmapIndex)" ]]; then
      # echo -e "\033[35mStart to build Genmap index...\033[0m"
      # rm -rf $GenmapIndex
      # genmap index -F $genome -I $GenmapIndex
      # echo -e "\033[32mComplete Genmap index building.\033[0m"

      # for kmer in "${kmers[@]}"; do
      #   echo "====== Make Genmap mappability file with kmer:$kmer  ======"
      #   mkdir -p $GenmapIndex/Mappability/${kmer}mer
      #   cd $GenmapIndex/Mappability/${kmer}mer
      #   genmap map --index $GenmapIndex --errors 2 --length ${kmer} --threads $threads --wig --output genome.${kmer}mer.genmap
      #   wigToBigWig genome.${kmer}mer.genmap.wig genome.${kmer}mer.genmap.chrom.sizes genome.${kmer}mer.genmap.bigwig

      #   for window in "${windows[@]}"; do
      #     echo "====== Count GC and mappability within a silding window:$window  ======"
      #     mkdir -p $GenmapIndex/windows/$window
      #     cd $GenmapIndex/windows/$window
      #     gcCounter -w $window --forgiving $genome >genome.w${window}.gc.wig
      #     mapCounter -w $window $GenmapIndex/Mappability/${kmer}mer/genome.${kmer}mer.genmap.bigwig >genome.w${window}.${kmer}mer.genmap.wig
      #   done
      # done
      # echo -e "\033[32mComplete Genmap mappability building.\033[0m"
      # fi

      status="completed"
      color_echo "blue" "+++++ ${id}: completed +++++"
    done

    if [[ "$status" == "completed" ]]; then
      echo "Completed: $id" >>"$tmpfile"
    else
      echo "Interrupted: $id" >>"$tmpfile"
      color_echo "red" "ERROR! ${id} interrupted! Please check the processing log."
    fi

    color_echo "green" "***** Completed:$(cat "$tmpfile" | grep "Completed" | uniq | wc -l) | Interrupted:$(cat "$tmpfile" | grep "Interrupted" | uniq | wc -l) | Total:$total_task *****"

    echo >&1000
  } &
done
wait

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "****************** Alignment Done ******************\n"
