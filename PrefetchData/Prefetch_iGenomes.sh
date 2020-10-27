#!/usr/bin/env bash
trap "trap - SIGTERM && kill -- -$$" SIGINT SIGTERM

############# Paramaters #############################################################
iGenomes_dir="/reference/iGenomes/"
Species=("Homo_sapiens" "Mus_musculus" "Macaca_mulatta")
Sources=("Ensembl" "NCBI" "UCSC")
kmers=(150 100 50)
windows=(2000000 1000000 500000 100000)
total_threads=384
ntask_per_run="ALL"

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

######## Download the iGenomes #####
for s in "${Species[@]}"; do
  for i in "${Sources[@]}"; do
    echo -e "\033[32mDownloading the iGenomes: $s/$i\033[0m"
    igenome="s3://ngi-igenomes/igenomes/$s/$i"
    aws s3 --no-sign-request sync $igenome $iGenomes_dir/$s/$i --exclude "*/genome.fa" --include "WholeGenomeFasta/genome.fa"
    if [[ ! "$(ls -A $iGenomes_dir/$s/$i)" ]]; then
      echo -e "\033[33miGenomes do not exist: $s/$i\033[0m"
      rm -rf $iGenomes_dir/$s/$i
    fi
    index_dir=($(find $iGenomes_dir -name "*Index" -type d))
    for index in "${index_dir[@]}"; do
      echo -e "NGSmodule finished the job [Index]" >$index/IndexStatus.log
    done
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

###### fifo ######
tempfifo=$$.fifo
trap_add "exec 1000>&-;exec 1000<&-;rm -f $tempfifo" SIGINT SIGTERM EXIT
mkfifo $tempfifo
exec 1000<>$tempfifo
rm -f $tempfifo
for ((i = 1; i <= ntask_per_run; i++)); do
  echo >&1000
done

threads=$(((total_threads + ntask_per_run) / ntask_per_run - 1))

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

for genome in "${arr[@]}"; do
  read -u1000
  {
    SequenceDir=${genome%%WholeGenomeFasta/genome.fa}
    genome_size=$(ls -lL $genome | awk '{print$5}')
    echo "SequenceDir:$SequenceDir"
    cd $SequenceDir

    ####### Create WholeGenomeFasta/genome.fa softlink for other index #######
    arr2=($(find $SequenceDir -name "genome.fa" | grep -v "WholeGenomeFasta"))
    for genomeln in "${arr2[@]}"; do
      if [[ -L $genomeln ]]; then
        echo -e "\033[32mok: $genomeln\033[0m"
        ln -fs $genome $genomeln
      else
        genomeln_size=$(ls -lL $genomeln | awk '{print$5}')
        if [[ "$genomeln_size" == "$genome_size" ]]; then
          echo -e "\033[32mSoft link ok: $genomeln\033[0m"
          ln -fs $genome $genomeln
        else
          echo -e "\033[31mSoft link fail: $genomeln genomeln_size:$genomeln_size genome_size:$genome_size\033[0m"
        fi
      fi
    done

    gtf="$SequenceDir/../Annotation/Genes/genes.gtf"
    BWAIndex="$SequenceDir/BWAIndex"
    BowtieIndex="$SequenceDir/BowtieIndex"
    Bowtie2Index="$SequenceDir/Bowtie2Index"
    STARIndex="$SequenceDir/STARIndex"
    BismarkIndex="$SequenceDir/BismarkIndex"
    Hisat2Index="$SequenceDir/Hisat2Index"
    GemIndex="$SequenceDir/GemIndex/"
    GenmapIndex="$SequenceDir/GenmapIndex/"

    ####### clear existed logs #######
    logfiles=("IndexStatus.log" "KmerStatus.log" "GCStatus.log" "MappabilityStatus.log")
    globalcheck_logfile "$dir" logfiles[@] "$force" "$error_pattern" "$complete_pattern" "$SequenceDir"

    ####### genome.fa index #####
    rm -f ${genome}.fai ${genome%%fa}dict
    samtools faidx ${genome} &>${genome%%genome.fa}IndexStatus.log
    picard CreateSequenceDictionary R=$genome &>${genome%%genome.fa}IndexStatus.log
    # faidx --regex "^(chr)*(([1-9][0-9]*)|([X,Y]))$" $genome >$SequenceDir/WholeGenomeFasta/genome_main.fa

    ###### BWA index #####
    if [[ ! -d $BWAIndex ]] || [[ ! "$(ls -A $BWAIndex)" ]]; then
      echo -e "\033[35mStart to build BWA index...\033[0m"
      mkdir -p $BWAIndex
      ln -fs $genome $BWAIndex/genome.fa
      bwa index $BWAIndex/genome.fa &>$BWAIndex/IndexStatus.log
      echo -e "\033[32mComplete BWA index building.\033[0m"
    fi

    ###### Bowtie index #####
    if [[ ! -d $BowtieIndex ]] || [[ ! "$(ls -A $BowtieIndex)" ]]; then
      echo -e "\033[35mStart to build Bowtie index...\033[0m"
      mkdir -p $BowtieIndex
      ln -fs $genome $BowtieIndex/genome.fa
      bowtie-build --quiet --threads $threads $BowtieIndex/genome.fa $BowtieIndex/genome &>$BowtieIndex/IndexStatus.log
      echo -e "\033[32mComplete Bowtie index building.\033[0m"
    fi

    ###### Bowtie2 index #####
    if [[ ! -d $Bowtie2Index ]] || [[ ! "$(ls -A $Bowtie2Index)" ]]; then
      echo -e "\033[35mStart to build Bowtie2 index...\033[0m"
      mkdir -p $Bowtie2Index
      ln -fs $genome $Bowtie2Index/genome.fa
      bowtie2-build --quiet --threads $threads $Bowtie2Index/genome.fa $Bowtie2Index/genome &>$Bowtie2Index/IndexStatus.log
      echo -e "\033[32mComplete Bowtie2 index building.\033[0m"
    fi

    ###### Hisat2 index #####  Segmentation fault when thread is too large (>120?)
    if [[ ! -d $Hisat2Index ]] || [[ ! "$(ls -A $Hisat2Index)" ]]; then
      echo -e "\033[35mStart to build Hisat2 index...\033[0m"
      mkdir -p $Hisat2Index
      ln -fs $genome $Hisat2Index/genome.fa
      hisat2-build --quiet -p $hisat2_threads $Hisat2Index/genome.fa $Hisat2Index/genome &>$Hisat2Index/IndexStatus.log
      echo -e "\033[32mComplete Hisat2 index building.\033[0m"
    fi

    ###### STAR index #####
    if [[ ! -d $STARIndex ]] || [[ ! "$(ls -A $STARIndex)" ]]; then
      echo -e "\033[35mStart to build STAR index...\033[0m"
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
      echo -e "\033[32mComplete STAR index building.\033[0m"
    fi

    ###### rebuild bismark index from iGenome ######
    if [[ -f $BismarkIndex/genome.fa ]] && [[ -d $BismarkIndex/Bisulfite_Genome ]]; then
      echo -e "\033[35mBuild bismark_hisat2 index...\033[0m"
      mkdir -p $BismarkIndex/bowtie2
      mv $BismarkIndex/genome.fa $BismarkIndex/bowtie2/
      mv $BismarkIndex/Bisulfite_Genome $BismarkIndex/bowtie2/
      mkdir -p $BismarkIndex/hisat2
      ln -fs $genome $BismarkIndex/hisat2/genome.fa
      bismark_genome_preparation --genomic_composition --hisat2 --parallel $hisat2_threads $BismarkIndex/hisat2 &>$BismarkIndex/hisat2/IndexStatus.log
      echo -e "\033[32mComplete bismark_hisat2 index building.\033[0m"
    fi

    ###### bismark index #####
    if [[ ! -d $BismarkIndex/bowtie2 ]] || [[ ! "$(ls -A $BismarkIndex/bowtie2)" ]]; then
      echo -e "\033[35mStart to build bismark_bowtie2 index...\033[0m"
      mkdir -p $BismarkIndex/bowtie2
      ln -fs $genome $BismarkIndex/bowtie2/genome.fa
      bismark_genome_preparation --genomic_composition --bowtie2 --parallel $threads $BismarkIndex/bowtie2 &>$BismarkIndex/bowtie2/IndexStatus.log
      echo -e "\033[32mComplete bismark_bowtie2 index building.\033[0m"
    fi

    if [[ ! -d $BismarkIndex/hisat2 ]] || [[ ! "$(ls -A $BismarkIndex/hisat2)" ]]; then
      echo -e "\033[35mStart to build bismark_hisat2 index...\033[0m"
      mkdir -p $BismarkIndex/hisat2
      ln -fs $genome $BismarkIndex/hisat2/genome.fa
      bismark_genome_preparation --genomic_composition --hisat2 --parallel $threads $BismarkIndex/hisat2 &>$BismarkIndex/hisat2/IndexStatus.log
      echo -e "\033[32mComplete bismark_hisat2 index building.\033[0m"
    fi

    #### Gem #####
    if [[ ! -d $GemIndex ]] || [[ ! "$(ls -A $GemIndex)" ]]; then
      echo -e "\033[35mStart to build Gem index...\033[0m"
      mkdir -p $GemIndex
      gem-indexer --threads $threads -i $genome -o $GemIndex/genome &>$GemIndex/IndexStatus.log
      echo -e "\033[32mComplete Gem index building.\033[0m"

      for kmer in "${kmers[@]}"; do
        {
          echo "====== Make Gem mappability file with kmer:$kmer  ======"
          mkdir -p $GemIndex/Mappability/${kmer}mer
          cd $GemIndex/Mappability/${kmer}mer
          gem-mappability -T $map_threads -I $GemIndex/genome.gem -l ${kmer} -o genome.${kmer}mer.gem &>KmerStatus.log
          gem-2-wig -I $GemIndex/genome.gem -i genome.${kmer}mer.gem.mappability -o genome.${kmer}mer.gem &>>KmerStatus.log
          wigToBigWig genome.${kmer}mer.gem.wig genome.${kmer}mer.gem.sizes genome.${kmer}mer.gem.bigwig &>>KmerStatus.log
          rm -f genome.${kmer}mer.gem.mappability

          for window in "${windows[@]}"; do
            echo "====== Count GC and mappability within a silding window:$window  ======"
            mkdir -p $GemIndex/windows/$window
            cd $GemIndex/windows/$window
            gcCounter -w $window --forgiving $genome >genome.w${window}.gc.wig &>GCStatus.log
            mapCounter -w $window $GemIndex/Mappability/${kmer}mer/genome.${kmer}mer.gem.bigwig >genome.w${window}.${kmer}mer.gem.wig &>MappabilityStatus.log
          done
        } &
      done
      wait
      echo -e "\033[32mComplete Gem mappability building.\033[0m"
    fi

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

    echo >&1000
  } &
done
wait
echo "Done"
