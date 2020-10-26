#!/usr/bin/env bash
trap "trap - SIGTERM && kill -- -$$" SIGINT SIGTERM

############# Paramaters #############################################################
iGenomes_dir="/data/database/iGenomes/"
Species=("Homo_sapiens" "Mus_musculus" "Macaca_mulatta")
kmers=(150 100 50)
windows=(2000000 1000000 500000 100000)
threads=100

######## check command available ##########
faidx --version &>/dev/null
[ $? -ne 0 ] && {
  echo -e "Cannot find the tool pyfaidx. User can install it with the command 'conda install -c bioconda pyfaidx'.\n"
  exit 1
}
gem-mappability --help &>/dev/null
[ $? -ne 0 ] && {
  echo -e "Cannot find the command gem-mappability. User can install it from 'https://sourceforge.net/projects/gemlibrary'.\n"
  exit 1
}
genmap --version &>/dev/null
[ $? -ne 0 ] && {
  echo -e "Cannot find the command genmap. User can install it with the command 'conda install -c bioconda genmap'.\n"
  exit 1
}
wigToBigWig &>/dev/null
[ $? -ne 255 ] && {
  echo -e "Cannot find the command wigToBigWig. User can install it with the command 'conda install -c bioconda ucsc-wigtobigwig'.\n"
  exit 1
}
mapCounter --help &>/dev/null
[ $? -ne 0 ] && {
  echo -e "Cannot find the command mapCounter. User can install it from 'https://github.com/shahcompbio/hmmcopy_utils'.\n"
  exit 1
}
picard &>/dev/null
[ $? -eq 127 ] && {
  echo -e "Cannot find the command picard.\n"
  exit 1
}

######## Download the iGenomes #####
igenomes_file_manifest=($(curl "https://raw.githubusercontent.com/ewels/AWS-iGenomes/master/ngi-igenomes_file_manifest.txt" |cat))
for s in "${Species[@]}";do
  echo "Downloading the iGenomes for Species: $s"
  for file in "${igenomes_file_manifest[@]}";do
    igenome=${file##s3://ngi-igenomes/igenomes/}
    igenome=${igenome%%/*}
    if [[ $igenome == $s  ]]; then
      aws s3 --no-sign-request sync $file $iGenomes_dir
    fi
  done
done

#aws s3 --no-sign-request sync s3://ngi-igenomes/igenomes $iGenomes_dir

######## Start buiding #########
arr=($(find $iGenomes_dir -name "genome.fa" | grep "WholeGenomeFasta"))
if ((threads > 120)); then
  threads=120
else
  threads=$threads
fi
for genome in "${arr[@]}"; do

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

  if [[ ! -f ${genome%%.fa}.dict ]];then
    picard CreateSequenceDictionary R=$genome
  fi

  ####### Extract main genome fasta #####
  # echo "====== Fetch main chromosome sequence into genome_main.fa ======"
  # faidx --regex "^(chr)*(([1-9][0-9]*)|([X,Y]))$" $genome >$SequenceDir/WholeGenomeFasta/genome_main.fa

  # ###### BWA index #####
  # echo -e "\033[35mStart to build BWA index...\033[0m"
  # mkdir -p $BWAIndex
  # ln -fs $genome $BWAIndex/genome.fa
  # cd $BWAIndex
  # bwa index $BWAIndex/genome.fa
  # echo -e "\033[32mComplete BWA index building.\033[0m"

  # ###### Bowtie index #####
  # echo -e "\033[35mStart to build Bowtie index...\033[0m"
  # mkdir -p $BowtieIndex
  # ln -fs $genome $BowtieIndex/genome.fa
  # bowtie-build --quiet --threads $threads $BowtieIndex/genome.fa $BowtieIndex/genome
  # echo -e "\033[32mComplete Bowtie index building.\033[0m"

  # ###### Bowtie2 index #####
  # echo -e "\033[35mStart to build Bowtie2 index...\033[0m"
  # mkdir -p $Bowtie2Index
  # ln -fs $genome $Bowtie2Index/genome.fa
  # bowtie2-build --quiet --threads $threads $Bowtie2Index/genome.fa $Bowtie2Index/genome
  # echo -e "\033[32mComplete Bowtie2 index building.\033[0m"

  ###### Hisat2 index #####  Segmentation fault when thread is too large (>120?)
  echo -e "\033[35mStart to build Hisat2 index...\033[0m"
  mkdir -p $Hisat2Index
  ln -fs $genome $Hisat2Index/genome.fa
  hisat2-build --quiet -p $threads $Hisat2Index/genome.fa $Hisat2Index/genome
  echo -e "\033[32mComplete Hisat2 index building.\033[0m"

  # ###### STAR index #####
  # echo -e "\033[35mStart to build STAR index...\033[0m"
  # mkdir -p $STARIndex
  # ln -fs $genome $STARIndex/genome.fa
  # STAR --runMode genomeGenerate --runThreadN $threads \
  # --genomeDir $STARIndex \
  # --genomeFastaFiles $STARIndex/genome.fa \
  # --sjdbGTFfile $gtf
  # echo -e "\033[32mComplete STAR index building.\033[0m"

  # ###### bismark index #####
  # echo -e "\033[35mStart to build bismark index...\033[0m"
  # mkdir -p $BismarkIndex/bowtie2
  # mkdir -p $BismarkIndex/hisat2
  # ln -fs $genome $BismarkIndex/bowtie2/genome.fa
  # ln -fs $genome $BismarkIndex/hisat2/genome.fa
  # bismark_genome_preparation --genomic_composition --bowtie2 --parallel $threads $BismarkIndex/bowtie2
  # bismark_genome_preparation --genomic_composition --hisat2 --parallel $threads $BismarkIndex/hisat2
  # echo -e "\033[32mComplete bismark index building.\033[0m"

  ###### rebuild bismark_hisat2 index ######
  echo -e "\033[35mBuild bismark_hisat2 index...\033[0m"
  mkdir -p $BismarkIndex/bowtie2 $BismarkIndex/hisat2
  if [[ -f $BismarkIndex/genome.fa ]] && [[ -d $BismarkIndex/Bisulfite_Genome ]]; then
    mv $BismarkIndex/genome.fa $BismarkIndex/bowtie2/
    mv $BismarkIndex/Bisulfite_Genome $BismarkIndex/bowtie2/
  fi
  ln -fs $genome $BismarkIndex/hisat2/genome.fa
  bismark_genome_preparation --genomic_composition --hisat2 --parallel $threads $BismarkIndex/hisat2
  echo -e "\033[32mComplete bismark_hisat2 index building.\033[0m"

  #### Gem #####
  echo -e "\033[35mStart to build Gem index...\033[0m"
  if [[ ! -d $GemIndex ]]; then
    mkdir -p $GemIndex
    gem-indexer --threads $threads -i $genome -o $GemIndex/genome
  fi
  echo -e "\033[32mComplete Gem index building.\033[0m"

  for kmer in "${kmers[@]}"; do
    echo "====== Make Gem mappability file  ======"
    mkdir -p $GemIndex/Mappability/${kmer}mer
    cd $GemIndex/Mappability/${kmer}mer
    gem-mappability -T $threads -I $GemIndex/genome.gem -l ${kmer} -o genome.${kmer}mer.gem
    gem-2-wig -I $GemIndex/genome.gem -i genome.${kmer}mer.gem.mappability -o genome.${kmer}mer.gem
    wigToBigWig genome.${kmer}mer.gem.wig genome.${kmer}mer.gem.sizes genome.${kmer}mer.gem.bigwig
    rm -f genome.${kmer}mer.gem.mappability

    echo "====== Count GC and mappability within a silding window  ======"
    for window in "${windows[@]}"; do
      mkdir -p $GemIndex/windows/$window
      cd $GemIndex/windows/$window
      gcCounter -w $window --forgiving $genome >genome.w${window}.gc.wig
      mapCounter -w $window $GemIndex/Mappability/${kmer}mer/genome.${kmer}mer.gem.bigwig >genome.w${window}.${kmer}mer.gem.wig
    done
  done

  # ##### Genmap #####
  # echo -e "\033[35mStart to build Genmap index...\033[0m"
  # rm -rf $GenmapIndex
  # genmap index -F $genome -I $GenmapIndex
  # echo -e "\033[32mComplete Genmap index building.\033[0m"

  # for kmer in "${kmers[@]}"; do
  #   echo "====== Make Genmap mappability file  ======"
  #   mkdir -p $GenmapIndex/Mappability/${kmer}mer
  #   cd $GenmapIndex/Mappability/${kmer}mer
  #   genmap map --index $GenmapIndex --errors 2 --length ${kmer} --threads $threads --wig --output genome.${kmer}mer.genmap
  #   wigToBigWig genome.${kmer}mer.genmap.wig genome.${kmer}mer.genmap.chrom.sizes genome.${kmer}mer.genmap.bigwig

  #   echo "====== Count GC and mappability within a silding window  ======"
  #   for window in "${windows[@]}"; do
  #     mkdir -p $GenmapIndex/windows/$window
  #     cd $GenmapIndex/windows/$window
  #     gcCounter -w $window --forgiving $genome >genome.w${window}.gc.wig
  #     mapCounter -w $window $GenmapIndex/Mappability/${kmer}mer/genome.${kmer}mer.genmap.bigwig >genome.w${window}.${kmer}mer.genmap.wig
  #   done
  # done

done

echo "Done"
