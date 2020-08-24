#!/usr/bin/env bash
trap "trap - SIGTERM && kill -- -$$" SIGINT SIGTERM

############# Paramaters #############################################################
rootdir="/data/database/iGenomes/"
threads=120
kmers=(150 140 100 90 50 40)
windows=(2000000 1000000 500000 100000)

######## check command available ##########
faidx --version &>/dev/null
[ $? -ne 0 ] && {
  color_echo "red" "Cannot find the tool pyfaidx. User can install it with the command 'conda install -c bioconda pyfaidx'.\n"
  exit 1
}
gem-mappability --help &>/dev/null
[ $? -ne 0 ] && {
  color_echo "red" "Cannot find the command gem-mappability. User can install it from 'https://sourceforge.net/projects/gemlibrary'.\n"
  exit 1
}
genmap --version &>/dev/null
[ $? -ne 0 ] && {
  color_echo "red" "Cannot find the command genmap. User can install it with the command 'conda install -c bioconda genmap'.\n"
  exit 1
}
wigToBigWig --version &>/dev/null
[ $? -ne 0 ] && {
  color_echo "red" "Cannot find the command wigToBigWig. User can install it with the command 'conda install -c bioconda ucsc-wigtobigwig'.\n"
  exit 1
}
mapCounter --help &>/dev/null
[ $? -ne 0 ] && {
  color_echo "red" "Cannot find the command mapCounter. User can install it from 'https://github.com/shahcompbio/hmmcopy_utils'.\n"
  exit 1
}


######## Start buiding #########
arr=($(find $rootdir -name "WholeGenomeFasta"))
# if ((threads > 120)); then
#   threads=120
# else
#   threads=$threads
# fi
for maindir in "${arr[@]}"; do
  {
    maindir=${maindir%%WholeGenomeFasta}
    echo "maindir:$maindir"

    genome="$maindir/WholeGenomeFasta/genome.fa"
    genome_main="$maindir/WholeGenomeFasta/genome_main.fa"

    BWAIndex="$maindir/BWAIndex"
    BowtieIndex="$maindir/BowtieIndex"
    Bowtie2Index="$maindir/Bowtie2Index"
    Hisat2Index="$maindir/Hisat2Index"
    STARIndex="$maindir/STARIndex"
    BismarkIndex="$maindir/BismarkIndex"
    gtf="$maindir/../Annotation/Genes/genes.gtf"

    GemIndex="$maindir/GemIndex/"
    GenmapIndex="$maindir/GenmapIndex/"

    ####### Extract main genome fasta #####
    if [[ ! -f $genome ]]; then
      echo "File:$genome do not exist. Please check the destination."
      break
    fi
    echo "====== Fetch main chromosome sequence into genome_main.fa ======"
    faidx --regex "^(chr)*(([1-9][0-9]*)|([X,Y]))$" $genome >$genome_main

    # ###### BWA index #####
    # mkdir -p $BWAIndex
    # ln -fs $genome $BWAIndex/genome.fa
    # ln -fs $genome_main $BWAIndex/genome_main.fa
    # cd $BWAIndex
    # bwa index $BWAIndex/genome.fa &
    # bwa index $BWAIndex/genome_main.fa &

    # ###### Bowtie index #####
    # mkdir -p $BowtieIndex
    # ln -fs $genome $BowtieIndex/genome.fa
    # ln -fs $genome_main $BowtieIndex/genome_main.fa
    # eval "$(conda shell.bash hook)"
    # conda activate py2
    # bowtie-build --threads $threads $BowtieIndex/genome.fa $BowtieIndex/genome
    # bowtie-build --threads $threads $BowtieIndex/genome_main.fa $BowtieIndex/genome_main
    # conda deactivate

    # ###### Bowtie2 index #####
    # mkdir -p $Bowtie2Index
    # ln -fs $genome $Bowtie2Index/genome.fa
    # ln -fs $genome_main $Bowtie2Index/genome_main.fa
    # bowtie2-build --threads $threads $Bowtie2Index/genome.fa $Bowtie2Index/genome
    # bowtie2-build --threads $threads $Bowtie2Index/genome_main.fa $Bowtie2Index/genome_main

    # ###### Hisat2 index #####  Segmentation fault when thread is too large (>120?)
    # mkdir -p $Hisat2Index
    # ln -fs $genome $Hisat2Index/genome.fa
    # ln -fs $genome_main $Hisat2Index/genome_main.fa
    # hisat2-build -p $threads $Hisat2Index/genome.fa $Hisat2Index/genome
    # hisat2-build -p $threads $Hisat2Index/genome_main.fa $Hisat2Index/genome_main
    # wait

    # ###### STAR index #####
    # mkdir -p $STARIndex/genome
    # mkdir -p $STARIndex/genome_main
    # ln -fs $genome $STARIndex/genome/genome.fa
    # ln -fs $genome $STARIndex/genome_main/genome_main.fa

    # STAR --runMode genomeGenerate --runThreadN $threads \
    # --genomeDir $STARIndex/genome \
    # --genomeFastaFiles $STARIndex/genome/genome.fa \
    # --sjdbGTFfile $gtf
    # STAR --runMode genomeGenerate --runThreadN $threads \
    # --genomeDir $STARIndex/genome_main \
    # --genomeFastaFiles $STARIndex/genome_main/genome_main.fa \
    # --sjdbGTFfile $gtf

    # ###### bismark index #####
    # mkdir -p $BismarkIndex/genome/bowtie2
    # mkdir -p $BismarkIndex/genome/hisat2
    # mkdir -p $BismarkIndex/genome_main/bowtie2
    # mkdir -p $BismarkIndex/genome_main/hisat2
    # ln -fs $genome $BismarkIndex/genome/bowtie2/genome.fa
    # ln -fs $genome $BismarkIndex/genome/hisat2/genome.fa
    # ln -fs $genome_main $BismarkIndex/genome_main/bowtie2/genome_main.fa
    # ln -fs $genome_main $BismarkIndex/genome_main/hisat2/genome_main.fa
    # bismark_genome_preparation --genomic_composition --bowtie2 --parallel $((($threads) / 2)) $BismarkIndex/genome/bowtie2
    # bismark_genome_preparation --genomic_composition --hisat2 --parallel $threads $BismarkIndex/genome/hisat2
    # bismark_genome_preparation --genomic_composition --bowtie2 --parallel $((($threads) / 2)) $BismarkIndex/genome_main/bowtie2
    # bismark_genome_preparation --genomic_composition --hisat2 --parallel $threads $BismarkIndex/genome_main/hisat2

    ##### Gem #####
    if [[ ! -d $GemIndex ]]; then
      mkdir $GemIndex
      echo "====== Make GemIndex  ======"
      gem-indexer --threads ${threads} -i ${genome_main} -o ${GemIndex}/genome_main
    fi

    for kmer in "${kmers[@]}"; do

      echo "====== Make gem mappability file  ======"
      mkdir -p $GemIndex/Mappability/${kmer}mer
      cd $GemIndex/Mappability/${kmer}mer
      gem-mappability -T ${threads} -I ${GemIndex}/genome_main.gem -l ${kmer} -o genome_main.gem
      gem-2-wig -I ${GemIndex}/genome_main.gem -i genome_main.gem.mappability -o genome_main.gem
      wigToBigWig genome_main.gem.wig genome_main.gem.sizes genome_main.gem.bigwig
      rm -f genome_main.gem.mappability

      echo "====== Count GC and mappability within a silding window  ======"
      for window in "${windows[@]}"; do
        mkdir -p $GemIndex/windows/$window
        cd $GemIndex/windows/$window
        gcCounter -w $window --forgiving $genome_main >genome_main.w$window.gc.wig
        mapCounter -w $window $GemIndex/Mappability/${kmer}mer/genome_main.gem.bigwig >genome_main.w$window.${kmer}mer.gem.wig
      done

    done

    ###### Genmap #####
    if [[ -d $GenmapIndex ]]; then
      rm -rf $GenmapIndex
    fi
    echo "====== Make GenmapIndex  ======"
    genmap index -F $genome_main -I $GenmapIndex

    for kmer in "${kmers[@]}"; do

      echo "====== Make genmap mappability file  ======"
      mkdir -p $GenmapIndex/Mappability/${kmer}mer
      cd $GenmapIndex/Mappability/${kmer}mer
      genmap map -I $GenmapIndex -O ./ -E 2 -K $kmer -w
      wigToBigWig genome_main.genmap.wig genome_main.genmap.chrom.sizes genome_main.genmap.bigwig

      echo "====== Count GC and mappability within a silding window  ======"
      for window in "${windows[@]}"; do
        mkdir -p $GenmapIndex/windows/$window
        cd $GenmapIndex/windows/$window
        gcCounter -w $window --forgiving $genome_main >genome_main.w$window.gc.wig
        mapCounter -w $window $GenmapIndex/Mappability/${kmer}mer/genome_main.genmap.bigwig >genome_main.w$window.${kmer}mer.genmap.wig
      done

    done

  } &
done

wait

echo "Done"
