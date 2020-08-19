#!/usr/bin/bash
maindir=($(pwd))
data_dir=$maindir/rawdata/
SampleGrepPattern=""
total_threads=380
ntask_per_run="ALL"

Aligner="bwa"
#gsnap_index="/data/database/iGenomes/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/gmap/genome_rmchrM.fa"
#gem_index="/data/database/iGenomes/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome_rmchrM.fa.gem"
#baseqCNV_config="/home/reprod/bin/singlecellseqcnv-code/baseqCNV/baseqCNV.ini"
#FREEC_config="/home/reprod/bin/FREEC/config_FREEC.hg19.ini"

#######################################################################################
arr=($(find $data_dir -mindepth 1 -maxdepth 1 \( -type l -o -type d \) -printf '%P\n' | grep -P "$SampleGrepPattern"))
total_task=${#arr[@]}
if [[ "$ntask_per_run" =~ ^[0-9]+$ ]]; then
  ntask_per_run=$ntask_per_run
elif [[ "$ntask_per_run" = "ALL" ]]; then
  ntask_per_run=$total_task
else
  echo "ntask_per_run should be 'ALL' or an interger "
  exit 1
fi
threads=$((($total_threads + $ntask_per_run - 1) / $ntask_per_run))
echo "threads=$threads"

###### fifo ######
tempfifo=$$.fifo
trap "exec 1000>&-;exec 1000<&-;exit 0" 2
mkfifo $tempfifo
exec 1000<>$tempfifo
rm -rf $tempfifo
for ((i = 1; i <= $ntask_per_run; i++)); do
  echo >&1000
done

###### processbar <current> <total> ######
processbar() {
  local current=$1
  local total=$2
  local maxlen=80
  local barlen=66
  local perclen=14
  local format="%-${barlen}s%$((maxlen - barlen))s"
  local perc="[$current/$total]"
  local progress=$((current * barlen / total))
  local prog=$(for i in $(seq 0 $progress); do printf '#'; done)
  printf "\r$format\n" $prog $perc
}
bar=0

#######################################################################################
for sample in "${arr[@]}"; do
  read -u1000
  {
    dir=$data_dir/$sample
    mkdir -p $dir/$Aligner
    cd $dir/$Aligner
    echo "===== $sample ====="

    #####
    ##### HMMcopy #####
    mkdir -p HMMcopy
    cd HMMcopy
    ####Use readCounter in HMMcopy to generate a wiggle file for each BAM file.
    #readCounter -w 1000000 ../${sample}.${Aligner}.rmdup.bam >${sample}.${Aligner}.1M_input.wig
    eval "$(conda shell.bash hook)"
    conda activate HmmCopy
    run_hmmcopy.r 0.995 ${sample}.${Aligner}.1M_input.wig "/data/database/iGenomes/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/1M_CNV/hg19_1M.gc.wig" "/data/database/iGenomes/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/1M_CNV/hg19_40mer_1M.map.wig" 1 ${sample}.${Aligner}.HMMcopy
    conda deactivate

    #####
    ##### BaseqCNV #####
    #mkdir -p $dir/baseqCNV
    #cd $dir/baseqCNV
    ###baseqCNV align -t $threads -1 $fq1 -2 $fq2 -g hg19 -c $baseqCNV_config -o ${sample}.hg19.bam
    #baseqCNV bincount -g hg19 -i $dir/${sample}.hg19.bwamem.rmdup.bam -c $baseqCNV_config -o ${sample}.bincounts
    #baseqCNV normalize -g hg19 -i ${sample}.bincounts -c $baseqCNV_config -o ${sample}.bincounts_norm.txt
    #baseqCNV cbs -i ${sample}.bincounts_norm.txt -o ${sample}.cbs.txt
    #baseqCNV plotgenome -i ${sample}.bincounts_norm.txt -c ${sample}.cbs.txt -o ${sample}.baseqCNV

    #####
    ##### Control-FREEC #####
    #mkdir -p $dir/FREEC
    #cd $dir/FREEC
    #sed "s#PE_sample.bam#$dir/${sample}.hg19.bowtie2.rmdup.bam#" $FREEC_config >FREEC_config.ini
    #sed -i "s#thread#$thread#" $FREEC_config
    #freec -conf FREEC_config.ini
    #makeGraph.R ${sample}.hg19.bowtie2.rmdup.bam_ratio.txt 2 ${sample}.FREEC

    #####
    ##### nQuire estimate ploidy #####
    #nQuire create -f 0 -b ${sample}.bowtie2.hg19.filter.rmdup.bam -x -o ${sample}.nQuire
    #nQuire denoise ${sample}.nQuire.bin -o ${sample}.nQuire.denoise
    #nQuire histo ${sample}.nQuire.denoise.bin > ${sample}.nQuire.histo.txt
    #nQuire lrdmodel -t 88 ${sample}.nQuire.denoise.bin
    #nQuire modeltest ${sample}.nQuire.denoise.bin
    #nQuire histotest ${sample}.nQuire.denoise.bin

    #####
    ##### ploidyNGS estimate ploidy ##### Need large memory
    #mkdir -p $dir/ploidyNGS
    #cd $dir/ploidyNGS
    #ploidyNGS.py --out ${sample}.ploidyNGS --bam $dir/${sample}.hg19.bwa.rmdup.bam

    ##### call SNP
    ### bcftools
    #bcftools mpileup --threads $threads -d 500 -Ou -f $genome ${sample}.bowtie2.hg19.filter.rmdup.bam | bcftools call --threads $threads -mv -Oz | bcftools view -i '%QUAL>=20' -Oz -o ${sample}.calls.vcf

    ### GATK3 #####
    rm -rf $dir/GATK3
    #mkdir -p GATK3
    ##picard CreateSequenceDictionary R=$genome
    ##gatk ValidateSamFile -M SUMMARY -I ${sample}.bowtie2.hg19.rmdup.bam
    #gatk3 -T HaplotypeCaller -Xmx30000m -nct $threads -R $genome -I $dir/${sample}.${version}.bwamem.rmdup.bam -o ${sample}.${version}.bwamem.vcf
    #bcftools filter -i 'TYPE="snp" && MIN(FORMAT/DP)>=4 && QUAL>=20' -Ov -o ${sample}.${version}.bwamem.filter.vcf ${sample}.${version}.bwamem.vcf
    #grep -Ev '^(chrY)' ${sample}.${version}.bwamem.filter.vcf > ${sample}.${version}.bwamem.filter.rmchrY.vcf
    #determining_ploidy.R ${sample}.${version}.bwamem.filter.rmchrY.vcf ${sample}.${version}

    ### Strelka #####
    #mkdir -p $dir/Strelka
    #cd $dir/Strelka
    #rm -rf ./*
    #/home/reprod/bin/pyenv/versions/anaconda2-5.1.0/bin/configureStrelkaGermlineWorkflow.py \
    #        --bam $dir/${sample}.hg19.bwamem.rmdup.bam \
    #        --referenceFasta $genome \
    #        --runDir ./
    #./runWorkflow.py -m local -j $thread
    #bcftools view results/variants/variants.vcf.gz | bcftools filter -i 'TYPE="snp" && MIN(FORMAT/DP)>=4 && QUAL>=20' -Ov -o results/variants/filter.variants.vcf
    #determining_ploidy.R results/variants/filter.variants.vcf ${sample}

    #Identify CNVs
    #1. Exclude cells for which the VS calculated in 4.4.6 exceeds 0.26.
    #2. Filter the segments called by HMMcopy in 4.4.6 to include only those for which the median log2 ratio is greater than 0.4 (putative gain)
    #or less than -0.35 (putative loss).
    #3. Filter the segments called by DNAcopy in 4.5.5 to include only those for which the segment mean is greater than 1.32 or less than 0.6.
    #4. Overlap the segments from 4.6.2 and 4.6.3, calling CNVs only in regions where both HMMcopy and DNAcopy identify a putative gain
    #or putative loss.

    echo >&1000
  } &
  ((bar++))
  processbar $bar $total_task
done
wait

echo "===== All tasks finished ====="

#### Preparation ####

###Use gcCounter in HMMcopy to generate a GC percentage reference file for the genome. Use option "-w 500000" to specify window size
#gcCounter -w 1000000 $genome >/data/database/iGenomes/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/1M_CNV/hg19_1M.gc.wig
#gcCounter -w 500000 $genome >/data/database/iGenomes/human/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/0.5M_CNV/mm10_0.5M.gc.wig
###Use generateMap.pl in HMMcopy to generate a wiggle file for mappability. Use option "-w 40" to specify read length.
#generateMap.pl -b $genome
##generateMap.pl -w 10 -i $genome $genome -o /data/database/iGenomes/human/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/2M_CNV/hg38_130mer.bigwig
##generateMap.pl -w 40 -i $genome $genome -o /data/database/iGenomes/mouse/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/0.5M_CNV/mm10_40mer.bigwig
####Use mapCounter in HMMcopy to generate a mappability reference file for the genome. Use option "-w 500000" to specify window size
#mapCounter -w 1000000 /data/database/iGenomes/human/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/1M_CNV/hg38_150mer.bigwig >$cnvdir/hg38_150mer_1M.map.wig
#mapCounter -w 2000000 /data/database/iGenomes/human/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/2M_CNV/hg38_150mer.bigwig >$cnvdir/hg38_150mer_2M.map.wig
