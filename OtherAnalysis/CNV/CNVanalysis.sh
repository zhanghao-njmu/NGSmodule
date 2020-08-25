#!/usr/bin/bash
trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

#######################################################################################

readCounter --help &>/dev/null
[ $? -ne 0 ] && {
  color_echo "red" "Cannot find the command readCounter. User can install it from 'https://github.com/shahcompbio/hmmcopy_utils'.\n"
  exit 1
}

GC_bin="/data/database/iGenomes/Macaca_fascicularis/UCSC/Macaca_fascicularis_5.0/Sequence-Plus_rhesus_chY/GemIndex/windows/1000000/genome_main.w1000000.gc.wig"
if [[ ! -f $GC_bin ]]; then
  color_echo "red" "ERROR! Cannot find the wig file containing GC content per bin: ${GC_bin}\n"
  exit 1
fi

Map_bin="/data/database/iGenomes/Macaca_fascicularis/UCSC/Macaca_fascicularis_5.0/Sequence-Plus_rhesus_chY/GemIndex/windows/1000000/genome_main.w1000000.130mer.gem.wig"
if [[ ! -f $Map_bin ]]; then
  color_echo "red" "ERROR! Cannot find the wig file containing average mappability per bin: ${Map_bin}\n"
  exit 1
fi

#######################################################################################
for sample in "${arr[@]}"; do
  read -u1000
  {
    dir=${work_dir}/${sample}
    mkdir -p $dir/$Aligner/CNV
    cd $dir/$Aligner/CNV

    echo "===== $sample ====="

    #####
    ##### HMMcopy #####
    mkdir -p $dir/$Aligner/CNV/HMMcopy
    cd $dir/$Aligner/CNV/HMMcopy
    readCounter -w $Window ${dir}/$Aligner/${sample}.${Aligner}.dedup.bam >${sample}.${Aligner}.w$Window.wig
    Rscript $1 0.995 ${sample}.${Aligner}.w$Window.wig $GC_bin $Map_bin $HypotheticalPloidy ${sample}.${Aligner}.HMMcopy

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
    #mkdir -p $dir/Control_FREEC
    #cd $dir/Control_FREEC
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
    #rm -rf $dir/GATK3
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

    echo >&1000
  } &
  ((bar++))
  processbar $bar $total_task
done
wait

echo "===== All tasks finished ====="
