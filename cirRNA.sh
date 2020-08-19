#!/usr/bin/bash

maindir=/data/lab/zhanghao/shz/analysis/
data_dir=$maindir/rawdata/
genome="/data/database/iGenomes/human/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
bowtieIndex="/data/database/iGenomes/human/Homo_sapiens/Ensembl/GRCh38/Sequence/BowtieIndex/Homo_sapiens.GRCh38.dna.primary_assembly"
bowtie2Index="/data/database/iGenomes/human/Homo_sapiens/Ensembl/GRCh38/Sequence/Bowtie2Index/Homo_sapiens.GRCh38.dna.primary_assembly"
gtf="/data/database/human/Homo_sapiens/Ensembl/GRCh38/Annotation/Homo_sapiens.GRCh38.94.gtf"
ref="/data/database/CIRCexplorer2/hg38/GRCh38_ref_all.txt"

total_threads=88
ntask_per_run=22

arr=` find $data_dir -mindepth 1 -maxdepth 1 -type d `
num_task=` find $data_dir -mindepth 1 -maxdepth 1 -type d |wc -l`
threads=$((($total_threads+$ntask_per_run-1)/$ntask_per_run))

#########
tempfifo=$$.fifo
trap "exec 1000>&-;exec 1000<&-;exit 0" 2
mkfifo $tempfifo
exec 1000<>$tempfifo
rm -rf $tempfifo
for ((i=1; i<=$ntask_per_run; i++))
do
    echo >&1000
done

# processbar <current> <total>  
processbar() {  
  local current=$1; local total=$2;  
  local maxlen=80; local barlen=66; local perclen=14;  
  local format="%-${barlen}s%$((maxlen-barlen))s"  
  local perc="[$current/$total]"  
  local progress=$((current*barlen/total))  
  local prog=$(for i in `seq 0 $progress`; do printf '#'; done)  
  printf "\r$format\n" $prog $perc  
}  
bar=0


eval "$(conda shell.bash hook)"
for dir in ${arr[@]}
do
  cd $dir
  read -u1000
#  {
      sample=${dir##*/}
      echo "======= Sample name: $sample ======="
      raw_fq1=$dir/${sample}_sequence_1.fastq.gz
      raw_fq2=$dir/${sample}_sequence_2.fastq.gz
      fq1=$dir/${sample}_1_trim.fq.gz
      fq2=$dir/${sample}_2_trim.fq.gz
      CIRI_trim_dir=$dir/CIRI-standard
      CIRI_dir=$dir/CIRI-full
      CIRCexplorer2_dir=$dir/CIRCexplorer2
      find_circ_dir=$dir/find_circ
      KNIFE_dir=$dir/KNIFE
      CircMaker_dir=$dir/CircMaker
      CircRNADBG_dir=$dir/CircRNADBG
      sailfishCir_dir=$dir/sailfishCir_dir
      
##### CIRI #####
### ciri-full
#      mkdir -p $CIRI_dir
#      cd $CIRI_dir
#      java -jar "/home/zhanghao/Program/circRNA/CIRI-full_v2.0/CIRI-full.jar" Pipeline \
#                -t $threads -1 $raw_fq1 -2 $raw_fq2 -a $gtf -r $genome -d $CIRI_dir -o ${sample}
#      unset DISPLAY
#      java -jar "/home/zhanghao/Program/circRNA/CIRI-full_v2.0/CIRI-vis.jar" \
#                -i $CIRI_dir/CIRI-full_output/${sample}_merge_circRNA_detail.anno \
#                -l $CIRI_dir/CIRI-AS_output/${sample}_library_length.list \
#                -r $genome -d $CIRI_dir/CIRI-vis_out -min 1 -iso 20 -o ${sample}

#### ciri-standard
#      mkdir -p $CIRI_trim_dir
#      cd $CIRI_trim_dir
#      bwa mem -T 19 -t $threads $genome $fq1 $fq2 > $CIRI_trim_dir/${sample}.sam
#      perl "/home/zhanghao/Program/circRNA/CIRI-full_v2.0/bin/CIRI2.pl" \
#           --thread_num $threads -I $CIRI_trim_dir/${sample}.sam -O $CIRI_trim_dir/${sample}.ciri -F $genome -A $gtf 
           

##### CIRCexplorer2 #####
## !!! Require tophat2.1.0 bowtie version 1.1.2 cufflinks v2.2.1 !!! ## 
#      conda activate py2
#      mkdir -p $CIRCexplorer2_dir $CIRCexplorer2_dir/align $CIRCexplorer2_dir/parse $CIRCexplorer2_dir/annotate $CIRCexplorer2_dir/assemble $CIRCexplorer2_dir/denovo 
#      cd $CIRCexplorer2_dir/align
#      unpigz -c $fq1 >tmp_R1.fq
#      unpigz -c $fq2 >tmp_R2.fq
#      tophat2 -o $CIRCexplorer2_dir/align -p $threads --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search $bowtieIndex tmp_R1.fq tmp_R2.fq
#      rm -f tmp_R1.fq tmp_R2.fq
#
#      cd $CIRCexplorer2_dir/parse
#      CIRCexplorer2 parse --pe -t TopHat-Fusion $CIRCexplorer2_dir/align/accepted_hits.bam > CIRCexplorer2_parse.log
#      cd $CIRCexplorer2_dir/annotate
#      CIRCexplorer2 annotate -r ${ref} -g ${genome} -b $CIRCexplorer2_dir/parse/back_spliced_junction.bed -o circularRNA_known.txt > CIRCexplorer2_annotate.log
#      
#      cd $CIRCexplorer2_dir
#      CIRCexplorer2 assemble -p $threads -r ${ref} -m $CIRCexplorer2_dir/align -o assemble > CIRCexplorer2_assemble.log
#      CIRCexplorer2 denovo -r ${ref} -g $genome \
#                           -b $CIRCexplorer2_dir/parse/back_spliced_junction.bed \
#                           -d $CIRCexplorer2_dir/assemble \
#                           -m $CIRCexplorer2_dir/align \
#                           --abs abs -o denovo > CIRCexplorer2_denovo.log
#      conda deactivate


##### find_circ #####
#      conda activate py2
#      mkdir -p $find_circ_dir
#      cd $find_circ_dir
#      bowtie2 -p $threads --very-sensitive --mm --score-min=C,-15,0 -x $bowtie2Index -q -1 ${fq1} -2 ${fq2} 2> bowtie2.log | samtools view -@ $threads -Shub -f4 - | samtools sort -@ $threads -o unmapped.bam
#      unmapped2anchors.py unmapped.bam | pigz > anchors.qfa.gz
#      bowtie2 -p $threads --reorder --mm --score-min=C,-15,0 -q -x $bowtie2Index -U anchors.qfa.gz \
#      | find_circ.py \
#                --genome=$genome \
#                --name=$sample \
#                --prefix=${sample}_ \
#                --stats=stats.txt \
#                --reads=spliced_reads.fa \
#                    > splice_sites.bed
#      grep CIRCULAR splice_sites.bed |grep -v chrM |awk '$5>=2'| grep UNAMBIGUOUS_BP | grep ANCHOR_UNIQUE | maxlength.py 100000 >circ_candidates.bed
#      conda deactivate


##### KNIFE #####
##### !!! cannot  parallel !!!
#      conda activate py2
#      rm -rf $KNIFE_dir
#      mkdir -p $KNIFE_dir
#      ln -fs $fq1 $KNIFE_dir/${sample}_trim_1.fq.gz
#      ln -fs $fq2 $KNIFE_dir/${sample}_trim_2.fq.gz
#      cd /home/zhanghao/Program/circRNA/KNIFE/circularRNApipeline_Standalone/
#      echo "THREADS=80" >"/home/zhanghao/Program/circRNA/KNIFE/circularRNApipeline_Standalone/config.sh"
#      bash /home/zhanghao/Program/circRNA/KNIFE/circularRNApipeline_Standalone/completeRun.sh $KNIFE_dir complete $KNIFE_dir output 13 GRCh38 
#      conda deactivate


##### CircMaker #####
##### !!! cannot  parallel !!!
#      mkdir -p $CircMaker_dir
#      cd $CircMaker_dir
#      cp "/home/zhanghao/Program/circRNA/CircMarker/CircRnaDetectDraft/CircRNA_Detect_Draft/CircRnaDetectDraft/config.ini" $CircMaker_dir/config.ini
#      config_file=$CircMaker_dir/config.ini
#      sed -i "/^Reference=/s|=.*|=$genome|" $config_file 
#      sed -i "/^GTF=/s|=.*|=$gtf|" $config_file 
#      sed -i "/^Reads1=/s|=.*|=$fq1|" $config_file 
#      sed -i "/^Reads2=/s|=.*|=$fq2|" $config_file 
#      sed -i "/^ReadsLen=/s|=.*|=140|" $config_file 
#      CircRnaDetectDraft $config_file
      
##### CircRNADBG #####
### !!! uncompleted !!!
#      mkdir -p $CircRNADBG_dir
#      cd $CircRNADBG_dir
#      cp "/home/zhanghao/Program/circRNA/CircDBG/CircDBG/CircDBG/CircRNADBG/config.ini" $CircRNADBG_dir/config.ini
#      config_file=$CircRNADBG_dir/config.ini
#      sed -i "/^Reference=/s|=.*|=$genome|" $config_file 
#      sed -i "/^GTF=/s|=.*|=$gtf|" $config_file 
#      sed -i "/^Reads1=/s|=.*|=$fq1|" $config_file 
#      sed -i "/^Reads2=/s|=.*|=$fq2|" $config_file 
#      sed -i "/^ReadLen=/s|=.*|=140|" $config_file 
#      sed -i "/^ThreadsNum=/s|=.*|=80|" $config_file 
#      CircRNADBG $config_file
      



#### sailfish-cir ####
### need make a .db file for gtf when first running
##      python 
##      import gffutils 
##      gffutils.create_db("/data/database/human/Homo_sapiens/Ensembl/GRCh38/Annotation/Homo_sapiens.GRCh38.94.gtf","/data/database/human/Homo_sapiens/Ensembl/GRCh38/Annotation/Homo_sapiens.GRCh38.94.db")
    
#      mkdir -p $sailfishCir_dir
#      unpigz -c $fq1 >tmp_R1.fq
#      unpigz -c $fq2 >tmp_R2.fq
#      python "/home/zhanghao/Program/circRNA/sailfish-cir/sailfish_cir.py" -g $genome -a $gtf -1 tmp_R1.fq -2 tmp_R2.fq -c $CIRI_trim_dir/${sample}.ciri -o $sailfishCir_dir 
#      rm -f tmp_R1.fq tmp_R2.fq





    echo >&1000
#  }&
      ((bar++))
      processbar $bar $num_task
done
wait

echo "TASK COMPLETE"

#for dir in ${arr[@]}
#do
#    sample=${dir##*/}
#    CIRI_dir=$dir/CIRI
#    cd $CIRI_dir
#    text=`grep "Number of circular RNAs found:" $sample.ciri.log`
#    echo $sample  $text
#done
    
    
