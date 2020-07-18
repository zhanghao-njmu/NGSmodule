#!/usr/bin/env bash

#######################################################################################
##### PreAlignment-QC&Trim.sh include four QC or trimming steps: 
##### 1. Do QC for the untrimmed reads with FastQC
##### 2. Trim and filter the reads with fastp
##### 3. Align trimmed reads on multi-genomes to detect the contaminantion with fastq-screen
##### 4. Remove the rRNA with SortmeRNA

trap 'j=`ps aux | grep -P "$work_dir" |grep -P "(fastqc)|(fastp)|(fastq_screen)|(reformat.sh)|(sortmerna)|(bgzip)|(pigz)|(bbmap)"| awk '"'"'{print $2}'"'"'`;kill -9 $j;kill -9 $(jobs -p);echo -e "\nKilling all background processes......\nExiting the script......\n";exit 1' SIGINT

fastqc --version &>/dev/null;[ $? -ne 0 ] && { echo -e "Cannot find the command fastqc.\n";exit 1; }
fastp --version &>/dev/null;[ $? -ne 0 ] && { echo -e "Cannot find the command fastp.\n";exit 1; }
fastq_screen --version &>/dev/null;[ $? -ne 0 ] && { echo -e "Cannot find the command fastq_screen.\n";exit 1; }
reformat.sh --version  &>/dev/null;[ $? -ne 0 ] && { echo -e "Cannot find the command reformat.sh which is a tool in BBmap.\n";exit 1; }
sortmerna --version &>/dev/null;[ $? -ne 0 ] && { echo -e "Cannot find the command sortmerna.\n";exit 1; }
pigz --version &>/dev/null;[ $? -ne 0 ] && { echo -e "Cannot find the command pigz.\n";exit 1; }
v=$(sortmerna --version 2>/dev/null |grep -oP "(?<=SortMeRNA version ).*"|cut -d. -f1)
(( v < 4 )) || [[ $v == "" ]] && { echo -e "sortmerna vertsion need to be updated to 4.x.x(https://github.com/biocore/sortmerna).\n";exit 1; }

if [[ ! -f $SortmeRNA_ref ]];then
  echo -e "Cannot find the SortmeRNA_ref file: ${SortmeRNA_ref}"
  exit 1
fi

if [[ ! -f $FastqScreen_config ]];then
  echo -e "Cannot find the FastqScreen_config file: ${FastqScreen_config}"
  exit 1
fi

echo -e "########################### preAlignmentQC Parameters ##########################\n"
echo -e "Fastp\n  trim_front1: ${trim_front1}\n  trim_tail1: ${trim_tail1}\n  trim_front2: ${trim_front2}\n  trim_tail2: ${trim_tail2}\n  qualified_quality_phred: ${qualified_quality_phred}\n  unqualified_percent_limit: ${unqualified_percent_limit}\n  read_cutting: ${read_cutting}\n  cut_window_size: ${cut_window_size}\n  cut_mean_quality: ${cut_mean_quality}\n  length_required: ${length_required}\n"
echo -e "FastqScreen\n  FastqScreen_config: ${FastqScreen_config}\n"
echo -e "SortmeRNA\n  SortmeRNA_ref: ${SortmeRNA_ref}\n"
echo -e "################################################################################\n\n\n"

echo -e "****************** Start preAlignmentQC ******************\n"
SECONDS=0

for sample in ${arr[@]}
do
  read -u1000
  {
  dir=$work_dir/$sample
	cd $dir
  mkdir -p $dir/PreAlignmentQC
  echo "+++++ $sample +++++"
  layout=${Layout_dict[$sample]}

  if [[ -f $dir/PreAlignmentQC/sortmerna/sortmerna.log ]];then
    echo -e "The last log file exist: $dir/PreAlignmentQC/sortmerna/sortmerna.log. Sample: ${sample} skipped."
  else
    if [[ $layout == "SE" ]]; then
      fq1=$dir/$(ls |grep -P "(.fastq.gz)|(.fq.gz)" | grep -Pv "(_R\d.fastq.gz)|(_R\d.fq.gz)|(_trim.fq.gz)")
      mkdir -p $dir/PreAlignmentQC/fastqc
      fastqc -o $dir/PreAlignmentQC/fastqc -t $threads ${fq1} >$dir/PreAlignmentQC/fastqc/fastqc.log 2>&1
      echo "+++++ $sample: FastQC done +++++"
      
      mkdir -p $dir/PreAlignmentQC/fastp
      fastp --thread $threads_fastp --trim_front1 $trim_front1 --trim_tail1 $trim_tail1 \
            --qualified_quality_phred $qualified_quality_phred --unqualified_percent_limit $unqualified_percent_limit \
            $read_cutting --cut_window_size $cut_window_size --cut_mean_quality $cut_mean_quality \
            --trim_poly_x --trim_poly_g --overrepresentation_analysis \
            --length_required $length_required \
            --in1 ${fq1} \
            --out1 ${sample}.fq  \
            -j $dir/PreAlignmentQC/fastp/${sample}.fastp.json \
            -h $dir/PreAlignmentQC/fastp/${sample}.fastp.html 2>$dir/PreAlignmentQC/fastp/fastp.log
      echo "+++++ $sample: Fastp done +++++"

      fq1=$dir/${sample}.fq
      
      mkdir -p $dir/PreAlignmentQC/fastq_screen
      fastq_screen  --force --aligner bowtie2 $FastqScreen_mode --conf $FastqScreen_config --threads $threads $fq1 \
                    --outdir $dir/PreAlignmentQC/fastq_screen 2>$dir/PreAlignmentQC/fastq_screen/fastq_screen.log
      echo "+++++ $sample: FastQ_Screen done +++++"
      

      rm -rf $dir/PreAlignmentQC/sortmerna_tmp
      mkdir -p $dir/PreAlignmentQC/sortmerna_tmp
      mkdir -p $dir/PreAlignmentQC/sortmerna
      sortmerna --ref ${SortmeRNA_ref} \
                --reads ${sample}.fq \
                --threads $threads \
                --workdir $dir/PreAlignmentQC/sortmerna_tmp \
                --fastx \
                --num_alignments 1 \
                --aligned aligned \
                --other other \
                -v &>$dir/PreAlignmentQC/sortmerna/sortmerna.process.log 
      size=$(du -sb other.fq | awk '{ print $1 }')
      if ! grep -i -q "error" $dir/PreAlignmentQC/sortmerna/sortmerna.process.log && ((size>1000)) ;then
        mv other.fq $dir/${sample}_trim.fq
        rm -rf ${sample}.fq aligned.fq $dir/PreAlignmentQC/sortmerna_tmp
        pigz -p $threads -f $dir/${sample}_trim.fq
        mv aligned.log $dir/PreAlignmentQC/sortmerna/sortmerna.log
        echo "+++++ $sample: SortMeRNA done +++++"
      else
        echo "+++++ !!! ${sample}: SortMeRNA corrupted !!! +++++"
      fi

      
    elif [[ $layout == "PE" ]]; then
      fq1=$dir/$(ls |grep -P "(_1.fastq.gz)|(_R1.fastq.gz)|(_1.fq.gz)|(_R1.fq.gz)" | grep -Pv "_trim.fq.gz")
      fq2=$dir/$(ls |grep -P "(_2.fastq.gz)|(_R2.fastq.gz)|(_2.fq.gz)|(_R2.fq.gz)" | grep -Pv "_trim.fq.gz")
      
      ##To verify that reads appear to be correctly paired 
      reformat.sh in1=$fq1 in2=$fq2 vpair >/dev/null 2>&1 
      [ $? -ne 0 ] && { echo -e "ERROR:$fq1 and $fq2 appear to have different numbers of reads!\n"; continue; }
      
      mkdir -p $dir/PreAlignmentQC/fastqc
      fastqc -o $dir/PreAlignmentQC/fastqc -t $threads ${fq1} ${fq2} >$dir/PreAlignmentQC/fastqc/fastqc.log 2>&1
      echo "+++++ $sample: FastQC done +++++"
      
      mkdir -p $dir/PreAlignmentQC/fastp
      fastp --thread $threads_fastp --trim_front1 $trim_front1 --trim_tail1 $trim_tail1 --trim_front2 $trim_front2 --trim_tail2 $trim_tail2 \
            --qualified_quality_phred $qualified_quality_phred --unqualified_percent_limit $unqualified_percent_limit \
            $read_cutting --cut_window_size $cut_window_size --cut_mean_quality $cut_mean_quality \
            --trim_poly_x --trim_poly_g --overrepresentation_analysis \
            --length_required $length_required --detect_adapter_for_pe --correction \
            --in1 ${fq1} --in2 ${fq2} \
            --out1 ${sample}_1_trim.fq --out2 ${sample}_2_trim.fq \
            -j $dir/PreAlignmentQC/fastp/${sample}.fastp.json \
            -h $dir/PreAlignmentQC/fastp/${sample}.fastp.html 2>$dir/PreAlignmentQC/fastp/fastp.log
      echo "+++++ $sample: Fastp done +++++"
      
      fq1=$dir/${sample}_1_trim.fq
      fq2=$dir/${sample}_2_trim.fq
      
      mkdir -p $dir/PreAlignmentQC/fastq_screen
      fastq_screen  --force --aligner bowtie2 $FastqScreen_mode --conf $FastqScreen_config --threads $threads $fq1 $fq2 \
                    --outdir $dir/PreAlignmentQC/fastq_screen 2>$dir/PreAlignmentQC/fastq_screen/fastq_screen.log
      echo "+++++ $sample: FastQ_Screen done +++++"
      

      rm -rf $dir/PreAlignmentQC/sortmerna_tmp
      mkdir -p $dir/PreAlignmentQC/sortmerna_tmp
      mkdir -p $dir/PreAlignmentQC/sortmerna
      reformat.sh in1=$fq1 in2=$fq2 out=$dir/${sample}.fq overwrite=true 2>$dir/PreAlignmentQC/sortmerna/reformat_merge.log
      sortmerna --ref ${SortmeRNA_ref} \
                --reads ${sample}.fq --paired_in \
                --threads $threads \
                --workdir $dir/PreAlignmentQC/sortmerna_tmp \
                --fastx \
                --num_alignments 1 \
                --aligned aligned \
                --other other \
                -v &>$dir/PreAlignmentQC/sortmerna/sortmerna.process.log 
      size=$(du -sb other.fq | awk '{ print $1 }')
      if ! grep -i -q "error" $dir/PreAlignmentQC/sortmerna/sortmerna.process.log && ((size>1000)) ;then
        reformat.sh in=other.fq out1=$fq1 out2=$fq2 overwrite=true 2>$dir/PreAlignmentQC/sortmerna/reformat_split.log
        rm -rf ${sample}.fq aligned.fq other.fq $dir/PreAlignmentQC/sortmerna_tmp 
        pigz -p $threads -f $fq1 $fq2 
        mv aligned.log $dir/PreAlignmentQC/sortmerna/sortmerna.log
        echo "+++++ $sample: SortMeRNA done +++++"
      else
        echo "+++++ !!! ${sample}: SortMeRNA corrupted !!! +++++"
      fi

    fi
  fi

  echo >&1000
  }&
  ((bar++))
  processbar $bar $total_task
done
wait 

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "****************** preAlignmentQC Done ******************\n"



