#!/usr/bin/env bash

work_dir=$maindir/NGSmodule_work/
if [[ ! -d $work_dir ]] && [[ $1 != "prepare" ]];then
  echo -e "Error! Can not find the work_dir: $work_dir\nPlease run 'NGSmodule PrepareWorkDir -c <Config_file>' first!"
  exit 1
fi
################################################################################################################
declare -A Species_arr=( ["human"]="Homo_sapiens" ["mouse"]="Mus_musculus" ["machin"]="Macaca_fascicularis" ["rhesus"]="Macaca_mulatta" ["fly"]="Drosophila_melanogaster" )
if [[ $SortmeRNA_ref_direct == "" ]];then
  SortmeRNA_ref="$SortmeRNA_Dir/$SortmeRNA_Type.${Species_arr[$Species]}.${SortmeRNA_DataVersion}.fa"
else
  SortmeRNA_ref=$SortmeRNA_ref_direct
fi
genome="$iGenomes_Dir/${Species_arr[$Species]}/$Database/$Genome_build/Sequence/WholeGenomeFasta/$Genome_name"
gtf="$iGenomes_Dir/${Species_arr[$Species]}/$Database/$Genome_build/Annotation/Genes/genes.gtf"
bwa_index="$iGenomes_Dir/${Species_arr[$Species]}/$Database/$Genome_build/Sequence/BWAIndex/$Genome_name"
bowtie_index="$iGenomes_Dir/${Species_arr[$Species]}/$Database/$Genome_build/Sequence/BowtieIndex/${Genome_name%%.fa}"
bowtie2_index="$iGenomes_Dir/${Species_arr[$Species]}/$Database/$Genome_build/Sequence/Bowtie2Index/${Genome_name%%.fa}"
hisat2_index="$iGenomes_Dir/${Species_arr[$Species]}/$Database/$Genome_build/Sequence/Hisat2Index/${Genome_name%%.fa}"
star_index="$iGenomes_Dir/${Species_arr[$Species]}/$Database/$Genome_build/Sequence/STARIndex/${Genome_name%%.fa}"
bismark_bowtie2_index="$iGenomes_Dir/${Species_arr[$Species]}/$Database/$Genome_build/Sequence/BismarkIndex/${Genome_name%%.fa}/bowtie2"
bismark_hisat2_index="$iGenomes_Dir/${Species_arr[$Species]}/$Database/$Genome_build/Sequence/BismarkIndex/${Genome_name%%.fa}/hisat2"
tophat2_index=$bowtie2_index
if [[ "$SequenceType" == "BSdna" ]] && [[ "$Aligner" =~ bismark_* ]];then
  FastqScreen_mode="--bisulfite"
else 
  FastqScreen_mode=""
fi


############# Load SampleInfoFile ###################################################################
declare -A Sample_dict
declare -A Layout_dict
if [[ -f $SampleInfoFile ]];then
  while IFS=',' read -r RunID SampleID Group Layout BatchID BatchInfo Other; do
      Sample_dict[$RunID]=$SampleID
      Layout_dict[$SampleID]=$Layout
  done < $SampleInfoFile
else
  echo -e "ERROR! Cannot find SampleInfoFile: $SampleInfoFile. Please check your config!\n"
  exit 1
fi


###### START ######
if [[ -d $work_dir ]];then
  arr=(` find $work_dir -mindepth 1 -maxdepth 1 -type l -o -type d -printf '%P\n' | grep -P "$SampleGrepPattern" `)
  total_task=${#arr[@]}
  if [[ "$ntask_per_run" =~ ^[0-9]+$ ]];then
    ntask_per_run=$ntask_per_run
  elif [[ "$ntask_per_run" = "ALL" ]]; then
    if (( total_task > total_threads ));then
      ntask_per_run=$total_threads
    else
      ntask_per_run=$total_task
    fi
  else
    echo "ERROR! ntask_per_run should be 'ALL' or an interger!"
    exit 1
  fi
  threads=$((($total_threads+$ntask_per_run)/$ntask_per_run-1))

  if (( threads > 120 ));then
    threads=120
  else
    threads=$threads
  fi

  if (( threads > 16 ));then
    threads_fastp=16
  else
    threads_fastp=$threads
  fi

  if (( threads > 64 ));then
    threads_featurecounts=64
  else
    threads_featurecounts=$threads
  fi

  ###### fifo ######
  tempfifo=$$.fifo
  trap "exec 1000>&-;exec 1000<&-;exit 0" 2
  mkfifo $tempfifo;exec 1000<>$tempfifo;rm -rf $tempfifo
  for ((i=1; i<=$ntask_per_run; i++));do
      echo >&1000
  done

else
  total_task="Waiting for the preparation of the workdir"
  ntask_per_run="Waiting for the preparation of the workdir"
  threads="1"
fi

###### processbar <current> <total> ###### 
processbar() {  
  local current=$1; local total=$2;  
  local maxlen=60; local barlen=50; local perclen=14;  
  local format="%-${barlen}s%$((maxlen-barlen))s"  
  local perc="[$current/$total]"  
  local progress=$((current*barlen/total))  
  local prog=$(for i in `seq 0 $progress`; do printf '='; done)  
  printf "\r$format\n" $prog $perc  
}  
bar=0

################################################################################################################
echo -e "########################### Global config patameters ###########################\n"
echo -e "  SequenceType: $SequenceType\n  maindir: ${maindir}\n  rawdata_dir: ${rawdata_dir}\n  work_dir: ${work_dir}\n  SampleInfoFile: ${SampleInfoFile}\n  SampleGrepPattern: ${SampleGrepPattern}\n\n  Total_tasks: ${total_task}\n  nTask_per_run: ${ntask_per_run}\n  Total_threads: ${total_threads}\n  Threads_per_task: ${threads} (max=120)\n\n"
echo -e "################################################################################\n\n\n"
