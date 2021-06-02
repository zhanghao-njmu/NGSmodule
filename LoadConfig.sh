#!/usr/bin/env bash
trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

################################################################################################################
work_dir=$maindir/NGSmodule_work/
if [[ ! -d $work_dir ]] && [[ $1 != "prepare" ]]; then
  color_echo "red" "Error! Can not find the work_dir: $work_dir\nPlease run 'NGSmodule CreateWorkDir -c <Config_file>' first!\n"
  exit 1
fi

###### START ######
if [[ -d $work_dir ]] && [[ $1 != "prepare" ]]; then

  arr=()
  while IFS='' read -r line; do
    arr+=("$line")
  done < <(find "$work_dir" -mindepth 1 -maxdepth 1 -type l -o -type d -printf '%P\n' | grep -P "$SampleGrepPattern" | sort)

  total_task=${#arr[@]}
  if [[ "$total_task" == 0 ]]; then
    color_echo "red" "ERROR! No sample sub-directory found in the work_dir:$work_dir\n"
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

  if ((threads > 48)); then
    threads=48
  else
    threads=$threads
  fi

  if ((threads > 16)); then
    threads_fastp=16
  else
    threads_fastp=$threads
  fi

  if ((threads > 48)); then
    threads_featurecounts=48
  else
    threads_featurecounts=$threads
  fi

  if ((((threads / 8)) == 0)); then
    threads_bismark=1
  else
    threads_bismark=$((threads / 8))
  fi

  types=("rna" "dna" "BSdna")
  if [[ " ${types[*]} " != *" $SequenceType "* ]]; then
    color_echo "red" "ERROR! SequenceType is wrong.\nPlease check the paramaters in your ConfigFile.\n"
    exit 1
  fi

  if [[ $SortmeRNA_ref_direct == "" ]]; then
    SortmeRNA_ref="${SortmeRNA_Dir}/${SortmeRNA_Type}.${Species}.${SortmeRNA_DataVersion}.fa"
  else
    SortmeRNA_ref=$SortmeRNA_ref_direct
  fi

  if [[ "$SequenceType" == "BSdna" ]]; then
    FastqScreen_mode="--bisulfite"
  else
    FastqScreen_mode=""
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

  ###### fifo ######
  fifo $ntask_per_run

else
  total_task="Waiting for creating the workdir"
  ntask_per_run="Waiting for creating the workdir"
  threads="Waiting for creating the workdir"
fi

################################################################################################################
echo -e "########################### Global config patameters ###########################\n"
echo -e "  SequenceType: $SequenceType\n  maindir: ${maindir}\n  rawdata_dir: ${rawdata_dir}\n  work_dir: ${work_dir}\n  SampleInfoFile: ${SampleInfoFile}\n  SampleGrepPattern: ${SampleGrepPattern}\n\n  Total_tasks: ${total_task}\n  nTask_per_run: ${ntask_per_run}\n  Total_threads: ${total_threads}\n  Threads_per_task: ${threads} (max=120)\n"
echo -e "################################################################################\n\n\n"
