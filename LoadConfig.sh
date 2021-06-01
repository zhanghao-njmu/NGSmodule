#!/usr/bin/env bash
trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

################################################################################################################
work_dir=$maindir/NGSmodule_work/
if [[ ! -d $work_dir ]] && [[ $1 != "prepare" ]]; then
  color_echo "red" "Error! Can not find the work_dir: $work_dir\nPlease run 'NGSmodule CreateWorkDir -c <Config_file>' first!\n"
  exit 1
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

############# Load SampleInfoFile ###################################################################
declare -A Sample_dict
declare -A Layout_dict
if [[ -f $SampleInfoFile ]]; then
  echo -e ">>> Find the SampleInfoFile: $SampleInfoFile\n"
  sed -i '/^$/d' $SampleInfoFile

  if [[ ! $(echo $SampleInfoFile | grep ".csv") ]]; then
    color_echo "red" "ERROR! SampleInfoFile name must end with '.csv'.\n"
    exit 1
  fi

  validation=$(awk 'BEGIN {FS=","; v = "TRUE" } NR == 1 { n = NF; next } NF != n || NF<2 { v = "FALSE"; exit }END{printf(v)}' $SampleInfoFile)
  if [[ $validation == "FALSE" ]]; then
    color_echo "red" "ERROR! Content in SampleInfoFile is not in a valid comma-separated format.\n.\n"
    exit 1
  fi

  dos2unix $SampleInfoFile &>/dev/null
  while IFS=',' read -r RunID SampleID Group Layout BatchID BatchInfo Other; do
    RunID="$(echo -e "${RunID}" | tr -d '[:space:]')"
    SampleID="$(echo -e "${SampleID}" | tr -d '[:space:]')"
    Sample_dict[$RunID]=$SampleID
    Layout_dict[$SampleID]=$Layout
  done <$SampleInfoFile
else
  color_echo "red" "ERROR! Cannot find SampleInfoFile: $SampleInfoFile. Please check your config!\n"
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

  if ((threads > 64)); then
    threads=64
  else
    threads=$threads
  fi

  if ((threads > 16)); then
    threads_fastp=16
  else
    threads_fastp=$threads
  fi

  if ((threads > 64)); then
    threads_featurecounts=64
  else
    threads_featurecounts=$threads
  fi

  if ((((threads / 8)) == 0)); then
    bismark_threads=1
  else
    bismark_threads=$((threads / 8))
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
