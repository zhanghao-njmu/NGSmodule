#!/usr/bin/env bash

###### trap_add <command> <trap_name> ######
###### e.g. trap_add 'echo "in trap SIGINT"' SIGINT ######
log() { printf '%s\n' "$*"; }
error() { log "ERROR: $*" >&2; }
fatal() {
    error "$@"
    exit 1
}

trap_add() {
    trap_add_cmd=$1
    shift || fatal "trap_add usage error"
    for trap_add_name in "$@"; do
        trap -- "$(
            # helper fn to get existing trap command from output
            # of trap -p
            extract_trap_cmd() { printf '%s\n' "$3"; }
            # print existing trap command with newline
            eval "extract_trap_cmd $(trap -p "${trap_add_name}")"
            # print the new trap command
            printf '%s;\n' "${trap_add_cmd}"
        )" "${trap_add_name}" || fatal "unable to add to trap ${trap_add_name}"
    done
}
declare -f -t trap_add

###### check_logfile <sample> <tool> <logfile> ######
color_echo() {
    local color=$1
    local text=$2
    if [[ $color == "red" ]]; then
        echo -e "\033[31m$text\033[0m"
    elif [[ $color == "green" ]]; then
        echo -e "\033[32m$text\033[0m"
    elif [[ $color == "yellow" ]]; then
        echo -e "\033[33m$text\033[0m"
    elif [[ $color == "blue" ]]; then
        echo -e "\033[34m$text\033[0m"
    elif [[ $color == "purple" ]]; then
        echo -e "\033[35m$text\033[0m"
    fi
}

###### processbar <current> <total> <label> ######
processbar() {
    local current=$1
    local total=$2
    local label=$3
    local maxlen=60
    local barlen=50
    local format="%-${barlen}s%$((maxlen - barlen))s %s"
    local perc="[$current/$total]"
    local progress=$((current * barlen / total))
    local prog=$(for i in $(seq 0 $progress); do printf '='; done)
    printf "\r$format\n" "$prog" "$perc" "$label"
}
bar=0

###### check_logfile <sample> <tool> <logfile> <error_pattern> <complete_pattern> <mode>######
error_pattern="(error)|(fatal)|(terrible)|(corrupt)|(interrupt)|(unexpected)|(denied)|(refused)|(Failed to process)|(java.io.EOFException)|(no such file or directory)"
complete_pattern="(FastqCheck passed)|(Analysis complete)|(fastp.json)|(Processing complete)|(Done Reports generation)|(Coverage by database)|(successfully)|(All done)|(Output created)"

check_logfile() {
    local sample=$1
    local tool=$2
    local logfile=$3
    local error_pattern=$4
    local complete_pattern=$5
    local mode=$6

    if [[ -f $logfile ]]; then
        error=$(grep -iP "${error_pattern}" "${logfile}")
        complete=$(grep -iP "${complete_pattern}" "${logfile}")
        if [[ $error ]]; then
            if [[ $mode == "precheck" ]]; then
                color_echo "yellow" "Warning! ${sample}: Detected problems in ${tool} logfile: ${logfile}. Restart ${tool}."
            elif [[ $mode == "postcheck" ]]; then
                color_echo "yellow" "Warning! ${sample}: Detected problems in ${tool} logfile: ${logfile}."
            fi
            return 1
        elif [[ $complete ]]; then
            if [[ $mode == "precheck" ]]; then
                color_echo "blue" "+++++ ${sample}: ${tool} skipped +++++"
            elif [[ $mode == "postcheck" ]]; then
                color_echo "blue" "+++++ ${sample}: ${tool} done +++++"
            fi
            return 0
        else
            if [[ $mode == "precheck" ]]; then
                color_echo "yellow" "Warning! ${sample}: Unable to determine ${tool} status. Restart ${tool}."
            elif [[ $mode == "postcheck" ]]; then
                color_echo "yellow" "Warning! ${sample}: Unable to determine ${tool} status."
            fi
            return 1
        fi
    else
        if [[ $mode == "precheck" ]]; then
            color_echo "blue" "+++++ ${sample}: Start ${tool} +++++"
        elif [[ $mode == "postcheck" ]]; then
            color_echo "yellow" "Warning! ${sample}: Cannot find the log file for the tool ${tool}."
        fi
        return 1
    fi
}

################################################################################################################

work_dir=$maindir/NGSmodule_SCP_work/
if [[ ! -d $work_dir ]] && [[ $1 != "prepare" ]]; then
    color_echo "red" "Error! Can not find the work_dir: $work_dir\nPlease run 'NGSmodule CreateSCPWorkDir -c <Config_file>' first!\n"
    exit 1
fi

############# Load SampleInfoFile ###################################################################
declare -A Sample_dict
if [[ -f $SampleInfoFile ]]; then
    dos2unix $SampleInfoFile &>/dev/null
    while IFS=',' read -r RunID SampleID; do
        Sample_dict[$RunID]=$SampleID
    done <$SampleInfoFile
else
    color_echo "red" "ERROR! Cannot find SampleInfoFile: $SampleInfoFile. Please check your config!\n"
    exit 1
fi

###### START ######
if [[ -d $work_dir ]]; then

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

    threads=$(((total_threads + ntask_per_run) / ntask_per_run - 1))
    memory=$(((total_memory + ntask_per_run) / ntask_per_run - 1))

    ###### fifo ######
    tempfifo=$$.fifo
    trap_add "exec 1000>&-;exec 1000<&-;rm -f $tempfifo" SIGINT SIGTERM EXIT
    mkfifo $tempfifo
    exec 1000<>$tempfifo
    rm -f $tempfifo
    for ((i = 1; i <= ntask_per_run; i++)); do
        echo >&1000
    done

    ###### temp file ######
    tmpfile=$(mktemp /tmp/NGSmodule_SCP.XXXXXXXXXXXXXX) || exit 1
    trap_add "rm -f $tmpfile" SIGINT SIGTERM EXIT

else
    total_task="Waiting for creating the workdir"
    ntask_per_run="Waiting for creating the workdir"
    threads="Waiting for creating the workdir"
    memory="Waiting for creating the workdir"
fi

################################################################################################################
echo -e "########################### Global config patameters ###########################\n"
echo -e "  maindir: ${maindir}\n  rawdata_dir: ${rawdata_dir}\n  work_dir: ${work_dir}\n  SampleInfoFile: ${SampleInfoFile}\n  SampleGrepPattern: ${SampleGrepPattern}\n\n  Total_tasks: ${total_task}\n  nTask_per_run: ${ntask_per_run}\n  Total_threads: ${total_threads}\n  Threads_per_task: ${threads}\n  Total_memory: ${total_memory}\n  Memory_per_task: ${memory}\n"
echo -e "################################################################################\n\n\n"
