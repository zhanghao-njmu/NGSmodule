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
error_pattern="(error)|(fatal)|(terrible)|(corrupt)|(interrupt)|(denied)|(refused)|(failed)|(EOFException)|(no such file or directory)"
complete_pattern="(NGSmodule finished the job)"

check_logfile() {
    local sample=$1
    local tool=$2
    local logfile=$3
    local error_pattern=$4
    local complete_pattern=$5
    local mode=$6

    if [[ -f $logfile ]]; then
        error=$(grep -ioP "${error_pattern}" "${logfile}" | sort | uniq | paste -sd "|")
        complete=$(grep -ioP "${complete_pattern}" "${logfile}" | sort | uniq | paste -sd "|")
        if [[ $error ]]; then
            if [[ $mode == "precheck" ]]; then
                color_echo "yellow" "Warning! ${sample}: Detected problems($error) in ${tool} logfile: ${logfile}. Restart ${tool}."
                return 1
            elif [[ $mode == "postcheck" ]]; then
                color_echo "yellow" "Warning! ${sample}: Detected problems($error) in ${tool} logfile: ${logfile}."
                return 1
            fi
        elif [[ $complete ]]; then
            if [[ $mode == "precheck" ]]; then
                color_echo "blue" "+++++ ${sample}: ${tool} skipped +++++"
                return 0
            elif [[ $mode == "postcheck" ]]; then
                color_echo "blue" "+++++ ${sample}: ${tool} done +++++"
                echo -e "NGSmodule finished the job [${tool}]" >>"${logfile}"
                return 0
            fi
        else
            if [[ $mode == "precheck" ]]; then
                color_echo "yellow" "Warning! ${sample}: Unable to determine ${tool} status. Restart ${tool}."
                return 1
            elif [[ $mode == "postcheck" ]]; then
                color_echo "blue" "+++++ ${sample}: ${tool} done with no problem +++++"
                echo -e "NGSmodule finished the job [${tool}]" >>"${logfile}"
                return 0
            fi
        fi
    else
        if [[ $mode == "precheck" ]]; then
            color_echo "blue" "+++++ ${sample}: Start ${tool} +++++"
            return 1
        elif [[ $mode == "postcheck" ]]; then
            color_echo "yellow" "Warning! ${sample}: Cannot find the log file for the tool ${tool}."
            return 1
        fi
    fi
}

###### global-check existed logs ######
globalcheck_logfile() {
    local dir=$1
    local logfiles=("${!2}")
    local force=$3
    local error_pattern=$4
    local complete_pattern=$5
    local sample=$6

    echo ${logfiles[@]}
    find_par=$(printf -- " -o -name %s" "${logfiles[@]}")
    find_par=${find_par:3}

    existlogs=()
    while IFS='' read -r line; do
        existlogs+=("$line")
    done < <(find "${dir}" "$find_par")

    if ((${#existlogs[*]} >= 1)); then
        for log in "${existlogs[@]}"; do
            if [[ $(grep -iP "${error_pattern}" "${log}") ]] || [[ ! $(grep -iP "${complete_pattern}" "${log}") ]]; then
                color_echo "yellow" "Warning! ${sample}: Detected problems in logfile: ${log}."
                rm -f "${log}"
            fi
            if [[ $force == "TRUE" ]]; then
                color_echo "yellow" "Warning! ${sample}: Force to perform a complete workflow."
                rm -f "${log}"
            fi
        done
    fi
}

###### fifo ######
fifo() {
    local ntask_per_run=$1
    tempfifo=$$.fifo
    trap_add "exec 1000>&-;exec 1000<&-;rm -f $tempfifo" SIGINT SIGTERM EXIT
    mkfifo $tempfifo
    exec 1000<>$tempfifo
    rm -f $tempfifo
    for ((i = 1; i <= ntask_per_run; i++)); do
        echo >&1000
    done
}

###### temp file ######
temp_file() {
    tmpfile=$(mktemp /tmp/NGSmodule.XXXXXXXXXXXXXX) || exit 1
    trap_add "rm -f $tmpfile" SIGINT SIGTERM EXIT
}
