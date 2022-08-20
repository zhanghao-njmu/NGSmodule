#!/usr/bin/env bash
trap_add 'trap - SIGTERM && kill -- -$$' SIGINT SIGTERM

Rscript &>/dev/null
[ $? -eq 127 ] && {
    echo -e "Cannot find the command Rscript.\n"
    exit 1
}

R_packages=("SCP" "Seurat" "SeuratDisk" "dplyr" "reshape2" "ggplot2" "ggridges" "ggrepel" "ggupset")
for package in "${R_packages[@]}"; do
    Rscript -e "installed.packages()" | awk '{print $1}' | grep "$package" &>/dev/null
    [ $? -ne 0 ] && {
        color_echo "red" "Cannot find the R package $package.\n"
        exit 1
    }
done

analysis_dir=$maindir/NGSmodule_SCP_analysis/

echo -e "########################### RunCellQC Parameters ##############################\n"
echo -e "*** base parameters ***"
echo -e "    analysis_dir: ${analysis_dir}\n    samples: ${samples}\n    simple_analysis: ${simple_analysis}\n"
echo -e "*** cell-filtering parameters ***"
echo -e "    db_method: ${db_method}\n    outlier_cutoff: ${outlier_cutoff}\n    outlier_n: ${outlier_n}"
echo -e "    gene_threshold: ${gene_threshold}\n    UMI_threshold: ${UMI_threshold}"
echo -e "    mito_threshold: ${mito_threshold}\n    ribo_threshold: ${ribo_threshold}\n    ribo_mito_ratio_min: ${ribo_mito_ratio_min}\n    ribo_mito_ratio_max: ${ribo_mito_ratio_max}" 
echo -e "    species: ${species}\n    species_gene_prefix: ${species_gene_prefix}\n    species_percent: ${species_percent}"
echo -e "    exogenous_genes: ${exogenous_genes}\n    features_inspect: ${features_inspect}\n"
echo -e "################################################################################\n"

echo -e "****************** Start RunCellQC ******************\n"
SECONDS=0
SCP_path=$1

if [[ ! -e "$maindir/NGSmodule_SCP_analysis/Prepare" ]]; then
    color_echo "red" "Cannot find the prepared data.\n"
    exit 1
fi

mkdir -p $analysis_dir/CellQC
cd $analysis_dir/CellQC

force=${force_complete}
status="uncompleted"

### clear existed logs
logfiles=("RunCellQCStatus.log")
globalcheck_logfile "$analysis_dir" logfiles[@] "$force" "$error_pattern" "$complete_pattern" "RunCellQC"

check_logfile "CellQC" "RunCellQC" "$analysis_dir"/CellQC/RunCellQCStatus.log "$error_pattern" "$complete_pattern" "precheck"
if [[ $? == 1 ]]; then
    script="$SCP_path/SCP-Analysis/CellQC/RunCellQC.R"
    Rscript "$script" "${analysis_dir}" "${samples}" "${simple_analysis}" \
        "${db_method}" "${outlier_cutoff}" "${outlier_n}"\
        "${gene_threshold}" "${UMI_threshold}"  "${mito_threshold}" "${ribo_threshold}" "${ribo_mito_ratio_min}" "${ribo_mito_ratio_max}"\
        "${species}" "${species_gene_prefix}" "${species_percent}" \
        "${exogenous_genes}" "${features_inspect}" | tee RunCellQCStatus.log

    check_logfile "CellQC" "RunCellQC" "$analysis_dir"/CellQC/RunCellQCStatus.log "$error_pattern" "$complete_pattern" "postcheck"
    if [[ $? == 1 ]]; then
        status="uncompleted"
    else
        status="completed"
    fi
else
    status="completed"
fi

if [[ "$status" == "completed" ]]; then
    echo "Completed: CellQC" >>"$tmpfile"
else
    echo "Interrupted: CellQC" >>"$tmpfile"
    color_echo "red" "ERROR in CellQC! Please check the processing log."
fi

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo -e "\n$ELAPSED"
echo -e "****************** RunCellQC Done ******************\n"
