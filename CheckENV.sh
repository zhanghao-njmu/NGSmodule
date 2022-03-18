#!/usr/bin/env bash

#######################################################################################
trap_add 'conda deactivate' SIGINT SIGTERM EXIT

conda &>/dev/null
[ $? -eq 127 ] && {
  color_echo "red" "Cannot find the conda. Please install conda and add conda to your PATH environment variable.\n"
  exit 1
}
eval "$(conda shell.bash hook)"

ENVS=($(conda env list | awk '{print $1}' ))
if [[ " ${ENVS[*]} " != *" NGSmodule "* ]]; then
  color_echo "yellow" ">>> Create NGSmodule environment...\n"
  conda create -y -q --name "NGSmodule"
else
  color_echo "green" ">>> NGSmodule environment found.\n"
fi

conda activate NGSmodule
