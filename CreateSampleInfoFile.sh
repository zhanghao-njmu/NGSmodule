#!/usr/bin/env bash

tmp=(`date +"%Y%m%d%H%M%S"`)
cat << EOF >temp_${tmp}.Sample_info.csv
SampleID(required),SampleName(required),Group(required),Layout(required,PE/SE)
R19051073,Hom1-80S,Hom-80S,PE
R19051077,Hom2-80S,Hom-80S,PE
R19051085,Hom3-80S,Hom-80S,PE
R19051061,WT1-80S,WT-80S,PE
R19051065,WT2-80S,WT-80S,PE
R19051069,WT3-80S,WT-80S,PE
R19051072,Hom1-Input,Hom-Input,PE
R19051076,Hom2-Input,Hom-Input,PE
R19051084,Hom3-Input,Hom-Input,PE
R19051060,WT1-Input,WT-Input,PE
R19051064,WT2-Input,WT-Input,PE
R19051068,WT3-Input,WT-Input,PE
EOF
echo -e "Task finished.\nSampleInfoFile: temp_${tmp}.Sample_info.csv\n"
