#!/usr/bin/env bash

tmp=(`date +"%Y%m%d%H%M%S"`)
cat << EOF >temp_${tmp}.Sample_info.csv
SampleID,SampleName,Group,Layout(PE/SE),Batch(optional)
R19051073,Hom1-enrich,Hom-enrich,PE,1st
R19051077,Hom2-enrich,Hom-enrich,PE,2nd
R19051085,Hom3-enrich,Hom-enrich,PE,2nd
R19051061,WT1-enrich,WT-enrich,PE,1st
R19051065,WT2-enrich,WT-enrich,PE,2nd
R19051069,WT3-enrich,WT-enrich,PE,2nd
R19051072,Hom1-Input,Hom-Input,PE,1st
R19051076,Hom2-Input,Hom-Input,PE,2nd
R19051084,Hom3-Input,Hom-Input,PE,2nd
R19051060,WT1-Input,WT-Input,PE,1st
R19051064,WT2-Input,WT-Input,PE,2nd
R19051068,WT3-Input,WT-Input,PE,2nd
EOF
echo -e "Task finished.\nSampleInfoFile: temp_${tmp}.Sample_info.csv\n"
