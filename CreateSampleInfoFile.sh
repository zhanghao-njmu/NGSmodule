#!/usr/bin/env bash

tmp=(`date +"%Y%m%d%H%M%S"`)
cat << EOF >temp_${tmp}.Sample_info.csv
RunID,SampleID,Group,Layout(PE/SE),BatchID,BatchInfo,Other
R19051073,Hom1-enrich,Hom-enrich,PE,1,1st-experiment,none
R19051074,Hom1-enrich,Hom-enrich,PE,1,1st-experiment,none
R19051075,Hom2-enrich,Hom-enrich,PE,2,2nd-experiment,none
R19051076,Hom2-enrich,Hom-enrich,PE,2,2nd-experiment,none
R19051077,Hom3-enrich,Hom-enrich,PE,2,2nd-experiment,none
R19051078,Hom3-enrich,Hom-enrich,PE,2,2nd-experiment,none
R19051079,WT1-enrich,WT-enrich,PE,1,1st-experiment,none
R19051080,WT1-enrich,WT-enrich,PE,1,1st-experiment,none
R19051081,WT2-enrich,WT-enrich,PE,2,2nd-experiment,none
R19051082,WT2-enrich,WT-enrich,PE,2,2nd-experiment,none
R19051083,WT3-enrich,WT-enrich,PE,2,2nd-experiment,none
R19051084,WT3-enrich,WT-enrich,PE,2,2nd-experiment,none
R19051085,Hom1-Input,Hom-Input,PE,1,1st-experiment,none
R19051086,Hom1-Input,Hom-Input,PE,1,1st-experiment,none
R19051087,Hom2-Input,Hom-Input,PE,2,2nd-experiment,none
R19051088,Hom2-Input,Hom-Input,PE,2,2nd-experiment,none
R19051089,Hom3-Input,Hom-Input,PE,2,2nd-experiment,none
R19051090,Hom3-Input,Hom-Input,PE,2,2nd-experiment,none
R19051091,WT1-Input,WT-Input,PE,1,1st-experiment,none
R19051092,WT1-Input,WT-Input,PE,1,1st-experiment,none
R19051093,WT2-Input,WT-Input,PE,2,2nd-experiment,none
R19051094,WT2-Input,WT-Input,PE,2,2nd-experiment,none
R19051095,WT3-Input,WT-Input,PE,2,2nd-experiment,none
R19051096,WT3-Input,WT-Input,PE,2,2nd-experiment,none

EOF
echo -e "Task finished.\nSampleInfoFile: temp_${tmp}.Sample_info.csv\n"
