#!/usr/bin/env bash

tmp=(`date +"%Y%m%d%H%M%S"`)
cat <<- EOF >temp_${tmp}.Sample_info.csv
RunID,SampleID,Group,Layout(PE/SE),BatchID,BatchInfo,Other
R19051073,KO1,KO,PE,1,1st-experiment,none
R19051074,KO1,KO,PE,1,1st-experiment,none
R19051075,KO2,KO,PE,2,2nd-experiment,none
R19051076,KO2,KO,PE,2,2nd-experiment,none
R19051077,KO3,KO,PE,2,2nd-experiment,none
R19051078,KO3,KO,PE,2,2nd-experiment,none
R19051079,WT1,WT,PE,1,1st-experiment,none
R19051080,WT1,WT,PE,1,1st-experiment,none
R19051081,WT2,WT,PE,2,2nd-experiment,none
R19051082,WT2,WT,PE,2,2nd-experiment,none
R19051083,WT3,WT,PE,2,2nd-experiment,none
R19051084,WT3,WT,PE,2,2nd-experiment,none


EOF

echo -e "Task finished.\nSampleInfoFile: temp_${tmp}.Sample_info.csv\n"
