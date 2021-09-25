#!/usr/bin/env bash

tmp=($(date +"%Y%m%d%H%M%S"))
cat <<-EOF >temp_${tmp}.Sample_info.csv
RunID,SampleID,Group
Sample1,Brain,Brain
Sample2,Testis,Testis
Sample3,Kidney-1,Kidney
Sample4,Kidney-2,Kidney



EOF

echo -e "SampleInfoFile: temp_${tmp}.Sample_info.csv\n"
echo -e "\n****************** CreateSampleInfoFile Finished ******************\n"
