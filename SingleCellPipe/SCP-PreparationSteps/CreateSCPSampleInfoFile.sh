#!/usr/bin/env bash

tmp=($(date +"%Y%m%d%H%M%S"))
cat <<-EOF >temp_${tmp}.Sample_info.csv
RunID,SampleID
Sample1,Brain
Sample2,Testis
Sample3,Kidney
Sample3,Pancreas



EOF

echo -e "SampleInfoFile: temp_${tmp}.Sample_info.csv\n"
echo -e "\n****************** CreateSampleInfoFile Finished ******************\n"
