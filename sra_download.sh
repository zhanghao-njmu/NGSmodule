#pysradb srp-to-srr --detailed --desc --expand --saveto ${SRP}.tsv ${SRP}

rawdata_dir=/data/database/SRR_collection/human/germ_cell_specification/rawdata/
SRPfile=$(find $rawdata_dir -maxdepth 1 -name "meta.csv")
IFS=','

for file in ${SRPfile[@]};
do
  while IFS=$IFS read line
  do
    srp=$(awk -F "$IFS" '{print $1}' <<<"$line")
    srr=$(awk -F "$IFS" '{print $2}' <<<"$line")
    
    if [[ $srr != "run_accession" ]];then
      {
      echo $srp $srr
      prefetch --output-directory ${srp} --max-size 1000000000000 {}
      cd $rawdata_dir/$srp/$srr
      rm -rf ./fasterq.tmp*
      fasterq-dump -f --split-3 ${srr}.sra -o $srr
      if [ -e ${srr} ]; then
        mv ${srr} ${srr}.fastq
      fi
      
      ls ./ | grep -E "(*.fastq$)|(*.fq$)" | xargs -i pigz -f {}
      
      if [ -e ${srr}_2.fastq.gz ] ;then 
        echo "$srp/$srr  pair-end"
        reformat.sh in1=${srr}_1.fastq.gz in2=${srr}_2.fastq.gz vpair >/dev/null 2>&1 
        [ $? -ne 0 ] && { echo -e "ERROR:${srr}_1.fastq.gz and ${srr}_2.fastq.gz appear to have different numbers of reads!\n"; continue; } 
      else
        echo "$srp/$srr  single-end"
      fi
      }
    fi
  done < $file
done

wait 
echo "done"