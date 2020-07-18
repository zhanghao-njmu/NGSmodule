#pysradb srp-to-srr --detailed --desc --expand --saveto ${SRP}.tsv ${SRP}

root_dir=/data/database/SRR_collection/human/early_embyro/rawdata/
SRPfile=$(find $root_dir -maxdepth 1 -name "SRP*.tsv")

for file in ${SRPfile[@]};
do
  srp=${file##*/}
  srp=${srp%%.tsv}
  srr_list=($(awk '{print $2}' $file))
  
  for srr in ${srr_list[@]}
  do
    if [[ $srr != "run_accession" ]];then
      {
      cd $root_dir
      prefetch --output-directory ${srp} --max-size 10000000000000 {}
      cd $root_dir/$srp/$srr
      rm -rf ./fasterq.tmp*
      fasterq-dump -f --split-3 ${srr}.sra -e $threads -o $srr
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
      }&
    fi
  done
done

wait

echo "Done"
